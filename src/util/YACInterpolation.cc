#include <vector>
#include <cmath>

#include "pism/util/projection.hh"
#include "pism/util/Grid.hh"
#include "pism/util/error_handling.hh"

#include "pism/util/Context.hh"
#include "pism/util/Logger.hh"
#include "pism/util/petscwrappers/Vec.hh"
#include "pism/util/array/Scalar.hh"

#include "YACInterpolation.hh"

extern "C" {
#include "yac_interface.h"
}

namespace pism {

//! Get projection info from a NetCDF file.
static pism::MappingInfo mapping(const pism::File &file, pism::units::System::Ptr sys) {

  pism::MappingInfo result("mapping", sys);
  result.proj = file.read_text_attribute("PISM_GLOBAL", "proj");

  return result;
}

/*!
 * Grid definition using coordinates in radians.
 */
struct LonLatGrid {
  std::vector<double> lon;
  std::vector<double> lat;

  /*!
   *
   * Converts a Cartesian grid in a `projection` that uses coordinates
   * `x` and `y` in meters into the form that can be used to define a
   * curvilinear grid in YAC.
   *
   * The `projection` string has to use the format compatible with PROJ.
   */
  LonLatGrid(const std::vector<double> &x, const std::vector<double> &y,
             const std::string &projection) {

    int nrow = y.size();
    int ncol = x.size();
    int N    = nrow * ncol;

    lon.resize(N);
    lat.resize(N);

    // convert from (row, col) to the linear index in "cell" arrays
    auto C = [ncol](int row, int col) { return row * ncol + col; };

    // convert from degrees to radians
    auto deg2rad = [](double degree) { return degree * M_PI / 180; };

    pism::LonLatCalculator mapping(projection);

    for (int row = 0; row < nrow; ++row) {
      for (int col = 0; col < ncol; ++col) {
        auto coords = mapping.lonlat(x[col], y[row]);

        lon[C(row, col)] = deg2rad(coords[0]);
        lat[C(row, col)] = deg2rad(coords[1]);
      }
    }
  }
};


/*!
 * Define the PISM grid. Each PE defines its own subdomain.
 *
 * Returns the point ID that can be used to define a "field".
 */
int YACInterpolation::define_grid(const pism::Grid &grid, const std::string &grid_name,
                                  const std::string &projection) {

  if (projection.empty()) {
    throw pism::RuntimeError::formatted(
        PISM_ERROR_LOCATION, "grid '%s' has no projection information", grid_name.c_str());
  }

  int xs = grid.xs();
  int ys = grid.ys();
  int xm = grid.xm();
  int ym = grid.ym();

  std::vector<double> x(xm);
  std::vector<double> y(ym);

  // Set x and y to coordinates of cell centers:
  {
    for (int k = 0; k < xm; ++k) {
      x[k] = grid.x(xs + k);
    }
    for (int k = 0; k < ym; ++k) {
      y[k] = grid.y(ys + k);
    }
  }

  // Compute lon,lat coordinates of cell centers:
  LonLatGrid cells(x, y, projection);

  // Shift x and y by half a grid spacing and add one more row and
  // column to get coordinates of cell corners:
  {
    double dx = x[1] - x[0];
    double dy = y[1] - y[0];

    double x_last = x.back() + 0.5 * dx;
    for (int k = 0; k < x.size(); ++k) {
      x[k] -= 0.5 * dx;
    }
    x.push_back(x_last);

    double y_last = y.back() + 0.5 * dy;
    for (int k = 0; k < y.size(); ++k) {
      y[k] -= 0.5 * dy;
    }
    y.push_back(y_last);
  }

  // Compute lon,lat coordinates of cell corners:
  LonLatGrid nodes(x, y, projection);
  int n_nodes[2] = { (int)x.size(), (int)y.size() };

  std::vector<int> cell_global_index(xm * ym);
  {
    int Mx = grid.Mx();
    int k  = 0;
    for (int j = ys; j < ys + ym; ++j) {
      for (int i = xs; i < xs + xm; ++i) {
        cell_global_index[k] = j * Mx + i;
        ++k;
      }
    }
  }

  int point_id = 0;
  {
    int cyclic[] = { 0, 0 };

    int grid_id = 0;

    yac_cdef_grid_curve2d(grid_name.c_str(), n_nodes, cyclic, nodes.lon.data(), nodes.lat.data(),
                          &grid_id);

    yac_cset_global_index(cell_global_index.data(), YAC_LOCATION_CELL, grid_id);

    int n_cells[2] = { xm, ym };
    yac_cdef_points_curve2d(grid_id, n_cells, YAC_LOCATION_CELL, cells.lon.data(), cells.lat.data(),
                            &point_id);
  }
  return point_id;
}

/*!
 * Return the YAC field_id corresponding to a given PISM grid and its
 * projection info.
 *
 * @param[in] component_id YAC component ID
 * @param[in] pism_grid PISM's grid
 * @param[in] name string describing this grid and field
 */
int YACInterpolation::define_field(int component_id, const pism::Grid &pism_grid,
                                   const std::string &name) {

  int point_id = define_grid(pism_grid, name, pism_grid.get_mapping_info().proj);

  const char *time_step_length = "1";
  const int point_set_size     = 1;
  const int collection_size    = 1;

  int field_id = 0;
  yac_cdef_field(name.c_str(), component_id, &point_id, point_set_size, collection_size,
                 time_step_length, YAC_TIME_UNIT_SECOND, &field_id);
  return field_id;
}

int YACInterpolation::interpolation_coarse_to_fine(double missing_value) {
  int id = 0;
  yac_cget_interp_stack_config(&id);

  int partial_coverage = 0;

  // average over source grid nodes containing a target point, weighted using
  // barycentric local coordinates
  yac_cadd_interp_stack_config_average(id, YAC_AVG_BARY, partial_coverage);

  // nearest neighbor
  int n_neighbors = 1;
  double scaling  = 1.0;
  yac_cadd_interp_stack_config_nnn(id, YAC_NNN_DIST, n_neighbors, scaling);

  // constant if all of the above failed
  yac_cadd_interp_stack_config_fixed(id, missing_value);

  return id;
}

int YACInterpolation::interpolation_fine_to_coarse(double missing_value) {
  int id = 0;
  yac_cget_interp_stack_config(&id);

  int order                = 1;
  int enforce_conservation = 1;
  int partial_coverage     = 0;

  // conservative
  yac_cadd_interp_stack_config_conservative(id, order, enforce_conservation, partial_coverage,
                                            YAC_CONSERV_DESTAREA);

  // average over source grid nodes containing a target point, weighted using
  // barycentric local coordinates
  yac_cadd_interp_stack_config_average(id, YAC_AVG_BARY, partial_coverage);

  // nearest neighbor
  int n_neighbors = 1;
  double scaling  = 1.0;
  yac_cadd_interp_stack_config_nnn(id, YAC_NNN_DIST, n_neighbors, scaling);

  // constant if all of the above failed
  yac_cadd_interp_stack_config_fixed(id, missing_value);

  return id;
}

/*!
 * Return the string that describes a 2D grid present in a NetCDF file.
 *
 * Here `variable_name` is the name of a 2D variable used to extract
 * grid information.
 *
 * We assume that a file may contain more than one grid, so the file
 * name alone is not sufficient.
 *
 * The output has the form "input_file.nc:y:x".
 */
std::string YACInterpolation::grid_name(const pism::File &file, const std::string &variable_name,
                                        pism::units::System::Ptr sys) {
  std::string result = file.name();
  for (const auto &d : file.dimensions(variable_name)) {
    auto type = file.dimension_type(d, sys);

    if (type == pism::X_AXIS or type == pism::Y_AXIS) {
      result += ":";
      result += d;
    }
  }
  return result;
}

static void pism_yac_error_handler(MPI_Comm comm, const char *msg, const char *source, int line) {
  throw pism::RuntimeError::formatted(pism::ErrorLocation(source, line), "YAC error: %s", msg);
}

YACInterpolation::YACInterpolation(const pism::Grid &target_grid, const pism::File &file,
                                   const std::string &variable_name) {
  auto ctx = target_grid.ctx();

  yac_set_abort_handler((yac_abort_func)pism_yac_error_handler);

  try {
    auto log = ctx->log();

    auto source_grid = pism::Grid::FromFile(ctx, file, { variable_name }, pism::grid::CELL_CENTER);

    std::string source_grid_name = grid_name(file, variable_name, ctx->unit_system());

    source_grid->set_mapping_info(mapping(file, ctx->unit_system()));

    log->message(2, "Input:\n");
    source_grid->report_parameters();

    m_buffer = std::make_shared<pism::array::Scalar>(source_grid, variable_name);

    std::string target_grid_name = "internal";
    double fill_value            = NAN;
    {
      // Initialize YAC:
      {
        yac_cinit_instance(&m_instance_id);
        yac_cdef_calendar(YAC_YEAR_OF_365_DAYS);
        // Note: zero-padding of months and days *is* required.
        yac_cdef_datetime_instance(m_instance_id, "-1-01-01", "+1-01-01");
      }

      // Define components: this has to be done using *one* call
      // (cannot call yac_cdef_comp?_instance() more than once)
      const int n_comps               = 2;
      const char *comp_names[n_comps] = { "source_component", "target_component" };
      int comp_ids[n_comps]           = { 0, 0 };
      yac_cdef_comps_instance(m_instance_id, comp_names, n_comps, comp_ids);

      m_source_field_id = define_field(comp_ids[0], *source_grid, source_grid_name);
      m_target_field_id = define_field(comp_ids[1], target_grid, target_grid_name);

      // Define the interpolation stack:
      {
        std::string direction;
        int interp_stack_id = 0;
        if (source_grid->dx() < target_grid.dx() or source_grid->dy() < target_grid.dy()) {
          interp_stack_id = interpolation_fine_to_coarse(fill_value);
          direction       = "fine to coarse";
        } else {
          interp_stack_id = interpolation_coarse_to_fine(fill_value);
          direction       = "coarse to fine";
        }
        log->message(2, "Interpolation direction: %s\n", direction.c_str());

        // Define the coupling between fields:
        const int src_lag = 0;
        const int tgt_lag = 0;
        yac_cdef_couple_instance(m_instance_id,
                                 "source_component",       // source component name
                                 source_grid_name.c_str(), // source grid name
                                 source_grid_name.c_str(), // source field name
                                 "target_component",       // target component name
                                 target_grid_name.c_str(), // target grid name
                                 target_grid_name.c_str(), // target field name
                                 "1",                      // time step length in units below
                                 YAC_TIME_UNIT_SECOND,     // time step length units
                                 YAC_REDUCTION_TIME_NONE,  // reduction in time (for
                                                           // asynchronous coupling)
                                 interp_stack_id, src_lag, tgt_lag);

        // free the interpolation stack config now that we defined the coupling
        yac_cfree_interp_stack_config(interp_stack_id);
      }

      double start = MPI_Wtime();
      yac_cenddef_instance(m_instance_id);
      double end = MPI_Wtime();
      log->message(2, "Initialized interpolation from %s in %f seconds.\n",
                   source_grid_name.c_str(), end - start);
    }
  } catch (pism::RuntimeError &e) {
    e.add_context("initializing interpolation from %s to the internal grid", file.name().c_str());
    throw;
  }
}

YACInterpolation::~YACInterpolation() {
  yac_ccleanup_instance(m_instance_id);
}

double YACInterpolation::interpolate(const pism::array::Scalar &source,
                                     pism::array::Scalar &target) const {

  pism::petsc::VecArray input_array(source.vec());
  pism::petsc::VecArray output_array(target.vec());

  double *send_field_    = input_array.get();
  double **send_field[1] = { &send_field_ };

  double *recv_field[1] = { output_array.get() };

  int ierror          = 0;
  int send_info       = 0;
  int recv_info       = 0;
  int collection_size = 1;
  double start        = MPI_Wtime();
  yac_cexchange(m_source_field_id, m_target_field_id, collection_size, send_field, recv_field,
                &send_info, &recv_info, &ierror);
  double end = MPI_Wtime();

  return end - start;
}

void YACInterpolation::regrid(const pism::File &file, pism::io::Default default_value,
                              pism::array::Scalar &target) const {

  m_buffer->metadata(0) = target.metadata(0);

  // FIXME: could be m_buffer->read(...)
  m_buffer->regrid(file, default_value);

  double time_spent = interpolate(*m_buffer, target);

  auto log = target.grid()->ctx()->log();

  log->message(2, "Interpolation took %f seconds.\n", time_spent);
}

} // namespace pism
