// Copyright (C) 2004-2021, 2023 Jed Brown, Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include <cassert>

#include <array>
#include <gsl/gsl_interp.h>
#include <map>
#include <numeric>
#include <petscsys.h>

#include "pism/util/ConfigInterface.hh"
#include "pism/util/Grid.hh"
#include "pism/util/error_handling.hh"
#include "pism/pism_config.hh"
#include "pism/util/Context.hh"
#include "pism/util/Logger.hh"
#include "pism/util/Vars.hh"
#include "pism/util/io/File.hh"
#include "pism/util/petscwrappers/DM.hh"
#include "pism/util/projection.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/io/IO_Flags.hh"

#if (Pism_USE_PIO == 1)
// Why do I need this???
#define _NETCDF
#include <pio.h>
#endif

namespace pism {

//! Internal structures of Grid.
struct Grid::Impl {
  Impl(std::shared_ptr<const Context> context);

  std::shared_ptr<petsc::DM> create_dm(int da_dof, int stencil_width) const;
  void set_ownership_ranges(const std::vector<unsigned int> &procs_x,
                            const std::vector<unsigned int> &procs_y);

  void compute_horizontal_coordinates();

  std::shared_ptr<const Context> ctx;

  MappingInfo mapping_info;

  // int to match types used by MPI
  int rank;
  int size;

  //! @brief array containing lenghts (in the x-direction) of processor sub-domains
  std::vector<PetscInt> procs_x;
  //! @brief array containing lenghts (in the y-direction) of processor sub-domains
  std::vector<PetscInt> procs_y;

  grid::Periodicity periodicity;

  grid::Registration registration;

  //! x-coordinates of grid points
  std::vector<double> x;
  //! y-coordinates of grid points
  std::vector<double> y;
  //! vertical grid levels in the ice; correspond to the storage grid
  std::vector<double> z;

  int xs, xm, ys, ym;
  //! horizontal grid spacing
  double dx;
  //! horizontal grid spacing
  double dy;
  //! cell area (meters^2)
  double cell_area;
  //! number of grid points in the x-direction
  unsigned int Mx;
  //! number of grid points in the y-direction
  unsigned int My;

  int max_patch_size;

  //! x-coordinate of the grid center
  double x0;
  //! y-coordinate of the grid center
  double y0;

  //! half width of the ice model grid in x-direction (m)
  double Lx;
  //! half width of the ice model grid in y-direction (m)
  double Ly;

  std::map<std::array<unsigned int, 2>, std::weak_ptr<petsc::DM> > dms;

  // This DM is used for I/O operations and is not owned by any
  // array::Array (so far, anyway). We keep a pointer to it here to
  // avoid re-allocating it many times.
  std::shared_ptr<petsc::DM> dm_scalar_global;

  //! @brief A dictionary with pointers to array::Arrays, for passing
  //! them from the one component to another (e.g. from IceModel to
  //! surface and ocean models).
  Vars variables;

  //! GSL binary search accelerator used to speed up kBelowHeight().
  gsl_interp_accel *bsearch_accel;

  //! ParallelIO I/O decompositions.
  std::map<std::array<int, 2>, int> io_decompositions;
};

Grid::Impl::Impl(std::shared_ptr<const Context> context)
    : ctx(context), mapping_info("mapping", ctx->unit_system()) {
  // empty
}

/*! @brief Initialize a uniform, shallow (3 z-levels) grid with half-widths (Lx,Ly) and Mx by My
 * nodes.
 */
std::shared_ptr<Grid> Grid::Shallow(std::shared_ptr<const Context> ctx, double Lx, double Ly,
                                    double x0, double y0, unsigned int Mx, unsigned int My,
                                    grid::Registration registration,
                                    grid::Periodicity periodicity) {
  try {
    grid::Parameters p(*ctx->config());
    p.Lx           = Lx;
    p.Ly           = Ly;
    p.x0           = x0;
    p.y0           = y0;
    p.Mx           = Mx;
    p.My           = My;
    p.registration = registration;
    p.periodicity  = periodicity;

    double Lz = ctx->config()->get_number("grid.Lz");
    p.z       = { 0.0, 0.5 * Lz, Lz };

    p.ownership_ranges_from_options(ctx->size());

    return std::make_shared<Grid>(ctx, p);
  } catch (RuntimeError &e) {
    e.add_context("initializing a shallow grid");
    throw;
  }
}

//! @brief Create a PISM distributed computational grid.
Grid::Grid(std::shared_ptr<const Context> context, const grid::Parameters &p)
    : com(context->com()), m_impl(new Impl(context)) {

  try {
    m_impl->bsearch_accel = gsl_interp_accel_alloc();
    if (m_impl->bsearch_accel == NULL) {
      throw RuntimeError(PISM_ERROR_LOCATION, "Failed to allocate a GSL interpolation accelerator");
    }

    MPI_Comm_rank(com, &m_impl->rank);
    MPI_Comm_size(com, &m_impl->size);

    p.validate();

    m_impl->Mx           = p.Mx;
    m_impl->My           = p.My;
    m_impl->x0           = p.x0;
    m_impl->y0           = p.y0;
    m_impl->Lx           = p.Lx;
    m_impl->Ly           = p.Ly;
    m_impl->registration = p.registration;
    m_impl->periodicity  = p.periodicity;
    m_impl->z            = p.z;
    m_impl->set_ownership_ranges(p.procs_x, p.procs_y);

    m_impl->compute_horizontal_coordinates();

    {
      int stencil_width = (int)context->config()->get_number("grid.max_stencil_width");

      try {
        std::shared_ptr<petsc::DM> tmp = this->get_dm(1, stencil_width);
      } catch (RuntimeError &e) {
        e.add_context("distributing a %d x %d grid across %d processors.", Mx(), My(), size());
        throw;
      }

      // hold on to a DM corresponding to dof=1, stencil_width=0 (it will
      // be needed for I/O operations)
      m_impl->dm_scalar_global = this->get_dm(1, 0);

      DMDALocalInfo info;
      PetscErrorCode ierr = DMDAGetLocalInfo(*m_impl->dm_scalar_global, &info);
      PISM_CHK(ierr, "DMDAGetLocalInfo");

      m_impl->xs = info.xs;
      m_impl->xm = info.xm;
      m_impl->ys = info.ys;
      m_impl->ym = info.ym;
    }

    int patch_size = m_impl->xm * m_impl->ym;
    GlobalMax(com, &patch_size, &m_impl->max_patch_size, 1);

  } catch (RuntimeError &e) {
    e.add_context("allocating Grid");
    throw;
  }
}

//! Create a grid from a file, get information from variable `var_name`.
static std::shared_ptr<Grid> Grid_FromFile(std::shared_ptr<const Context> ctx, const File &file,
                                           const std::string &var_name, grid::Registration r) {
  try {
    const Logger &log = *ctx->log();

    // The following call may fail because var_name does not exist. (And this is fatal!)
    // Note that this sets defaults using configuration parameters, too.
    grid::Parameters p(*ctx, file, var_name, r);

    // if we have no vertical grid information, create a fake 2-level vertical grid.
    if (p.z.size() < 2) {
      double Lz = ctx->config()->get_number("grid.Lz");
      log.message(3,
                  "WARNING: Can't determine vertical grid information using '%s' in %s'\n"
                  "         Using 2 levels and Lz of %3.3fm\n",
                  var_name.c_str(), file.filename().c_str(), Lz);

      p.z = { 0.0, Lz };
    }


    p.ownership_ranges_from_options(ctx->size());

    return std::make_shared<Grid>(ctx, p);
  } catch (RuntimeError &e) {
    e.add_context("initializing computational grid from variable \"%s\" in \"%s\"",
                  var_name.c_str(), file.filename().c_str());
    throw;
  }
}

//! Create a grid using one of variables in `var_names` in `file`.
std::shared_ptr<Grid> Grid::FromFile(std::shared_ptr<const Context> ctx,
                                     const std::string &filename,
                                     const std::vector<std::string> &var_names,
                                     grid::Registration r) {

  File file(ctx->com(), filename, io::PISM_NETCDF3, io::PISM_READONLY);

  for (const auto &name : var_names) {
    if (file.find_variable(name)) {
      return Grid_FromFile(ctx, file, name, r);
    }
  }

  throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                "file %s does not have any of %s."
                                " Cannot initialize the grid.",
                                filename.c_str(), join(var_names, ",").c_str());
}

Grid::~Grid() {
  gsl_interp_accel_free(m_impl->bsearch_accel);

#if (Pism_USE_PIO == 1)
  for (auto p : m_impl->io_decompositions) {
    int ierr = PIOc_freedecomp(m_impl->ctx->pio_iosys_id(), p.second);
    if (ierr != PIO_NOERR) {
      m_impl->ctx->log()->message(1, "Failed to de-allocate a ParallelIO decomposition");
    }
  }
#endif

  delete m_impl;
}


//! Return the index `k` into `zlevels[]` so that `zlevels[k] <= height < zlevels[k+1]` and `k < Mz`.
unsigned int Grid::kBelowHeight(double height) const {

  const double eps = 1.0e-6;
  if (height < 0.0 - eps) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "height = %5.4f is below base of ice"
                                  " (height must be non-negative)\n",
                                  height);
  }

  if (height > Lz() + eps) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "height = %5.4f is above top of computational"
                                  " grid Lz = %5.4f\n",
                                  height, Lz());
  }

  return gsl_interp_accel_find(m_impl->bsearch_accel, m_impl->z.data(), m_impl->z.size(), height);
}

//! \brief Computes the number of processors in the X- and Y-directions.
static void compute_nprocs(unsigned int Mx, unsigned int My, unsigned int size, unsigned int &Nx,
                           unsigned int &Ny) {

  if (My <= 0) {
    throw RuntimeError(PISM_ERROR_LOCATION, "'My' is invalid.");
  }

  Nx = (unsigned int)(0.5 + sqrt(((double)Mx) * ((double)size) / ((double)My)));
  Ny = 0;

  if (Nx == 0) {
    Nx = 1;
  }

  while (Nx > 0) {
    Ny = size / Nx;
    if (Nx * Ny == (unsigned int)size) {
      break;
    }
    Nx--;
  }

  if (Mx > My and Nx < Ny) {
    // Swap Nx and Ny
    auto tmp = Nx;
    Nx       = Ny;
    Ny       = tmp;
  }

  if ((Mx / Nx) < 2) { // note: integer division
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "Can't split %d grid points into %d parts (X-direction).", Mx,
                                  (int)Nx);
  }

  if ((My / Ny) < 2) { // note: integer division
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "Can't split %d grid points into %d parts (Y-direction).", My,
                                  (int)Ny);
  }
}


//! \brief Computes processor ownership ranges corresponding to equal area
//! distribution among processors.
static std::vector<unsigned int> ownership_ranges(unsigned int Mx, unsigned int Nx) {

  std::vector<unsigned int> result(Nx);

  for (unsigned int i = 0; i < Nx; i++) {
    result[i] = Mx / Nx + static_cast<unsigned int>((Mx % Nx) > i);
  }
  return result;
}

//! Set processor ownership ranges. Takes care of type conversion (`unsigned int` -> `PetscInt`).
void Grid::Impl::set_ownership_ranges(const std::vector<unsigned int> &input_procs_x,
                                      const std::vector<unsigned int> &input_procs_y) {
  if (input_procs_x.size() * input_procs_y.size() != (size_t)size) {
    throw RuntimeError(PISM_ERROR_LOCATION, "length(procs_x) * length(procs_y) != MPI size");
  }

  procs_x.resize(input_procs_x.size());
  for (unsigned int k = 0; k < input_procs_x.size(); ++k) {
    procs_x[k] = static_cast<PetscInt>(input_procs_x[k]);
  }

  procs_y.resize(input_procs_y.size());
  for (unsigned int k = 0; k < input_procs_y.size(); ++k) {
    procs_y[k] = static_cast<PetscInt>(input_procs_y[k]);
  }
}

struct OwnershipRanges {
  std::vector<unsigned int> x, y;
};

//! Compute processor ownership ranges using the grid size, MPI communicator size, and command-line
//! options `-Nx`, `-Ny`, `-procs_x`, `-procs_y`.
static OwnershipRanges compute_ownership_ranges(unsigned int Mx, unsigned int My,
                                                unsigned int size) {
  OwnershipRanges result;

  unsigned int Nx_default, Ny_default;
  compute_nprocs(Mx, My, size, Nx_default, Ny_default);

  // check -Nx and -Ny
  options::Integer Nx("-Nx", "Number of processors in the x direction", Nx_default);
  options::Integer Ny("-Ny", "Number of processors in the y direction", Ny_default);

  // validate results (compute_nprocs checks its results, but we also need to validate command-line
  // options)
  if ((Mx / Nx) < 2) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "Can't split %d grid points into %d parts (X-direction).", Mx,
                                  (int)Nx);
  }

  if ((My / Ny) < 2) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "Can't split %d grid points into %d parts (Y-direction).", My,
                                  (int)Ny);
  }

  if (Nx * Ny != (int)size) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Nx * Ny has to be equal to %d.", size);
  }

  // check -procs_x and -procs_y
  options::IntegerList procs_x("-procs_x", "Processor ownership ranges (x direction)", {});
  options::IntegerList procs_y("-procs_y", "Processor ownership ranges (y direction)", {});

  if (procs_x.is_set()) {
    if (procs_x->size() != (unsigned int)Nx) {
      throw RuntimeError(PISM_ERROR_LOCATION, "-Nx has to be equal to the -procs_x size.");
    }

    result.x.resize(procs_x->size());
    for (unsigned int k = 0; k < procs_x->size(); ++k) {
      result.x[k] = procs_x[k];
    }

  } else {
    result.x = ownership_ranges(Mx, Nx);
  }

  if (procs_y.is_set()) {
    if (procs_y->size() != (unsigned int)Ny) {
      throw RuntimeError(PISM_ERROR_LOCATION, "-Ny has to be equal to the -procs_y size.");
    }

    result.y.resize(procs_y->size());
    for (unsigned int k = 0; k < procs_y->size(); ++k) {
      result.y[k] = procs_y[k];
    }
  } else {
    result.y = ownership_ranges(My, Ny);
  }

  if (result.x.size() * result.y.size() != size) {
    throw RuntimeError(PISM_ERROR_LOCATION, "length(procs_x) * length(procs_y) != MPI size");
  }

  return result;
}

//! Compute horizontal grid spacing. See compute_horizontal_coordinates() for more.
static double compute_horizontal_spacing(double half_width, unsigned int M, bool cell_centered) {
  if (cell_centered) {
    return 2.0 * half_width / M;
  }

  return 2.0 * half_width / (M - 1);
}

//! Compute grid coordinates for one direction (X or Y).
static std::vector<double> compute_coordinates(unsigned int M, double delta, double v_min,
                                               double v_max, bool cell_centered) {
  std::vector<double> result(M);

  // Here v_min, v_max define the extent of the computational domain,
  // which is not necessarily the same thing as the smallest and
  // largest values of grid coordinates.
  if (cell_centered) {
    for (unsigned int i = 0; i < M; ++i) {
      result[i] = v_min + (i + 0.5) * delta;
    }
    result[M - 1] = v_max - 0.5 * delta;
  } else {
    for (unsigned int i = 0; i < M; ++i) {
      result[i] = v_min + i * delta;
    }
    result[M - 1] = v_max;
  }
  return result;
}

//! Compute horizontal spacing parameters `dx` and `dy` and grid coordinates using `Mx`, `My`, `Lx`, `Ly` and periodicity.
/*!
The grid used in PISM, in particular the PETSc DAs used here, are periodic in x and y.
This means that the ghosted values ` foo[i+1][j], foo[i-1][j], foo[i][j+1], foo[i][j-1]`
for all 2D Vecs, and similarly in the x and y directions for 3D Vecs, are always available.
That is, they are available even if i,j is a point at the edge of the grid.  On the other
hand, by default, `dx`  is the full width  `2 * Lx`  divided by  `Mx - 1`.
This means that we conceive of the computational domain as starting at the `i = 0`
grid location and ending at the  `i = Mx - 1`  grid location, in particular.
This idea is not quite compatible with the periodic nature of the grid.

The upshot is that if one computes in a truly periodic way then the gap between the
`i = 0`  and  `i = Mx - 1`  grid points should \em also have width  `dx`.
Thus we compute  `dx = 2 * Lx / Mx`.
 */
void Grid::Impl::compute_horizontal_coordinates() {

  bool cell_centered = registration == grid::CELL_CENTER;

  dx = compute_horizontal_spacing(Lx, Mx, cell_centered);

  dy = compute_horizontal_spacing(Ly, My, cell_centered);

  cell_area = dx * dy;

  double x_min = x0 - Lx, x_max = x0 + Lx;

  x = compute_coordinates(Mx, dx, x_min, x_max, cell_centered);

  double y_min = y0 - Ly, y_max = y0 + Ly;

  y = compute_coordinates(My, dy, y_min, y_max, cell_centered);
}

//! \brief Report grid parameters.
void Grid::report_parameters() const {

  const Logger &log      = *this->ctx()->log();
  units::System::Ptr sys = this->ctx()->unit_system();

  log.message(2, "computational domain and grid:\n");

  units::Converter km(sys, "m", "km");

  // report on grid
  log.message(2, "                grid size   %d x %d x %d\n", Mx(), My(), Mz());

  // report on computational box
  log.message(2, "           spatial domain   %.2f km x %.2f km x %.2f m\n", km(2 * Lx()),
              km(2 * Ly()), Lz());

  // report on grid cell dims
  const double one_km = 1000.0;
  if (std::min(dx(), dy()) > one_km) {
    log.message(2, "     horizontal grid cell   %.2f km x %.2f km\n", km(dx()), km(dy()));
  } else {
    log.message(2, "     horizontal grid cell   %.0f m x %.0f m\n", dx(), dy());
  }
  if (fabs(dz_max() - dz_min()) <= 1.0e-8) {
    log.message(2, "  vertical spacing in ice   dz = %.3f m (equal spacing)\n", dz_min());
  } else {
    log.message(2, "  vertical spacing in ice   uneven, %d levels, %.3f m < dz < %.3f m\n", Mz(),
                dz_min(), dz_max());
  }

  // if -verbose (=-verbose 3) then (somewhat redundantly) list parameters of grid
  {
    log.message(3, "  Grid parameters:\n");
    log.message(3, "            Lx = %6.2f km, Ly = %6.2f km, Lz = %6.2f m, \n", km(Lx()), km(Ly()),
                Lz());
    log.message(3, "            x0 = %6.2f km, y0 = %6.2f km, (coordinates of center)\n", km(x0()),
                km(y0()));
    log.message(3, "            Mx = %d, My = %d, Mz = %d, \n", Mx(), My(), Mz());
    log.message(3, "            dx = %6.3f km, dy = %6.3f km, \n", km(dx()), km(dy()));
    log.message(3, "            Nx = %d, Ny = %d\n", (int)m_impl->procs_x.size(),
                (int)m_impl->procs_y.size());

    log.message(3, "            Registration: %s\n",
                registration_to_string(m_impl->registration).c_str());
    log.message(3, "            Periodicity: %s\n",
                periodicity_to_string(m_impl->periodicity).c_str());
  }

  {
    log.message(5, "  REALLY verbose output on Grid:\n");
    log.message(5, "    vertical levels in ice (Mz=%d, Lz=%5.4f): ", Mz(), Lz());
    for (unsigned int k = 0; k < Mz(); k++) {
      log.message(5, " %5.4f, ", z(k));
    }
    log.message(5, "\n");
  }
}


//! \brief Computes indices of grid points to the lower left and upper right from (X,Y).
/*!
 * \code
 * 3       2
 * o-------o
 * |       |
 * |    +  |
 * o-------o
 * 0       1
 * \endcode
 *
 * If "+" is the point (X,Y), then (i_left, j_bottom) corresponds to
 * point "0" and (i_right, j_top) corresponds to point "2".
 *
 * Does not check if the resulting indexes are in the current
 * processor's domain. Ensures that computed indexes are within the
 * grid.
 */
void Grid::compute_point_neighbors(double X, double Y, int &i_left, int &i_right, int &j_bottom,
                                   int &j_top) const {
  i_left   = (int)floor((X - m_impl->x[0]) / m_impl->dx);
  j_bottom = (int)floor((Y - m_impl->y[0]) / m_impl->dy);

  i_right = i_left + 1;
  j_top   = j_bottom + 1;

  i_left  = std::max(i_left, 0);
  i_right = std::max(i_right, 0);

  i_left  = std::min(i_left, (int)m_impl->Mx - 1);
  i_right = std::min(i_right, (int)m_impl->Mx - 1);

  j_bottom = std::max(j_bottom, 0);
  j_top    = std::max(j_top, 0);

  j_bottom = std::min(j_bottom, (int)m_impl->My - 1);
  j_top    = std::min(j_top, (int)m_impl->My - 1);
}

std::vector<int> Grid::point_neighbors(double X, double Y) const {
  int i_left, i_right, j_bottom, j_top;
  this->compute_point_neighbors(X, Y, i_left, i_right, j_bottom, j_top);
  return { i_left, i_right, j_bottom, j_top };
}

//! \brief Compute 4 interpolation weights necessary for linear interpolation
//! from the current grid. See compute_point_neighbors for the ordering of
//! neighbors.
std::vector<double> Grid::interpolation_weights(double X, double Y) const {
  int i_left = 0, i_right = 0, j_bottom = 0, j_top = 0;
  // these values (zeros) are used when interpolation is impossible
  double alpha = 0.0, beta = 0.0;

  compute_point_neighbors(X, Y, i_left, i_right, j_bottom, j_top);

  if (i_left != i_right) {
    assert(m_impl->x[i_right] - m_impl->x[i_left] != 0.0);
    alpha = (X - m_impl->x[i_left]) / (m_impl->x[i_right] - m_impl->x[i_left]);
  }

  if (j_bottom != j_top) {
    assert(m_impl->y[j_top] - m_impl->y[j_bottom] != 0.0);
    beta = (Y - m_impl->y[j_bottom]) / (m_impl->y[j_top] - m_impl->y[j_bottom]);
  }

  return { (1.0 - alpha) * (1.0 - beta), alpha * (1.0 - beta), alpha * beta, (1.0 - alpha) * beta };
}

//! @brief Get a PETSc DM ("distributed array manager") object for given `dof` (number of degrees of
//! freedom per grid point) and stencil width.
std::shared_ptr<petsc::DM> Grid::get_dm(unsigned int dm_dof, unsigned int stencil_width) const {

  std::array<unsigned int, 2> key{ dm_dof, stencil_width };

  if (m_impl->dms[key].expired()) {
    // note: here "result" is needed because m_impl->dms is a std::map of weak_ptr
    //
    // m_impl->dms[j] = m_impl->create_dm(dm_dof, stencil_width);
    //
    // would create a shared_ptr, then assign it to a weak_ptr. At this point the
    // shared_ptr (the right hand side) will be destroyed and the corresponding weak_ptr
    // will be a nullptr.
    auto result      = m_impl->create_dm(dm_dof, stencil_width);
    m_impl->dms[key] = result;
    return result;
  }

  return m_impl->dms[key].lock();
}

//! Return grid periodicity.
grid::Periodicity Grid::periodicity() const {
  return m_impl->periodicity;
}

grid::Registration Grid::registration() const {
  return m_impl->registration;
}

//! Return execution context this grid corresponds to.
std::shared_ptr<const Context> Grid::ctx() const {
  return m_impl->ctx;
}

//! @brief Create a DM with the given number of `dof` (degrees of freedom per grid point) and
//! stencil width.
std::shared_ptr<petsc::DM> Grid::Impl::create_dm(int da_dof, int stencil_width) const {

  ctx->log()->message(3, "* Creating a DM with dof=%d and stencil_width=%d...\n", da_dof,
                      stencil_width);

  DM result;
  PetscErrorCode ierr =
      DMDACreate2d(ctx->com(), DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DMDA_STENCIL_BOX, Mx, My,
                   (PetscInt)procs_x.size(), (PetscInt)procs_y.size(), da_dof, stencil_width,
                   procs_x.data(), procs_y.data(), // lx, ly
                   &result);
  PISM_CHK(ierr, "DMDACreate2d");

#if PETSC_VERSION_GE(3, 8, 0)
  ierr = DMSetUp(result);
  PISM_CHK(ierr, "DMSetUp");
#endif

  return std::make_shared<petsc::DM>(result);
}

//! MPI rank.
int Grid::rank() const {
  return m_impl->rank;
}

//! MPI communicator size.
unsigned int Grid::size() const {
  return m_impl->size;
}

//! Dictionary of variables (2D and 3D fields) associated with this grid.
Vars &Grid::variables() {
  return m_impl->variables;
}

//! Dictionary of variables (2D and 3D fields) associated with this grid.
const Vars &Grid::variables() const {
  return m_impl->variables;
}

//! Global starting index of this processor's subset.
int Grid::xs() const {
  return m_impl->xs;
}

//! Global starting index of this processor's subset.
int Grid::ys() const {
  return m_impl->ys;
}

//! Width of this processor's sub-domain.
int Grid::xm() const {
  return m_impl->xm;
}

//! Width of this processor's sub-domain.
int Grid::ym() const {
  return m_impl->ym;
}

//! Total grid size in the X direction.
unsigned int Grid::Mx() const {
  return m_impl->Mx;
}

//! Total grid size in the Y direction.
unsigned int Grid::My() const {
  return m_impl->My;
}

//! Number of vertical levels.
unsigned int Grid::Mz() const {
  return m_impl->z.size();
}

//! X-coordinates.
const std::vector<double> &Grid::x() const {
  return m_impl->x;
}

//! Get a particular x coordinate.
double Grid::x(size_t i) const {
  return m_impl->x[i];
}

//! Y-coordinates.
const std::vector<double> &Grid::y() const {
  return m_impl->y;
}

//! Get a particular y coordinate.
double Grid::y(size_t i) const {
  return m_impl->y[i];
}

//! Z-coordinates within the ice.
const std::vector<double> &Grid::z() const {
  return m_impl->z;
}

//! Get a particular z coordinate.
double Grid::z(size_t i) const {
  return m_impl->z[i];
}

//! Horizontal grid spacing.
double Grid::dx() const {
  return m_impl->dx;
}

//! Horizontal grid spacing.
double Grid::dy() const {
  return m_impl->dy;
}

double Grid::cell_area() const {
  return m_impl->cell_area;
}

//! Minimum vertical spacing.
double Grid::dz_min() const {
  double result = m_impl->z.back();
  for (unsigned int k = 0; k < m_impl->z.size() - 1; ++k) {
    const double dz = m_impl->z[k + 1] - m_impl->z[k];
    result          = std::min(dz, result);
  }
  return result;
}

//! Maximum vertical spacing.
double Grid::dz_max() const {
  double result = 0.0;
  for (unsigned int k = 0; k < m_impl->z.size() - 1; ++k) {
    const double dz = m_impl->z[k + 1] - m_impl->z[k];
    result          = std::max(dz, result);
  }
  return result;
}

//! Half-width of the computational domain.
double Grid::Lx() const {
  return m_impl->Lx;
}

//! Half-width of the computational domain.
double Grid::Ly() const {
  return m_impl->Ly;
}

//! Height of the computational domain.
double Grid::Lz() const {
  return m_impl->z.back();
}

//! X-coordinate of the center of the domain.
double Grid::x0() const {
  return m_impl->x0;
}

//! Y-coordinate of the center of the domain.
double Grid::y0() const {
  return m_impl->y0;
}

/*!
 * Return the size of the biggest sub-domain (part owned by a MPI process)
 */
int Grid::max_patch_size() const {
  return m_impl->max_patch_size;
}


namespace grid {
//! \brief Set the vertical levels in the ice according to values in `Mz` (number of levels), `Lz`
//! (domain height), `spacing` (quadratic or equal) and `lambda` (quadratic spacing parameter).
/*!
  - When `vertical_spacing == EQUAL`, the vertical grid in the ice is equally spaced:
    `zlevels[k] = k dz` where `dz = Lz / (Mz - 1)`.
  - When `vertical_spacing == QUADRATIC`, the spacing is a quadratic function.  The intent
    is that the spacing is smaller near the base than near the top.

    In particular, if
    \f$\zeta_k = k / (\mathtt{Mz} - 1)\f$ then `zlevels[k] = Lz *`
    ((\f$\zeta_k\f$ / \f$\lambda\f$) * (1.0 + (\f$\lambda\f$ - 1.0)
    * \f$\zeta_k\f$)) where \f$\lambda\f$ = 4.  The value \f$\lambda\f$
    indicates the slope of the quadratic function as it leaves the base.
    Thus a value of \f$\lambda\f$ = 4 makes the spacing about four times finer
    at the base than equal spacing would be.
 */
std::vector<double> compute_vertical_levels(double new_Lz, unsigned int new_Mz,
                                            grid::VerticalSpacing spacing, double lambda) {

  if (new_Mz < 2) {
    throw RuntimeError(PISM_ERROR_LOCATION, "Mz must be at least 2");
  }

  if (new_Lz <= 0) {
    throw RuntimeError(PISM_ERROR_LOCATION, "Lz must be positive");
  }

  if (spacing == grid::QUADRATIC and lambda <= 0) {
    throw RuntimeError(PISM_ERROR_LOCATION, "lambda must be positive");
  }

  std::vector<double> result(new_Mz);

  // Fill the levels in the ice:
  switch (spacing) {
  case grid::EQUAL: {
    double dz = new_Lz / ((double)new_Mz - 1);

    // Equal spacing
    for (unsigned int k = 0; k < new_Mz - 1; k++) {
      result[k] = dz * ((double)k);
    }
    result[new_Mz - 1] = new_Lz; // make sure it is exactly equal
    break;
  }
  case grid::QUADRATIC: {
    // this quadratic scheme is an attempt to be less extreme in the fineness near the base.
    for (unsigned int k = 0; k < new_Mz - 1; k++) {
      const double zeta = ((double)k) / ((double)new_Mz - 1);
      result[k]         = new_Lz * ((zeta / lambda) * (1.0 + (lambda - 1.0) * zeta));
    }
    result[new_Mz - 1] = new_Lz; // make sure it is exactly equal
    break;
  }
  default:
    throw RuntimeError(PISM_ERROR_LOCATION, "spacing can not be UNKNOWN");
  }

  return result;
}

//! Convert a string to Periodicity.
Periodicity string_to_periodicity(const std::string &keyword) {
  if (keyword == "none") {
    return NOT_PERIODIC;
  }

  if (keyword == "x") {
    return X_PERIODIC;
  }

  if (keyword == "y") {
    return Y_PERIODIC;
  }

  if (keyword == "xy") {
    return XY_PERIODIC;
  }

  throw RuntimeError::formatted(PISM_ERROR_LOCATION, "grid periodicity type '%s' is invalid.",
                                keyword.c_str());
}

//! Convert Periodicity to a STL string.
std::string periodicity_to_string(Periodicity p) {
  switch (p) {
  case NOT_PERIODIC:
    return "none";
  case X_PERIODIC:
    return "x";
  case Y_PERIODIC:
    return "y";
  default:
  case XY_PERIODIC:
    return "xy";
  }
}

//! Convert an STL string to SpacingType.
VerticalSpacing string_to_spacing(const std::string &keyword) {
  if (keyword == "quadratic") {
    return QUADRATIC;
  }

  if (keyword == "equal") {
    return EQUAL;
  }

  throw RuntimeError::formatted(PISM_ERROR_LOCATION, "ice vertical spacing type '%s' is invalid.",
                                keyword.c_str());
}

//! Convert SpacingType to an STL string.
std::string spacing_to_string(VerticalSpacing s) {
  switch (s) {
  case EQUAL:
    return "equal";
  default:
  case QUADRATIC:
    return "quadratic";
  }
}

Registration string_to_registration(const std::string &keyword) {
  if (keyword == "center") {
    return CELL_CENTER;
  }

  if (keyword == "corner") {
    return CELL_CORNER;
  }

  throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid grid registration: %s",
                                keyword.c_str());
}

std::string registration_to_string(Registration registration) {
  switch (registration) {
  case CELL_CORNER:
    return "corner";
  default:
  case CELL_CENTER:
    return "center";
  }
}

void InputGridInfo::reset() {

  filename = "";

  t_len = 0;

  x0    = 0;
  Lx    = 0;

  y0    = 0;
  Ly    = 0;

  z_min = 0;
  z_max = 0;
}

void InputGridInfo::report(const Logger &log, int threshold, units::System::Ptr s) const {
  units::Converter km(s, "m", "km");

  log.message(threshold, "  x:  %5d points, [%10.3f, %10.3f] km, x0 = %10.3f km, Lx = %10.3f km\n",
              (int)this->x.size(), km(this->x0 - this->Lx), km(this->x0 + this->Lx), km(this->x0),
              km(this->Lx));

  log.message(threshold, "  y:  %5d points, [%10.3f, %10.3f] km, y0 = %10.3f km, Ly = %10.3f km\n",
              (int)this->y.size(), km(this->y0 - this->Ly), km(this->y0 + this->Ly), km(this->y0),
              km(this->Ly));

  log.message(threshold, "  z:  %5d points, [%10.3f, %10.3f] m\n", (int)this->z.size(), this->z_min,
              this->z_max);

  log.message(threshold, "  t:  %5d records\n\n", this->t_len);
}

InputGridInfo::InputGridInfo(const File &file, const std::string &variable,
                             units::System::Ptr unit_system, Registration r) {
  try {
    reset();

    filename      = file.filename();
    variable_name = variable;

    // try "variable" as the standard_name first, then as the short name:
    auto var = file.find_variable(variable, variable);

    if (not var.exists) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "variable \"%s\" is missing",
                                    variable.c_str());
    }

    auto dimensions = file.dimensions(var.name);

    bool time_dimension_processed = false;
    for (const auto &dimension_name : dimensions) {

      AxisType dimtype = file.dimension_type(dimension_name, unit_system);

      this->dimension_types[dimension_name] = dimtype;

      switch (dimtype) {
      case X_AXIS: {
        this->x      = file.read_dimension(dimension_name);
        double x_min = vector_min(this->x), x_max = vector_max(this->x);
        this->x0 = 0.5 * (x_min + x_max);
        this->Lx = 0.5 * (x_max - x_min);
        if (r == CELL_CENTER) {
          const double dx = this->x[1] - this->x[0];
          this->Lx += 0.5 * dx;
        }
        break;
      }
      case Y_AXIS: {
        this->y      = file.read_dimension(dimension_name);
        double y_min = vector_min(this->y), y_max = vector_max(this->y);
        this->y0 = 0.5 * (y_min + y_max);
        this->Ly = 0.5 * (y_max - y_min);
        if (r == CELL_CENTER) {
          const double dy = this->y[1] - this->y[0];
          this->Ly += 0.5 * dy;
        }
        break;
      }
      case Z_AXIS: {
        this->z     = file.read_dimension(dimension_name);
        this->z_min = vector_min(this->z);
        this->z_max = vector_max(this->z);
        break;
      }
      case T_AXIS: {
        if (time_dimension_processed) {
          // ignore the second, third, etc dimension interpreted as "time" and override
          // the dimension type: it is not "time" in the sense of "record dimension"
          this->dimension_types[dimension_name] = UNKNOWN_AXIS;
        } else {
          this->t_len          = file.dimension_length(dimension_name);
          time_dimension_processed = false;
        }
        break;
      }
      case UNKNOWN_AXIS:
      default: {
        // ignore unknown axes
        break;
      }
      } // switch
    }   // for loop
  } catch (RuntimeError &e) {
    e.add_context("getting grid information using variable '%s' in '%s'", variable.c_str(),
                  file.filename().c_str());
    throw;
  }
}

Parameters::Parameters(const Config &config) {
  Lx = config.get_number("grid.Lx");
  Ly = config.get_number("grid.Ly");

  x0 = 0.0;
  y0 = 0.0;

  Mx = static_cast<unsigned int>(config.get_number("grid.Mx"));
  My = static_cast<unsigned int>(config.get_number("grid.My"));

  periodicity  = string_to_periodicity(config.get_string("grid.periodicity"));
  registration = string_to_registration(config.get_string("grid.registration"));

  double Lz         = config.get_number("grid.Lz");
  unsigned int Mz   = config.get_number("grid.Mz");
  double lambda     = config.get_number("grid.lambda");
  VerticalSpacing s = string_to_spacing(config.get_string("grid.ice_vertical_spacing"));
  z                 = compute_vertical_levels(Lz, Mz, s, lambda);
  // does not set ownership ranges because we don't know if these settings are final
}

void Parameters::ownership_ranges_from_options(unsigned int size) {
  OwnershipRanges procs = compute_ownership_ranges(Mx, My, size);
  procs_x               = procs.x;
  procs_y               = procs.y;
}

void Parameters::init_from_file(const Context &ctx, const File &file,
                                const std::string &variable_name, Registration r) {
  int size = 0;
  MPI_Comm_size(ctx.com(), &size);

  InputGridInfo input_grid(file, variable_name, ctx.unit_system(), r);

  Lx           = input_grid.Lx;
  Ly           = input_grid.Ly;
  x0           = input_grid.x0;
  y0           = input_grid.y0;
  Mx           = input_grid.x.size();
  My           = input_grid.y.size();
  registration = r;
  z            = input_grid.z;
}

Parameters::Parameters(const Context &ctx, const File &file, const std::string &variable_name,
                       Registration r) {
  init_from_file(ctx, file, variable_name, r);
}

Parameters::Parameters(const Context &ctx, const std::string &filename,
                       const std::string &variable_name, Registration r) {
  File file(ctx.com(), filename, io::PISM_NETCDF3, io::PISM_READONLY);
  init_from_file(ctx, file, variable_name, r);
}


void Parameters::horizontal_size_from_options() {
  Mx = options::Integer("-Mx", "grid size in X direction", Mx);
  My = options::Integer("-My", "grid size in Y direction", My);
}

void Parameters::horizontal_extent_from_options(std::shared_ptr<units::System> unit_system) {
  // Domain size
  {
    const double km = 1000.0;
    Lx = km * options::Real(unit_system, "-Lx", "Half of the grid extent in the Y direction, in km",
                            "km", Lx / km);
    Ly = km * options::Real(unit_system, "-Ly", "Half of the grid extent in the X direction, in km",
                            "km", Ly / km);
  }

  // Alternatively: domain size and extent
  {
    options::RealList x_range("-x_range", "min,max x coordinate values", {});
    options::RealList y_range("-y_range", "min,max y coordinate values", {});

    if (x_range.is_set() and y_range.is_set()) {
      if (x_range->size() != 2 or y_range->size() != 2) {
        throw RuntimeError(PISM_ERROR_LOCATION, "-x_range and/or -y_range argument is invalid.");
      }
      x0 = (x_range[0] + x_range[1]) / 2.0;
      y0 = (y_range[0] + y_range[1]) / 2.0;
      Lx = (x_range[1] - x_range[0]) / 2.0;
      Ly = (y_range[1] - y_range[0]) / 2.0;
    }
  }
}

void Parameters::vertical_grid_from_options(Config::ConstPtr config) {
  double Lz = (not z.empty()) ? z.back() : config->get_number("grid.Lz");
  int Mz    = (not z.empty()) ? z.size() : config->get_number("grid.Mz");

  z = compute_vertical_levels(Lz, Mz,
                              string_to_spacing(config->get_string("grid.ice_vertical_spacing")),
                              config->get_number("grid.lambda"));
}

void Parameters::validate() const {
  if (Mx < 3) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "Mx = %d is invalid (has to be 3 or greater)", Mx);
  }

  if (My < 3) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "My = %d is invalid (has to be 3 or greater)", My);
  }

  if (Lx <= 0.0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Lx = %f is invalid (negative)", Lx);
  }

  if (Ly <= 0.0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Ly = %f is invalid (negative)", Ly);
  }

  if (not is_increasing(z)) {
    throw RuntimeError(PISM_ERROR_LOCATION, "z levels are not increasing");
  }

  if (z[0] > 1e-6) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "first z level is not zero: %f", z[0]);
  }

  if (z.back() < 0.0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "last z level is negative: %f", z.back());
  }

  if (std::accumulate(procs_x.begin(), procs_x.end(), 0.0) != Mx) {
    throw RuntimeError(PISM_ERROR_LOCATION, "procs_x don't sum up to Mx");
  }

  if (std::accumulate(procs_y.begin(), procs_y.end(), 0.0) != My) {
    throw RuntimeError(PISM_ERROR_LOCATION, "procs_y don't sum up to My");
  }
}

} // namespace grid

//! Create a grid using command-line options and (possibly) an input file.
/** Processes options -i, -bootstrap, -Mx, -My, -Mz, -Lx, -Ly, -Lz, -x_range, -y_range.
 */
std::shared_ptr<Grid> Grid::FromOptions(std::shared_ptr<const Context> ctx) {
  auto config = ctx->config();

  auto input_file = config->get_string("input.file");
  bool bootstrap  = config->get_flag("input.bootstrap");

  auto r = grid::string_to_registration(config->get_string("grid.registration"));

  auto log = ctx->log();

  if (not input_file.empty() and (not bootstrap)) {
    // get grid from a PISM input file
    return Grid::FromFile(ctx, input_file, { "enthalpy", "temp" }, r);
  }

  if (not input_file.empty() and bootstrap) {
    // bootstrapping; get domain size defaults from an input file, allow overriding all grid
    // parameters using command-line options

    grid::Parameters input_grid(*config);

    bool grid_info_found = false;

    File file(ctx->com(), input_file, io::PISM_NETCDF3, io::PISM_READONLY);

    for (const auto *name : { "land_ice_thickness", "bedrock_altitude", "thk", "topg" }) {

      grid_info_found = file.find_variable(name);
      if (not grid_info_found) {
        // Failed to find using a short name. Try using name as a
        // standard name...
        grid_info_found = file.find_variable("unlikely_name", name).exists;
      }

      if (grid_info_found) {
        input_grid = grid::Parameters(*ctx, file, name, r);
        break;
      }
    }

    if (not grid_info_found) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no geometry information found in '%s'",
                                    input_file.c_str());
    }

    // process all possible options controlling grid parameters, overriding values read
    // from a file
    input_grid.horizontal_size_from_options();
    input_grid.horizontal_extent_from_options(ctx->unit_system());
    input_grid.vertical_grid_from_options(config);
    input_grid.ownership_ranges_from_options(ctx->size());

    auto result = std::make_shared<Grid>(ctx, input_grid);

    units::System::Ptr sys = ctx->unit_system();
    units::Converter km(sys, "m", "km");

    // report on resulting computational box
    log->message(2,
                 "  setting computational box for ice from '%s' and\n"
                 "    user options: [%6.2f km, %6.2f km] x [%6.2f km, %6.2f km] x [0 m, %6.2f m]\n",
                 input_file.c_str(), km(result->x0() - result->Lx()),
                 km(result->x0() + result->Lx()), km(result->y0() - result->Ly()),
                 km(result->y0() + result->Ly()), result->Lz());

    return result;
  }

  {
    // This covers the two remaining cases "-i is not set, -bootstrap is set" and "-i is
    // not set, -bootstrap is not set either".

    // Use defaults from the configuration database
    grid::Parameters P(*ctx->config());
    P.horizontal_size_from_options();
    P.horizontal_extent_from_options(ctx->unit_system());
    P.vertical_grid_from_options(ctx->config());
    P.ownership_ranges_from_options(ctx->size());

    return std::make_shared<Grid>(ctx, P);
  }
}

const MappingInfo &Grid::get_mapping_info() const {
  return m_impl->mapping_info;
}

void Grid::set_mapping_info(const MappingInfo &info) {
  m_impl->mapping_info = info;
  // FIXME: re-compute lat/lon coordinates
}

/*!
 * initialize an I/O decomposition
 *
 * @param[in] dof size of the last dimension (usually z)
 * @param[in] output_datatype an integer specifying a data type (`PIO_DOUBLE`, etc)
 */
int Grid::pio_io_decomposition(int dof, int output_datatype) const {
  int result = 0;
#if (Pism_USE_PIO == 1)
  {
    std::array<int, 2> key{ dof, output_datatype };

    result = m_impl->io_decompositions[key];

    if (result == 0) {

      int ndims = dof < 2 ? 2 : 3;

      // the last element is not used if ndims == 2
      std::vector<int> gdimlen{(int)My(), (int)Mx(), dof};
      std::vector<long int> start{ys(), xs(), 0}, count{ym(), xm(), dof};

      int stat = PIOc_InitDecomp_bc(m_impl->ctx->pio_iosys_id(),
                                    output_datatype, ndims, gdimlen.data(),
                                    start.data(), count.data(), &result);
      m_impl->io_decompositions[key] = result;
      if (stat != PIO_NOERR) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "Failed to create a ParallelIO I/O decomposition");
      }
    }
  }
#else
  (void) dof;
  (void) output_datatype;
#endif
  return result;
}

PointsWithGhosts::PointsWithGhosts(const Grid &grid, unsigned int stencil_width) {
  m_i_first = grid.xs() - stencil_width;
  m_i_last  = grid.xs() + grid.xm() + stencil_width - 1;
  m_j_first = grid.ys() - stencil_width;
  m_j_last  = grid.ys() + grid.ym() + stencil_width - 1;

  m_i    = m_i_first;
  m_j    = m_j_first;
  m_done = false;
}

namespace grid {

double radius(const Grid &grid, int i, int j) {
  return sqrt(grid.x(i) * grid.x(i) + grid.y(j) * grid.y(j));
}

} // namespace grid

} // end of namespace pism
