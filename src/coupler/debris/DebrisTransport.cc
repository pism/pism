// Copyright (C) 2026 Andy Aschwanden and Constantine Khroulev
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

#include <algorithm>
#include <cmath>

#include "pism/coupler/debris/DebrisTransport.hh"

#include "pism/geometry/Geometry.hh"
#include "pism/geometry/TransportScheme.hh"
#include "pism/util/Config.hh"
#include "pism/util/Context.hh"
#include "pism/util/Grid.hh"
#include "pism/util/Logger.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/Time.hh"
#include "pism/util/array/Vector.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/File.hh"
#include "pism/util/io/IO_Flags.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/stressbalance/timestepping.hh"

namespace pism {
namespace debris {

// Helpers used to initialize sub-models from the configuration in the initializer list.
namespace {

const Config &config_of(const Grid &grid) {
  return *grid.ctx()->config();
}

double config_solid_density(const Config &config) {
  const double rho = config.get_number("debris.transport.density"),
               phi = config.get_number("debris.transport.porosity");
  if (phi < 0.0 or phi >= 1.0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "debris.transport.porosity has to be in [0, 1) (got %f)", phi);
  }
  return (1.0 - phi) * rho;
}

double config_mobility(const Config &config) {
  return config.get_number("debris.transport.mobility", "Pa-1 m2 s-1");
}

std::shared_ptr<TransportScheme3D> scheme_3d(std::shared_ptr<const Grid> grid, const Config &config) {
  return TransportScheme3D::create(grid, config.get_string("debris.transport.englacial.scheme"),
                                   (int)config.get_number("debris.transport.mpdata.iterations"),
                                   config.get_flag("debris.transport.mpdata.nonoscillatory"));
}

std::shared_ptr<TransportScheme2D> scheme_2d(std::shared_ptr<const Grid> grid, const Config &config) {
  return TransportScheme2D::create(grid, config.get_string("debris.transport.supraglacial.scheme"),
                                   (int)config.get_number("debris.transport.mpdata.iterations"),
                                   config.get_flag("debris.transport.mpdata.nonoscillatory"));
}

} // end of anonymous namespace

DebrisTransport::DebrisTransport(std::shared_ptr<const Grid> grid)
  : DebrisModel(grid),
    m_solid_density(config_solid_density(config_of(*grid))),
    m_vertical_cfl_ratio(config_of(*grid).get_number("debris.transport.vertical_cfl_ratio")),
    m_cover_coefficient(config_of(*grid).get_number("debris.transport.cover_fraction_coefficient")),
    m_do_gravity(config_of(*grid).get_flag("debris.transport.gravitational_transport")),
    m_do_terminus(config_of(*grid).get_flag("debris.transport.terminus_removal")),
    m_debris_thickness(grid, "debris_thickness"),
    m_concentration(grid, "englacial_debris_concentration", array::WITHOUT_GHOSTS, grid->z()),
    m_input(grid),
    m_englacial(grid, scheme_3d(grid, config_of(*grid)), m_solid_density),
    m_supraglacial(grid, scheme_2d(grid, config_of(*grid))),
    m_gravity(grid, config_mobility(config_of(*grid)), m_solid_density,
              config_of(*grid).get_number("debris.transport.diffusion_cfl_ratio")),
    m_terminus(grid, config_of(*grid).get_number("debris.transport.marginal_length_scale"),
               config_mobility(config_of(*grid)), m_solid_density),
    m_H_old(grid, "debris_ice_thickness_old"),
    m_cell_type_old(grid, "debris_cell_type_old"),
    m_surface_input(grid, "debris_surface_input"),
    m_surface_loss(grid, "debris_surface_loss"),
    m_cumulative_input(grid, "debris_cumulative_input"),
    m_cumulative_melt_out(grid, "debris_cumulative_melt_out"),
    m_cumulative_removal(grid, "debris_cumulative_removal"),
    m_mass_englacial(0.0),
    m_mass_supraglacial(0.0),
    m_mass_initial(0.0),
    m_step_input(0.0),
    m_step_melt_out(0.0),
    m_step_output(0.0),
    m_step_lost(0.0),
    m_cumulative_input_mass(0.0),
    m_cumulative_output_mass(0.0),
    m_cumulative_lost_mass(0.0) {

  m_debris_thickness.metadata(0)
      .long_name("thickness of the supraglacial debris layer")
      .units("m");
  m_debris_thickness.metadata(0)["valid_min"] = { 0.0 };
  m_debris_thickness.set(0.0);

  m_concentration.metadata(0)
      .long_name("englacial debris mass concentration")
      .units("kg m^-3");
  m_concentration.metadata(0)["valid_min"] = { 0.0 };
  m_concentration.set(0.0);

  m_surface_input.metadata(0)
      .long_name("debris added to the supraglacial layer from the terrain during the last step")
      .units("m");
  m_surface_loss.metadata(0)
      .long_name("supraglacial debris lost in ice-free cells during the last step")
      .units("m");

  m_cumulative_input.metadata(0)
      .long_name("cumulative debris input from the surrounding terrain (solid debris thickness)")
      .units("m");
  m_cumulative_melt_out.metadata(0)
      .long_name("cumulative melt-out of englacial debris (solid debris thickness)")
      .units("m");
  m_cumulative_removal.metadata(0)
      .long_name("cumulative removal of supraglacial debris into the foreland")
      .units("m");

  for (auto *v : { &m_surface_input, &m_surface_loss, &m_cumulative_input,
                   &m_cumulative_melt_out, &m_cumulative_removal }) {
    v->set(0.0);
  }
  m_H_old.set(0.0);
}

const EnglacialTransport &DebrisTransport::englacial() const {
  return m_englacial;
}

const SupraglacialTransport &DebrisTransport::supraglacial() const {
  return m_supraglacial;
}

const GravitationalTransport &DebrisTransport::gravitational() const {
  return m_gravity;
}

const TerminusRemoval &DebrisTransport::terminus() const {
  return m_terminus;
}

const DebrisInput &DebrisTransport::input() const {
  return m_input;
}

const array::Array3D &DebrisTransport::concentration() const {
  return m_concentration;
}

const array::Scalar &DebrisTransport::cumulative_input() const {
  return m_cumulative_input;
}

const array::Scalar &DebrisTransport::cumulative_melt_out() const {
  return m_cumulative_melt_out;
}

const array::Scalar &DebrisTransport::cumulative_removal() const {
  return m_cumulative_removal;
}

double DebrisTransport::solid_density() const {
  return m_solid_density;
}

double DebrisTransport::cover_fraction_coefficient() const {
  return m_cover_coefficient;
}

double DebrisTransport::englacial_mass() const {
  return m_mass_englacial;
}

double DebrisTransport::supraglacial_mass() const {
  return m_mass_supraglacial;
}

double DebrisTransport::total_mass() const {
  return m_mass_englacial + m_mass_supraglacial;
}

double DebrisTransport::initial_mass() const {
  return m_mass_initial;
}

double DebrisTransport::input_last_step() const {
  return m_step_input;
}

double DebrisTransport::melt_out_last_step() const {
  return m_step_melt_out;
}

double DebrisTransport::output_last_step() const {
  return m_step_output;
}

double DebrisTransport::lost_last_step() const {
  return m_step_lost;
}

double DebrisTransport::cumulative_input_mass() const {
  return m_cumulative_input_mass;
}

double DebrisTransport::cumulative_output_mass() const {
  return m_cumulative_output_mass;
}

double DebrisTransport::cumulative_lost_mass() const {
  return m_cumulative_lost_mass;
}

double DebrisTransport::conservation_error() const {
  return (m_mass_englacial + m_mass_supraglacial) - m_mass_initial - m_cumulative_input_mass +
         m_cumulative_output_mass + m_cumulative_lost_mass;
}

void DebrisTransport::init_impl(const Geometry &geometry) {

  m_log->message(2,
                 "* Initializing the debris transport model (englacial and supraglacial\n"
                 "  debris; Verhaegen and Huybrechts, 2026)...\n");

  auto opts = process_input_options(m_grid->com, m_config);

  if (opts.type == INIT_RESTART) {
    m_log->message(2, "  - Reading the debris state from '%s'...\n", opts.filename.c_str());

    File file(m_grid->com, opts.filename, io::PISM_GUESS, io::PISM_READONLY);

    m_debris_thickness.read(file, opts.record);
    m_concentration.read(file, opts.record);
  } else {
    m_debris_thickness.set(0.0);
    m_concentration.set(0.0);

    if (opts.type == INIT_BOOTSTRAP) {
      m_log->message(2, "  - Bootstrapping: using 'debris_thickness' and"
                        " 'englacial_debris_concentration' from '%s' if present\n"
                        "    (zero otherwise)...\n", opts.filename.c_str());

      File file(m_grid->com, opts.filename, io::PISM_GUESS, io::PISM_READONLY);

      m_debris_thickness.regrid(file, io::Default(0.0));
      m_concentration.regrid(file, io::Default(0.0));
    }
  }

  // -regrid_file
  regrid("debris transport", m_debris_thickness, REGRID_WITHOUT_REGRID_VARS);
  regrid("debris transport", m_concentration, REGRID_WITHOUT_REGRID_VARS);

  m_englacial.set_concentration(m_concentration, geometry.ice_thickness);
  m_debris_thickness.update_ghosts();

  m_input.init();

  m_H_old.copy_from(geometry.ice_thickness);
  m_cell_type_old.copy_from(geometry.cell_type);

  for (auto *v : { &m_surface_input, &m_surface_loss, &m_cumulative_input,
                   &m_cumulative_melt_out, &m_cumulative_removal }) {
    v->set(0.0);
  }

  m_log->message(2,
                 "  - Schemes: englacial '%s', supraglacial '%s' (MPDATA passes: %d, FCT: %s)\n"
                 "  - Gravitational transport: %s; removal at the margin: %s"
                 " (%d cells up-glacier)\n",
                 m_config->get_string("debris.transport.englacial.scheme").c_str(),
                 m_config->get_string("debris.transport.supraglacial.scheme").c_str(),
                 (int)m_config->get_number("debris.transport.mpdata.iterations"),
                 m_config->get_flag("debris.transport.mpdata.nonoscillatory") ? "on" : "off",
                 m_do_gravity ? "on" : "off", m_do_terminus ? "on" : "off",
                 m_terminus.upstream_cells());

  // mass budget
  {
    const double cell_area = m_grid->cell_area();
    double local = 0.0;
    {
      array::AccessScope list{ &m_debris_thickness };
      for (auto p : m_grid->points()) {
        local += m_debris_thickness(p.i(), p.j());
      }
    }
    m_mass_supraglacial = GlobalSum(m_grid->com, local) * cell_area * m_solid_density;
    m_mass_englacial    = GlobalSum(m_grid->com, m_englacial.local_mass());
    m_mass_initial      = m_mass_englacial + m_mass_supraglacial;

    m_step_input = m_step_melt_out = m_step_output = m_step_lost = 0.0;
    m_cumulative_input_mass = m_cumulative_output_mass = m_cumulative_lost_mass = 0.0;
  }
}

void DebrisTransport::update_impl(const Inputs &inputs, double t, double dt) {
  inputs.check_transport();

  const Geometry &geometry = *inputs.geometry;
  const array::Scalar &H = geometry.ice_thickness;
  const array::CellType2 &cell_type = geometry.cell_type;

  m_input.update(t, dt);
  const array::Scalar &rate = m_input.rate();

  // 1. englacial transport, burial and melt-out
  m_englacial.step(dt, m_H_old, m_cell_type_old, H, cell_type,
                   *inputs.top_surface_mass_balance, *inputs.bottom_surface_mass_balance, rate,
                   *inputs.u3, *inputs.v3, *inputs.w3);

  // 2. sources of the supraglacial layer: melt-out (equation 19) and direct input in the
  // ablation zone (equation 18)
  {
    const array::Scalar &melt_out = m_englacial.melt_out(), &burial = m_englacial.burial(),
                        &top_smb = *inputs.top_surface_mass_balance;

    array::AccessScope list{ &H, &cell_type, &top_smb, &rate, &melt_out, &burial,
                             &m_debris_thickness, &m_surface_input, &m_cumulative_input,
                             &m_cumulative_melt_out };

    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();

      m_surface_input(i, j) = 0.0;

      if (not cell_type.icy(i, j)) {
        continue;
      }

      const double melted = melt_out(i, j) / m_solid_density; // m of solid debris
      m_debris_thickness(i, j) += melted;
      m_cumulative_melt_out(i, j) += melted;

      // debris arriving where there is no accumulation stays on the surface
      const double input = std::max(rate(i, j), 0.0) * dt;
      if (top_smb(i, j) <= 0.0 and input > 0.0) {
        m_debris_thickness(i, j) += input;
        m_surface_input(i, j) = input;
      }

      m_cumulative_input(i, j) += m_surface_input(i, j) + burial(i, j) / m_solid_density;
    }
  }
  m_debris_thickness.update_ghosts();

  // 3. advection by the surface velocity
  m_supraglacial.update_surface_velocity(*inputs.u3, *inputs.v3, H);
  m_supraglacial.step(dt, cell_type, m_debris_thickness);
  m_debris_thickness.update_ghosts();

  // 4. gravitational redistribution
  if (m_do_gravity) {
    m_gravity.step(dt, cell_type, geometry.ice_surface_elevation, m_debris_thickness);
  }

  // 5. removal at the margin
  if (m_do_terminus) {
    m_terminus.step(dt, cell_type, geometry.ice_surface_elevation, m_debris_thickness);

    array::AccessScope list{ &m_cumulative_removal, &m_terminus.removal_rate() };
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();
      m_cumulative_removal(i, j) += m_terminus.removal_rate()(i, j) * dt;
    }
  }

  // 6. debris in ice-free cells is lost
  {
    array::AccessScope list{ &cell_type, &m_debris_thickness, &m_surface_loss };
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();
      if (cell_type.icy(i, j)) {
        m_surface_loss(i, j) = 0.0;
      } else {
        m_surface_loss(i, j) = m_debris_thickness(i, j);
        m_debris_thickness(i, j) = 0.0;
      }
    }
  }
  m_debris_thickness.update_ghosts();

  update_budget(geometry, dt);

  // remember the geometry this velocity field corresponds to
  m_H_old.copy_from(H);
  m_cell_type_old.copy_from(cell_type);

  m_englacial.concentration(H, m_concentration);
}

void DebrisTransport::update_budget(const Geometry &geometry, double dt) {
  (void) dt;

  const double cell_area = m_grid->cell_area();

  double input = 0.0, melt_out = 0.0, output = 0.0, lost = 0.0, supraglacial = 0.0;
  {
    const array::Scalar &burial = m_englacial.burial(), &melted = m_englacial.melt_out(),
                        &basal = m_englacial.basal_loss(), &ice_free = m_englacial.ice_free_loss(),
                        &removal = m_terminus.removal_rate();

    array::AccessScope list{ &burial, &melted, &basal, &ice_free, &removal, &m_surface_input,
                             &m_surface_loss, &m_debris_thickness };

    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();

      input    += burial(i, j) + m_surface_input(i, j) * m_solid_density;
      melt_out += melted(i, j);
      output   += m_do_terminus ? removal(i, j) * dt * m_solid_density : 0.0;
      lost     += basal(i, j) + ice_free(i, j) + m_surface_loss(i, j) * m_solid_density;
      supraglacial += m_debris_thickness(i, j) * m_solid_density;
    }
  }

  double local[5] = { input, melt_out, output, lost, supraglacial }, global[5];
  GlobalSum(m_grid->com, local, global, 5);

  m_step_input    = global[0] * cell_area;
  m_step_melt_out = global[1] * cell_area;
  m_step_output   = global[2] * cell_area;
  m_step_lost     = global[3] * cell_area;

  m_mass_supraglacial = global[4] * cell_area;
  m_mass_englacial    = GlobalSum(m_grid->com, m_englacial.local_mass());

  m_cumulative_input_mass  += m_step_input;
  m_cumulative_output_mass += m_step_output;
  m_cumulative_lost_mass   += m_step_lost;

  (void) geometry;
}

MaxTimestep DebrisTransport::max_timestep_impl(double t, const CFLData *cfl_data) const {
  (void) t;

  if (cfl_data == nullptr) {
    return {};
  }

  double dt = cfl_data->dt_max.value();

  // explicit vertical CFL: the thinnest control volume is half the smallest spacing
  if (cfl_data->w_max > 0.0) {
    dt = std::min(dt, m_vertical_cfl_ratio * 0.5 * m_grid->dz_min() / cfl_data->w_max);
  }

  return MaxTimestep(dt, "debris transport");
}

const array::Scalar &DebrisTransport::debris_impl() const {
  return m_debris_thickness;
}

std::set<VariableMetadata> DebrisTransport::state_impl() const {
  return array::metadata({ &m_debris_thickness, &m_concentration });
}

void DebrisTransport::write_state_impl(const OutputFile &output) const {
  m_debris_thickness.write(output);
  m_concentration.write(output);
}

void DebrisTransport::begin_pointwise_access_impl() const {
  m_debris_thickness.begin_access();
}

void DebrisTransport::end_pointwise_access_impl() const {
  m_debris_thickness.end_access();
}

void DebrisTransport::init_timeseries_impl(const std::vector<double> &ts) const {
  m_ts_times = ts;
}

void DebrisTransport::debris_time_series_impl(int i, int j, std::vector<double> &result) const {
  for (size_t k = 0; k < m_ts_times.size(); ++k) {
    result[k] = m_debris_thickness(i, j);
  }
}

// diagnostics

namespace diagnostics {

/*! @brief Englacial debris mass integrated over the column. */
class ColumnMass : public Diag<DebrisTransport> {
public:
  ColumnMass(const DebrisTransport *m) : Diag<DebrisTransport>(m) {
    m_vars = { { m_sys, "englacial_debris_column_mass", *m_grid } };
    m_vars[0].long_name("englacial debris mass per unit area").units("kg m^-2");
  }

protected:
  std::shared_ptr<array::Array> compute_impl(const Geometry & /*geometry*/) const {
    auto result = allocate<array::Scalar>("englacial_debris_column_mass");
    model->englacial().column_mass(*result);
    return result;
  }
};

/*! @brief Fraction of a cell covered by debris (Reid and Brock, 2010). */
class CoverFraction : public Diag<DebrisTransport> {
public:
  CoverFraction(const DebrisTransport *m) : Diag<DebrisTransport>(m) {
    m_vars = { { m_sys, "debris_cover_fraction", *m_grid } };
    m_vars[0].long_name("fraction of the ice surface covered by debris").units("1");
  }

protected:
  std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const {
    auto result = allocate<array::Scalar>("debris_cover_fraction");

    const double C = model->cover_fraction_coefficient();
    const array::Scalar &h = model->debris();

    array::AccessScope list{ &h, &geometry.cell_type, result.get() };
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();
      (*result)(i, j) = geometry.cell_type.icy(i, j) ? 1.0 - std::exp(-C * h(i, j)) : 0.0;
    }
    return result;
  }
};

/*! @brief Gravitational debris flux (cell-centered). */
class GravitationalFlux : public Diag<DebrisTransport> {
public:
  GravitationalFlux(const DebrisTransport *m) : Diag<DebrisTransport>(m) {
    m_vars = { { m_sys, "debris_gravitational_flux_x", *m_grid },
               { m_sys, "debris_gravitational_flux_y", *m_grid } };
    m_vars[0].long_name("x-component of the gravitational supraglacial debris flux")
        .units("m^2 s^-1").output_units("m^2 year^-1");
    m_vars[1].long_name("y-component of the gravitational supraglacial debris flux")
        .units("m^2 s^-1").output_units("m^2 year^-1");
  }

protected:
  std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const {
    auto result = allocate<array::Vector>("debris_gravitational_flux");
    array::staggered_to_regular(geometry.cell_type, model->gravitational().flux(), true, *result);
    return result;
  }
};

/*! @brief Time-averaged rate from a cumulative field. */
class AverageRate : public DiagAverageRate<DebrisTransport> {
public:
  AverageRate(const DebrisTransport *m, const std::string &name, const std::string &long_name,
              const array::Scalar &cumulative)
    : DiagAverageRate<DebrisTransport>(m, name, TOTAL_CHANGE), m_cumulative(cumulative) {

    m_accumulator.metadata()["units"] = "m";

    m_vars = { { m_sys, name, *m_grid } };
    m_vars[0].long_name(long_name).units("m s^-1").output_units("m year^-1");
    m_vars[0]["cell_methods"] = "time: mean";
    m_vars[0]["_FillValue"]   = { fill_value() };
  }

protected:
  const array::Scalar &model_input() {
    return m_cumulative;
  }

  const array::Scalar &m_cumulative;
};

/*! @brief Scalar snapshot of a mass (kg). */
class Mass : public TSDiag<TSSnapshotDiagnostic, DebrisTransport> {
public:
  typedef double (DebrisTransport::*Getter)() const;

  Mass(const DebrisTransport *m, const std::string &name, const std::string &long_name, Getter getter)
    : TSDiag<TSSnapshotDiagnostic, DebrisTransport>(m, name), m_getter(getter) {
    set_units("kg", "kg");
    m_variable["long_name"] = long_name;
  }

  double compute() {
    return (model->*m_getter)();
  }

private:
  Getter m_getter;
};

/*! @brief Scalar mass flux (kg s^-1) from the amount transferred during the last step. */
class MassFlux : public TSDiag<TSFluxDiagnostic, DebrisTransport> {
public:
  typedef double (DebrisTransport::*Getter)() const;

  MassFlux(const DebrisTransport *m, const std::string &name, const std::string &long_name, Getter getter)
    : TSDiag<TSFluxDiagnostic, DebrisTransport>(m, name), m_getter(getter) {
    set_units("kg s^-1", "kg year^-1");
    m_variable["long_name"] = long_name;
  }

  double compute() {
    return (model->*m_getter)();
  }

private:
  Getter m_getter;
};

} // end of namespace diagnostics

DiagnosticList DebrisTransport::spatial_diagnostics_impl() const {
  using namespace diagnostics;

  DiagnosticList result = DebrisModel::spatial_diagnostics_impl();

  result["debris_thickness"]               = Diagnostic::wrap(m_debris_thickness);
  result["englacial_debris_concentration"] = Diagnostic::wrap(m_concentration);
  result["englacial_debris_column_mass"]   = Diagnostic::Ptr(new ColumnMass(this));
  result["debris_cover_fraction"]          = Diagnostic::Ptr(new CoverFraction(this));
  result["debris_surface_velocity"]        = Diagnostic::wrap(m_supraglacial.surface_velocity());
  result["debris_gravitational_flux"]      = Diagnostic::Ptr(new GravitationalFlux(this));
  result["debris_input_rate"]              = Diagnostic::Ptr(new AverageRate(
      this, "debris_input_rate", "rate of debris input from the surrounding terrain (solid debris thickness)",
      m_cumulative_input));
  result["debris_melt_out_rate"] = Diagnostic::Ptr(new AverageRate(
      this, "debris_melt_out_rate", "rate of melt-out of englacial debris (solid debris thickness)",
      m_cumulative_melt_out));
  result["debris_removal_rate"] = Diagnostic::Ptr(new AverageRate(
      this, "debris_removal_rate", "rate of removal of supraglacial debris into the foreland",
      m_cumulative_removal));

  return result;
}

TSDiagnosticList DebrisTransport::scalar_diagnostics_impl() const {
  using namespace diagnostics;

  TSDiagnosticList result = DebrisModel::scalar_diagnostics_impl();

  result["englacial_debris_mass"] = TSDiagnostic::Ptr(new Mass(
      this, "englacial_debris_mass", "total mass of englacial debris", &DebrisTransport::englacial_mass));
  result["supraglacial_debris_mass"] = TSDiagnostic::Ptr(new Mass(
      this, "supraglacial_debris_mass", "total mass of supraglacial debris", &DebrisTransport::supraglacial_mass));
  result["total_debris_mass"] = TSDiagnostic::Ptr(new Mass(
      this, "total_debris_mass", "total mass of englacial and supraglacial debris",
      &DebrisTransport::total_mass));
  result["debris_mass_conservation_error"] = TSDiagnostic::Ptr(new Mass(
      this, "debris_mass_conservation_error",
      "debris mass minus its initial value, inputs, outputs and losses since the start of this run",
      &DebrisTransport::conservation_error));

  result["debris_input_mass_flux"] = TSDiagnostic::Ptr(new MassFlux(
      this, "debris_input_mass_flux", "rate of debris input from the surrounding terrain",
      &DebrisTransport::input_last_step));
  result["debris_melt_out_mass_flux"] = TSDiagnostic::Ptr(new MassFlux(
      this, "debris_melt_out_mass_flux", "rate of transfer of englacial debris to the surface by melt",
      &DebrisTransport::melt_out_last_step));
  result["debris_output_mass_flux"] = TSDiagnostic::Ptr(new MassFlux(
      this, "debris_output_mass_flux", "rate of removal of supraglacial debris into the foreland",
      &DebrisTransport::output_last_step));
  result["debris_lost_mass_flux"] = TSDiagnostic::Ptr(new MassFlux(
      this, "debris_lost_mass_flux",
      "rate of loss of debris with ice removed at the base or in cells that became ice-free",
      &DebrisTransport::lost_last_step));

  return result;
}

} // end of namespace debris
} // end of namespace pism
