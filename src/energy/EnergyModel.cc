/* Copyright (C) 2016, 2017, 2018 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "EnergyModel.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/stressbalance/StressBalance.hh"
#include "pism/util/io/PIO.hh"
#include "pism/util/Vars.hh"
#include "utilities.hh"
#include "pism/util/EnthalpyConverter.hh"
#include "bootstrapping.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/IceModelVec2CellType.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/Profiling.hh"

namespace pism {
namespace energy {

static void check_input(const IceModelVec *ptr, const char *name) {
  if (ptr == NULL) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "energy balance model input %s was not provided", name);
  }
}

Inputs::Inputs() {
  basal_frictional_heating = NULL;
  basal_heat_flux          = NULL;
  cell_type                = NULL;
  ice_thickness            = NULL;
  surface_liquid_fraction  = NULL;
  shelf_base_temp          = NULL;
  surface_temp             = NULL;
  till_water_thickness     = NULL;

  volumetric_heating_rate  = NULL;
  u3                       = NULL;
  v3                       = NULL;
  w3                       = NULL;

  no_model_mask = NULL;
}

void Inputs::check() const {
  check_input(cell_type,                "cell_type");
  check_input(basal_frictional_heating, "basal_frictional_heating");
  check_input(basal_heat_flux,          "basal_heat_flux");
  check_input(ice_thickness,            "ice_thickness");
  check_input(surface_liquid_fraction,  "surface_liquid_fraction");
  check_input(shelf_base_temp,          "shelf_base_temp");
  check_input(surface_temp,             "surface_temp");
  check_input(till_water_thickness,     "till_water_thickness");

  check_input(volumetric_heating_rate, "volumetric_heating_rate");
  check_input(u3, "u3");
  check_input(v3, "v3");
  check_input(w3, "w3");
}

EnergyModelStats::EnergyModelStats() {
  bulge_counter            = 0;
  reduced_accuracy_counter = 0;
  low_temperature_counter  = 0;
  liquified_ice_volume     = 0.0;
}

EnergyModelStats& EnergyModelStats::operator+=(const EnergyModelStats &other) {
  bulge_counter            += other.bulge_counter;
  reduced_accuracy_counter += other.reduced_accuracy_counter;
  low_temperature_counter  += other.low_temperature_counter;
  liquified_ice_volume     += other.liquified_ice_volume;
  return *this;
}


bool marginal(const IceModelVec2S &thickness, int i, int j, double threshold) {
  int
    n = j + 1,
    e = i + 1,
    s = j - 1,
    w = i - 1;

  const double
    N  = thickness(i, n),
    E  = thickness(e, j),
    S  = thickness(i, s),
    W  = thickness(w, j),
    NW = thickness(w, n),
    SW = thickness(w, s),
    NE = thickness(e, n),
    SE = thickness(e, s);

  return ((E  < threshold) or
          (NE < threshold) or
          (N  < threshold) or
          (NW < threshold) or
          (W  < threshold) or
          (SW < threshold) or
          (S  < threshold) or
          (SE < threshold));
}


void EnergyModelStats::sum(MPI_Comm com) {
  bulge_counter            = GlobalSum(com, bulge_counter);
  reduced_accuracy_counter = GlobalSum(com, reduced_accuracy_counter);
  low_temperature_counter  = GlobalSum(com, low_temperature_counter);
  liquified_ice_volume     = GlobalSum(com, liquified_ice_volume);
}


EnergyModel::EnergyModel(IceGrid::ConstPtr grid,
                         stressbalance::StressBalance *stress_balance)
  : Component(grid), m_stress_balance(stress_balance) {

  const unsigned int WIDE_STENCIL = m_config->get_double("grid.max_stencil_width");

  {
    m_ice_enthalpy.create(m_grid, "enthalpy", WITH_GHOSTS, WIDE_STENCIL);
    // POSSIBLE standard name = land_ice_enthalpy
    m_ice_enthalpy.set_attrs("model_state",
                             "ice enthalpy (includes sensible heat, latent heat, pressure)",
                             "J kg-1", "");
  }

  {
    m_basal_melt_rate.create(m_grid, "basal_melt_rate_grounded", WITHOUT_GHOSTS);
    // ghosted to allow the "redundant" computation of tauc
    m_basal_melt_rate.set_attrs("model_state",
                                "ice basal melt rate from energy conservation, in ice thickness per time (valid in grounded areas)",
                                "m s-1", "");
    // We could use land_ice_basal_melt_rate, but that way both basal_melt_rate_grounded and bmelt
    // have this standard name.
    m_basal_melt_rate.metadata().set_string("glaciological_units", "m year-1");
    m_basal_melt_rate.metadata().set_string("comment", "positive basal melt rate corresponds to ice loss");
  }

  // a 3d work vector
  {
    m_work.create(m_grid, "work_vector", WITHOUT_GHOSTS);
    m_work.set_attrs("internal",
                     "usually new values of temperature or enthalpy during time step",
                     "", "");
  }
}

void EnergyModel::init_enthalpy(const PIO &input_file, bool do_regrid, int record) {

  if (input_file.inq_var("enthalpy")) {
    if (do_regrid) {
      m_ice_enthalpy.regrid(input_file, CRITICAL);
    } else {
      m_ice_enthalpy.read(input_file, record);
    }
  } else if (input_file.inq_var("temp")) {
    IceModelVec3
      &temp    = m_work,
      &liqfrac = m_ice_enthalpy;

    {
      temp.set_name("temp");
      temp.metadata(0).set_name("temp");
      temp.set_attrs("temporary", "ice temperature", "Kelvin", "land_ice_temperature");

      if (do_regrid) {
        temp.regrid(input_file, CRITICAL);
      } else {
        temp.read(input_file, record);
      }
    }

    const IceModelVec2S & ice_thickness = *m_grid->variables().get_2d_scalar("land_ice_thickness");

    if (input_file.inq_var("liqfrac")) {
      SpatialVariableMetadata enthalpy_metadata = m_ice_enthalpy.metadata();

      liqfrac.set_name("liqfrac");
      liqfrac.metadata(0).set_name("liqfrac");
      liqfrac.set_attrs("temporary", "ice liquid water fraction",
                        "1", "");

      if (do_regrid) {
        liqfrac.regrid(input_file, CRITICAL);
      } else {
        liqfrac.read(input_file, record);
      }

      m_ice_enthalpy.set_name(enthalpy_metadata.get_name());
      m_ice_enthalpy.metadata() = enthalpy_metadata;

      m_log->message(2,
                     " - Computing enthalpy using ice temperature,"
                     "  liquid water fraction and thickness...\n");
      compute_enthalpy(temp, liqfrac, ice_thickness, m_ice_enthalpy);
    } else {
      m_log->message(2,
                     " - Computing enthalpy using ice temperature and thickness "
                     "and assuming zero liquid water fraction...\n");
      compute_enthalpy_cold(temp, ice_thickness, m_ice_enthalpy);
    }
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "neither enthalpy nor temperature was found in '%s'.\n",
                                  input_file.inq_filename().c_str());
  }
}

/*!
 * The `-regrid_file` may contain enthalpy, temperature, or *both* temperature and liquid water
 * fraction.
 */
void EnergyModel::regrid_enthalpy() {

  auto regrid_filename = m_config->get_string("input.regrid.file");
  auto regrid_vars     = set_split(m_config->get_string("input.regrid.vars"), ',');


  if (regrid_filename.empty()) {
    return;
  }

  std::string enthalpy_name = m_ice_enthalpy.metadata().get_name();

  if (regrid_vars.empty() or member(enthalpy_name, regrid_vars)) {
    PIO regrid_file(m_grid->com, "guess_mode", regrid_filename, PISM_READONLY);
    init_enthalpy(regrid_file, true, 0);
  }
}


void EnergyModel::restart(const PIO &input_file, int record) {
  this->restart_impl(input_file, record);
}

void EnergyModel::bootstrap(const PIO &input_file,
                            const IceModelVec2S &ice_thickness,
                            const IceModelVec2S &surface_temperature,
                            const IceModelVec2S &climatic_mass_balance,
                            const IceModelVec2S &basal_heat_flux) {
  this->bootstrap_impl(input_file,
                       ice_thickness, surface_temperature,
                       climatic_mass_balance, basal_heat_flux);
}

void EnergyModel::initialize(const IceModelVec2S &basal_melt_rate,
                             const IceModelVec2S &ice_thickness,
                             const IceModelVec2S &surface_temperature,
                             const IceModelVec2S &climatic_mass_balance,
                             const IceModelVec2S &basal_heat_flux) {
  this->initialize_impl(basal_melt_rate,
                        ice_thickness,
                        surface_temperature,
                        climatic_mass_balance,
                        basal_heat_flux);
}

void EnergyModel::update(double t, double dt, const Inputs &inputs) {
  // reset standard out flags at the beginning of every time step
  m_stdout_flags = "";
  m_stats = EnergyModelStats();

  const Profiling &profiling = m_grid->ctx()->profiling();

  profiling.begin("ice_energy");
  {
    // this call should fill m_work with new values of enthalpy
    this->update_impl(t, dt, inputs);

    m_work.update_ghosts(m_ice_enthalpy);
  }
  profiling.end("ice_energy");

  // globalize m_stats and update m_stdout_flags
  {
    char buffer[50] = "";

    m_stats.sum(m_grid->com);

    if (m_stats.reduced_accuracy_counter > 0.0) { // count of when BOMBPROOF switches to lower accuracy
      const double reduced_accuracy_percentage = 100.0 * (m_stats.reduced_accuracy_counter / (m_grid->Mx() * m_grid->My()));
      const double reporting_threshold = 5.0; // only report if above 5%

      if (reduced_accuracy_percentage > reporting_threshold and m_log->get_threshold() > 2) {
        snprintf(buffer, 50, "  [BPsacr=%.4f%%] ", reduced_accuracy_percentage);
        m_stdout_flags = buffer + m_stdout_flags;
      }
    }

    if (m_stats.bulge_counter > 0) {
      // count of when advection bulges are limited; frequently it is identically zero
      snprintf(buffer, 50, " BULGE=%d ", m_stats.bulge_counter);
      m_stdout_flags = buffer + m_stdout_flags;
    }
  }
}

MaxTimestep EnergyModel::max_timestep_impl(double t) const {
  // silence a compiler warning
  (void) t;

  if (m_stress_balance == NULL) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "EnergyModel: no stress balance provided."
                                  " Cannot compute max. time step.");
  }

  return MaxTimestep(m_stress_balance->max_timestep_cfl_3d().dt_max.value(), "energy");
}

const std::string& EnergyModel::stdout_flags() const {
  return m_stdout_flags;
}

const EnergyModelStats& EnergyModel::stats() const {
  return m_stats;
}

const IceModelVec3 & EnergyModel::enthalpy() const {
  return m_ice_enthalpy;
}

/*! @brief Basal melt rate in grounded areas. (It is set to zero elsewhere.) */
const IceModelVec2S & EnergyModel::basal_melt_rate() const {
  return m_basal_melt_rate;
}

/*! @brief Ice loss "flux" due to ice liquefaction. */
class LiquifiedIceFlux : public TSDiag<TSFluxDiagnostic,EnergyModel> {
public:
  LiquifiedIceFlux(const EnergyModel *m)
    : TSDiag<TSFluxDiagnostic, EnergyModel>(m, "liquified_ice_flux") {

    m_ts.variable().set_string("units", "m3 / second");
    m_ts.variable().set_string("glaciological_units", "m3 / year");
    m_ts.variable().set_string("long_name",
                               "rate of ice loss due to liquefaction,"
                               " averaged over the reporting interval");
    m_ts.variable().set_string("comment", "positive means ice loss");
    m_ts.variable().set_string("cell_methods", "time: mean");
  }
protected:
  double compute() {
    // liquified ice volume during the last time step
    return model->stats().liquified_ice_volume;
  }
};

namespace diagnostics {
/*! @brief Report ice enthalpy. */
class Enthalpy : public Diag<EnergyModel>
{
public:
  Enthalpy(const EnergyModel *m)
    : Diag<EnergyModel>(m) {
    m_vars = {model->enthalpy().metadata()};
  }

protected:
  IceModelVec::Ptr compute_impl() const {

    IceModelVec3::Ptr result(new IceModelVec3(m_grid, "enthalpy", WITHOUT_GHOSTS));
    result->metadata(0) = m_vars[0];

    const IceModelVec3 &input = model->enthalpy();

    // FIXME: implement IceModelVec3::copy_from()

    IceModelVec::AccessList list {result.get(), &input};
    ParallelSection loop(m_grid->com);
    try {
      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        result->set_column(i, j, input.get_column(i, j));
      }
    } catch (...) {
      loop.failed();
    }
    loop.check();


    return result;
  }
};

} // end of namespace diagnostics

DiagnosticList EnergyModel::diagnostics_impl() const {
  DiagnosticList result;
  result = {
    {"enthalpy",                 Diagnostic::Ptr(new diagnostics::Enthalpy(this))},
    {"basal_melt_rate_grounded", Diagnostic::wrap(m_basal_melt_rate)}
  };
  return result;
}

TSDiagnosticList EnergyModel::ts_diagnostics_impl() const {
  return {
    {"liquified_ice_flux", TSDiagnostic::Ptr(new LiquifiedIceFlux(this))}
  };
}

} // end of namespace energy

} // end of namespace pism
