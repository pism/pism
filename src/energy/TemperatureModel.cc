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

#include "TemperatureModel.hh"
#include "pism/energy/tempSystem.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/energy/utilities.hh"
#include "pism/util/IceModelVec2CellType.hh"
#include "pism/util/Vars.hh"
#include "pism/util/io/PIO.hh"

namespace pism {
namespace energy {

TemperatureModel::TemperatureModel(IceGrid::ConstPtr grid,
                                   stressbalance::StressBalance *stress_balance)
  : EnergyModel(grid, stress_balance) {

  m_ice_temperature.create(m_grid, "temp", WITH_GHOSTS);
  m_ice_temperature.set_attrs("model_state",
                              "ice temperature", "K", "land_ice_temperature");
  m_ice_temperature.metadata().set_double("valid_min", 0.0);
}

const IceModelVec3 & TemperatureModel::temperature() const {
  return m_ice_temperature;
}

void TemperatureModel::restart_impl(const PIO &input_file, int record) {

  m_log->message(2, "* Restarting the temperature-based energy balance model from %s...\n",
                 input_file.inq_filename().c_str());

  m_basal_melt_rate.read(input_file, record);

  const IceModelVec2S &ice_thickness = *m_grid->variables().get_2d_scalar("land_ice_thickness");

  if (input_file.inq_var(m_ice_temperature.metadata().get_name())) {
    m_ice_temperature.read(input_file, record);
  } else {
    init_enthalpy(input_file, false, record);
    compute_temperature(m_ice_enthalpy, ice_thickness, m_ice_temperature);
  }

  regrid("Temperature-based energy balance model", m_basal_melt_rate, REGRID_WITHOUT_REGRID_VARS);
  regrid("Temperature-based energy balance model", m_ice_temperature, REGRID_WITHOUT_REGRID_VARS);

  compute_enthalpy_cold(m_ice_temperature, ice_thickness, m_ice_enthalpy);
}

void TemperatureModel::bootstrap_impl(const PIO &input_file,
                                      const IceModelVec2S &ice_thickness,
                                      const IceModelVec2S &surface_temperature,
                                      const IceModelVec2S &climatic_mass_balance,
                                      const IceModelVec2S &basal_heat_flux) {

  m_log->message(2, "* Bootstrapping the temperature-based energy balance model from %s...\n",
                 input_file.inq_filename().c_str());

  m_basal_melt_rate.regrid(input_file, OPTIONAL,
                           m_config->get_double("bootstrapping.defaults.bmelt"));
  regrid("Temperature-based energy balance model", m_basal_melt_rate, REGRID_WITHOUT_REGRID_VARS);

  int temp_revision = m_ice_temperature.state_counter();
  regrid("Temperature-based energy balance model", m_ice_temperature, REGRID_WITHOUT_REGRID_VARS);

  if (temp_revision == m_ice_temperature.state_counter()) {
    bootstrap_ice_temperature(ice_thickness, surface_temperature, climatic_mass_balance,
                              basal_heat_flux, m_ice_temperature);
  }

  compute_enthalpy_cold(m_ice_temperature, ice_thickness, m_ice_enthalpy);
}

void TemperatureModel::initialize_impl(const IceModelVec2S &basal_melt_rate,
                                       const IceModelVec2S &ice_thickness,
                                       const IceModelVec2S &surface_temperature,
                                       const IceModelVec2S &climatic_mass_balance,
                                       const IceModelVec2S &basal_heat_flux) {

  m_log->message(2, "* Bootstrapping the temperature-based energy balance model...\n");

  m_basal_melt_rate.copy_from(basal_melt_rate);
  regrid("Temperature-based energy balance model", m_basal_melt_rate, REGRID_WITHOUT_REGRID_VARS);

  int temp_revision = m_ice_temperature.state_counter();
  regrid("Temperature-based energy balance model", m_ice_temperature, REGRID_WITHOUT_REGRID_VARS);

  if (temp_revision == m_ice_temperature.state_counter()) {
    bootstrap_ice_temperature(ice_thickness, surface_temperature, climatic_mass_balance,
                              basal_heat_flux, m_ice_temperature);
  }

  compute_enthalpy_cold(m_ice_temperature, ice_thickness, m_ice_enthalpy);
}

//! Takes a semi-implicit time-step for the temperature equation.
/*!
This method should be kept because it is worth having alternative physics, and
  so that older results can be reproduced.  In particular, this method is
  documented by papers [\ref BBL,\ref BBssasliding].   The new browser page
  \ref bombproofenth essentially documents the cold-ice-BOMBPROOF method here, but
  the newer enthalpy-based method is slightly different and (we hope) a superior
  implementation of the conservation of energy principle.

  The conservation of energy equation written in terms of temperature is
  \f[ \rho c_p(T) \frac{dT}{dt} = k \frac{\partial^2 T}{\partial z^2} + \Sigma,\f]
  where \f$T(t,x,y,z)\f$ is the temperature of the ice.  This equation is the shallow approximation
  of the full 3D conservation of energy.  Note \f$dT/dt\f$ stands for the material
  derivative, so advection is included.  Here \f$\rho\f$ is the density of ice,
  \f$c_p\f$ is its specific heat, and \f$k\f$ is its conductivity.  Also \f$\Sigma\f$ is the volume
  strain heating (with SI units of \f$J/(\text{s} \text{m}^3) = \text{W}\,\text{m}^{-3}\f$).

  We handle horizontal advection explicitly by first-order upwinding.  We handle
  vertical advection implicitly by centered differencing when possible, and retreat to
  implicit first-order upwinding when necessary.  There is a CFL condition
  for the horizontal explicit upwinding [\ref MortonMayers].  We report
  any CFL violations, but they are designed to not occur.

  The vertical conduction term is always handled implicitly (%i.e. by backward Euler).

    We work from the bottom of the column upward in building the system to solve
    (in the semi-implicit time-stepping scheme).  The excess energy above pressure melting
    is converted to melt-water, and that a fraction of this melt water is transported to
    the ice base according to the scheme in excessToFromBasalMeltLayer().

    The method uses equally-spaced calculation but the columnSystemCtx
    methods coarse_to_fine(), fine_to_coarse() interpolate
    back-and-forth from this equally-spaced computational grid to the
    (usually) non-equally spaced storage grid.

    An instance of tempSystemCtx is used to solve the tridiagonal system set-up here.

    In this procedure two scalar fields are modified: basal_melt_rate and m_work.
    But basal_melt_rate will never need to communicate ghosted values (horizontal stencil
    neighbors).  The ghosted values for m_ice_temperature are updated from the values in m_work in the
    communication done by energy_step().

  The (older) scheme cold-ice-BOMBPROOF, implemented here, is very reliable, but there is
  still an extreme and rare fjord situation which causes trouble.  For example, it
  occurs in one column of ice in one fjord perhaps only once
  in a 200ka simulation of the whole sheet, in my (ELB) experience modeling the Greenland
  ice sheet.  It causes the discretized advection bulge to give temperatures below that
  of the coldest ice anywhere, a continuum impossibility.  So as a final protection
  there is a "bulge limiter" which sets the temperature to the surface temperature of
  the column minus the bulge maximum (15 K) if it is below that level.  The number of
  times this occurs is reported as a "BPbulge" percentage.
  */
void TemperatureModel::update_impl(double t, double dt, const Inputs &inputs) {
  // current time does not matter here
  (void) t;

  using mask::ocean;

  Logger log(MPI_COMM_SELF, m_log->get_threshold());

  const double
    ice_density        = m_config->get_double("constants.ice.density"),
    ice_c              = m_config->get_double("constants.ice.specific_heat_capacity"),
    L                  = m_config->get_double("constants.fresh_water.latent_heat_of_fusion"),
    melting_point_temp = m_config->get_double("constants.fresh_water.melting_point_temperature"),
    beta_CC_grad       = m_config->get_double("constants.ice.beta_Clausius_Clapeyron") * ice_density * m_config->get_double("constants.standard_gravity");

  const bool allow_above_melting = m_config->get_boolean("energy.allow_temperature_above_melting");

  // this is bulge limit constant in K; is max amount by which ice
  //   or bedrock can be lower than surface temperature
  const double bulge_max  = m_config->get_double("energy.enthalpy.cold_bulge_max") / ice_c;

  inputs.check();
  const IceModelVec3
    &strain_heating3 = *inputs.volumetric_heating_rate,
    &u3              = *inputs.u3,
    &v3              = *inputs.v3,
    &w3              = *inputs.w3;

  const IceModelVec2CellType &cell_type = *inputs.cell_type;

  const IceModelVec2S
    &basal_frictional_heating = *inputs.basal_frictional_heating,
    &basal_heat_flux          = *inputs.basal_heat_flux,
    &ice_thickness            = *inputs.ice_thickness,
    &shelf_base_temp          = *inputs.shelf_base_temp,
    &ice_surface_temp         = *inputs.surface_temp,
    &till_water_thickness     = *inputs.till_water_thickness;

  IceModelVec::AccessList list{&ice_surface_temp, &shelf_base_temp, &ice_thickness,
      &cell_type, &basal_heat_flux, &till_water_thickness, &basal_frictional_heating,
      &u3, &v3, &w3, &strain_heating3, &m_basal_melt_rate, &m_ice_temperature, &m_work};

  energy::tempSystemCtx system(m_grid->z(), "temperature",
                               m_grid->dx(), m_grid->dy(), dt,
                               *m_config,
                               m_ice_temperature, u3, v3, w3, strain_heating3);

  double dz = system.dz();
  const std::vector<double>& z_fine = system.z();
  size_t Mz_fine = z_fine.size();
  std::vector<double> x(Mz_fine);// space for solution of system
  std::vector<double> Tnew(Mz_fine); // post-processed solution

  // counts unreasonably low temperature values; deprecated?
  unsigned int maxLowTempCount = m_config->get_double("energy.max_low_temperature_count");
  const double T_minimum = m_config->get_double("energy.minimum_allowed_temperature");

  double margin_threshold = m_config->get_double("energy.margin_ice_thickness_limit");

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      MaskValue mask = static_cast<MaskValue>(cell_type.as_int(i,j));

      const double H = ice_thickness(i, j);
      const double T_surface = ice_surface_temp(i, j);

      system.initThisColumn(i, j,
                            marginal(ice_thickness, i, j, margin_threshold),
                            mask, H);

      const int ks = system.ks();

      if (ks > 0) { // if there are enough points in ice to bother ...

        if (system.lambda() < 1.0) {
          m_stats.reduced_accuracy_counter += 1; // count columns with lambda < 1
        }

        // set boundary values for tridiagonal system
        system.setSurfaceBoundaryValuesThisColumn(T_surface);
        system.setBasalBoundaryValuesThisColumn(basal_heat_flux(i,j),
                                                shelf_base_temp(i,j),
                                                basal_frictional_heating(i,j));

        // solve the system for this column; melting not addressed yet
        system.solveThisColumn(x);
      }       // end of "if there are enough points in ice to bother ..."

      // prepare for melting/refreezing
      double bwatnew = till_water_thickness(i,j);

      // insert solution for generic ice segments
      for (int k=1; k <= ks; k++) {
        if (allow_above_melting) { // in the ice
          Tnew[k] = x[k];
        } else {
          const double
            Tpmp = melting_point_temp - beta_CC_grad * (H - z_fine[k]); // FIXME issue #15
          if (x[k] > Tpmp) {
            Tnew[k] = Tpmp;
            double Texcess = x[k] - Tpmp; // always positive
            column_drainage(ice_density, ice_c, L, z_fine[k], dz, &Texcess, &bwatnew);
            // Texcess  will always come back zero here; ignore it
          } else {
            Tnew[k] = x[k];
          }
        }
        if (Tnew[k] < T_minimum) {
          log.message(1,
                      "  [[too low (<200) ice segment temp T = %f at %d, %d, %d;"
                      " proc %d; mask=%d; w=%f m year-1]]\n",
                      Tnew[k], i, j, k, m_grid->rank(), mask,
                      units::convert(m_sys, system.w(k), "m second-1", "m year-1"));

          m_stats.low_temperature_counter++;
        }
        if (Tnew[k] < T_surface - bulge_max) {
          Tnew[k] = T_surface - bulge_max;
          m_stats.bulge_counter += 1;
        }
      }

      // insert solution for ice base segment
      if (ks > 0) {
        if (allow_above_melting == true) { // ice/rock interface
          Tnew[0] = x[0];
        } else {  // compute diff between x[k0] and Tpmp; melt or refreeze as appropriate
          const double Tpmp = melting_point_temp - beta_CC_grad * H; // FIXME issue #15
          double Texcess = x[0] - Tpmp; // positive or negative
          if (ocean(mask)) {
            // when floating, only half a segment has had its temperature raised
            // above Tpmp
            column_drainage(ice_density, ice_c, L, 0.0, dz/2.0, &Texcess, &bwatnew);
          } else {
            column_drainage(ice_density, ice_c, L, 0.0, dz, &Texcess, &bwatnew);
          }
          Tnew[0] = Tpmp + Texcess;
          if (Tnew[0] > (Tpmp + 0.00001)) {
            throw RuntimeError(PISM_ERROR_LOCATION, "updated temperature came out above Tpmp");
          }
        }
        if (Tnew[0] < T_minimum) {
          log.message(1,
                      "  [[too low (<200) ice/bedrock segment temp T = %f at %d,%d;"
                      " proc %d; mask=%d; w=%f]]\n",
                      Tnew[0],i,j,m_grid->rank(), mask,
                      units::convert(m_sys, system.w(0), "m second-1", "m year-1"));

          m_stats.low_temperature_counter++;
        }
        if (Tnew[0] < T_surface - bulge_max) {
          Tnew[0] = T_surface - bulge_max;
          m_stats.bulge_counter += 1;
        }
      }

      // set to air temp above ice
      for (unsigned int k = ks; k < Mz_fine; k++) {
        Tnew[k] = T_surface;
      }

      // transfer column into m_work; communication later
      system.fine_to_coarse(Tnew, i, j, m_work);

      // basal_melt_rate(i,j) is rate of mass loss at bottom of ice
      if (ocean(mask)) {
        m_basal_melt_rate(i,j) = 0.0;
      } else {
        // basalMeltRate is rate of change of bwat;  can be negative
        //   (subglacial water freezes-on); note this rate is calculated
        //   *before* limiting or other nontrivial modelling of bwat,
        //   which is Hydrology's job
        m_basal_melt_rate(i,j) = (bwatnew - till_water_thickness(i,j)) / dt;
      } // end of the grounded case
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  m_stats.low_temperature_counter = GlobalSum(m_grid->com, m_stats.low_temperature_counter);
  if (m_stats.low_temperature_counter > maxLowTempCount) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "too many low temps: %d",
                                  m_stats.low_temperature_counter);
  }

  // copy to m_ice_temperature, updating ghosts
  m_work.update_ghosts(m_ice_temperature);

  // Set ice enthalpy in place. EnergyModel::update will scatter ghosts
  compute_enthalpy_cold(m_work, ice_thickness, m_work);
}

void TemperatureModel::define_model_state_impl(const PIO &output) const {
  m_ice_temperature.define(output);
  m_basal_melt_rate.define(output);
  // ice enthalpy is not a part of the model state
}

void TemperatureModel::write_model_state_impl(const PIO &output) const {
  m_ice_temperature.write(output);
  m_basal_melt_rate.write(output);
  // ice enthalpy is not a part of the model state
}

//! Compute the melt water which should go to the base if \f$T\f$ is above pressure-melting.
void TemperatureModel::column_drainage(const double rho, const double c, const double L,
                                       const double z, const double dz,
                                       double *Texcess, double *bwat) const {

  const double
    darea      = m_grid->cell_area(),
    dvol       = darea * dz,
    dE         = rho * c * (*Texcess) * dvol,
    massmelted = dE / L;

  if (*Texcess >= 0.0) {
    // T is at or above pressure-melting temp, so temp needs to be set to
    // pressure-melting; a fraction of excess energy is turned into melt water at base
    // note massmelted is POSITIVE!
    const double FRACTION_TO_BASE
                         = (z < 100.0) ? 0.2 * (100.0 - z) / 100.0 : 0.0;
    // note: ice-equiv thickness:
    *bwat += (FRACTION_TO_BASE * massmelted) / (rho * darea);
    *Texcess = 0.0;
  } else {  // neither Texcess nor bwat needs to change if Texcess < 0.0
    // Texcess negative; only refreeze (i.e. reduce bwat) if at base and bwat > 0.0
    // note ONLY CALLED IF AT BASE!   note massmelted is NEGATIVE!
    if (z > 0.00001) {
      throw RuntimeError(PISM_ERROR_LOCATION, "excessToBasalMeltLayer() called with z not at base and negative Texcess");
    }
    if (*bwat > 0.0) {
      const double thicknessToFreezeOn = - massmelted / (rho * darea);
      if (thicknessToFreezeOn <= *bwat) { // the water *is* available to freeze on
        *bwat -= thicknessToFreezeOn;
        *Texcess = 0.0;
      } else { // only refreeze bwat thickness of water; update Texcess
        *bwat = 0.0;
        const double dTemp = L * (*bwat) / (c * dz);
        *Texcess += dTemp;
      }
    }
    // note: if *bwat == 0 and Texcess < 0.0 then Texcess unmolested; temp will go down
  }
}

} // end of namespace energy
} // end of namespace
