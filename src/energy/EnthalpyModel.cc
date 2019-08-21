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

#include "EnthalpyModel.hh"

#include "DrainageCalculator.hh"
#include "pism/util/EnthalpyConverter.hh"
#include "pism/energy/enthSystem.hh"
#include "pism/util/IceModelVec2CellType.hh"
#include "pism/util/io/PIO.hh"
#include "utilities.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {
namespace energy {

EnthalpyModel::EnthalpyModel(IceGrid::ConstPtr grid,
                             stressbalance::StressBalance *stress_balance)
  : EnergyModel(grid, stress_balance) {
  // empty
}

void EnthalpyModel::restart_impl(const PIO &input_file, int record) {

  m_log->message(2, "* Restarting the enthalpy-based energy balance model from %s...\n",
                 input_file.inq_filename().c_str());

  m_basal_melt_rate.read(input_file, record);
  init_enthalpy(input_file, false, record);

  regrid("Energy balance model", m_basal_melt_rate, REGRID_WITHOUT_REGRID_VARS);
  regrid_enthalpy();
}

void EnthalpyModel::bootstrap_impl(const PIO &input_file,
                                   const IceModelVec2S &ice_thickness,
                                   const IceModelVec2S &surface_temperature,
                                   const IceModelVec2S &climatic_mass_balance,
                                   const IceModelVec2S &basal_heat_flux) {

  m_log->message(2, "* Bootstrapping the enthalpy-based energy balance model from %s...\n",
                 input_file.inq_filename().c_str());

  m_basal_melt_rate.regrid(input_file, OPTIONAL,
                           m_config->get_double("bootstrapping.defaults.bmelt"));

  regrid("Energy balance model", m_basal_melt_rate, REGRID_WITHOUT_REGRID_VARS);

  int enthalpy_revision = m_ice_enthalpy.state_counter();
  regrid_enthalpy();

  if (enthalpy_revision == m_ice_enthalpy.state_counter()) {
    bootstrap_ice_enthalpy(ice_thickness, surface_temperature, climatic_mass_balance,
                           basal_heat_flux, m_ice_enthalpy);
  }
}

void EnthalpyModel::initialize_impl(const IceModelVec2S &basal_melt_rate,
                                    const IceModelVec2S &ice_thickness,
                                    const IceModelVec2S &surface_temperature,
                                    const IceModelVec2S &climatic_mass_balance,
                                    const IceModelVec2S &basal_heat_flux) {

  m_log->message(2, "* Bootstrapping the enthalpy-based energy balance model...\n");

  m_basal_melt_rate.copy_from(basal_melt_rate);

  regrid("Energy balance model", m_basal_melt_rate, REGRID_WITHOUT_REGRID_VARS);

  int enthalpy_revision = m_ice_enthalpy.state_counter();
  regrid_enthalpy();

  if (enthalpy_revision == m_ice_enthalpy.state_counter()) {
    bootstrap_ice_enthalpy(ice_thickness, surface_temperature, climatic_mass_balance,
                           basal_heat_flux, m_ice_enthalpy);
  }
}

//! Update ice enthalpy field based on conservation of energy.
/*!
This method is documented by the page \ref bombproofenth and by [\ref
AschwandenBuelerKhroulevBlatter].

This method updates IceModelVec3 m_work and IceModelVec2S basal_melt_rate.
No communication of ghosts is done for any of these fields.

We use an instance of enthSystemCtx.

Regarding drainage, see [\ref AschwandenBuelerKhroulevBlatter] and references therein.
 */

void EnthalpyModel::update_impl(double t, double dt, const Inputs &inputs) {
  // current time does not matter here
  (void) t;

  EnthalpyConverter::Ptr EC = m_grid->ctx()->enthalpy_converter();

  const double
    ice_density           = m_config->get_double("constants.ice.density"), // kg m-3
    bulgeEnthMax          = m_config->get_double("energy.enthalpy.cold_bulge_max"), // J kg-1
    target_water_fraction = m_config->get_double("energy.drainage_target_water_fraction");

  energy::DrainageCalculator dc(*m_config);

  inputs.check();

  // give them names that are a bit shorter...
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
    &surface_liquid_fraction  = *inputs.surface_liquid_fraction,
    &shelf_base_temp          = *inputs.shelf_base_temp,
    &ice_surface_temp         = *inputs.surface_temp,
    &till_water_thickness     = *inputs.till_water_thickness;

  energy::enthSystemCtx system(m_grid->z(), "energy.enthalpy", m_grid->dx(), m_grid->dy(), dt,
                               *m_config, m_ice_enthalpy, u3, v3, w3, strain_heating3, EC);

  const size_t Mz_fine = system.z().size();
  const double dz = system.dz();
  std::vector<double> Enthnew(Mz_fine); // new enthalpy in column

  IceModelVec::AccessList list{&ice_surface_temp, &shelf_base_temp, &surface_liquid_fraction,
      &ice_thickness, &basal_frictional_heating, &basal_heat_flux, &till_water_thickness,
      &cell_type, &u3, &v3, &w3, &strain_heating3, &m_basal_melt_rate, &m_ice_enthalpy,
      &m_work};

  double margin_threshold = m_config->get_double("energy.margin_ice_thickness_limit");

  unsigned int liquifiedCount = 0;

  ParallelSection loop(m_grid->com);
  try {
    for (Points pt(*m_grid); pt; pt.next()) {
      const int i = pt.i(), j = pt.j();

      const double H = ice_thickness(i, j);

      system.init(i, j,
                  marginal(ice_thickness, i, j, margin_threshold),
                  H);

      // enthalpy and pressures at top of ice
      const double
        depth_ks = H - system.ks() * dz,
        p_ks     = EC->pressure(depth_ks); // FIXME issue #15

      const double Enth_ks = EC->enthalpy_permissive(ice_surface_temp(i, j),
                                                     surface_liquid_fraction(i, j), p_ks);

      const bool ice_free_column = (system.ks() == 0);

      // deal completely with columns with no ice; enthalpy and basal_melt_rate need setting
      if (ice_free_column) {
        m_work.set_column(i, j, Enth_ks);
        // The floating basal melt rate will be set later; cover this
        // case and set to zero for now. Also, there is no basal melt
        // rate on ice free land and ice free ocean
        m_basal_melt_rate(i, j) = 0.0;
        continue;
      } // end of if (ice_free_column)

      if (system.lambda() < 1.0) {
        m_stats.reduced_accuracy_counter += 1; // count columns with lambda < 1
      }

      const bool
        is_floating        = cell_type.ocean(i, j),
        base_is_warm       = system.Enth(0) >= system.Enth_s(0),
        above_base_is_warm = system.Enth(1) >= system.Enth_s(1);

      // set boundary conditions and update enthalpy
      {
        system.set_surface_dirichlet_bc(Enth_ks);

        // determine lowest-level equation at bottom of ice; see
        // decision chart in the source code browser and page
        // documenting BOMBPROOF
        if (is_floating) {
          // floating base: Dirichlet application of known temperature from ocean
          //   coupler; assumes base of ice shelf has zero liquid fraction
          double Enth0 = EC->enthalpy_permissive(shelf_base_temp(i, j), 0.0, EC->pressure(H));

          system.set_basal_dirichlet_bc(Enth0);
        } else {
          // grounded ice warm and wet
          if (base_is_warm && (till_water_thickness(i, j) > 0.0)) {
            if (above_base_is_warm) {
              // temperate layer at base (Neumann) case:  q . n = 0  (K0 grad E . n = 0)
              system.set_basal_heat_flux(0.0);
            } else {
              // only the base is warm: E = E_s(p) (Dirichlet)
              // ( Assumes ice has zero liquid fraction. Is this a valid assumption here?
              system.set_basal_dirichlet_bc(system.Enth_s(0));
            }
          } else {
            // (Neumann) case:  q . n = q_lith . n + F_b
            // a) cold and dry base, or
            // b) base that is still warm from the last time step, but without basal water
            system.set_basal_heat_flux(basal_heat_flux(i, j) + basal_frictional_heating(i, j));
          }
        }

        // solve the system
        system.solve(Enthnew);

      }

      // post-process (drainage and bulge-limiting)
      double Hdrainedtotal = 0.0;
      double Hfrozen = 0.0;
      {
        // drain ice segments by mechanism in [\ref AschwandenBuelerKhroulevBlatter],
        //   using DrainageCalculator dc
        for (unsigned int k=0; k < system.ks(); k++) {
          if (Enthnew[k] > system.Enth_s(k)) { // avoid doing any more work if cold

            const double
              depth = H - k * dz,
              p     = EC->pressure(depth), // FIXME issue #15
              T_m   = EC->melting_temperature(p),
              L     = EC->L(T_m);

            if (Enthnew[k] >= system.Enth_s(k) + 0.5 * L) {
              liquifiedCount++; // count these rare events...
              Enthnew[k] = system.Enth_s(k) + 0.5 * L; //  but lose the energy
            }

            double omega = EC->water_fraction(Enthnew[k], p);

            if (omega > target_water_fraction) {
              double fractiondrained = dc.get_drainage_rate(omega) * dt; // pure number

              fractiondrained  = std::min(fractiondrained,
                                          omega - target_water_fraction);
              Hdrainedtotal   += fractiondrained * dz; // always a positive contribution
              Enthnew[k]      -= fractiondrained * L;
            }
          }
        }

        // apply bulge limiter
        const double lowerEnthLimit = Enth_ks - bulgeEnthMax;
        for (unsigned int k=0; k < system.ks(); k++) {
          if (Enthnew[k] < lowerEnthLimit) {
            // Count grid points which have very large cold limit advection bulge... enthalpy not
            // too low.
            m_stats.bulge_counter += 1;
            Enthnew[k] = lowerEnthLimit;
          }
        }

        // if there is subglacial water, don't allow ice base enthalpy to be below
        // pressure-melting; that is, assume subglacial water is at the pressure-
        // melting temperature and enforce continuity of temperature
        {
          if (Enthnew[0] < system.Enth_s(0) && till_water_thickness(i,j) > 0.0) {
            const double E_difference = system.Enth_s(0) - Enthnew[0];

            const double depth = H,
              pressure         = EC->pressure(depth),
              T_m              = EC->melting_temperature(pressure);

            Enthnew[0] = system.Enth_s(0);
            // This adjustment creates energy out of nothing. We will
            // freeze some basal water, subtracting an equal amount of
            // energy, to make up for it.
            //
            // Note that [E_difference] = J/kg, so
            //
            // U_difference = E_difference * ice_density * dx * dy * (0.5*dz)
            //
            // is the amount of energy created (we changed enthalpy of
            // a block of ice with the volume equal to
            // dx*dy*(0.5*dz); note that the control volume
            // corresponding to the grid point at the base of the
            // column has thickness 0.5*dz, not dz).
            //
            // Also, [L] = J/kg, so
            //
            // U_freeze_on = L * ice_density * dx * dy * Hfrozen,
            //
            // is the amount of energy created by freezing a water
            // layer of thickness Hfrozen (using units of ice
            // equivalent thickness).
            //
            // Setting U_difference = U_freeze_on and solving for
            // Hfrozen, we find the thickness of the basal water layer
            // we need to freeze co restore energy conservation.

            Hfrozen = E_difference * (0.5*dz) / EC->L(T_m);
          }
        }

      } // end of post-processing

      // compute basal melt rate
      {
        bool base_is_cold = (Enthnew[0] < system.Enth_s(0)) && (till_water_thickness(i,j) == 0.0);
        // Determine melt rate, but only preliminarily because of
        // drainage, from heat flux out of bedrock, heat flux into
        // ice, and frictional heating
        if (is_floating) {
          // The floating basal melt rate will be set later; cover
          // this case and set to zero for now. Note that
          // Hdrainedtotal is discarded (the ocean model determines
          // the basal melt).
          m_basal_melt_rate(i, j) = 0.0;
        } else {
          if (base_is_cold) {
            m_basal_melt_rate(i, j) = 0.0;  // zero melt rate if cold base
          } else {
            const double
              p_0 = EC->pressure(H),
              p_1 = EC->pressure(H - dz), // FIXME issue #15
              Tpmp_0 = EC->melting_temperature(p_0);

            const bool k1_istemperate = EC->is_temperate(Enthnew[1], p_1); // level  z = + \Delta z
            double hf_up = 0.0;
            if (k1_istemperate) {
              const double
                Tpmp_1 = EC->melting_temperature(p_1);

              hf_up = -system.k_from_T(Tpmp_0) * (Tpmp_1 - Tpmp_0) / dz;
            } else {
              double T_0 = EC->temperature(Enthnew[0], p_0);
              const double K_0 = system.k_from_T(T_0) / EC->c();

              hf_up = -K_0 * (Enthnew[1] - Enthnew[0]) / dz;
            }

            // compute basal melt rate from flux balance:
            //
            // basal_melt_rate = - Mb / rho in [\ref AschwandenBuelerKhroulevBlatter];
            //
            // after we compute it we make sure there is no refreeze if
            // there is no available basal water
            m_basal_melt_rate(i, j) = (basal_frictional_heating(i, j) + basal_heat_flux(i, j) - hf_up) / (ice_density * EC->L(Tpmp_0));

            if (till_water_thickness(i, j) <= 0 && m_basal_melt_rate(i, j) < 0) {
              m_basal_melt_rate(i, j) = 0.0;
            }
          }

          // Add drained water from the column to basal melt rate.
          m_basal_melt_rate(i, j) += (Hdrainedtotal - Hfrozen) / dt;
        } // end of the grounded case
      } // end of the basal melt rate computation

      system.fine_to_coarse(Enthnew, i, j, m_work);
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  m_stats.liquified_ice_volume = ((double) liquifiedCount) * dz * m_grid->cell_area();
}

void EnthalpyModel::define_model_state_impl(const PIO &output) const {
  m_ice_enthalpy.define(output);
  m_basal_melt_rate.define(output);
}

void EnthalpyModel::write_model_state_impl(const PIO &output) const {
  m_ice_enthalpy.write(output);
  m_basal_melt_rate.write(output);
}

} // end of namespace energy
} // end of namespace pism
