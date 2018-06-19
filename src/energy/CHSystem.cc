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
#include <algorithm>            // std::max

#include "CHSystem.hh"

#include "DrainageCalculator.hh"
#include "pism/util/EnthalpyConverter.hh"
#include "pism/energy/enthSystem.hh"
#include "pism/util/IceModelVec2CellType.hh"
#include "pism/util/io/PIO.hh"
#include "utilities.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {
namespace energy {

/*!
 * Model of the energy content in the cryo-hydrologic system.
 *
 * This model can be used to include the influence of cryo-hydrologic warming on the
 * internal energy of the ice.
 *
 * This model is based on
 *
 * @article{PhillipsRajaramSteffan2010,
 * author = {Phillips, T. and Rajaram, H. and Steffen, K.},
 * doi = {10.1029/2010GL044397},
 * journal = {Geophy. Res. Lett.},
 * number = {20},
 * pages = {1--5},
 * title = {{Cryo-hydrologic warming: A potential mechanism for rapid thermal response of ice sheets}},
 * volume = {37},
 * year = {2010}
 * }
 *
 * We assume that during the melt season the enthalpy in the cryo-hydrologic system is
 * equal to the enthalpy of the ice with a fixed water fraction corresponding to the
 * amount of water remaining in the system once all of the moving water drains out (0.5%
 * by default). During the winter the CH system is allowed to cool.
*/

CHSystem::CHSystem(IceGrid::ConstPtr grid,
                   stressbalance::StressBalance *stress_balance)
  : EnergyModel(grid, stress_balance) {

  m_ice_enthalpy.set_name("ch_enthalpy");
  m_ice_enthalpy.metadata().set_name("ch_enthalpy");
  m_ice_enthalpy.metadata().set_string("long_name",
                                       "enthalpy of the cryo-hydrologic system");
}

CHSystem::~CHSystem() {
  // empty
}

void CHSystem::restart_impl(const PIO &input_file, int record) {

  m_log->message(2, "* Restarting the cryo-hydrologic system from %s...\n",
                 input_file.inq_filename().c_str());

  init_enthalpy(input_file, false, record);

  regrid_enthalpy();
}

void CHSystem::bootstrap_impl(const PIO &input_file,
                              const IceModelVec2S &ice_thickness,
                              const IceModelVec2S &surface_temperature,
                              const IceModelVec2S &climatic_mass_balance,
                              const IceModelVec2S &basal_heat_flux) {

  m_log->message(2, "* Bootstrapping the cryo-hydrologic warming model from %s...\n",
                 input_file.inq_filename().c_str());

  int enthalpy_revision = m_ice_enthalpy.state_counter();
  regrid_enthalpy();

  if (enthalpy_revision == m_ice_enthalpy.state_counter()) {
    bootstrap_ice_enthalpy(ice_thickness, surface_temperature, climatic_mass_balance,
                           basal_heat_flux, m_ice_enthalpy);
  }
}

void CHSystem::initialize_impl(const IceModelVec2S &basal_melt_rate,
                               const IceModelVec2S &ice_thickness,
                               const IceModelVec2S &surface_temperature,
                               const IceModelVec2S &climatic_mass_balance,
                               const IceModelVec2S &basal_heat_flux) {
  (void) basal_melt_rate;

  m_log->message(2, "* Bootstrapping the cryo-hydrologic warming model...\n");

  int enthalpy_revision = m_ice_enthalpy.state_counter();
  regrid_enthalpy();

  if (enthalpy_revision == m_ice_enthalpy.state_counter()) {
    bootstrap_ice_enthalpy(ice_thickness, surface_temperature, climatic_mass_balance,
                           basal_heat_flux, m_ice_enthalpy);
  }
}

//! Update the enthalpy of the cryo-hydrologic system.
/*!
  This method updates IceModelVec3 m_work. No communication of ghosts is done.
*/
void CHSystem::update_impl(double t, double dt, const Inputs &inputs) {
  // current time does not matter here
  (void) t;

  EnthalpyConverter::Ptr EC = m_grid->ctx()->enthalpy_converter();

  inputs.check();

  // give them names that are a bit shorter...
  const IceModelVec3
    &volumetric_heat = *inputs.volumetric_heating_rate,
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
    &ice_surface_temp         = *inputs.surface_temp;

  energy::enthSystemCtx system(m_grid->z(), "energy.ch_warming", m_grid->dx(), m_grid->dy(), dt,
                               *m_config, m_ice_enthalpy, u3, v3, w3, volumetric_heat, EC);

  const size_t Mz_fine = system.z().size();
  const double dz = system.dz();
  std::vector<double> Enthnew(Mz_fine); // new enthalpy in column

  IceModelVec::AccessList list{&ice_surface_temp, &shelf_base_temp, &surface_liquid_fraction,
      &ice_thickness, &basal_frictional_heating, &basal_heat_flux,
      &cell_type, &u3, &v3, &w3, &volumetric_heat, &m_ice_enthalpy,
      &m_work};

  double
    margin_threshold = m_config->get_double("energy.margin_ice_thickness_limit"),
    T_pm = m_config->get_double("constants.fresh_water.melting_point_temperature"),
    residual_water_fraction = m_config->get_double("energy.ch_warming.residual_water_fraction");

  const std::vector<double> &z = m_grid->z();
  const unsigned int Mz = m_grid->Mz();

  ParallelSection loop(m_grid->com);
  try {
    for (Points pt(*m_grid); pt; pt.next()) {
      const int i = pt.i(), j = pt.j();

      const double H = ice_thickness(i, j);

      if (ice_surface_temp(i, j) >= T_pm) {
        // We use surface temperature to determine if we're in a melt season or not. It
        // probably makes sense to use the surface mass balance instead.

        double *column = m_work.get_column(i, j);
        for (unsigned int k = 0; k < Mz; ++k) {
          double
            depth = std::max(H - z[k], 0.0),
            P     = EC->pressure(depth);
            column[k] = EC->enthalpy(EC->melting_temperature(P),
                                     residual_water_fraction,
                                     P);
        }
        continue;
      }

      // enthalpy and pressures at top of ice
      const double
        depth_ks = H - system.ks() * dz,
        p_ks     = EC->pressure(depth_ks); // FIXME issue #15

      const double Enth_ks = EC->enthalpy_permissive(ice_surface_temp(i, j),
                                                     surface_liquid_fraction(i, j), p_ks);

      system.init(i, j,
                  marginal(ice_thickness, i, j, margin_threshold),
                  H);

      const bool ice_free_column = (system.ks() == 0);

      // deal completely with columns with no ice
      if (ice_free_column) {
        m_work.set_column(i, j, Enth_ks);
        continue;
      } // end of if (ice_free_column)

      if (system.lambda() < 1.0) {
        m_stats.reduced_accuracy_counter += 1; // count columns with lambda < 1
      }

      // set boundary conditions and update enthalpy
      {
        system.set_surface_dirichlet_bc(Enth_ks);

        if (cell_type.ocean(i, j)) {
          // floating base: Dirichlet application of known temperature from ocean coupler;
          //   assumes base of ice shelf has zero liquid fraction
          double Enth0 = EC->enthalpy_permissive(shelf_base_temp(i, j), 0.0, EC->pressure(H));

          system.set_basal_dirichlet_bc(Enth0);
        } else {
          // grounded
          system.set_basal_heat_flux(basal_heat_flux(i, j) + basal_frictional_heating(i, j));
        }
        // solve the system
        system.solve(Enthnew);
      }

      system.fine_to_coarse(Enthnew, i, j, m_work);
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

void CHSystem::define_model_state_impl(const PIO &output) const {
  m_ice_enthalpy.define(output);
}

void CHSystem::write_model_state_impl(const PIO &output) const {
  m_ice_enthalpy.write(output);
}

/*!
 * Compute the heat flux corresponding to the cryo-hydrologic warming.
 *
 * `Q = (k / R**2) * (T_ch - T_ice),`
 *
 * where `k` is the thermal conductivity of ice and `R` us the average spacing of
 * channels in the cryo-hydrologic system.
 */
void cryo_hydrologic_warming_flux(double k,
                                  double R,
                                  const IceModelVec2S &ice_thickness,
                                  const IceModelVec3 &ice_enthalpy,
                                  const IceModelVec3 &ch_enthalpy,
                                  IceModelVec3 &result) {

  auto grid = result.grid();

  const auto &z = grid->z();
  auto Mz = grid->Mz();

  auto EC = grid->ctx()->enthalpy_converter();

  IceModelVec::AccessList access{&ice_thickness, &ice_enthalpy, &ch_enthalpy, &result};

  double C = k / (R * R);

  ParallelSection loop(grid->com);
  try {
    for (Points p(*grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double
        *E_ch   = ch_enthalpy.get_column(i, j),
        *E_ice  = ice_enthalpy.get_column(i, j);
      double *Q = result.get_column(i, j);

      for (unsigned int m = 0; m < Mz; ++m) {
        double
          depth = ice_thickness(i, j) - z[m];

        if (depth > 0.0) {
          double P = EC->pressure(depth);
          Q[m] = std::max(C * (EC->temperature(E_ch[m], P) - EC->temperature(E_ice[m], P)),
                          0.0);
        } else {
          Q[m] = 0.0;
        }
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}


} // end of namespace energy
} // end of namespace pism
