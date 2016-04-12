// Copyright (C) 2009-2016 Andreas Aschwanden and Ed Bueler and Constantine Khroulev
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

#include "iceModel.hh"
#include "DrainageCalculator.hh"
#include "base/energy/bedrockThermalUnit.hh"
#include "base/energy/enthSystem.hh"
#include "base/hydrology/PISMHydrology.hh"
#include "base/stressbalance/PISMStressBalance.hh"
#include "base/util/IceGrid.hh"
#include "base/util/Mask.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/error_handling.hh"
#include "base/util/pism_options.hh"
#include "coupler/PISMOcean.hh"
#include "coupler/PISMSurface.hh"
#include "enthalpyConverter.hh"

namespace pism {

//! \file iMenthalpy.cc Methods of IceModel which implement the enthalpy formulation of conservation of energy.


//! Compute m_ice_enthalpy from temperature m_ice_temperature by assuming the ice has zero liquid fraction.
/*!
First this method makes sure the temperatures is at most the pressure-melting
value, before computing the enthalpy for that temperature, using zero liquid
fraction.

Because of how EnthalpyConverter::pressure() works, the energy
content in the air is set to the value that ice would have if it a chunk of it
occupied the air; the atmosphere actually has much lower energy content.  It is
done this way for regularity (i.e. dEnth/dz computations).

Because m_ice_enthalpy gets set, does ghost communication to finish.
 */
void IceModel::compute_enthalpy_cold(const IceModelVec3 &temperature, IceModelVec3 &result) {

  EnthalpyConverter::Ptr EC = m_ctx->enthalpy_converter();

  IceModelVec::AccessList list;
  list.add(temperature);
  list.add(result);
  list.add(m_ice_thickness);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double *Tij = temperature.get_column(i,j);
    double *Enthij = result.get_column(i,j);

    for (unsigned int k = 0; k < m_grid->Mz(); ++k) {
      const double depth = m_ice_thickness(i, j) - m_grid->z(k); // FIXME issue #15
      Enthij[k] = EC->enthalpy_permissive(Tij[k], 0.0, EC->pressure(depth));
    }
  }

  result.inc_state_counter();

  result.update_ghosts();
}


//! Compute m_ice_enthalpy from temperature m_ice_temperature and liquid fraction.
/*!
Because m_ice_enthalpy gets set, does ghost communication to finish.
 */
void IceModel::compute_enthalpy(const IceModelVec3 &temperature,
                                const IceModelVec3 &liquid_water_fraction,
                                IceModelVec3 &result) {

  EnthalpyConverter::Ptr EC = m_ctx->enthalpy_converter();

  IceModelVec::AccessList list;
  list.add(temperature);
  list.add(liquid_water_fraction);
  list.add(result);
  list.add(m_ice_thickness);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double *Tij = temperature.get_column(i,j);
    const double *Liqfracij = liquid_water_fraction.get_column(i,j);
    double *Enthij = result.get_column(i,j);

    for (unsigned int k=0; k<m_grid->Mz(); ++k) {
      const double depth = m_ice_thickness(i,j) - m_grid->z(k); // FIXME issue #15
      Enthij[k] = EC->enthalpy_permissive(Tij[k], Liqfracij[k],
                                        EC->pressure(depth));
    }
  }

  result.update_ghosts();

  result.inc_state_counter();
}

//! Compute the liquid fraction corresponding to m_ice_enthalpy, and put in a global IceModelVec3 provided by user.
/*!
Does not communicate ghosts for IceModelVec3 result
 */
void IceModel::compute_liquid_water_fraction(const IceModelVec3 &enthalpy,
                                             IceModelVec3 &result) {

  EnthalpyConverter::Ptr EC = m_ctx->enthalpy_converter();

  result.set_name("liqfrac");
  result.metadata(0).set_name("liqfrac");
  result.set_attrs("diagnostic",
                   "liquid water fraction in ice (between 0 and 1)",
                   "1", "", 0);

  IceModelVec::AccessList list;
  list.add(result);
  list.add(enthalpy);
  list.add(m_ice_thickness);

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double *Enthij = enthalpy.get_column(i,j);
      double *omegaij = result.get_column(i,j);

      for (unsigned int k=0; k<m_grid->Mz(); ++k) {
        const double depth = m_ice_thickness(i,j) - m_grid->z(k); // FIXME issue #15
        omegaij[k] = EC->water_fraction(Enthij[k],EC->pressure(depth));
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  result.inc_state_counter();
}

//! @brief Compute the CTS field, CTS = E/E_s(p), from m_ice_enthalpy, and put
//! in a global IceModelVec3 provided by user.
/*!
 * The actual cold-temperate transition surface (CTS) is the level set CTS = 1.
 *
 * Does not communicate ghosts for IceModelVec3 result.
 */
void IceModel::setCTSFromEnthalpy(IceModelVec3 &result) {

  EnthalpyConverter::Ptr EC = m_ctx->enthalpy_converter();

  result.set_name("cts");
  result.metadata(0).set_name("cts");
  result.set_attrs("diagnostic",
                   "cts = E/E_s(p), so cold-temperate transition surface is at cts = 1",
                   "", "", 0);

  IceModelVec::AccessList list;
  list.add(result);
  list.add(m_ice_enthalpy);
  list.add(m_ice_thickness);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double *CTS  = result.get_column(i,j);
    const double *enthalpy = m_ice_enthalpy.get_column(i,j);

    for (unsigned int k=0; k<m_grid->Mz(); ++k) {
      const double depth = m_ice_thickness(i,j) - m_grid->z(k); // FIXME issue #15
      CTS[k] = enthalpy[k] / EC->enthalpy_cts(EC->pressure(depth));
    }
  }
}

//! Update ice enthalpy field based on conservation of energy.
/*!
This method is documented by the page \ref bombproofenth and by [\ref
AschwandenBuelerKhroulevBlatter].

This method updates IceModelVec3 vWork3d = vEnthnew and IceModelVec2S basal_melt_rate.
No communication of ghosts is done for any of these fields.

We use an instance of enthSystemCtx.

Regarding drainage, see [\ref AschwandenBuelerKhroulevBlatter] and references therein.

\image html BC-decision-chart.png "Setting the basal boundary condition"
 */
void IceModel::enthalpyAndDrainageStep(unsigned int *vertSacrCount,
                                       double* liquifiedVol,
                                       unsigned int *bulgeCount) {

  EnthalpyConverter::Ptr EC = m_ctx->enthalpy_converter();

  assert(not m_config->get_boolean("do_cold_ice_methods"));

  // essentially physical constants:
  const double
    ice_density = m_config->get_double("ice_density"),              // kg m-3
    // constants controlling the numerical method:
    bulgeEnthMax = m_config->get_double("enthalpy_cold_bulge_max"); // J kg-1

  bool viewOneColumn = options::Bool("-view_sys",
                                     "save column system information to a file");

  energy::DrainageCalculator dc(*m_config);

  const IceModelVec2S &Rb = m_stress_balance->basal_frictional_heating();
  const IceModelVec3
    &u3 = m_stress_balance->velocity_u(),
    &v3 = m_stress_balance->velocity_v(),
    &w3 = m_stress_balance->velocity_w();

  const IceModelVec3 &strain_heating3 = m_stress_balance->volumetric_strain_heating();

  energy::enthSystemCtx system(m_grid->z(), "enth", m_grid->dx(), m_grid->dy(), dt_TempAge,
                               *m_config, m_ice_enthalpy, u3, v3, w3, strain_heating3, EC);

  size_t Mz_fine = system.z().size();
  double dz = system.dz();
  std::vector<double> Enthnew(Mz_fine); // new enthalpy in column

  // Now get map-plane coupler fields: Dirichlet upper surface
  // boundary and mass balance lower boundary under shelves
  assert(m_surface != NULL);
  m_surface->ice_surface_temperature(m_ice_surface_temp);
  m_surface->ice_surface_liquid_water_fraction(m_liqfrac_surface);

  assert(m_ocean != NULL);
  m_ocean->shelf_base_temperature(m_shelfbtemp);

  assert(btu != NULL);
  const IceModelVec2S &basal_heat_flux = btu->upward_geothermal_flux();

  IceModelVec2S &till_water_thickness = vWork2d[0];
  till_water_thickness.set_attrs("internal", "current amount of basal water in the till",
                                 "m", "");
  assert(subglacial_hydrology != NULL);
  subglacial_hydrology->till_water_thickness(till_water_thickness);

  IceModelVec::AccessList list;
  list.add(m_ice_surface_temp);

  list.add(m_shelfbtemp);

  // get other map-plane fields
  list.add(m_liqfrac_surface);
  list.add(m_ice_thickness);
  list.add(m_basal_melt_rate);
  list.add(Rb);
  list.add(basal_heat_flux);
  list.add(till_water_thickness);
  list.add(m_cell_type);

  // these are accessed a column at a time
  list.add(u3);
  list.add(v3);
  list.add(w3);
  list.add(strain_heating3);
  list.add(m_ice_enthalpy);
  list.add(vWork3d);

  unsigned int liquifiedCount = 0;

  ParallelSection loop(m_grid->com);
  try {
    for (Points pt(*m_grid); pt; pt.next()) {
      const int i = pt.i(), j = pt.j();

      system.init(i, j, m_ice_thickness(i, j));

      // enthalpy and pressures at top of ice
      const double
        depth_ks = m_ice_thickness(i, j) - system.ks() * dz,
        p_ks     = EC->pressure(depth_ks); // FIXME issue #15

      double Enth_ks = EC->enthalpy_permissive(m_ice_surface_temp(i, j), m_liqfrac_surface(i, j),
                                               p_ks);

      const bool ice_free_column = (system.ks() == 0);

      // deal completely with columns with no ice; enthalpy and basal_melt_rate need setting
      if (ice_free_column) {
        vWork3d.set_column(i, j, Enth_ks);
        // The floating basal melt rate will be set later; cover this
        // case and set to zero for now. Also, there is no basal melt
        // rate on ice free land and ice free ocean
        m_basal_melt_rate(i, j) = 0.0;
        continue;
      } // end of if (ice_free_column)

      if (system.lambda() < 1.0) {
        *vertSacrCount += 1; // count columns with lambda < 1
      }

      const bool is_floating = m_cell_type.ocean(i, j);
      bool base_is_warm = system.Enth(0) >= system.Enth_s(0);
      bool above_base_is_warm = system.Enth(1) >= system.Enth_s(1);

      // set boundary conditions and update enthalpy
      {
        system.set_surface_dirichlet_bc(Enth_ks);

        // determine lowest-level equation at bottom of ice; see
        // decision chart in the source code browser and page
        // documenting BOMBPROOF
        if (is_floating) {
          // floating base: Dirichlet application of known temperature from ocean
          //   coupler; assumes base of ice shelf has zero liquid fraction
          double Enth0 = EC->enthalpy_permissive(m_shelfbtemp(i, j), 0.0,
                                               EC->pressure(m_ice_thickness(i, j)));

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
            system.set_basal_heat_flux(basal_heat_flux(i, j) + Rb(i, j));
          }
        }

        // solve the system
        system.solve(Enthnew);

        if (viewOneColumn && (i == m_id && j == m_jd)) {
          system.save_to_file(Enthnew);
        }
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
              depth = m_ice_thickness(i, j) - k * dz,
              p     = EC->pressure(depth), // FIXME issue #15
              T_m   = EC->melting_temperature(p),
              L     = EC->L(T_m),
              omega = EC->water_fraction(Enthnew[k], p);

            if (Enthnew[k] >= system.Enth_s(k) + 0.5 * L) {
              liquifiedCount++; // count these rare events...
              Enthnew[k] = system.Enth_s(k) + 0.5 * L; //  but lose the energy
            }

            if (omega > 0.01) {                          // FIXME: make "0.01" configurable here
              double fractiondrained = dc.get_drainage_rate(omega) * dt_TempAge; // pure number

              fractiondrained  = std::min(fractiondrained, omega - 0.01); // only drain down to 0.01
              Hdrainedtotal   += fractiondrained * dz; // always a positive contribution
              Enthnew[k]      -= fractiondrained * L;
            }
          }
        }

        // apply bulge limiter
        const double lowerEnthLimit = Enth_ks - bulgeEnthMax;
        for (unsigned int k=0; k < system.ks(); k++) {
          if (Enthnew[k] < lowerEnthLimit) {
            *bulgeCount += 1;      // count the columns which have very large cold
            Enthnew[k] = lowerEnthLimit;  // limit advection bulge ... enthalpy not too low
          }
        }

        // if there is subglacial water, don't allow ice base enthalpy to be below
        // pressure-melting; that is, assume subglacial water is at the pressure-
        // melting temperature and enforce continuity of temperature
        {
          if (Enthnew[0] < system.Enth_s(0) && till_water_thickness(i,j) > 0.0) {
            const double E_difference = system.Enth_s(0) - Enthnew[0];

            const double depth = m_ice_thickness(i, j),
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
        if (is_floating == true) {
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
              p_0 = EC->pressure(m_ice_thickness(i, j)),
              p_1 = EC->pressure(m_ice_thickness(i, j) - dz), // FIXME issue #15
              Tpmp_0 = EC->melting_temperature(p_0);

            const bool k1_istemperate = EC->is_temperate(Enthnew[1], p_1); // level  z = + \Delta z
            double hf_up;
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
            m_basal_melt_rate(i, j) = (Rb(i, j) + basal_heat_flux(i, j) - hf_up) / (ice_density * EC->L(Tpmp_0));

            if (till_water_thickness(i, j) <= 0 && m_basal_melt_rate(i, j) < 0) {
              m_basal_melt_rate(i, j) = 0.0;
            }
          }

          // Add drained water from the column to basal melt rate.
          m_basal_melt_rate(i, j) += (Hdrainedtotal - Hfrozen) / dt_TempAge;
        } // end of the grounded case
      } // end of the basal melt rate computation

      system.fine_to_coarse(Enthnew, i, j, vWork3d);
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();


  // FIXME: use cell areas
  *liquifiedVol = ((double) liquifiedCount) * dz * m_grid->dx() * m_grid->dy();
}

/*
  The decision above was produced using this TikZ source code:

% Define block styles
\tikzstyle{decision} = [ellipse, draw, text width=7em, text badly centered, inner sep=2pt]
\tikzstyle{block} = [rectangle, draw, text width=5em, text badly centered, rounded corners, minimum height=4em]
\tikzstyle{line} = [draw, -latex']

\begin{tikzpicture}[node distance = 3cm, auto]
    % Place nodes
    \node (invisiblestart) {};

    \node [decision, below of=invisiblestart, text height=0.2cm] (coldvstemp) {$H<H_{\text s}(p)$ ?};
    \node [decision, left of=coldvstemp, xshift=-4em] (excludebad) {$\eta_{\text b}>0$ ?};
    \node [block, below of=excludebad, text width=6em] (fixbad) {$H := H_{\text s}(p)$};

    % edges
    \path [line] (invisiblestart) -- (coldvstemp);
    \path [line] (excludebad) -- node [text width=6em] {yes (consider base to be temperate)} (fixbad);

    % cold branch:
    \node [block, left of=fixbad, text width=7.5em] (coldmodeltype) {Eqn (49) is Neumann b.c.~for Eqn (67); $M_b=0$};
    % edges
    \path [line] (coldvstemp) -- node {yes} (excludebad);
    \path [line] (excludebad) -- node {no} (coldmodeltype);

    % temperate branch
    \node [block, below of=coldvstemp, text width=12em] (qtemperate) {$\nabla H \cdot \bn=0$ is Neumann b.c.~for Eqn (67)};
    \node [decision, below left of=qtemperate, text width=8em] (tempthick) {positive thickness of temperate ice at base?};
    \node [block, below right of=tempthick, text width=10em] (Mbforqtemperate) {$\bq = - k(H,p)\nabla T_{\text m}(p)$ \\ at ice base};
    \node [block, below left of=tempthick, text width=9em, xshift=-4em] (Mbforqcold) {$\bq = - K_{\text i}(H) \nabla H$ \\ at ice base};
    \node [block, below left of=Mbforqtemperate, text width=9em] (getMbtemp) {compute $M_b$ from Eqn (50) or (66)};

    % edges
    \path [line] (fixbad) -- (qtemperate);
    \path [line] (coldvstemp) -- node {no} (qtemperate);
    \path [line] (tempthick) -- node [above] {no} (Mbforqcold);
    \path [line] (tempthick) -- node {yes} (Mbforqtemperate);
    \path [line] (qtemperate) -- (tempthick);
    \path [line] (Mbforqcold) -- node {} (getMbtemp);
    \path [line] (Mbforqtemperate) -- node {} (getMbtemp);
\end{tikzpicture}
 */

} // end of namespace pism
