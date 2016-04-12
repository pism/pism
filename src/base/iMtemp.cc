// Copyright (C) 2004-2016 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "iceModel.hh"

#include "base/energy/bedrockThermalUnit.hh"
#include "base/energy/tempSystem.hh"
#include "base/hydrology/PISMHydrology.hh"
#include "base/stressbalance/PISMStressBalance.hh"
#include "base/util/IceGrid.hh"
#include "base/util/Mask.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/error_handling.hh"
#include "base/util/pism_options.hh"
#include "coupler/PISMOcean.hh"
#include "coupler/PISMSurface.hh"

namespace pism {

//! \file iMtemp.cc Methods of IceModel which implement the cold-ice, temperature-based formulation of conservation of energy.


//! Compute the melt water which should go to the base if \f$T\f$ is above pressure-melting.
void IceModel::excessToFromBasalMeltLayer(const double rho, const double c, const double L,
                                          const double z, const double dz,
                                          double *Texcess, double *bwat) {

  const double
    darea      = m_grid->dx() * m_grid->dy(),
    dvol       = darea * dz,
    dE         = rho * c * (*Texcess) * dvol,
    massmelted = dE / L;

  assert(not m_config->get_boolean("temperature_allow_above_melting"));

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
      throw RuntimeError("excessToBasalMeltLayer() called with z not at base and negative Texcess");
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

    In this procedure two scalar fields are modified: basal_melt_rate and vWork3d.
    But basal_melt_rate will never need to communicate ghosted values (horizontal stencil
    neighbors).  The ghosted values for m_ice_temperature are updated from the values in vWork3d in the
    communication done by energyStep().

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
void IceModel::temperatureStep(unsigned int *vertSacrCount, unsigned int *bulgeCount) {
  PetscErrorCode  ierr;

  bool viewOneColumn = options::Bool("-view_sys", "save column system information to a file");

  const double
    ice_density        = m_config->get_double("ice_density"),
    ice_c              = m_config->get_double("ice_specific_heat_capacity"),
    L                  = m_config->get_double("water_latent_heat_fusion"),
    melting_point_temp = m_config->get_double("water_melting_point_temperature"),
    beta_CC_grad       = m_config->get_double("beta_CC") * ice_density * m_config->get_double("standard_gravity");

  const bool allow_above_melting = m_config->get_boolean("temperature_allow_above_melting");


  // this is bulge limit constant in K; is max amount by which ice
  //   or bedrock can be lower than surface temperature
  const double bulgeMax  = m_config->get_double("enthalpy_cold_bulge_max") / ice_c;

  // now get map-plane fields, starting with coupler fields
  assert(m_surface != NULL);
  m_surface->ice_surface_temperature(m_ice_surface_temp);

  assert(m_ocean != NULL);
  m_ocean->shelf_base_temperature(m_shelfbtemp);

  assert(btu != NULL);
  const IceModelVec2S &G0 = btu->upward_geothermal_flux();

  IceModelVec2S &bwatcurr = vWork2d[0];
  bwatcurr.set_attrs("internal", "current amount of basal water", "m", "");
  bwatcurr.metadata().set_string("glaciological_units", "m");

  assert(subglacial_hydrology != NULL);
  subglacial_hydrology->subglacial_water_thickness(bwatcurr);

  IceModelVec::AccessList list;
  list.add(m_ice_surface_temp);
  list.add(m_shelfbtemp);

  list.add(m_ice_thickness);
  list.add(m_basal_melt_rate);
  list.add(m_cell_type);
  list.add(G0);
  list.add(bwatcurr);

  assert(m_stress_balance != NULL);
  // basal frictional heating
  const IceModelVec2S &Rb = m_stress_balance->basal_frictional_heating();
  const IceModelVec3 &strain_heating3 = m_stress_balance->volumetric_strain_heating();
  const IceModelVec3
    &u3 = m_stress_balance->velocity_u(),
    &v3 = m_stress_balance->velocity_v(),
    &w3 = m_stress_balance->velocity_w();

  energy::tempSystemCtx system(m_grid->z(), "temperature",
                               m_grid->dx(), m_grid->dy(), dt_TempAge,
                               *m_config,
                               m_ice_temperature, u3, v3, w3, strain_heating3);

  double dz = system.dz();
  const std::vector<double>& z_fine = system.z();
  size_t Mz_fine = z_fine.size();
  std::vector<double> x(Mz_fine);// space for solution of system
  std::vector<double> Tnew(Mz_fine);

  list.add(Rb);

  list.add(u3);
  list.add(v3);
  list.add(w3);
  list.add(strain_heating3);
  list.add(m_ice_temperature);
  list.add(vWork3d);

  // counts unreasonably low temperature values; deprecated?
  int myLowTempCount = 0;
  int maxLowTempCount = static_cast<int>(m_config->get_double("max_low_temp_count"));
  double globalMinAllowedTemp = m_config->get_double("global_min_allowed_temp");

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      MaskValue mask_value = static_cast<MaskValue>(m_cell_type.as_int(i,j));

      // false means "don't ignore horizontal advection and strain heating near margins"
      system.initThisColumn(i, j, false, mask_value, m_ice_thickness(i,j));

      const int ks = system.ks();

      if (ks > 0) { // if there are enough points in ice to bother ...

        if (system.lambda() < 1.0) {
          *vertSacrCount += 1; // count columns with lambda < 1
        }

        // set boundary values for tridiagonal system
        system.setSurfaceBoundaryValuesThisColumn(m_ice_surface_temp(i,j));
        system.setBasalBoundaryValuesThisColumn(G0(i,j), m_shelfbtemp(i,j), Rb(i,j));

        // solve the system for this column; melting not addressed yet
        system.solveThisColumn(x);

        if (viewOneColumn && (i == m_id && j == m_jd)) {
          ierr = PetscPrintf(m_grid->com,
                             "\n"
                             "in temperatureStep(): viewing tempSystemCtx at (i,j)=(%d,%d) to m-file... \n",
                             i, j);
          PISM_CHK(ierr, "PetscPrintf");

          system.save_to_file(x);
        }

      }       // end of "if there are enough points in ice to bother ..."

      // prepare for melting/refreezing
      double bwatnew = bwatcurr(i,j);

      // insert solution for generic ice segments
      for (int k=1; k <= ks; k++) {
        if (allow_above_melting == true) { // in the ice
          Tnew[k] = x[k];
        } else {
          const double
            Tpmp = melting_point_temp - beta_CC_grad * (m_ice_thickness(i,j) - z_fine[k]); // FIXME issue #15
          if (x[k] > Tpmp) {
            Tnew[k] = Tpmp;
            double Texcess = x[k] - Tpmp; // always positive
            excessToFromBasalMeltLayer(ice_density, ice_c, L, z_fine[k], dz, &Texcess, &bwatnew);
            // Texcess  will always come back zero here; ignore it
          } else {
            Tnew[k] = x[k];
          }
        }
        if (Tnew[k] < globalMinAllowedTemp) {
          ierr = PetscPrintf(PETSC_COMM_SELF,
                             "  [[too low (<200) ice segment temp T = %f at %d, %d, %d;"
                             " proc %d; mask=%d; w=%f m year-1]]\n",
                             Tnew[k], i, j, k, m_grid->rank(), m_cell_type.as_int(i, j),
                             units::convert(m_sys, system.w(k), "m second-1", "m year-1"));
          PISM_CHK(ierr, "PetscPrintf");

          myLowTempCount++;
        }
        if (Tnew[k] < m_ice_surface_temp(i,j) - bulgeMax) {
          Tnew[k] = m_ice_surface_temp(i,j) - bulgeMax;
          *bulgeCount += 1;
        }
      }

      // insert solution for ice base segment
      if (ks > 0) {
        if (allow_above_melting == true) { // ice/rock interface
          Tnew[0] = x[0];
        } else {  // compute diff between x[k0] and Tpmp; melt or refreeze as appropriate
          const double Tpmp = melting_point_temp - beta_CC_grad * m_ice_thickness(i,j); // FIXME issue #15
          double Texcess = x[0] - Tpmp; // positive or negative
          if (m_cell_type.ocean(i,j)) {
            // when floating, only half a segment has had its temperature raised
            // above Tpmp
            excessToFromBasalMeltLayer(ice_density, ice_c, L, 0.0, dz/2.0, &Texcess, &bwatnew);
          } else {
            excessToFromBasalMeltLayer(ice_density, ice_c, L, 0.0, dz, &Texcess, &bwatnew);
          }
          Tnew[0] = Tpmp + Texcess;
          if (Tnew[0] > (Tpmp + 0.00001)) {
            throw RuntimeError("updated temperature came out above Tpmp");
          }
        }
        if (Tnew[0] < globalMinAllowedTemp) {
          ierr = PetscPrintf(PETSC_COMM_SELF,
                             "  [[too low (<200) ice/bedrock segment temp T = %f at %d,%d;"
                             " proc %d; mask=%d; w=%f]]\n",
                             Tnew[0],i,j,m_grid->rank(),m_cell_type.as_int(i,j),
                             units::convert(m_sys, system.w(0), "m second-1", "m year-1"));
          PISM_CHK(ierr, "PetscPrintf");

          myLowTempCount++;
        }
        if (Tnew[0] < m_ice_surface_temp(i,j) - bulgeMax) {
          Tnew[0] = m_ice_surface_temp(i,j) - bulgeMax;
          *bulgeCount += 1;
        }
      }

      // set to air temp above ice
      for (unsigned int k = ks; k < Mz_fine; k++) {
        Tnew[k] = m_ice_surface_temp(i,j);
      }

      // transfer column into vWork3d; communication later
      system.fine_to_coarse(Tnew, i, j, vWork3d);

      // basal_melt_rate(i,j) is rate of mass loss at bottom of ice
      if (m_cell_type.ocean(i,j)) {
        m_basal_melt_rate(i,j) = 0.0;
      } else {
        // basalMeltRate is rate of change of bwat;  can be negative
        //   (subglacial water freezes-on); note this rate is calculated
        //   *before* limiting or other nontrivial modelling of bwat,
        //   which is Hydrology's job
        m_basal_melt_rate(i,j) = (bwatnew - bwatcurr(i,j)) / dt_TempAge;
      } // end of the grounded case
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  if (myLowTempCount > maxLowTempCount) {
    throw RuntimeError::formatted("too many low temps: %d", myLowTempCount);
  }
}


} // end of namespace pism
