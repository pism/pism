// Copyright (C) 2004-2014 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include "tempSystem.hh"
#include "Mask.hh"
#include "PISMSurface.hh"
#include "PISMOcean.hh"
#include "PISMStressBalance.hh"
#include "PISMHydrology.hh"
#include "bedrockThermalUnit.hh"
#include "pism_options.hh"

#include <cassert>

#include "error_handling.hh"

namespace pism {

//! \file iMtemp.cc Methods of IceModel which implement the cold-ice, temperature-based formulation of conservation of energy.


//! Compute the melt water which should go to the base if \f$T\f$ is above pressure-melting.
void IceModel::excessToFromBasalMeltLayer(const double rho, const double c, const double L,
                                          const double z, const double dz,
                                          double *Texcess, double *bwat) {

  const double
    darea      = grid.dx() * grid.dy(),
    dvol       = darea * dz,
    dE         = rho * c * (*Texcess) * dvol,
    massmelted = dE / L;

  assert(config.get_flag("temperature_allow_above_melting") == false);

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

    The method uses equally-spaced calculation but the methods getValColumn(),
    setValColumn() interpolate back-and-forth from this equally-spaced calculational
    grid to the (usually) non-equally spaced storage grid.

    An instance of tempSystemCtx is used to solve the tridiagonal system set-up here.

    In this procedure two scalar fields are modified: basal_melt_rate and vWork3d.
    But basal_melt_rate will never need to communicate ghosted values (horizontal stencil
    neighbors).  The ghosted values for T3 are updated from the values in vWork3d in the
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
void IceModel::temperatureStep(double* vertSacrCount, double* bulgeCount) {
  PetscErrorCode  ierr;

  bool viewOneColumn;
  OptionsIsSet("-view_sys", viewOneColumn);

  const double
    ice_density        = config.get("ice_density"),
    ice_c              = config.get("ice_specific_heat_capacity"),
    L                  = config.get("water_latent_heat_fusion"),
    melting_point_temp = config.get("water_melting_point_temperature"),
    beta_CC_grad       = config.get("beta_CC") * ice_density * config.get("standard_gravity");

  const bool allow_above_melting = config.get_flag("temperature_allow_above_melting");


  // this is bulge limit constant in K; is max amount by which ice
  //   or bedrock can be lower than surface temperature
  const double bulgeMax  = config.get("enthalpy_cold_bulge_max") / ice_c;

  // now get map-plane fields, starting with coupler fields
  assert(surface != NULL);
  surface->ice_surface_temperature(ice_surface_temp);

  assert(ocean != NULL);
  ocean->shelf_base_temperature(shelfbtemp);

  IceModelVec2S &G0 = vWork2d[0];
  G0.set_attrs("internal", "upward geothermal flux at z=0", "W m-2", "");
  G0.set_glaciological_units("mW m-2");

  assert(btu != NULL);
  btu->get_upward_geothermal_flux(G0);

  IceModelVec2S &bwatcurr = vWork2d[1];
  bwatcurr.set_attrs("internal", "current amount of basal water", "m", "");
  bwatcurr.set_glaciological_units("m");

  assert(subglacial_hydrology != NULL);
  subglacial_hydrology->subglacial_water_thickness(bwatcurr);

  IceModelVec::AccessList list;
  list.add(ice_surface_temp);
  list.add(shelfbtemp);

  list.add(ice_thickness);
  list.add(basal_melt_rate);
  list.add(vMask);
  list.add(G0);
  list.add(bwatcurr);

  IceModelVec2S *Rb;            // basal frictional heating
  assert(stress_balance != NULL);
  stress_balance->get_basal_frictional_heating(Rb);

  IceModelVec3 *u3, *v3, *w3, *strain_heating3;
  stress_balance->get_3d_velocity(u3, v3, w3);
  stress_balance->get_volumetric_strain_heating(strain_heating3);

  tempSystemCtx system(grid.z(), "temperature",
                       grid.dx(), grid.dy(), dt_TempAge,
                       config,
                       &T3, u3, v3, w3, strain_heating3);

  double dz = system.dz();
  const std::vector<double>& z_fine = system.z();
  size_t Mz_fine = z_fine.size();
  std::vector<double> x(Mz_fine);// space for solution of system
  std::vector<double> Tnew(Mz_fine);

  list.add(*Rb);

  list.add(*u3);
  list.add(*v3);
  list.add(*w3);
  list.add(*strain_heating3);
  list.add(T3);
  list.add(vWork3d);

  // counts unreasonably low temperature values; deprecated?
  int myLowTempCount = 0;
  int maxLowTempCount = static_cast<int>(config.get("max_low_temp_count"));
  double globalMinAllowedTemp = config.get("global_min_allowed_temp");

  MaskQuery mask(vMask);

  const double thickness_threshold = config.get("energy_advection_ice_thickness_threshold");

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // if isMarginal then only do vertical conduction for ice; ignore advection
    //   and strain heating if isMarginal
    const bool isMarginal = checkThinNeigh(ice_thickness, i, j, thickness_threshold);
    MaskValue mask_value = static_cast<MaskValue>(vMask.as_int(i,j));

    system.initThisColumn(i, j, isMarginal, mask_value, ice_thickness(i,j));

    const int ks = system.ks();

    if (ks > 0) { // if there are enough points in ice to bother ...

      if (system.lambda() < 1.0) {
        *vertSacrCount += 1; // count columns with lambda < 1
      }

      // set boundary values for tridiagonal system
      system.setSurfaceBoundaryValuesThisColumn(ice_surface_temp(i,j));
      system.setBasalBoundaryValuesThisColumn(G0(i,j),shelfbtemp(i,j),(*Rb)(i,j));

      // solve the system for this column; melting not addressed yet
      system.solveThisColumn(x);

      if (viewOneColumn && (i == id && j == jd)) {
        ierr = PetscPrintf(grid.com,
                           "\n"
                           "in temperatureStep(): viewing tempSystemCtx at (i,j)=(%d,%d) to m-file... \n",
                           i, j);
        PISM_PETSC_CHK(ierr, "PetscPrintf");
        system.viewColumnInfoMFile(x);
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
          Tpmp = melting_point_temp - beta_CC_grad * (ice_thickness(i,j) - z_fine[k]); // FIXME issue #15
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
                           " proc %d; mask=%d; w=%f m/year]]\n",
                           Tnew[k], i, j, k, grid.rank(), vMask.as_int(i, j),
                           grid.convert(system.w(k), "m/s", "m/year"));
        PISM_PETSC_CHK(ierr, "PetscPrintf");
        myLowTempCount++;
      }
      if (Tnew[k] < ice_surface_temp(i,j) - bulgeMax) {
        Tnew[k] = ice_surface_temp(i,j) - bulgeMax;  bulgeCount++;   }
    }

    // insert solution for ice base segment
    if (ks > 0) {
      if (allow_above_melting == true) { // ice/rock interface
        Tnew[0] = x[0];
      } else {  // compute diff between x[k0] and Tpmp; melt or refreeze as appropriate
        const double Tpmp = melting_point_temp - beta_CC_grad * ice_thickness(i,j); // FIXME issue #15
        double Texcess = x[0] - Tpmp; // positive or negative
        if (mask.ocean(i,j)) {
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
                           Tnew[0],i,j,grid.rank(),vMask.as_int(i,j),
                           grid.convert(system.w(0), "m/s", "m/year"));
        PISM_PETSC_CHK(ierr, "PetscPrintf");
        myLowTempCount++;
      }
      if (Tnew[0] < ice_surface_temp(i,j) - bulgeMax) {
        Tnew[0] = ice_surface_temp(i,j) - bulgeMax;   bulgeCount++;   }
    }

    // set to air temp above ice
    for (unsigned int k = ks; k < Mz_fine; k++) {
      Tnew[k] = ice_surface_temp(i,j);
    }

    // transfer column into vWork3d; communication later
    vWork3d.setValColumnPL(i,j,Tnew);

    // basal_melt_rate(i,j) is rate of mass loss at bottom of ice
    if (mask.ocean(i,j)) {
      basal_melt_rate(i,j) = 0.0;
    } else {
      // basalMeltRate is rate of change of bwat;  can be negative
      //   (subglacial water freezes-on); note this rate is calculated
      //   *before* limiting or other nontrivial modelling of bwat,
      //   which is Hydrology's job
      basal_melt_rate(i,j) = (bwatnew - bwatcurr(i,j)) / dt_TempAge;
    } // end of the grounded case
  }

  if (myLowTempCount > maxLowTempCount) {
    throw RuntimeError::formatted("too many low temps: %d", myLowTempCount);
  }
}


} // end of namespace pism
