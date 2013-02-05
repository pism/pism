// Copyright (C) 2004-2012 Jed Brown, Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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


//! \file iMtemp.cc Methods of IceModel which implement the cold-ice, temperature-based formulation of conservation of energy.


//! Compute the melt water which should go to the base if \f$T\f$ is above pressure-melting.
PetscErrorCode IceModel::excessToFromBasalMeltLayer(
                const PetscScalar rho, const PetscScalar c, const PetscScalar L,
                const PetscScalar z, const PetscScalar dz,
                PetscScalar *Texcess, PetscScalar *bwat) {

  const PetscScalar darea = grid.dx * grid.dy,
                    dvol = darea * dz,
                    dE = rho * c * (*Texcess) * dvol,
                    massmelted = dE / L;

  if (allowAboveMelting == PETSC_TRUE) {
    SETERRQ(grid.com, 1,"IceModel::excessToBasalMeltLayer() called but allowAboveMelting==TRUE");
  }
  if (*Texcess >= 0.0) {
    // T is at or above pressure-melting temp, so temp needs to be set to
    // pressure-melting; a fraction of excess energy is turned into melt water at base
    // note massmelted is POSITIVE!
    const PetscScalar FRACTION_TO_BASE
                         = (z < 100.0) ? 0.2 * (100.0 - z) / 100.0 : 0.0;
    // note: ice-equiv thickness:
    *bwat += (FRACTION_TO_BASE * massmelted) / (rho * darea);
    *Texcess = 0.0;
  } else {  // neither Texcess nor bwat needs to change if Texcess < 0.0
    // Texcess negative; only refreeze (i.e. reduce bwat) if at base and bwat > 0.0
    // note ONLY CALLED IF AT BASE!   note massmelted is NEGATIVE!
    if (z > 0.00001) {
      SETERRQ(grid.com, 1, "excessToBasalMeltLayer() called with z not at base and negative Texcess");
    }
    if (*bwat > 0.0) {
      const PetscScalar thicknessToFreezeOn = - massmelted / (rho * darea);
      if (thicknessToFreezeOn <= *bwat) { // the water *is* available to freeze on
        *bwat -= thicknessToFreezeOn;
        *Texcess = 0.0;
      } else { // only refreeze bwat thickness of water; update Texcess
        *bwat = 0.0;
        const PetscScalar dTemp = L * (*bwat) / (c * dz);
        *Texcess += dTemp;
      }
    }
    // note: if *bwat == 0 and Texcess < 0.0 then Texcess unmolested; temp will go down
  }
  return 0;
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

    In this procedure two scalar fields are modified: vbmr and vWork3d.
    But vbmr will never need to communicate ghosted values (horizontal stencil
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
PetscErrorCode IceModel::temperatureStep(PetscScalar* vertSacrCount, PetscScalar* bulgeCount) {
    PetscErrorCode  ierr;

    // set up fine grid in ice
    PetscInt    fMz    = grid.Mz_fine;
    PetscScalar fdz    = grid.dz_fine;
    vector<double> &fzlev = grid.zlevels_fine;

    ierr = verbPrintf(5,grid.com,
                      "\n  [entering temperatureStep(); fMz = %d, fdz = %5.3f]",
                      fMz, fdz); CHKERRQ(ierr);

    bool viewOneColumn;
    ierr = PISMOptionsIsSet("-view_sys", viewOneColumn); CHKERRQ(ierr);

    const PetscScalar
      ice_rho   = config.get("ice_density"),
      ice_k     = config.get("ice_thermal_conductivity"),
      ice_c     = config.get("ice_specific_heat_capacity"),
      L         = config.get("water_latent_heat_fusion"),
      melting_point_temp = config.get("water_melting_point_temperature"),
      beta_CC_grad = config.get("beta_CC") * ice_rho * config.get("standard_gravity");

    tempSystemCtx system(fMz, "temperature");
    system.dx              = grid.dx;
    system.dy              = grid.dy;
    system.dtTemp          = dt_TempAge; // same time step for temp and age, currently
    system.dzEQ            = fdz;
    system.ice_rho         = ice_rho;
    system.ice_k           = ice_k;
    system.ice_c_p         = ice_c;

    PetscScalar *x;
    x = new PetscScalar[fMz]; // space for solution of system

    // this is bulge limit constant in K; is max amount by which ice
    //   or bedrock can be lower than surface temperature
    const PetscScalar bulgeMax  = config.get("enthalpy_cold_bulge_max") / ice_c;

    PetscScalar *Tnew;
    // pointers to values in current column
    system.u     = new PetscScalar[fMz];
    system.v     = new PetscScalar[fMz];
    system.w     = new PetscScalar[fMz];
    system.Sigma = new PetscScalar[fMz];
    system.T     = new PetscScalar[fMz];
    Tnew         = new PetscScalar[fMz];

    // system needs access to T3 for T3.getPlaneStar_fine()
    system.T3 = &T3;

    // checks that all needed constants and pointers got set:
    ierr = system.initAllColumns(); CHKERRQ(ierr);

    // now get map-plane fields, starting with coupler fields
    PetscScalar  **basalMeltRate;

    if (surface != PETSC_NULL) {
      ierr = surface->ice_surface_temperature(artm); CHKERRQ(ierr);
    } else {
      SETERRQ(grid.com, 1,"PISM ERROR: surface == PETSC_NULL");
    }
    if (ocean != PETSC_NULL) {
      ierr = ocean->shelf_base_mass_flux(shelfbmassflux); CHKERRQ(ierr);
      ierr = ocean->shelf_base_temperature(shelfbtemp); CHKERRQ(ierr);
    } else {
      SETERRQ(grid.com, 1,"PISM ERROR: ocean == PETSC_NULL");
    }

    IceModelVec2S G0 = vWork2d[0];
    ierr = G0.set_attrs("internal", "upward geothermal flux at z=0", "W m-2", ""); CHKERRQ(ierr);
    ierr = G0.set_glaciological_units("mW m-2");
    if (btu) {
      ierr = btu->get_upward_geothermal_flux(G0); CHKERRQ(ierr);
    } else {
      SETERRQ(grid.com, 3,"PISM ERROR: PISMBedThermalUnit* btu == PETSC_NULL in temperatureStep()");
    }

    IceModelVec2S bwatcurr = vWork2d[1];
    ierr = bwatcurr.set_attrs("internal", "current amount of basal water", "m", ""); CHKERRQ(ierr);
    ierr = bwatcurr.set_glaciological_units("m");
    if (subglacial_hydrology) {
      ierr = subglacial_hydrology->subglacial_water_thickness(bwatcurr); CHKERRQ(ierr);
    } else {
      SETERRQ(grid.com, 3,"PISM ERROR: PISMHydrology* subglacial_hydrology is NULL in temperatureStep()");
    }

    ierr = artm.begin_access(); CHKERRQ(ierr);
    ierr = shelfbmassflux.begin_access(); CHKERRQ(ierr);
    ierr = shelfbtemp.begin_access(); CHKERRQ(ierr);

    ierr = vH.begin_access(); CHKERRQ(ierr);
    ierr = vbmr.get_array(basalMeltRate); CHKERRQ(ierr);
    ierr = vMask.begin_access(); CHKERRQ(ierr);
    ierr = G0.begin_access(); CHKERRQ(ierr);
    ierr = bwatcurr.begin_access(); CHKERRQ(ierr);

    IceModelVec2S *Rb;            // basal frictional heating
    ierr = stress_balance->get_basal_frictional_heating(Rb); CHKERRQ(ierr);

    IceModelVec3 *u3, *v3, *w3, *Sigma3;
    ierr = stress_balance->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr);
    ierr = stress_balance->get_volumetric_strain_heating(Sigma3); CHKERRQ(ierr);

    ierr = Rb->begin_access(); CHKERRQ(ierr);

    ierr = u3->begin_access(); CHKERRQ(ierr);
    ierr = v3->begin_access(); CHKERRQ(ierr);
    ierr = w3->begin_access(); CHKERRQ(ierr);
    ierr = Sigma3->begin_access(); CHKERRQ(ierr);
    ierr = T3.begin_access(); CHKERRQ(ierr);
    ierr = vWork3d.begin_access(); CHKERRQ(ierr);

    // counts unreasonably low temperature values; deprecated?
    PetscInt myLowTempCount = 0;
    PetscInt maxLowTempCount = static_cast<PetscInt>(config.get("max_low_temp_count"));
    PetscReal globalMinAllowedTemp = config.get("global_min_allowed_temp");

    MaskQuery mask(vMask);

    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

        // this should *not* be replaced by call to grid.kBelowHeight():
        const PetscInt  ks = static_cast<PetscInt>(floor(vH(i,j)/fdz));

        if (ks>0) { // if there are enough points in ice to bother ...
          ierr = system.setIndicesAndClearThisColumn(i,j,ks); CHKERRQ(ierr);

          ierr = u3->getValColumn(i,j,ks,system.u); CHKERRQ(ierr);
          ierr = v3->getValColumn(i,j,ks,system.v); CHKERRQ(ierr);
          ierr = w3->getValColumn(i,j,ks,system.w); CHKERRQ(ierr);
          ierr = Sigma3->getValColumn(i,j,ks,system.Sigma); CHKERRQ(ierr);
          ierr = T3.getValColumn(i,j,ks,system.T); CHKERRQ(ierr);

          // go through column and find appropriate lambda for BOMBPROOF
          PetscScalar lambda = 1.0;  // start with centered implicit for more accuracy
          for (PetscInt k = 1; k < ks; k++) {
            const PetscScalar denom = (PetscAbs(system.w[k]) + 0.000001/secpera)
              * ice_rho * ice_c * fdz;
            lambda = PetscMin(lambda, 2.0 * ice_k / denom);
          }
          if (lambda < 1.0)  *vertSacrCount += 1; // count columns with lambda < 1
          // if isMarginal then only do vertical conduction for ice; ignore advection
          //   and strain heating if isMarginal
          const bool isMarginal = checkThinNeigh(vH(i+1,j),vH(i+1,j+1),vH(i,j+1),vH(i-1,j+1),
                                                 vH(i-1,j),vH(i-1,j-1),vH(i,j-1),vH(i+1,j-1));
          PismMask mask_value = static_cast<PismMask>(vMask.as_int(i,j));
          ierr = system.setSchemeParamsThisColumn(mask_value, isMarginal, lambda);
          CHKERRQ(ierr);

          // set boundary values for tridiagonal system
          ierr = system.setSurfaceBoundaryValuesThisColumn(artm(i,j)); CHKERRQ(ierr);
          ierr = system.setBasalBoundaryValuesThisColumn(G0(i,j),shelfbtemp(i,j),(*Rb)(i,j)); CHKERRQ(ierr);

          // solve the system for this column; melting not addressed yet
          PetscErrorCode pivoterr;
          ierr = system.solveThisColumn(&x, pivoterr); CHKERRQ(ierr);

          if (pivoterr != 0) {
            ierr = PetscPrintf(PETSC_COMM_SELF,
              "\n\ntridiagonal solve of tempSystemCtx in temperatureStep() FAILED at (%d,%d)\n"
                  " with zero pivot position %d; viewing system to m-file ... \n",
              i, j, pivoterr); CHKERRQ(ierr);
            ierr = system.reportColumnZeroPivotErrorMFile(pivoterr); CHKERRQ(ierr);
            SETERRQ(grid.com, 1,"PISM ERROR in temperatureStep()\n");
          }
          if (viewOneColumn && issounding(i,j)) {
            ierr = PetscPrintf(grid.com,
              "\n\nin temperatureStep(): viewing tempSystemCtx at (i,j)=(%d,%d) to m-file ... \n\n",
              i, j); CHKERRQ(ierr);
            ierr = system.viewColumnInfoMFile(x, fMz); CHKERRQ(ierr);
          }

        }	// end of "if there are enough points in ice to bother ..."

        // prepare for melting/refreezing
        PetscScalar bwatnew = bwatcurr(i,j);

        // insert solution for generic ice segments
        for (PetscInt k=1; k <= ks; k++) {
          if (allowAboveMelting == PETSC_TRUE) { // in the ice
            Tnew[k] = x[k];
          } else {
            const PetscScalar
              Tpmp = melting_point_temp - beta_CC_grad * (vH(i,j) - fzlev[k]); // FIXME issue #15
            if (x[k] > Tpmp) {
              Tnew[k] = Tpmp;
              PetscScalar Texcess = x[k] - Tpmp; // always positive
              excessToFromBasalMeltLayer(ice_rho, ice_c, L, fzlev[k], fdz, &Texcess, &bwatnew);
              // Texcess  will always come back zero here; ignore it
            } else {
              Tnew[k] = x[k];
            }
          }
          if (Tnew[k] < globalMinAllowedTemp) {
            ierr = PetscPrintf(PETSC_COMM_SELF,
                               "  [[too low (<200) ice segment temp T = %f at %d,%d,%d;"
                               " proc %d; mask=%d; w=%f m/a]]\n",
                               Tnew[k],i,j,k,grid.rank,vMask.as_int(i,j),
                               convert(system.w[k], "m/s", "m/year")); CHKERRQ(ierr);
            myLowTempCount++;
          }
          if (Tnew[k] < artm(i,j) - bulgeMax) {
            Tnew[k] = artm(i,j) - bulgeMax;  bulgeCount++;   }
        }

        // insert solution for ice base segment
        if (ks > 0) {
          if (allowAboveMelting == PETSC_TRUE) { // ice/rock interface
            Tnew[0] = x[0];
          } else {  // compute diff between x[k0] and Tpmp; melt or refreeze as appropriate
            const PetscScalar Tpmp = melting_point_temp - beta_CC_grad * vH(i,j); // FIXME issue #15
            PetscScalar Texcess = x[0] - Tpmp; // positive or negative
            if (mask.ocean(i,j)) {
              // when floating, only half a segment has had its temperature raised
              // above Tpmp
              excessToFromBasalMeltLayer(ice_rho, ice_c, L, 0.0, fdz/2.0, &Texcess, &bwatnew);
            } else {
              excessToFromBasalMeltLayer(ice_rho, ice_c, L, 0.0, fdz, &Texcess, &bwatnew);
            }
            Tnew[0] = Tpmp + Texcess;
            if (Tnew[0] > (Tpmp + 0.00001)) {
              SETERRQ(grid.com, 1,"updated temperature came out above Tpmp");
            }
          }
          if (Tnew[0] < globalMinAllowedTemp) {
            ierr = PetscPrintf(PETSC_COMM_SELF,
                               "  [[too low (<200) ice/bedrock segment temp T = %f at %d,%d;"
                               " proc %d; mask=%d; w=%f]]\n",
                               Tnew[0],i,j,grid.rank,vMask.as_int(i,j),
                               convert(system.w[0], "m/s", "m/year")); CHKERRQ(ierr);
            myLowTempCount++;
          }
          if (Tnew[0] < artm(i,j) - bulgeMax) {
            Tnew[0] = artm(i,j) - bulgeMax;   bulgeCount++;   }
        }

        // set to air temp above ice
        for (PetscInt k=ks; k<fMz; k++) {
          Tnew[k] = artm(i,j);
        }

        // transfer column into vWork3d; communication later
        ierr = vWork3d.setValColumnPL(i,j,Tnew); CHKERRQ(ierr);

        // basalMeltRate[][] is rate of mass loss at bottom of ice; finalize it
        //   note massContExplicitStep() calls PISMOceanCoupler; FIXME: does there
        //   need to be a check that shelfbmassflux(i,j) is up to date?
        if (mask.ocean(i,j)) {
          if (mask.icy(i,j)) {
            // rate of mass loss at bottom of ice shelf;  can be negative (marine freeze-on)
            basalMeltRate[i][j] = shelfbmassflux(i,j); // set by PISMOceanCoupler
          } else {
            basalMeltRate[i][j] = 0.0;
          }
        } else {
          // basalMeltRate is rate of change of bwat;  can be negative
          //   (subglacial water freezes-on); note this rate is calculated
          //   *before* limiting or other nontrivial modelling of bwat,
          //   which is PISMHydrology's job
          basalMeltRate[i][j] = (bwatnew - bwatcurr(i,j)) / dt_TempAge;
        }

    }
  }

  if (myLowTempCount > maxLowTempCount) { SETERRQ(grid.com, 1,"too many low temps"); }

  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = Rb->end_access(); CHKERRQ(ierr);
  ierr = G0.end_access(); CHKERRQ(ierr);
  ierr = bwatcurr.end_access(); CHKERRQ(ierr);
  ierr = vbmr.end_access(); CHKERRQ(ierr);

  ierr = artm.end_access(); CHKERRQ(ierr);

  ierr = shelfbmassflux.end_access(); CHKERRQ(ierr);
  ierr = shelfbtemp.end_access(); CHKERRQ(ierr);

  ierr = u3->end_access(); CHKERRQ(ierr);
  ierr = v3->end_access(); CHKERRQ(ierr);
  ierr = w3->end_access(); CHKERRQ(ierr);
  ierr = Sigma3->end_access(); CHKERRQ(ierr);
  ierr = T3.end_access(); CHKERRQ(ierr);
  ierr = vWork3d.end_access(); CHKERRQ(ierr);

  delete [] x;
  delete [] system.T;  delete [] system.Sigma;
  delete [] system.u;  delete [] system.v;  delete [] system.w;
  delete [] Tnew;
  return 0;
}

