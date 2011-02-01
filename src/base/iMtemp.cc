// Copyright (C) 2004-2011 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <petscda.h>
#include "iceModelVec.hh"
#include "columnSystem.hh"
#include "iceModel.hh"


// it would be reasonable for code clarity to split this into three
// small source files, "iMenergy.cc", "iMtemperature.cc", "iMage.cc"


//! Manage the solution of the energy equation, and related parallel communication.
/*!
Normally calls the method enthalpyAndDrainageStep().  Calls temperatureStep() if
do_cold_ice_methods == true.

This method (energyStep()) \e must update these four fields:
  - IceModelVec3 Enth3
  - IceModelVec3Bedrock Tb3
  - IceModelVec2 vbasalMeltRate
  - IceModelVec2 vHmelt
That is, energyStep() is in charge of calling other methods that update, and
then it is in charge of doing the ghost communication as needed.  If
do_cold_ice_methods == true, then energyStep() must also update this field
  - IceModelVec3 T3
 */
PetscErrorCode IceModel::energyStep() {
  PetscErrorCode  ierr;

  PetscScalar  myCFLviolcount = 0.0,   // these are counts but they are type "PetscScalar"
               myVertSacrCount = 0.0,  //   because that type works with PetscGlobalSum()
               myBulgeCount = 0.0;
  PetscScalar gVertSacrCount, gBulgeCount;

  // always count CFL violations for sanity check (but can occur only if -skip N with N>1)
  ierr = countCFLViolations(&myCFLviolcount); CHKERRQ(ierr);

  if (config.get_flag("do_cold_ice_methods")) {
    // new temperature values go in vTnew; also updates Hmelt:
    ierr = temperatureStep(&myVertSacrCount,&myBulgeCount); CHKERRQ(ierr);  

    ierr = T3.beginGhostCommTransfer(vWork3d); CHKERRQ(ierr);
    ierr = T3.endGhostCommTransfer(vWork3d); CHKERRQ(ierr);

    // compute_enthalpy_cold() updates ghosts of Enth3 using a
    // begin/endGhostComm pair. Is not optimized because this
    // (do_cold_ice_methods) is a rare case.
    ierr = compute_enthalpy_cold(T3, Enth3);  CHKERRQ(ierr);

  } else {
    // new enthalpy values go in vWork3d; also updates (and communicates) Hmelt
    PetscScalar myLiquifiedVol = 0.0, gLiquifiedVol;

    ierr = enthalpyAndDrainageStep(&myVertSacrCount,&myLiquifiedVol,&myBulgeCount);
       CHKERRQ(ierr);

    ierr = Enth3.beginGhostCommTransfer(vWork3d); CHKERRQ(ierr);
    ierr = Enth3.endGhostCommTransfer(vWork3d); CHKERRQ(ierr);

    ierr = PetscGlobalSum(&myLiquifiedVol, &gLiquifiedVol, grid.com); CHKERRQ(ierr);
    if (gLiquifiedVol > 0.0) {
      ierr = verbPrintf(1,grid.com,
        "\n PISM WARNING: fully-liquified cells detected: volume liquified = %.3f km^3\n\n",
        gLiquifiedVol / 1.0e9); CHKERRQ(ierr);
    }
  }

  ierr = PetscGlobalSum(&myCFLviolcount, &CFLviolcount, grid.com); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&myVertSacrCount, &gVertSacrCount, grid.com); CHKERRQ(ierr);
  if (gVertSacrCount > 0.0) { // count of when BOMBPROOF switches to lower accuracy
    const PetscScalar bfsacrPRCNT = 100.0 * (gVertSacrCount / (grid.Mx * grid.My));
    const PetscScalar BPSACR_REPORT_VERB2_PERCENT = 5.0; // only report if above 5%
    if (   (bfsacrPRCNT > BPSACR_REPORT_VERB2_PERCENT) 
        && (getVerbosityLevel() > 2)                    ) {
      char tempstr[50] = "";
      snprintf(tempstr,50, "  [BPsacr=%.4f%%] ", bfsacrPRCNT);
      stdout_flags = tempstr + stdout_flags;
    }
  }

  ierr = PetscGlobalSum(&myBulgeCount, &gBulgeCount, grid.com); CHKERRQ(ierr);
  if (gBulgeCount > 0.0) {   // count of when advection bulges are limited;
                             //    frequently it is identically zero
    char tempstr[50] = "";
    snprintf(tempstr,50, " BULGE=%d ", static_cast<int>(ceil(gBulgeCount)) );
    stdout_flags = tempstr + stdout_flags;
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

    In this procedure four scalar fields are modified: vHmelt, vbmr, Tb3, and vWork3d.
    But vbmr and Tb3 will never need to communicate ghosted values (horizontal 
                                                                    stencil neighbors).  The ghosted values for T3 are updated from the values in vWork3d in the
  communication done by energyStep().  There is a diffusion model for vHmelt in 
  diffuseHmelt() which does communication for vHmelt.

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

    // set up fine grid in ice and bedrock
    PetscInt    fMz = grid.Mz_fine,
      fMbz = grid.Mbz_fine;
    PetscScalar fdz = grid.dz_fine,
      fdzb = fdz,
      *fzlev = grid.zlevels_fine;

    ierr = verbPrintf(5,grid.com,
                      "\n  [entering temperatureStep(); fMz = %d, fdz = %5.3f, fMbz = %d, fdzb = %5.3f]",
                      fMz, fdz, fMbz, fdzb); CHKERRQ(ierr);

    // diagnostic/DEBUG; added for comparison to IceEnthalpyModel
    bool viewOneColumn;
    ierr = PISMOptionsIsSet("-view_sys", viewOneColumn); CHKERRQ(ierr);

    tempSystemCtx system(fMz,fMbz);
    system.dx              = grid.dx;
    system.dy              = grid.dy;
    system.dtTemp          = dtTempAge; // same time step for temp and age, currently
    system.dzEQ            = fdz;
    system.dzbEQ           = fdzb;
    system.ice_rho         = ice->rho;
    system.ice_c_p         = ice->c_p;
    system.ice_k           = ice->k;
    system.bed_thermal_rho = config.get("bedrock_thermal_density");
    system.bed_thermal_c_p = config.get("bedrock_thermal_specific_heat_capacity");
    system.bed_thermal_k   = config.get("bedrock_thermal_conductivity");

    const PetscInt k0 = fMbz - 1;
    PetscScalar *x;  
    x = new PetscScalar[fMz + k0]; // space for solution of system; length = fMz + fMbz - 1 

    // constants needed after solution of system, in insertion phase
    const PetscScalar rho_c_I = ice->rho * ice->c_p,
      rho_c_br = config.get("bedrock_thermal_density") * config.get("bedrock_thermal_specific_heat_capacity"),
      rho_c_av = (fdz * rho_c_I + fdzb * rho_c_br) / (fdz + fdzb);
    // this is bulge limit constant in K; is max amount by which ice
    //   or bedrock can be lower than surface temperature
    const PetscScalar bulgeMax  = config.get("enthalpy_cold_bulge_max") / ice->c_p;

    PetscScalar *Tnew, *Tbnew;
    // pointers to values in current column
    system.u     = new PetscScalar[fMz];
    system.v     = new PetscScalar[fMz];
    system.w     = new PetscScalar[fMz];
    system.Sigma = new PetscScalar[fMz];
    system.T     = new PetscScalar[fMz];
    Tnew         = new PetscScalar[fMz];

    system.Tb    = new PetscScalar[fMbz];
    Tbnew        = new PetscScalar[fMbz];
  
    // system needs access to T3 for T3.getPlaneStar_fine()
    system.T3 = &T3;

    // checks that all needed constants and pointers got set:
    ierr = system.initAllColumns(); CHKERRQ(ierr);

    // now get map-plane fields, starting with coupler fields
    PetscScalar  **Hmelt, **basalMeltRate;
  
    if (surface != PETSC_NULL) {
      ierr = surface->ice_surface_temperature(grid.year, dtTempAge / secpera, artm); CHKERRQ(ierr);
    } else {
      SETERRQ(1,"PISM ERROR: surface == PETSC_NULL");
    }
    if (ocean != PETSC_NULL) {
      ierr = ocean->shelf_base_mass_flux(grid.year, dtTempAge / secpera, shelfbmassflux); CHKERRQ(ierr);
      ierr = ocean->shelf_base_temperature(grid.year, dtTempAge / secpera, shelfbtemp); CHKERRQ(ierr);
    } else {
      SETERRQ(1,"PISM ERROR: ocean == PETSC_NULL");
    }

    ierr = artm.begin_access(); CHKERRQ(ierr);
    ierr = shelfbmassflux.begin_access(); CHKERRQ(ierr);
    ierr = shelfbtemp.begin_access(); CHKERRQ(ierr);

    ierr = vH.begin_access(); CHKERRQ(ierr);
    ierr = vHmelt.get_array(Hmelt); CHKERRQ(ierr);
    ierr = vbmr.get_array(basalMeltRate); CHKERRQ(ierr);
    ierr = vMask.begin_access(); CHKERRQ(ierr);
    ierr = vGhf.begin_access(); CHKERRQ(ierr);

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

    ierr = Tb3.begin_access(); CHKERRQ(ierr);

    // counts unreasonably low temperature values; deprecated?
    PetscInt myLowTempCount = 0;
    PetscInt maxLowTempCount = static_cast<PetscInt>(config.get("max_low_temp_count"));
    PetscReal globalMinAllowedTemp = config.get("global_min_allowed_temp");

    PetscReal hmelt_max = config.get("hmelt_max");

    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

        // this should *not* be replaced by call to grid.kBelowHeight():
        const PetscInt  ks = static_cast<PetscInt>(floor(vH(i,j)/fdz));
      
        if (k0+ks>0) { // if there are enough points in bedrock&ice to bother ...
          ierr = system.setIndicesAndClearThisColumn(i,j,ks); CHKERRQ(ierr);

          ierr = Tb3.getValColumnPL(i,j,system.Tb); CHKERRQ(ierr);

          ierr = u3->getValColumn(i,j,ks,system.u); CHKERRQ(ierr);
          ierr = v3->getValColumn(i,j,ks,system.v); CHKERRQ(ierr);
          ierr = w3->getValColumn(i,j,ks,system.w); CHKERRQ(ierr);
          ierr = Sigma3->getValColumn(i,j,ks,system.Sigma); CHKERRQ(ierr);
          ierr = T3.getValColumn(i,j,ks,system.T); CHKERRQ(ierr);

          // go through column and find appropriate lambda for BOMBPROOF
          PetscScalar lambda = 1.0;  // start with centered implicit for more accuracy
          for (PetscInt k = 1; k < ks; k++) {   
            const PetscScalar denom = (PetscAbs(system.w[k]) + 0.000001/secpera)
              * ice->rho * ice->c_p * fdz;  
            lambda = PetscMin(lambda, 2.0 * ice->k / denom);
          }
          if (lambda < 1.0)  *vertSacrCount += 1; // count columns with lambda < 1
          // if isMarginal then only do vertical conduction for ice; ignore advection
          //   and strain heating if isMarginal
          const bool isMarginal = checkThinNeigh(vH(i+1,j),vH(i+1,j+1),vH(i,j+1),vH(i-1,j+1),
                                                 vH(i-1,j),vH(i-1,j-1),vH(i,j-1),vH(i+1,j-1));
          ierr = system.setSchemeParamsThisColumn(vMask.value(i,j), isMarginal, lambda);
          CHKERRQ(ierr);  

          // set boundary values for tridiagonal system
          ierr = system.setSurfaceBoundaryValuesThisColumn(artm(i,j)); CHKERRQ(ierr);
          ierr = system.setBasalBoundaryValuesThisColumn(vGhf(i,j),shelfbtemp(i,j),(*Rb)(i,j)); CHKERRQ(ierr);

          // solve the system for this column; melting not addressed yet
          ierr = system.solveThisColumn(&x); // no CHKERRQ(ierr) immediately because:
          if (ierr > 0) {
            PetscPrintf(grid.com,
                        "PISM ERROR in IceModel::temperatureStep():\n"
                        "   tridiagonal solve failed at (%d,%d) with zero pivot position %d\n"
                        "   1-norm = %.3e  and  diagonal-dominance ratio = %.5f\n"
                        "   ENDING! ...\n",
                        i, j, ierr, system.norm1(k0+ks+1), system.ddratio(k0+ks+1));
            PISMEnd();
          } else { CHKERRQ(ierr); }

          // diagnostic/DEBUG; added for comparison to IceEnthalpyModel
          if (viewOneColumn) {
            if ((i==id) && (j==jd)) {
              ierr = verbPrintf(1,grid.com,
                                "\nin IceModel::temperatureStep();\n"
                                "   fMz = %d, fdz = %5.3f, fMbz = %d, fdzb = %5.3f, k0 = %d, ks = %d\n\n",
                                fMz, fdz, fMbz, fdzb, k0, ks); CHKERRQ(ierr);
              ierr = verbPrintf(1,grid.com,
                                "viewing system and solution at (i,j)=(%d,%d):\n", i, j); CHKERRQ(ierr);
              ierr = verbPrintf(1,grid.com,
                                "   1-norm = %.3e  and  diagonal-dominance ratio = %.5f\n",
                                system.norm1(k0+ks+1), system.ddratio(k0+ks+1)); CHKERRQ(ierr);
              ierr = system.viewSystem(NULL,"system"); CHKERRQ(ierr);
              ierr = system.viewColumnValues(NULL, x, fMz+k0, "solution x"); CHKERRQ(ierr);
            }
          }

        }	// end of "// if there are enough points in bedrock&ice to bother ..."

        // insert bedrock solution; check for too low below
        for (PetscInt k=0; k < k0; k++) {
          Tbnew[k] = x[k];
        }

        // prepare for melting/refreezing
        PetscScalar Hmeltnew = Hmelt[i][j];
      
        // insert solution for generic ice segments
        for (PetscInt k=1; k <= ks; k++) {
          if (allowAboveMelting == PETSC_TRUE) { // in the ice
            Tnew[k] = x[k0 + k];
          } else {
            const PetscScalar 
              Tpmp = ice->triple_point_temp - ice->beta_CC_grad * (vH(i,j) - fzlev[k]); // FIXME task #7297
            if (x[k0 + k] > Tpmp) {
              Tnew[k] = Tpmp;
              PetscScalar Texcess = x[k0 + k] - Tpmp; // always positive
              excessToFromBasalMeltLayer(rho_c_I, fzlev[k], fdz, &Texcess, &Hmeltnew);
              // Texcess  will always come back zero here; ignore it
            } else {
              Tnew[k] = x[k0 + k];
            }
          }
          if (Tnew[k] < globalMinAllowedTemp) {
            ierr = PetscPrintf(PETSC_COMM_SELF,
                               "  [[too low (<200) ice segment temp T = %f at %d,%d,%d;"
                               " proc %d; mask=%d; w=%f]]\n",
                               Tnew[k],i,j,k,grid.rank,vMask.value(i,j),system.w[k]*secpera); CHKERRQ(ierr);
            myLowTempCount++;
          }
          if (Tnew[k] < artm(i,j) - bulgeMax) {
            Tnew[k] = artm(i,j) - bulgeMax;  bulgeCount++;   }
        }
      
        // insert solution for ice/rock interface (or base of ice shelf) segment
        if (ks > 0) {
          if (allowAboveMelting == PETSC_TRUE) { // ice/rock interface
            Tnew[0] = x[k0];
          } else {  // compute diff between x[k0] and Tpmp; melt or refreeze as appropriate
            const PetscScalar Tpmp = ice->triple_point_temp - ice->beta_CC_grad * vH(i,j); // FIXME task #7297
            PetscScalar Texcess = x[k0] - Tpmp; // positive or negative
            if (vMask.is_floating(i,j)) {
              // when floating, only half a segment has had its temperature raised
              // above Tpmp
              excessToFromBasalMeltLayer(rho_c_I/2, 0.0, fdz, &Texcess, &Hmeltnew);
            } else {
              excessToFromBasalMeltLayer(rho_c_av, 0.0, fdz, &Texcess, &Hmeltnew);
            }
            Tnew[0] = Tpmp + Texcess;
            if (Tnew[0] > (Tpmp + 0.00001)) {
              SETERRQ(1,"updated temperature came out above Tpmp");
            }
          }
          if (Tnew[0] < globalMinAllowedTemp) {
            ierr = PetscPrintf(PETSC_COMM_SELF,
                               "  [[too low (<200) ice/bedrock segment temp T = %f at %d,%d;"
                               " proc %d; mask=%d; w=%f]]\n",
                               Tnew[0],i,j,grid.rank,vMask.value(i,j),system.w[0]*secpera); CHKERRQ(ierr);
            myLowTempCount++;
          }
          if (Tnew[0] < artm(i,j) - bulgeMax) {
            Tnew[0] = artm(i,j) - bulgeMax;   bulgeCount++;   }
        } else {
          Hmeltnew = 0.0;
        }
      
        // we must agree on redundant values T(z=0) at top of bedrock and at bottom of ice
        if (ks > 0) {
          Tbnew[k0] = Tnew[0];
        } else {
          if (vMask.is_floating(i,j)) { // top of bedrock sees ocean
            Tbnew[k0] = shelfbtemp(i,j); // set by PISMOceanCoupler
          } else { // top of bedrock sees atmosphere
            Tbnew[k0] = artm(i,j);
          }
        }
        // check bedrock solution        
        for (PetscInt k=0; k <= k0; k++) {
          if (Tbnew[k] < globalMinAllowedTemp) {
            ierr = PetscPrintf(PETSC_COMM_SELF,
                               "  [[too low (<200) bedrock segment temp T = %f at %d,%d,%d;"
                               " proc %d; mask=%d]]\n",
                               Tbnew[k],i,j,k,grid.rank,vMask.value(i,j)); CHKERRQ(ierr);
            myLowTempCount++;
          }
          if (Tbnew[k] < artm(i,j) - bulgeMax) {
            Tbnew[k] = artm(i,j) - bulgeMax;   bulgeCount++;   }
        }

        // transfer column into Tb3; neighboring columns will not reference!
        ierr = Tb3.setValColumnPL(i,j,Tbnew); CHKERRQ(ierr);

        // set to air temp above ice
        for (PetscInt k=ks; k<fMz; k++) {
          Tnew[k] = artm(i,j);
        }

        // transfer column into vWork3d; communication later
        ierr = vWork3d.setValColumnPL(i,j,Tnew); CHKERRQ(ierr);

        // basalMeltRate[][] is rate of mass loss at bottom of ice everywhere;
        //   note massContExplicitStep() calls PISMOceanCoupler separately
        if (vMask.is_floating(i,j)) {
          // rate of mass loss at bottom of ice shelf;  can be negative (marine freeze-on)
          basalMeltRate[i][j] = shelfbmassflux(i,j); // set by PISMOceanCoupler
        } else {
          // rate of change of Hmelt[][];  can be negative (till water freeze-on)
          // Also note that this rate is calculated *before* limiting Hmelt.
          basalMeltRate[i][j] = (Hmeltnew - Hmelt[i][j]) / dtTempAge;
        }

        if (vMask.is_floating(i,j)) {
          // if floating assume maximally saturated till
          Hmelt[i][j] = hmelt_max;
        } else {
          // limit Hmelt by default max and store
          Hmelt[i][j] = PetscMin(hmelt_max, Hmeltnew);
      }

    } 
  }
  
  if (myLowTempCount > maxLowTempCount) { SETERRQ(1,"too many low temps"); }

  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vHmelt.end_access(); CHKERRQ(ierr);
  ierr = Rb->end_access(); CHKERRQ(ierr);
  ierr = vGhf.end_access(); CHKERRQ(ierr);
  ierr = vbmr.end_access(); CHKERRQ(ierr);

  ierr = artm.end_access(); CHKERRQ(ierr);

  ierr = shelfbmassflux.end_access(); CHKERRQ(ierr);
  ierr = shelfbtemp.end_access(); CHKERRQ(ierr);

  ierr = Tb3.end_access(); CHKERRQ(ierr);
  ierr = u3->end_access(); CHKERRQ(ierr);
  ierr = v3->end_access(); CHKERRQ(ierr);
  ierr = w3->end_access(); CHKERRQ(ierr);
  ierr = Sigma3->end_access(); CHKERRQ(ierr);
  ierr = T3.end_access(); CHKERRQ(ierr);
  ierr = vWork3d.end_access(); CHKERRQ(ierr);
  
  delete [] x;
  delete [] system.T;  delete [] system.Tb;  
  delete [] system.u;  delete [] system.v;  delete [] system.w;
  delete [] system.Sigma;
  
  delete [] Tbnew;  delete [] Tnew;

  return 0;
}


//! Compute the melt water which should go to the base if \f$T\f$ is above pressure-melting.
PetscErrorCode IceModel::excessToFromBasalMeltLayer(
                const PetscScalar rho_c, const PetscScalar z, const PetscScalar dz,
                PetscScalar *Texcess, PetscScalar *Hmelt) {

  const PetscScalar darea = grid.dx * grid.dy,
                    dvol = darea * dz,
                    dE = rho_c * (*Texcess) * dvol,
                    massmelted = dE / ice->latentHeat;

  if (allowAboveMelting == PETSC_TRUE) {
    SETERRQ(1,"IceModel::excessToBasalMeltLayer() called but allowAboveMelting==TRUE");
  }
  if (*Texcess >= 0.0) {
    if (updateHmelt == PETSC_TRUE) {
      // T is at or above pressure-melting temp, so temp needs to be set to 
      // pressure-melting; a fraction of excess energy is turned into melt water at base
      // note massmelted is POSITIVE!
      const PetscScalar FRACTION_TO_BASE
                           = (z < 100.0) ? 0.2 * (100.0 - z) / 100.0 : 0.0;
      // note: ice-equiv thickness:
      *Hmelt += (FRACTION_TO_BASE * massmelted) / (ice->rho * darea);  
    }
    *Texcess = 0.0;
  } else if (updateHmelt == PETSC_TRUE) {  // neither Texcess nor Hmelt need to change 
                                           // if Texcess < 0.0
    // Texcess negative; only refreeze (i.e. reduce Hmelt) if at base and Hmelt > 0.0
    // note ONLY CALLED IF AT BASE!   note massmelted is NEGATIVE!
    if (z > 0.00001) {
      SETERRQ(1, "excessToBasalMeltLayer() called with z not at base and negative Texcess");
    }
    if (*Hmelt > 0.0) {
      const PetscScalar thicknessToFreezeOn = - massmelted / (ice->rho * darea);
      if (thicknessToFreezeOn <= *Hmelt) { // the water *is* available to freeze on
        *Hmelt -= thicknessToFreezeOn;
        *Texcess = 0.0;
      } else { // only refreeze Hmelt thickness of water; update Texcess
        *Hmelt = 0.0;
        const PetscScalar dTemp = ice->latentHeat * ice->rho * (*Hmelt) / (rho_c * dz);
        *Texcess += dTemp;
      }
    } 
    // note: if *Hmelt == 0 and Texcess < 0.0 then Texcess unmolested; temp will go down
  }
  return 0;
}                           


//! Take a semi-implicit time-step for the age equation.
/*!
The age equation is\f$d\tau/dt = 1\f$, that is,
    \f[ \frac{\partial \tau}{\partial t} + u \frac{\partial \tau}{\partial x}
        + v \frac{\partial \tau}{\partial y} + w \frac{\partial \tau}{\partial z} = 1\f]
where \f$\tau(t,x,y,z)\f$ is the age of the ice and \f$(u,v,w)\f$  is the three dimensional
velocity field.  This equation is purely advective.  And it is hyperbolic.

The boundary condition is that when the ice falls as snow it has age zero.  
That is, \f$\tau(t,x,y,h(t,x,y)) = 0\f$ in accumulation areas, while there is no 
boundary condition elsewhere, as the characteristics go outward in the ablation zone.
(Some more numerical-analytic attention to this is worthwhile.)

If the velocity in the bottom cell of ice is upward (\code (w[i][j][0] > 0 \endcode)
then we also apply an age = 0 boundary condition.  This is the case where ice freezes
on at the base, either grounded basal ice freezing on stored water in till, or marine basal ice.

The numerical method is first-order upwind but the vertical advection term is computed
implicitly.  (Thus there is no CFL-type stability condition for that part.)

We use a finely-spaced, equally-spaced vertical grid in the calculation.  Note that the IceModelVec3 
methods getValColumn...() and setValColumn..() interpolate back and forth between the grid 
on which calculation is done and the storage grid.  Thus the storage grid can be either 
equally spaced or not.
 */
PetscErrorCode IceModel::ageStep() {
  PetscErrorCode  ierr;

  // set up fine grid in ice
  PetscInt    fMz = grid.Mz_fine;
  PetscScalar fdz = grid.dz_fine;

  PetscScalar *x;  
  x = new PetscScalar[fMz]; // space for solution

  ageSystemCtx system(fMz); // linear system to solve in each column
  system.dx    = grid.dx;
  system.dy    = grid.dy;
  system.dtAge = dtTempAge;
  system.dzEQ  = fdz;
  // pointers to values in current column
  system.u     = new PetscScalar[fMz];
  system.v     = new PetscScalar[fMz];
  system.w     = new PetscScalar[fMz];
  // system needs access to tau3 for planeStar()
  system.tau3  = &tau3;
  // this checks that all needed constants and pointers got set
  ierr = system.initAllColumns(); CHKERRQ(ierr);

  IceModelVec3 *u3, *v3, *w3;
  ierr = stress_balance->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr); 

  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = tau3.begin_access(); CHKERRQ(ierr);
  ierr = u3->begin_access(); CHKERRQ(ierr);
  ierr = v3->begin_access(); CHKERRQ(ierr);
  ierr = w3->begin_access(); CHKERRQ(ierr);
  ierr = vWork3d.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // this should *not* be replaced by a call to grid.kBelowHeight()
      const PetscInt  fks = static_cast<PetscInt>(floor(vH(i,j)/fdz));

      if (fks == 0) { // if no ice, set the entire column to zero age
        ierr = vWork3d.setColumn(i,j,0.0); CHKERRQ(ierr);
      } else { // general case: solve advection PDE; start by getting 3D velocity ...

	ierr = u3->getValColumn(i,j,fks,system.u); CHKERRQ(ierr);
	ierr = v3->getValColumn(i,j,fks,system.v); CHKERRQ(ierr);
	ierr = w3->getValColumn(i,j,fks,system.w); CHKERRQ(ierr);

        ierr = system.setIndicesAndClearThisColumn(i,j,fks); CHKERRQ(ierr);

        // solve the system for this column; call checks that params set
        ierr = system.solveThisColumn(&x); // no "CHKERRQ(ierr)" because:
        if (ierr > 0) {
          PetscPrintf(grid.com,
            "PISM ERROR in IceModel::ageStep():\n"
            "  tridiagonal solve failed at (%d,%d) with zero pivot position %d\n"
            "  1-norm = %.3e  and  diagonal-dominance ratio = %.5f\n"
            "  ENDING! ...\n\n",
            i, j, ierr, system.norm1(fks+1), system.ddratio(fks+1));
          PISMEnd();
        } else { CHKERRQ(ierr); }

        // x[k] contains age for k=0,...,ks, but set age of ice above (and at) surface to zero years
        for (PetscInt k=fks+1; k<fMz; k++) {
          x[k] = 0.0;
        }
        
        // put solution in IceModelVec3
        ierr = vWork3d.setValColumnPL(i,j,x); CHKERRQ(ierr);
      }
    }
  }

  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = tau3.end_access();  CHKERRQ(ierr);
  ierr = u3->end_access();  CHKERRQ(ierr);
  ierr = v3->end_access();  CHKERRQ(ierr);
  ierr = w3->end_access();  CHKERRQ(ierr);
  ierr = vWork3d.end_access();  CHKERRQ(ierr);

  delete [] x;  
  delete [] system.u;  delete [] system.v;  delete [] system.w;

  ierr = tau3.beginGhostCommTransfer(vWork3d); CHKERRQ(ierr);
  ierr = tau3.endGhostCommTransfer(vWork3d); CHKERRQ(ierr);

  return 0;
}


bool IceModel::checkThinNeigh(PetscScalar E, PetscScalar NE, PetscScalar N, PetscScalar NW, 
                              PetscScalar W, PetscScalar SW, PetscScalar S, PetscScalar SE) {
  const PetscScalar THIN = 100.0;  // thin = (at most 100m thick)
  return (   (E < THIN) || (NE < THIN) || (N < THIN) || (NW < THIN)
          || (W < THIN) || (SW < THIN) || (S < THIN) || (SE < THIN) );
}

