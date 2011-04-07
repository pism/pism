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

#include "iceModel.hh"
#include "tempSystem.hh"


// FIXME:  it might be reasonable for code clarity to split this in two:
//   "iMenergy.cc", "iMtemperature.cc"


//! Manage the solution of the energy equation, and related parallel communication.
/*!
This method updates three fields:
  - IceModelVec3 Enth3
  - IceModelVec2 vbmr
  - IceModelVec2 vHmelt
That is, energyStep() is in charge of calling other methods that actually update, and
then it is in charge of doing the ghost communication as needed.  If
do_cold_ice_methods == true, then energyStep() must also update this field
  - IceModelVec3 T3

Normally calls the method enthalpyAndDrainageStep().  Calls temperatureStep() if
do_cold_ice_methods == true.
 */
PetscErrorCode IceModel::energyStep() {
  PetscErrorCode  ierr;

  PetscScalar  myCFLviolcount = 0.0,   // these are counts but they are type "PetscScalar"
               myVertSacrCount = 0.0,  //   because that type works with PetscGlobalSum()
               myBulgeCount = 0.0;
  PetscScalar gVertSacrCount, gBulgeCount;

  // always count CFL violations for sanity check (but can occur only if -skip N with N>1)
  ierr = countCFLViolations(&myCFLviolcount); CHKERRQ(ierr);

  // operator-splitting occurs here (ice and bedrock energy updates are split):
  //   tell PISMBedThermalUnit* btu that we have an ice base temp; it will return
  //   the z=0 value of geothermal flux when called inside temperatureStep() or
  //   enthalpyAndDrainageStep()
  ierr = get_bed_top_temp(bedtoptemp); CHKERRQ(ierr);
  ierr = btu->update(t_years_TempAge, dt_years_TempAge); CHKERRQ(ierr);  // has ptr to bedtoptemp

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

  // Both cases above update the basal melt rate field; here we update its
  // ghosts, which are needed to compute tauc locally
  ierr = vbmr.beginGhostComm(); CHKERRQ(ierr);
  ierr = vbmr.endGhostComm(); CHKERRQ(ierr);

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


PetscErrorCode IceModel::get_bed_top_temp(IceModelVec2S &result) {
  PetscErrorCode  ierr;

  // will need coupler fields in ice-free land and 
  if (surface != PETSC_NULL) {
    ierr = surface->ice_surface_temperature(artm); CHKERRQ(ierr);
  } else {
    SETERRQ(1,"PISM ERROR: surface == PETSC_NULL");
  }
  if (ocean != PETSC_NULL) {
    ierr = ocean->shelf_base_temperature(shelfbtemp); CHKERRQ(ierr);
  } else {
    SETERRQ(5,"PISM ERROR: ocean == PETSC_NULL");
  }

  // start by grabbing 2D basal enthalpy field at z=0; converted to temperature if needed, below
  ierr = Enth3.getHorSlice(result, 0.0); CHKERRQ(ierr);

  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = artm.begin_access(); CHKERRQ(ierr);
  ierr = shelfbtemp.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if (vMask.is_grounded(i,j)) {
        if (vMask.value(i,j) == MASK_ICE_FREE_BEDROCK) { // no ice: sees air temp
          result(i,j) = artm(i,j);
        } else { // ice: sees temp of base of ice
          const PetscReal pressure = EC->getPressureFromDepth(vH(i,j));
          ierr = EC->getAbsTemp(result(i,j), pressure, result(i,j)); CHKERRQ(ierr);  
        }
      } else { // floating: apply shelf base temp as top of bedrock temp
        result(i,j) = shelfbtemp(i,j);
      }
    }
  }
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = artm.end_access(); CHKERRQ(ierr);
  ierr = shelfbtemp.end_access(); CHKERRQ(ierr);

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


//! \brief Is one of my neighbors below a critical thickness?
/*!
Used in ice energy methods, both enthalpy and temperature, to determine whether
advection should be included right at edge of ice sheet.

FIXME: silly hard-wired critical level, but we want to avoid config.get() in loops.
 */
bool IceModel::checkThinNeigh(PetscScalar E, PetscScalar NE, PetscScalar N, PetscScalar NW, 
                              PetscScalar W, PetscScalar SW, PetscScalar S, PetscScalar SE) {
  const PetscScalar THIN = 100.0;  // thin = (at most 100m thick)
  return (   (E < THIN) || (NE < THIN) || (N < THIN) || (NW < THIN)
          || (W < THIN) || (SW < THIN) || (S < THIN) || (SE < THIN) );
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

    In this procedure four scalar fields are modified: vHmelt, vbmr, and vWork3d.
    But vbmr will never need to communicate ghosted values (horizontal stencil
    neighbors).  The ghosted values for T3 are updated from the values in vWork3d in the
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

    // set up fine grid in ice
    PetscInt    fMz    = grid.Mz_fine;
    PetscScalar fdz    = grid.dz_fine;
    vector<double> &fzlev = grid.zlevels_fine;

    ierr = verbPrintf(5,grid.com,
                      "\n  [entering temperatureStep(); fMz = %d, fdz = %5.3f]",
                      fMz, fdz); CHKERRQ(ierr);

    bool viewOneColumn;
    ierr = PISMOptionsIsSet("-view_sys", viewOneColumn); CHKERRQ(ierr);

    tempSystemCtx system(fMz);
    system.dx              = grid.dx;
    system.dy              = grid.dy;
    system.dtTemp          = dt_years_TempAge * secpera; // same time step for temp and age, currently
    system.dzEQ            = fdz;
    system.ice_rho         = ice->rho;
    system.ice_c_p         = ice->c_p;
    system.ice_k           = ice->k;

    PetscScalar *x;  
    x = new PetscScalar[fMz]; // space for solution of system

    // constants needed after solution of system, in insertion phase
    const PetscScalar rho_c_I = ice->rho * ice->c_p;

    // this is bulge limit constant in K; is max amount by which ice
    //   or bedrock can be lower than surface temperature
    const PetscScalar bulgeMax  = config.get("enthalpy_cold_bulge_max") / ice->c_p;

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
    PetscScalar  **Hmelt, **basalMeltRate;
  
    if (surface != PETSC_NULL) {
      ierr = surface->ice_surface_temperature(artm); CHKERRQ(ierr);
    } else {
      SETERRQ(1,"PISM ERROR: surface == PETSC_NULL");
    }
    if (ocean != PETSC_NULL) {
      ierr = ocean->shelf_base_mass_flux(shelfbmassflux); CHKERRQ(ierr);
      ierr = ocean->shelf_base_temperature(shelfbtemp); CHKERRQ(ierr);
    } else {
      SETERRQ(1,"PISM ERROR: ocean == PETSC_NULL");
    }

    IceModelVec2S G0 = vWork2d[0];
    ierr = G0.set_attrs("internal", "upward geothermal flux at z=0", "W m-2", ""); CHKERRQ(ierr);
    ierr = G0.set_glaciological_units("mW m-2");
    if (btu) {
      ierr = btu->get_upward_geothermal_flux(G0); CHKERRQ(ierr);
    } else {
      SETERRQ(3,"PISM ERROR: PISMBedThermalUnit* btu == PETSC_NULL in temperatureStep()");
    }

    ierr = artm.begin_access(); CHKERRQ(ierr);
    ierr = shelfbmassflux.begin_access(); CHKERRQ(ierr);
    ierr = shelfbtemp.begin_access(); CHKERRQ(ierr);

    ierr = vH.begin_access(); CHKERRQ(ierr);
    ierr = vHmelt.get_array(Hmelt); CHKERRQ(ierr);
    ierr = vbmr.get_array(basalMeltRate); CHKERRQ(ierr);
    ierr = vMask.begin_access(); CHKERRQ(ierr);
    ierr = G0.begin_access(); CHKERRQ(ierr);

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

    PetscReal hmelt_max = config.get("hmelt_max");

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
              * ice->rho * ice->c_p * fdz;  
            lambda = PetscMin(lambda, 2.0 * ice->k / denom);
          }
          if (lambda < 1.0)  *vertSacrCount += 1; // count columns with lambda < 1
          // if isMarginal then only do vertical conduction for ice; ignore advection
          //   and strain heating if isMarginal
          const bool isMarginal = checkThinNeigh(vH(i+1,j),vH(i+1,j+1),vH(i,j+1),vH(i-1,j+1),
                                                 vH(i-1,j),vH(i-1,j-1),vH(i,j-1),vH(i+1,j-1));
          PismMask mask_value = static_cast<PismMask>(vMask.value(i,j));
          ierr = system.setSchemeParamsThisColumn(mask_value, isMarginal, lambda);
          CHKERRQ(ierr);  

          // set boundary values for tridiagonal system
          ierr = system.setSurfaceBoundaryValuesThisColumn(artm(i,j)); CHKERRQ(ierr);
          ierr = system.setBasalBoundaryValuesThisColumn(G0(i,j),shelfbtemp(i,j),(*Rb)(i,j)); CHKERRQ(ierr);

          // solve the system for this column; melting not addressed yet
          ierr = system.solveThisColumn(&x); // no CHKERRQ(ierr) immediately because:
          if (ierr > 0) {
            PetscPrintf(grid.com,
                        "PISM ERROR in IceModel::temperatureStep():\n"
                        "   tridiagonal solve failed at (%d,%d) with zero pivot position %d\n"
                        "   1-norm = %.3e  and  diagonal-dominance ratio = %.5f\n"
                        "   ENDING! ...\n",
                        i, j, ierr, system.norm1(ks+1), system.ddratio(ks+1));
            PISMEnd();
          } else { CHKERRQ(ierr); }

          // diagnostic/DEBUG; added for comparison to IceEnthalpyModel
          if (viewOneColumn) {
            if ((i==id) && (j==jd)) {
              ierr = verbPrintf(1,grid.com,
                                "\nin IceModel::temperatureStep();\n"
                                "   fMz = %d, fdz = %5.3f, ks = %d\n\n",
                                fMz, fdz, ks); CHKERRQ(ierr);
              ierr = verbPrintf(1,grid.com,
                                "viewing system and solution at (i,j)=(%d,%d):\n", i, j); CHKERRQ(ierr);
              ierr = verbPrintf(1,grid.com,
                                "   1-norm = %.3e  and  diagonal-dominance ratio = %.5f\n",
                                system.norm1(ks+1), system.ddratio(ks+1)); CHKERRQ(ierr);
              ierr = system.viewSystem(NULL,"system"); CHKERRQ(ierr);
              ierr = system.viewColumnValues(NULL, x, fMz, "solution x"); CHKERRQ(ierr);
            }
          }

        }	// end of "if there are enough points in ice to bother ..."

        // prepare for melting/refreezing
        PetscScalar Hmeltnew = Hmelt[i][j];
      
        // insert solution for generic ice segments
        for (PetscInt k=1; k <= ks; k++) {
          if (allowAboveMelting == PETSC_TRUE) { // in the ice
            Tnew[k] = x[k];
          } else {
            const PetscScalar 
              Tpmp = ice->triple_point_temp - ice->beta_CC_grad * (vH(i,j) - fzlev[k]); // FIXME task #7297
            if (x[k] > Tpmp) {
              Tnew[k] = Tpmp;
              PetscScalar Texcess = x[k] - Tpmp; // always positive
              excessToFromBasalMeltLayer(rho_c_I, fzlev[k], fdz, &Texcess, &Hmeltnew);
              // Texcess  will always come back zero here; ignore it
            } else {
              Tnew[k] = x[k];
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
      
        // insert solution for ice base segment
        if (ks > 0) {
          if (allowAboveMelting == PETSC_TRUE) { // ice/rock interface
            Tnew[0] = x[0];
          } else {  // compute diff between x[k0] and Tpmp; melt or refreeze as appropriate
            const PetscScalar Tpmp = ice->triple_point_temp - ice->beta_CC_grad * vH(i,j); // FIXME task #7297
            PetscScalar Texcess = x[0] - Tpmp; // positive or negative
            if (vMask.is_floating(i,j)) {
              // when floating, only half a segment has had its temperature raised
              // above Tpmp
              excessToFromBasalMeltLayer(rho_c_I/2, 0.0, fdz, &Texcess, &Hmeltnew);
            } else {
              excessToFromBasalMeltLayer(rho_c_I, 0.0, fdz, &Texcess, &Hmeltnew);
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
          basalMeltRate[i][j] = (Hmeltnew - Hmelt[i][j]) / (dt_years_TempAge * secpera);
        }

        if (vMask.is_floating(i,j)) {
          if (vH(i,j) < 0.1) {
            // truely no ice, so zero-out subglacial fields
            Hmelt[i][j] = 0.0;
          } else {
            // if floating assume maximally saturated till to avoid "shock" if grounding line advances
            Hmelt[i][j] = hmelt_max;
          }
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
  ierr = G0.end_access(); CHKERRQ(ierr);
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

