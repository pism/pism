// Copyright (C) 2009-2011 Andreas Aschwanden and Ed Bueler and Constantine Khroulev
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
#include "enthSystem.hh"

//! Compute Enth3 from temperature T3 by assuming the ice has zero liquid fraction.
/*!
First this method makes sure the temperatures is at most the pressure-melting
value, before computing the enthalpy for that temperature, using zero liquid
fraction.

Because of how EnthalpyConverter::getPressureFromDepth() works, the energy 
content in the air is set to the value that ice would have if it a chunk of it
occupied the air; the atmosphere actually has much lower energy content.  It is
done this way for regularity (i.e. dEnth/dz computations).

Because Enth3 gets set, does ghost communication to finish.
 */
PetscErrorCode IceModel::compute_enthalpy_cold(IceModelVec3 &temperature, IceModelVec3 &result) {
  PetscErrorCode ierr;
  
  ierr = temperature.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);

  PetscScalar *Tij, *Enthij; // columns of these values
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = temperature.getInternalColumn(i,j,&Tij); CHKERRQ(ierr);
      ierr = result.getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
      for (PetscInt k=0; k<grid.Mz; ++k) {
        const PetscScalar depth = vH(i,j) - grid.zlevels[k]; // FIXME task #7297
        ierr = EC->getEnthPermissive(Tij[k],0.0,EC->getPressureFromDepth(depth),
                                    Enthij[k]); CHKERRQ(ierr);
      }
    }
  }

  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = temperature.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = result.beginGhostComm(); CHKERRQ(ierr);
  ierr = result.endGhostComm(); CHKERRQ(ierr);
  return 0;
}


//! Compute Enth3 from temperature T3 and liquid fraction.
/*!
Because Enth3 gets set, does ghost communication to finish.
 */
PetscErrorCode IceModel::compute_enthalpy(IceModelVec3 &temperature,
                                          IceModelVec3 &liquid_water_fraction,
                                          IceModelVec3 &result) {
  PetscErrorCode ierr;
  
  ierr = temperature.begin_access(); CHKERRQ(ierr);
  ierr = liquid_water_fraction.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);

  PetscScalar *Tij, *Liqfracij, *Enthij; // columns of these values
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = temperature.getInternalColumn(i,j,&Tij); CHKERRQ(ierr);
      ierr = liquid_water_fraction.getInternalColumn(i,j,&Liqfracij); CHKERRQ(ierr);
      ierr = result.getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
      for (PetscInt k=0; k<grid.Mz; ++k) {
        const PetscScalar depth = vH(i,j) - grid.zlevels[k]; // FIXME task #7297
        ierr = EC->getEnthPermissive(Tij[k],Liqfracij[k],
                      EC->getPressureFromDepth(depth), Enthij[k]); CHKERRQ(ierr);
      }
    }
  }

  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = temperature.end_access(); CHKERRQ(ierr);
  ierr = liquid_water_fraction.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = result.beginGhostComm(); CHKERRQ(ierr);
  ierr = result.endGhostComm(); CHKERRQ(ierr);
  return 0;
}

//! Compute the liquid fraction corresponding to Enth3, and put in a global IceModelVec3 provided by user.
/*!
Does not communicate ghosts for IceModelVec3 useForLiquidFrac.
 */
PetscErrorCode IceModel::compute_liquid_water_fraction(IceModelVec3 &enthalpy,
                                                       IceModelVec3 &result) {
  PetscErrorCode ierr;

  ierr = result.set_name("liqfrac"); CHKERRQ(ierr);
  ierr = result.set_attrs(
     "diagnostic",
     "liquid water fraction in ice (between 0 and 1)",
     "", ""); CHKERRQ(ierr);

  PetscScalar *omegaij, *Enthij; // columns of these values
  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = enthalpy.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = result.getInternalColumn(i,j,&omegaij); CHKERRQ(ierr);
      ierr = enthalpy.getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
      for (PetscInt k=0; k<grid.Mz; ++k) {
        const PetscScalar depth = vH(i,j) - grid.zlevels[k]; // FIXME task #7297
        ierr = EC->getWaterFraction(Enthij[k],EC->getPressureFromDepth(depth),
                                   omegaij[k]); CHKERRQ(ierr);
      }
    }
  }
  ierr = enthalpy.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  return 0;
}

//! Compute the CTS field, CTS = E/E_s(p), from Enth3, and put in a global IceModelVec3 provided by user.
/*!
The actual cold-temperate transition surface (CTS) is the level set CTS = 1.

Does not communicate ghosts for IceModelVec3 useForCTS.
 */
PetscErrorCode IceModel::setCTSFromEnthalpy(IceModelVec3 &useForCTS) {
  PetscErrorCode ierr;

  ierr = useForCTS.set_name("cts"); CHKERRQ(ierr);
  ierr = useForCTS.set_attrs(
     "diagnostic",
     "cts = E/E_s(p), so cold-temperate transition surface is at cts = 1",
     "", ""); CHKERRQ(ierr);

  PetscScalar *CTSij, *Enthij; // columns of these values
  ierr = useForCTS.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = useForCTS.getInternalColumn(i,j,&CTSij); CHKERRQ(ierr);
      ierr = Enth3.getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
      for (PetscInt k=0; k<grid.Mz; ++k) {
        const PetscScalar depth = vH(i,j) - grid.zlevels[k]; // FIXME task #7297
        CTSij[k] = EC->getCTS(Enthij[k], EC->getPressureFromDepth(depth));
      }
    }
  }
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = useForCTS.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  return 0;
}


//! Compute the CTS value of enthalpy in an ice column, and the lambda for BOMBPROOF.
/*!
Return argument Enth_s[Mz] has the enthalpy value for the pressure-melting 
temperature at the corresponding z level.
 */
PetscErrorCode IceModel::getEnthalpyCTSColumn(PetscScalar p_air,
					      PetscScalar thk,
					      PetscInt ks,
					      const PetscScalar *Enth,
					      const PetscScalar *w,
					      PetscScalar *lambda,
					      PetscScalar **Enth_s) {

  *lambda = 1.0;  // start with centered implicit for more accuracy
  const PetscScalar
      ice_rho_c = ice->rho * ice->c_p,
      ice_k     = ice->k;
  for (PetscInt k = 0; k <= ks; k++) {
    (*Enth_s)[k] = EC->getEnthalpyCTS(EC->getPressureFromDepth(thk - grid.zlevels_fine[k]));

    if (Enth[k] > (*Enth_s)[k]) { // lambda = 0 if temperate ice present in column
      *lambda = 0.0;
    } else {
      const PetscScalar 
          denom = (PetscAbs(w[k]) + 0.000001/secpera) * ice_rho_c * grid.dz_fine;
      *lambda = PetscMin(*lambda, 2.0 * ice_k / denom);
    }
  }

  for (PetscInt k = ks+1; k < grid.Mz_fine; k++) {
    (*Enth_s)[k] = EC->getEnthalpyCTS(p_air);
  }

  return 0;
}


/******  next 3 are helper functions for enthalpyDrainageStep() ******/

PetscErrorCode reportColumnSolveError(
    const PetscErrorCode solve_ierr, columnSystemCtx &sys, 
    const char *prefix, const PetscInt i, const PetscInt j) {

  char fname[PETSC_MAX_PATH_LEN];
  snprintf(fname, PETSC_MAX_PATH_LEN, "%s_i%d_j%d_zeropivot%d.m",
           prefix,i,j,solve_ierr);
  PetscErrorCode ierr;
  ierr = PetscPrintf(PETSC_COMM_SELF,
    "\n\ntridiagonal solve in enthalpyAndDrainageStep(), for %sSystemCtx,\n"
        "   failed at (%d,%d) with zero pivot position %d\n"
        "   viewing system to file %s ... \n",
        prefix, i, j, solve_ierr, fname); CHKERRQ(ierr);
  PetscViewer viewer;
  ierr = PetscViewerCreate(PETSC_COMM_SELF, &viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer, PETSC_VIEWER_ASCII);CHKERRQ(ierr);
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, fname);CHKERRQ(ierr);
  ierr = sys.viewSystem(viewer,"system"); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode reportColumn(
    MPI_Comm com, columnSystemCtx &sys,
    const char *prefix, const PetscInt i, const PetscInt j,
    PetscScalar *x, PetscInt n) {

  char fname[PETSC_MAX_PATH_LEN];
  snprintf(fname, PETSC_MAX_PATH_LEN, "%s_i%d_j%d.m", prefix,i,j);
  PetscErrorCode ierr;
  ierr = PetscPrintf(com,
    "\n\nviewing %s system and solution at (i,j)=(%d,%d):\n"
        "   viewing system to file %s ... \n\n\n",
        prefix, i, j, fname); CHKERRQ(ierr);
  PetscViewer viewer;
  ierr = PetscViewerCreate(com, &viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer, PETSC_VIEWER_ASCII);CHKERRQ(ierr);
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, fname);CHKERRQ(ierr);

  ierr = PetscViewerASCIIPrintf(viewer,
        "   1-norm = %.3e  and  diagonal-dominance ratio = %.5f\n",
        sys.norm1(n), sys.ddratio(n)); CHKERRQ(ierr);
  ierr = sys.viewSystem(viewer,"system"); CHKERRQ(ierr);
  ierr = sys.viewColumnValues(viewer, x, n, "solution x"); CHKERRQ(ierr);

  ierr = PetscViewerDestroy(viewer); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode copyColumn(PetscScalar *src, PetscScalar *dest, const PetscInt n) {
  for (PetscInt k = 0; k < n; k++) {
    dest[k] = src[k];
  }
  return 0;
}


// do this define to turn on messages
#define DEBUG_SHOW_BMELT 0


//! Update ice enthalpy field based on conservation of energy.
/*!
This method is documented by the page \ref bombproofenth and by [\ref
AschwandenBuelerBlatter].


This method updates IceModelVec3 vWork3d, IceModelVec2S vbmr, and 
IceModelVec2S vHmelt.  No communication of ghosts is done for any of these fields.
We use an instance of enthSystemCtx.

Regarding drainage, see [\ref AschwandenBuelerBlatter] and references therein.
 */

PetscErrorCode IceModel::enthalpyAndDrainageStep(
                      PetscScalar* vertSacrCount, PetscScalar* liquifiedVol,
                      PetscScalar* bulgeCount) {
  PetscErrorCode  ierr;

  if (config.get_flag("do_cold_ice_methods")) {
    SETERRQ(1,
      "PISM ERROR:  enthalpyAndDrainageStep() called but do_cold_ice_methods==true\n");
  }

  // get fine grid levels in ice
  PetscInt    fMz = grid.Mz_fine;  
  PetscScalar fdz = grid.dz_fine;
  vector<double> &fzlev = grid.zlevels_fine;

  const PetscScalar
    p_air     = config.get("surface_pressure"),
    ice_K     = ice->k / ice->rho, // enthalpy-conductivity for cold ice
    L         = config.get("water_latent_heat_fusion"),  // J kg-1
    omega_max = config.get("liquid_water_fraction_max"), // pure
    bulgeEnthMax  = config.get("enthalpy_cold_bulge_max"), // J kg-1
    hmelt_max = config.get("hmelt_max");                 // m

  IceModelVec2S *Rb;            // basal frictional heating
  ierr = stress_balance->get_basal_frictional_heating(Rb); CHKERRQ(ierr);

  IceModelVec3 *u3, *v3, *w3, *Sigma3;
  ierr = stress_balance->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr);
  ierr = stress_balance->get_volumetric_strain_heating(Sigma3); CHKERRQ(ierr); 

  PetscScalar *Enthnew;
  Enthnew = new PetscScalar[fMz];  // new enthalpy in column

  enthSystemCtx esys(config, Enth3, fMz);
  ierr = esys.initAllColumns(grid.dx, grid.dy, (dt_years_TempAge * secpera), fdz); CHKERRQ(ierr);

  bool viewOneColumn;
  ierr = PISMOptionsIsSet("-view_sys", viewOneColumn); CHKERRQ(ierr);

  if (getVerbosityLevel() >= 4) {  // view: all column-independent constants correct?
    ierr = EC->viewConstants(NULL); CHKERRQ(ierr);
    ierr = esys.viewConstants(NULL, false); CHKERRQ(ierr);
  }

  // now get map-plane coupler fields: Dirichlet upper surface boundary and
  //    mass balance lower boundary under shelves
  if (surface != PETSC_NULL) {
    ierr = surface->ice_surface_temperature(t_years_TempAge, dt_years_TempAge, artm);
    ierr = surface->ice_surface_liquid_water_fraction(t_years_TempAge, dt_years_TempAge,
                                                      liqfrac_surface); CHKERRQ(ierr);
    CHKERRQ(ierr);
  } else {
    SETERRQ(4,"PISM ERROR: surface == PETSC_NULL");
  }
  if (ocean != PETSC_NULL) {
    ierr = ocean->shelf_base_mass_flux(t_years_TempAge, dt_years_TempAge, shelfbmassflux);
        CHKERRQ(ierr);
    ierr = ocean->shelf_base_temperature(t_years_TempAge, dt_years_TempAge, shelfbtemp);
        CHKERRQ(ierr);
  } else {
    SETERRQ(5,"PISM ERROR: ocean == PETSC_NULL");
  }

  // FIXME: is it inefficient to create an IceModelVec2S here?
  IceModelVec2S G0;
  ierr = G0.create(grid, "bheatflx0", false); CHKERRQ(ierr);
  ierr = G0.set_attrs("internal",
                      "upward geothermal flux at z=0", 
                      "W m-2", ""); CHKERRQ(ierr);
  ierr = G0.set_glaciological_units("mW m-2");
  if (btu) {
    ierr = btu->get_upward_geothermal_flux(G0); CHKERRQ(ierr);
  } else {
    SETERRQ(3,"PISM ERROR: PISMBedThermalUnit* btu == PETSC_NULL in enthalpyAndDrainageStep()");
  }

  ierr = artm.begin_access(); CHKERRQ(ierr);
  ierr = shelfbmassflux.begin_access(); CHKERRQ(ierr);
  ierr = shelfbtemp.begin_access(); CHKERRQ(ierr);

  // get other map-plane fields
  ierr = liqfrac_surface.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vHmelt.begin_access(); CHKERRQ(ierr);
  ierr = vbmr.begin_access(); CHKERRQ(ierr);
  ierr = Rb->begin_access(); CHKERRQ(ierr);
  ierr = G0.begin_access(); CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);

  // these are accessed a column at a time
  ierr = u3->begin_access(); CHKERRQ(ierr);
  ierr = v3->begin_access(); CHKERRQ(ierr);
  ierr = w3->begin_access(); CHKERRQ(ierr);
  ierr = Sigma3->begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = vWork3d.begin_access(); CHKERRQ(ierr);

  PetscInt liquifiedCount = 0;

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      // for fine grid; this should *not* be replaced by call to grid.kBelowHeight()
      const PetscInt ks = static_cast<PetscInt>(floor(vH(i,j)/fdz));

#ifdef PISM_DEBUG
      // check if ks is valid
      if ((ks < 0) || (ks >= grid.Mz_fine)) {
        PetscPrintf(grid.com,
                    "ERROR: ks = %d computed at i = %d, j = %d is invalid,"
                    " possibly because of invalid ice thickness.\n",
                    ks, i, j);
        SETERRQ(1, "invalid ks");
      }
#endif

      const bool ice_free_column = (ks == 0),
                 is_floating     = vMask.is_floating(i,j);

      // enthalpy and pressures at top of ice
      const PetscScalar p_ks = EC->getPressureFromDepth(vH(i,j) - fzlev[ks]); // FIXME task #7297
      PetscScalar Enth_ks;
      ierr = EC->getEnthPermissive(artm(i,j), liqfrac_surface(i,j), p_ks, Enth_ks); CHKERRQ(ierr);

      // deal completely with columns with no ice; enthalpy, vHmelt, vbmr all need setting
      if (ice_free_column) {
        ierr = vWork3d.setColumn(i,j,Enth_ks); CHKERRQ(ierr);
        if ((vH(i,j) > 0.1) &&  (is_floating)) {
          // ks == 0, but we have ice (though not enough to do an ice column solve)
          // if floating then assume-maximally saturated till to avoid "shock" if grounding line advances
          vHmelt(i,j) = hmelt_max;
          vbmr(i,j) = shelfbmassflux(i,j);
        } else {
          // either truely no ice or grounded or both; either way zero-out subglacial fields
          vHmelt(i,j) = 0.0;  // no stored water on ice free land
          vbmr(i,j) = 0.0;    // no basal melt rate; melting is a surface process
                              //   on ice free land
        }

        goto donewithcolumn;
      } // end of if (ice_free_column)

      { // explicit scoping to deal with goto and initializers

        // ignore advection and strain heating in ice if isMarginal
        const bool isMarginal = checkThinNeigh(
                                 vH(i+1,j),vH(i+1,j+1),vH(i,j+1),vH(i-1,j+1),
                                 vH(i-1,j),vH(i-1,j-1),vH(i,j-1),vH(i+1,j-1)  );

        ierr = Enth3.getValColumn(i,j,ks,esys.Enth); CHKERRQ(ierr);
        ierr = w3->getValColumn(i,j,ks,esys.w); CHKERRQ(ierr);

        PetscScalar lambda;
        ierr = getEnthalpyCTSColumn(p_air, vH(i,j), ks, esys.Enth, esys.w, // FIXME task #7297
                                    &lambda, &esys.Enth_s); CHKERRQ(ierr);
        if (lambda < 1.0)  *vertSacrCount += 1; // count columns with lambda < 1

        if (vHmelt(i,j) > 0.0) { // if there is basal water, force base enthalpy up to pressure-melting
          esys.Enth[0] = esys.Enth_s[0];
        }

        bool base_is_cold = (esys.Enth[0] < esys.Enth_s[0]);
        const PetscScalar p1 = EC->getPressureFromDepth(vH(i,j) - fdz); // FIXME task #7297
        const bool k1_istemperate = EC->isTemperate(esys.Enth[1], p1); // level  z = + \Delta z

        // can now determine melt, but only preliminarily because of drainage,
        //   from heat flux out of bedrock, heat flux into ice, and frictional heating
        if (is_floating) {
          vbmr(i,j) = shelfbmassflux(i,j);
        } else {
          if (base_is_cold) {
              vbmr(i,j) = 0.0;  // zero melt rate if cold base
          } else {
            PetscScalar hf_up;
            if (k1_istemperate) {
              const PetscScalar pbasal = EC->getPressureFromDepth(vH(i,j)); // FIXME task #7297
              hf_up = - ice->k * (EC->getMeltingTemp(p1) - EC->getMeltingTemp(pbasal)) / fdz;
            } else {
              hf_up = - ice_K * (esys.Enth[1] - esys.Enth[0]) / fdz;
            }
            vbmr(i,j) = ( (*Rb)(i,j) + G0(i,j) - hf_up ) / (ice->rho * L); // = - Mb / rho in paper
#if DEBUG_SHOW_BMELT == 1
            verbPrintf(3,grid.com,
               "\n [stage 1; i,j=%d,%d has warm base, is grounded;\n"
               "             k1_istemperate=%d, G0=%.4f, hf_up=%.4f, Rb=%.4f, vbmr=%.6f(m/a)]",
               i,j,k1_istemperate,G0(i,j),hf_up,(*Rb)(i,j),vbmr(i,j)*secpera);
#endif
          }
        }

        // now set-up for solve in ice; note esys.Enth[], esys.w[],
        //   esys.Enth_s[] are already filled
        ierr = esys.setIndicesAndClearThisColumn(i,j,ks); CHKERRQ(ierr);

        ierr = u3->getValColumn(i,j,ks,esys.u); CHKERRQ(ierr);
        ierr = v3->getValColumn(i,j,ks,esys.v); CHKERRQ(ierr);
        ierr = Sigma3->getValColumn(i,j,ks,esys.Sigma); CHKERRQ(ierr);

        ierr = esys.setSchemeParamsThisColumn(isMarginal, lambda); CHKERRQ(ierr);
        ierr = esys.setBoundaryValuesThisColumn(Enth_ks); CHKERRQ(ierr);

        // ***** determine lowest-level equation at bottom of ice
        //       see figure in efgis paper, and page documenting BOMBPROOF
        if (is_floating) {
          // floating base: Dirichlet application of known temperature from ocean coupler
          // FIXME: this is assuming base of ice shelf has zero liquid fraction
          //        and is at cryostatic pressure; is this right?
          PetscScalar Enth0;
          ierr = EC->getEnthPermissive(shelfbtemp(i,j), 0.0, EC->getPressureFromDepth(vH(i,j)),
                                       Enth0); CHKERRQ(ierr);
          ierr = esys.setLevel0EqnThisColumn(1.0,0.0,Enth0); CHKERRQ(ierr);
        } else if (base_is_cold) {
          // cold, grounded base case:  Neumann q . n = q_lith . n + F_b   and   q = - K_i \nabla H
          ierr = esys.setLevel0EqnThisColumn(
                   1.0,-1.0,(fdz / ice_K) * (G0(i,j) + (*Rb)(i,j))); CHKERRQ(ierr);
        } else {
          // warm, grounded base case
          if (k1_istemperate) {
            // positive thickness of temperate ice:  Neumann q . n = 0 and q = - K_0 \nabla H
            //   so H(k=1)-H(k=0) = 0
            ierr = esys.setLevel0EqnThisColumn(1.0,-1.0,0.0); CHKERRQ(ierr);
          } else {
            // no thickness of temperate ice:  Dirichlet  H = H_s(pbasal)
            ierr = esys.setLevel0EqnThisColumn(1.0,0.0,esys.Enth_s[0]); CHKERRQ(ierr);
          }
        }

        ierr = esys.solveThisColumn(&Enthnew);
        if (ierr) reportColumnSolveError(ierr, esys, "enth", i, j);
        CHKERRQ(ierr);
        if (viewOneColumn && issounding(i,j)) {
          ierr = reportColumn(grid.com, esys, "enth", 
                              i, j, Enthnew, fMz); CHKERRQ(ierr);
        }

      // thermodynamic basal melt rate causes water to be added to layer
      PetscScalar Hmeltnew = vHmelt(i,j);
      if (vMask.is_grounded(i,j)) {
        Hmeltnew += vbmr(i,j) * (dt_years_TempAge * secpera);
      }

      // FIXME:  make drainage match paper
      // drain ice segments; has result that Enthnew[] is ice with at most
      //   omega_max liquid
      PetscScalar Hdrainedtotal = 0.0;
      for (PetscInt k=0; k < ks; k++) {
        if (EC->isLiquified(Enthnew[k],EC->getPressureFromDepth(vH(i,j) - fzlev[k]))) { // FIXME task #7297
          liquifiedCount++;
        }
        // if there is liquid water already, thus temperate, consider whether there
        //   is enough to cause drainage;  FIXME: UNACCOUNTED ENERGY LOSS IF E>E_l
        const PetscScalar p     = EC->getPressureFromDepth(vH(i,j) - fzlev[k]); // FIXME task #7297
        PetscScalar omega;
        EC->getWaterFraction(Enthnew[k], p, omega);  // return code not checked;
                                                     // we ignore E>E_l situation here
        PetscScalar dHdrained;
        if (omega > omega_max) {
          // drain water:
          dHdrained = (omega - omega_max) * fdz;
          // update enthalpy because omega == omega_max now:
          ierr = EC->getEnthAtWaterFraction(omega_max, p, Enthnew[k]); CHKERRQ(ierr);
        } else {
          dHdrained = 0.0;
        }                                       
        Hdrainedtotal += dHdrained;  // always a positive contribution
      }

      // in grounded case, add to both basal melt rate and Hmelt; if floating,
      // Hdrainedtotal is discarded because ocean determines basal melt rate
      if (vMask.is_grounded(i,j)) {
        vbmr(i,j) += Hdrainedtotal / (dt_years_TempAge * secpera);
        Hmeltnew += Hdrainedtotal;
#if DEBUG_SHOW_BMELT == 1
        verbPrintf(3,grid.com,
               "\n [stage 2; i,j=%d,%d has vbmr=%.6f(m/a) from drainage]",
               i,j,vbmr(i,j)*secpera);
#endif
      }

      // Enthnew[] is finalized!:  apply bulge limiter and transfer column
      //   into vWork3d; communication will occur later
      const PetscReal lowerEnthLimit = Enth_ks - bulgeEnthMax;
      for (PetscInt k=0; k < ks; k++) {
        if (Enthnew[k] < lowerEnthLimit) {
          *bulgeCount += 1;      // count the columns which have very large cold 
          Enthnew[k] = lowerEnthLimit;  // advection bulge ... and then actually
                                        // limit how low the enthalpy
        }
      }
      ierr = vWork3d.setValColumnPL(i,j,Enthnew); CHKERRQ(ierr);

      // finalize Hmelt value
      if (is_floating) {
        // if floating assume maximally saturated till to avoid "shock" if grounding line advances
        // UNACCOUNTED MASS & ENERGY (LATENT) LOSS/GAIN (TO/FROM OCEAN)!!
        vHmelt(i,j) = hmelt_max;
      } else {
        // limit Hmelt to be in [0.0, hmelt_max]
        // UNACCOUNTED MASS & ENERGY (LATENT) LOSS (TO INFINITY AND BEYOND)!!
        vHmelt(i,j) = PetscMax(0.0, PetscMin(hmelt_max, Hmeltnew) );
      }

#if DEBUG_SHOW_BMELT == 1
      verbPrintf(3,grid.com,
               "\n [stage 3; i,j=%d,%d has vbmr=%.6f(m/a) from drainage]",
               i,j,vbmr(i,j)*secpera);
#endif

      } // end explicit scoping
      
      donewithcolumn: 
      { }  // odd thing: something needs to follow goto target to get compilation

    }
  }

  ierr = artm.end_access(); CHKERRQ(ierr);
  ierr = shelfbmassflux.end_access(); CHKERRQ(ierr);
  ierr = shelfbtemp.end_access(); CHKERRQ(ierr);

  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vHmelt.end_access(); CHKERRQ(ierr);
  ierr = Rb->end_access(); CHKERRQ(ierr);
  ierr = G0.end_access(); CHKERRQ(ierr);
  ierr = vbmr.end_access(); CHKERRQ(ierr);
  ierr = liqfrac_surface.end_access(); CHKERRQ(ierr);

  ierr = u3->end_access(); CHKERRQ(ierr);
  ierr = v3->end_access(); CHKERRQ(ierr);
  ierr = w3->end_access(); CHKERRQ(ierr);
  ierr = Sigma3->end_access(); CHKERRQ(ierr);
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = vWork3d.end_access(); CHKERRQ(ierr);

  delete [] Enthnew;

  *liquifiedVol = ((double) liquifiedCount) * fdz * grid.dx * grid.dy;
  return 0;
}

