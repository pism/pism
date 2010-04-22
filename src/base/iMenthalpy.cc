// Copyright (C) 2009-2010 Andreas Aschwanden and Ed Bueler and Constantine Khroulev
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
#include "iceModel.hh"
#include "iceModelVec.hh"
#include "enthalpyConverter.hh"
#include "bedrockOnlySystem.hh"
#include "iceenthOnlySystem.hh"
#include "combinedSystem.hh"


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
PetscErrorCode IceModel::setEnth3FromT3_ColdIce() {
  PetscErrorCode ierr;
  
  ierr = T3.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);

  PetscScalar *Tij, *Enthij; // columns of these values
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = T3.getInternalColumn(i,j,&Tij); CHKERRQ(ierr);
      ierr = Enth3.getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
      for (PetscInt k=0; k<grid.Mz; ++k) {
        const PetscScalar depth = vH(i,j) - grid.zlevels[k];
        ierr = EC->getEnthPermissive(Tij[k],0.0,EC->getPressureFromDepth(depth),
                                    Enthij[k]); CHKERRQ(ierr);
      }
    }
  }

  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = T3.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = Enth3.beginGhostComm(); CHKERRQ(ierr);
  ierr = Enth3.endGhostComm(); CHKERRQ(ierr);
  return 0;
}


//! Compute Enth3 from temperature T3 and liquid fraction.
/*!
Because Enth3 gets set, does ghost communication to finish.
 */
PetscErrorCode IceModel::setEnth3FromT3AndLiqfrac3(
                                          IceModelVec3 &Liqfrac3) {
  PetscErrorCode ierr;
  
  ierr = T3.begin_access(); CHKERRQ(ierr);
  ierr = Liqfrac3.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);

  PetscScalar *Tij, *Liqfracij, *Enthij; // columns of these values
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = T3.getInternalColumn(i,j,&Tij); CHKERRQ(ierr);
      ierr = Liqfrac3.getInternalColumn(i,j,&Liqfracij); CHKERRQ(ierr);
      ierr = Enth3.getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
      for (PetscInt k=0; k<grid.Mz; ++k) {
        const PetscScalar depth = vH(i,j) - grid.zlevels[k];
        ierr = EC->getEnthPermissive(Tij[k],Liqfracij[k],
                      EC->getPressureFromDepth(depth), Enthij[k]); CHKERRQ(ierr);
      }
    }
  }

  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = T3.end_access(); CHKERRQ(ierr);
  ierr = Liqfrac3.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = Enth3.beginGhostComm(); CHKERRQ(ierr);
  ierr = Enth3.endGhostComm(); CHKERRQ(ierr);
  return 0;
}


//! Compute the ice temperature corresponding to Enth3, and put in a given IceModelVec3.
/*!
  Computes ice temperature from enthalpy (as a diagnostic quantity).
 */
PetscErrorCode IceModel::compute_temp(IceModelVec3 &result) {
  PetscErrorCode ierr;

  PetscScalar *Tij, *Enthij; // columns of these values
  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = result.getInternalColumn(i,j,&Tij); CHKERRQ(ierr);
      ierr = Enth3.getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
      for (PetscInt k=0; k<grid.Mz; ++k) {
        const PetscScalar depth = vH(i,j) - grid.zlevels[k];
        ierr = EC->getAbsTemp(Enthij[k],EC->getPressureFromDepth(depth), Tij[k]); 
        if (ierr) {
          PetscPrintf(grid.com,
            "\n\nEnthalpyConverter.getAbsTemp() error at i=%d,j=%d,k=%d\n\n",
            i,j,k);
        }
        CHKERRQ(ierr);
      }
    }
  }
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  ierr = result.set_name("temp"); CHKERRQ(ierr);
  ierr = result.set_attrs("diagnostic","ice temperature",
			  "K", "land_ice_temperature"); CHKERRQ(ierr);
  ierr = result.set_attr("valid_min", 0.0); CHKERRQ(ierr);

  return 0;
}


//! Compute the liquid fraction corresponding to Enth3, and put in a global IceModelVec3 provided by user.
/*!
Does not communicate ghosts for IceModelVec3 useForLiquidFrac.
 */
PetscErrorCode IceModel::setLiquidFracFromEnthalpy(IceModelVec3 &useForLiquidFrac) {
  PetscErrorCode ierr;

  ierr = useForLiquidFrac.set_name("liqfrac"); CHKERRQ(ierr);
  ierr = useForLiquidFrac.set_attrs(
     "diagnostic",
     "liquid water fraction in ice (between 0 and 1)",
     "", ""); CHKERRQ(ierr);

  PetscScalar *omegaij, *Enthij; // columns of these values
  ierr = useForLiquidFrac.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = useForLiquidFrac.getInternalColumn(i,j,&omegaij); CHKERRQ(ierr);
      ierr = Enth3.getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
      for (PetscInt k=0; k<grid.Mz; ++k) {
        const PetscScalar depth = vH(i,j) - grid.zlevels[k];
        ierr = EC->getWaterFraction(Enthij[k],EC->getPressureFromDepth(depth),
                                   omegaij[k]); CHKERRQ(ierr);
      }
    }
  }
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = useForLiquidFrac.end_access(); CHKERRQ(ierr);
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
        const PetscScalar depth = vH(i,j) - grid.zlevels[k];
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


//! Update enthalpy field based on conservation of energy in ice and bedrock.
/*!
This method is documented by the page \ref bombproofenth.

This method uses instances of combinedSystemCtx, bedrockOnlySystemCtx, and
iceenthOnlySystemCtx.

This method modifies IceModelVec3 Enthnew3, IceModelVec3Bedrock Tb3,
IceModelVec2S vBasalMeltRate, and IceModelVec2S vHmelt.  No communication of
ghosts is done for any of these fields.

Regarding drainage, we move the liquid water fraction which is in excess of
the fixed fraction \c omega_max in a column segment [z,z+dz] to the base.
Heuristic: Once liquid water fraction exceeds a cap, all of it goes to the base.
Follows [\ref Greve97Greenland] and references therein.
 */
PetscErrorCode IceModel::enthalpyAndDrainageStep(
                      PetscScalar* vertSacrCount, PetscScalar* liquifiedVol,
                      PetscScalar* bulgeCount) {
  PetscErrorCode  ierr;

  if (config.get_flag("do_cold_ice_methods")) {
    SETERRQ(1,
      "PISM ERROR:  enthalpyAndDrainageStep() called but do_cold_ice_methods==true\n");
  }

  // get fine grid levels in ice and bedrock; guaranteed dz=dzb
  PetscInt    fMz = grid.Mz_fine,
    fMbz = grid.Mbz_fine;  
  PetscScalar fdz = grid.dz_fine,
    *fzlev = grid.zlevels_fine;

  if (fMbz == 2) {
    SETERRQ(2,
      "PISM ERROR:  enthalpyAndDrainageStep() does not currently allow fMbz == 2;\n"
      "   fMbz==1 and fMbz>2 are allowed\n");
  }

  const PetscScalar
    p_air     = config.get("surface_pressure"),
    ice_rho   = ice->rho,
    ice_c     = ice->c_p,
    ice_k     = ice->k,
    L         = config.get("water_latent_heat_fusion"),  // J kg-1
    omega_max = config.get("liquid_water_fraction_max"), // pure
    warm_dE   = config.get("warm_base_flux_enthalpy_fraction") * L,
    refreeze_rate = config.get("cold_base_refreeze_rate"), // m s-1
    bulgeEnthMax  = config.get("enthalpy_cold_bulge_max"), // J kg-1
    hmelt_max = config.get("hmelt_max");                 // m

  PetscScalar *Enthnew, *Tbnew;
  Enthnew = new PetscScalar[fMz];  // new enthalpy in column
  Tbnew   = new PetscScalar[fMbz]; // new bedrock temperature in column

  combinedSystemCtx    cbsys(config, Enth3, fMz, fMbz);
  ierr = cbsys.initAllColumns(grid.dx, grid.dy, dtTempAge, fdz, fdz); CHKERRQ(ierr);
  // space for solution when ice and bedrock are combined in one system
  PetscScalar *xcombined;
  xcombined = new PetscScalar[fMbz + fMz - 1];

  bedrockOnlySystemCtx bosys(config, fMbz);
  ierr = bosys.initAllColumns(dtTempAge, fdz); CHKERRQ(ierr);

  iceenthOnlySystemCtx iosys(config, Enth3, fMz);
  ierr = iosys.initAllColumns(grid.dx, grid.dy, dtTempAge, fdz); CHKERRQ(ierr);

  bool viewOneColumn;
  ierr = PISMOptionsIsSet("-view_sys", viewOneColumn); CHKERRQ(ierr);

  // FIXME: verbosity failure?: option "-verbose 4" does not generate true here?
  if (getVerbosityLevel() >= 4) {  // view: all column-independent constants correct?
    ierr = EC->viewConstants(NULL); CHKERRQ(ierr);
    ierr = cbsys.viewConstants(NULL, false); CHKERRQ(ierr);
    ierr = bosys.viewConstants(NULL, false); CHKERRQ(ierr);
    ierr = iosys.viewConstants(NULL, false); CHKERRQ(ierr);
  }

  // now get map-plane coupler fields: Dirichlet upper surface boundary and
  //    mass balance lower boundary under shelves
  if (surface != PETSC_NULL) {
    ierr = surface->ice_surface_temperature(grid.year, dtTempAge / secpera, artm);
        CHKERRQ(ierr);
  } else {
    SETERRQ(4,"PISM ERROR: surface == PETSC_NULL");
  }
  if (ocean != PETSC_NULL) {
    ierr = ocean->shelf_base_mass_flux(grid.year, dtTempAge / secpera, shelfbmassflux);
        CHKERRQ(ierr);
    ierr = ocean->shelf_base_temperature(grid.year, dtTempAge / secpera, shelfbtemp);
        CHKERRQ(ierr);
  } else {
    SETERRQ(5,"PISM ERROR: ocean == PETSC_NULL");
  }
  ierr = artm.begin_access(); CHKERRQ(ierr);
  ierr = shelfbmassflux.begin_access(); CHKERRQ(ierr);
  ierr = shelfbtemp.begin_access(); CHKERRQ(ierr);

  // get other map-plane fields
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vHmelt.begin_access(); CHKERRQ(ierr);
  ierr = vbmr.begin_access(); CHKERRQ(ierr);
  ierr = vRb.begin_access(); CHKERRQ(ierr);
  ierr = vGhf.begin_access(); CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);

  // these are accessed a column at a time
  ierr = u3.begin_access(); CHKERRQ(ierr);
  ierr = v3.begin_access(); CHKERRQ(ierr);
  ierr = w3.begin_access(); CHKERRQ(ierr);
  ierr = Sigma3.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = Enthnew3.begin_access(); CHKERRQ(ierr);
  ierr = Tb3.begin_access(); CHKERRQ(ierr);

  PetscInt liquifiedCount = 0;

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      // for fine grid; this should *not* be replaced by call to grid.kBelowHeight()
      const PetscInt ks = static_cast<PetscInt>(floor(vH(i,j)/fdz));

      // enthalpy and pressures at top of ice
      const PetscScalar p_ks = EC->getPressureFromDepth(vH(i,j) - fzlev[ks]);
      PetscScalar Enth_ks;
      // in theory we could have a water fraction at k=ks level, but for
      //   now there is no case where we have that:
      ierr = EC->getEnthPermissive(artm(i,j), 0.0, p_ks,  Enth_ks); CHKERRQ(ierr);

      const bool is_floating = vMask.is_floating(i,j);

      // deal completely with columns with no ice; note bedrock does need actions
      if (ks == 0) {
        ierr = Enthnew3.setColumn(i,j,Enth_ks); CHKERRQ(ierr);
        PetscScalar Tbtop;
        if (is_floating) {
          Tbtop = shelfbtemp(i,j);
        } else {
          ierr = EC->getAbsTemp(Enth_ks, p_ks, Tbtop); CHKERRQ(ierr);
        }
        if (fMbz > 1) { // bedrock layer if present
          ierr = bosys.setIndicesAndClearThisColumn(i,j,-1); CHKERRQ(ierr);  
          ierr = Tb3.getValColumnPL(i,j,bosys.Tb); CHKERRQ(ierr);
          ierr = bosys.setBoundaryValuesThisColumn(Tbtop, vGhf(i,j)); CHKERRQ(ierr);
          ierr = bosys.solveThisColumn(&Tbnew);
          if (ierr)  reportColumnSolveError(ierr, bosys, "bedrockOnly", i, j);
          CHKERRQ(ierr);
          if (viewOneColumn && issounding(i,j)) {
            ierr = reportColumn(grid.com, bosys, "bedrockOnly",
                                i, j, Tbnew, fMbz); CHKERRQ(ierr);
          }
          ierr = Tb3.setValColumnPL(i,j,Tbnew); CHKERRQ(ierr);
        } else { // no bedrock layer: we are actually just setting one value
          ierr = Tb3.setColumn(i,j,Tbtop); CHKERRQ(ierr);
        }
        if (is_floating) {
          vHmelt(i,j) = hmelt_max;
          vbmr(i,j) = shelfbmassflux(i,j);
        } else {
          vHmelt(i,j) = 0.0;  // no stored water on ice free land
          vbmr(i,j) = 0.0;    // no basal melt rate; melting is a surface process
                              //   on ice free land
        }
       
        goto donewithcolumn;
      }

      { // explicit scoping to deal with goto and initializers

      // ignore advection and strain heating in ice if isMarginal
      const bool isMarginal = checkThinNeigh(
                                 vH(i+1,j),vH(i+1,j+1),vH(i,j+1),vH(i-1,j+1),
                                 vH(i-1,j),vH(i-1,j-1),vH(i,j-1),vH(i+1,j-1)  );
      const PetscScalar p_basal = EC->getPressureFromDepth(vH(i,j));

      ierr = Enth3.getValColumn(i,j,ks,iosys.Enth); CHKERRQ(ierr);
      ierr = w3.getValColumn(i,j,ks,iosys.w); CHKERRQ(ierr);

      PetscScalar lambda;
      ierr = getEnthalpyCTSColumn(p_air, vH(i,j), ks, iosys.Enth, iosys.w,
                                  &lambda, &iosys.Enth_s); CHKERRQ(ierr);

      if (lambda < 1.0)  *vertSacrCount += 1; // count columns with lambda < 1

      // major decision: is cold base and grounded and has bedrock layer?:
      if ( (iosys.Enth[0] < iosys.Enth_s[0]) && (fMbz > 1) && (!is_floating) ) {

        // ***** COLD BASE, GROUNDED CASE WITH BEDROCK *****
        ierr = cbsys.setIndicesAndClearThisColumn(i,j,ks); CHKERRQ(ierr);

        ierr = copyColumn(iosys.Enth,cbsys.Enth,fMz); CHKERRQ(ierr);
        ierr = copyColumn(iosys.Enth_s,cbsys.Enth_s,fMz); CHKERRQ(ierr);
        ierr = u3.getValColumn(i,j,ks,cbsys.u); CHKERRQ(ierr);
        ierr = v3.getValColumn(i,j,ks,cbsys.v); CHKERRQ(ierr);
        ierr = copyColumn(iosys.w,cbsys.w,fMz); CHKERRQ(ierr);
        ierr = Sigma3.getValColumn(i,j,ks,cbsys.Sigma); CHKERRQ(ierr);
        ierr = Tb3.getValColumnPL(i,j,cbsys.Tb); CHKERRQ(ierr);

        ierr = cbsys.setSchemeParamsThisColumn(isMarginal, lambda);
            CHKERRQ(ierr);
        ierr = cbsys.setBoundaryValuesThisColumn(Enth_ks, vGhf(i,j), vRb(i,j));
            CHKERRQ(ierr);

        ierr = cbsys.solveThisColumn(&xcombined);
        if (ierr) reportColumnSolveError(ierr, cbsys, "combined", i, j);
        CHKERRQ(ierr);
        if (viewOneColumn && issounding(i,j)) {
          ierr = reportColumn(grid.com, cbsys, "combined", i, j,
                              xcombined, fMbz+fMz-1); CHKERRQ(ierr);
        }
        // break result x[] of combined system between Enthnew[fMz]
        //   and Tbnew[fMbz]
        for (PetscInt k = 0; k < fMbz-1; k++) {
          Tbnew[k] = xcombined[k];
        }
        // at this point we need a temperature from ice that could in extreme
        //   situations *be fully melted*; thus we catch the return code and
        //   count this phenomenon
        ierr = EC->getAbsTemp(xcombined[fMbz-1], p_basal, Tbnew[fMbz-1]);
        if (ierr==1) { // return code of 1 means block of ice melted completely
          liquifiedCount++;
        } else CHKERRQ(ierr);
        for (PetscInt k = 0; k < fMz; k++) {
          Enthnew[k] = xcombined[k + fMbz-1];
        }

        vbmr(i,j) = 0.0;  // zero melt rate if cold base

      } else {

        // ***** ALL OTHER CASES *****

        // ***** BEDROCK ONLY SOLVE *****
        PetscScalar hf_base;
        if (fMbz > 1) { // deal with bedrock layer first, if present
          // case of temperate bed and a bedrock layer
          ierr = bosys.setIndicesAndClearThisColumn(i,j,-1); CHKERRQ(ierr);  

          ierr = Tb3.getValColumnPL(i,j,bosys.Tb); CHKERRQ(ierr);

          const PetscScalar Tbtop = (is_floating ? shelfbtemp(i,j)
				     : EC->getMeltingTemp(p_basal));
          ierr = bosys.setBoundaryValuesThisColumn(Tbtop, vGhf(i,j));
              CHKERRQ(ierr);

          ierr = bosys.solveThisColumn(&Tbnew);
          if (ierr) reportColumnSolveError(ierr, bosys, "bedrockOnly", i, j);
          CHKERRQ(ierr);
          if (viewOneColumn && issounding(i,j)) {
            ierr = reportColumn(grid.com, bosys, "bedrockOnly",
                                i, j, Tbnew, fMbz); CHKERRQ(ierr);
          }

          hf_base = bosys.extractHeatFluxFromSoln(Tbnew);
        } else {
          hf_base = vGhf(i,j);
        }

        // can determine melt now from heat flux out of base, etc.
        if (is_floating) {
          vbmr(i,j) = shelfbmassflux(i,j);
        } else {
          if (iosys.Enth[0] < iosys.Enth_s[0]) {
            // this case occurs only if no bedrock thermal layer
            vbmr(i,j) = 0.0;  // zero melt rate if cold base
          } else {
            vbmr(i,j) = ( hf_base + vRb(i,j) ) / (ice_rho * L);
          }
        }

        // ***** ICE ONLY SOLVE *****
        // now set-up for solve in ice; note iosys.Enth[], iosys.w[],
        //   iosys.Enth_s[] are already filled
        ierr = iosys.setIndicesAndClearThisColumn(i,j,ks); CHKERRQ(ierr);

        ierr = u3.getValColumn(i,j,ks,iosys.u); CHKERRQ(ierr);
        ierr = v3.getValColumn(i,j,ks,iosys.v); CHKERRQ(ierr);
        ierr = Sigma3.getValColumn(i,j,ks,iosys.Sigma); CHKERRQ(ierr);

        ierr = iosys.setSchemeParamsThisColumn(isMarginal, lambda); CHKERRQ(ierr);
        ierr = iosys.setBoundaryValuesThisColumn(Enth_ks); CHKERRQ(ierr);

        // ***** determine lowest-level equation at bottom of ice
        //       see page documenting BOMBPROOF
        if (iosys.Enth[0] < iosys.Enth_s[0]) {
          // cold base case with fMbz==1: ice base equation says heat flux is known
          // this case only if no bedrock thermal layer
          const PetscScalar C = ice_c * fdz / ice_k;
          ierr = iosys.setLevel0EqnThisColumn(
                   1.0,-1.0,C * (hf_base + vRb(i,j))); CHKERRQ(ierr);
        } else {
          // we are in the warm base case, so velocity at bottom of ice in the
          //   last time step determines type of boundary condition, either
          //   (i) if w(0)<0 then outflow b.c. or (ii) if w(0)>=0 then Dirichlet
          // *but*
          // for basal ice only slightly above the pressure-melting temperature,
          //   we combine the boundary condition (either (i) or (ii)) with
          //   an amount of heat flux into the base; alpha is the amount of that flux
          PetscScalar a0, a1, rhs;
          if (iosys.w[0] < 0.0) {
            // outflow "boundary condition": apply diffusion-free, upwinded form
            //   of enthalpy equation (bombtwo)
            rhs  = iosys.Enth[0];
            if (!isMarginal) {
              planeStar ss;
              Enth3.getPlaneStar(i,j,0,&ss);
              const PetscScalar
                 UpEnthu = (iosys.u[0] < 0) ? iosys.u[0] * (ss.ip1 -  ss.ij) / grid.dx
                                            : iosys.u[0] * (ss.ij  - ss.im1) / grid.dx,
                 UpEnthv = (iosys.v[0] < 0) ? iosys.v[0] * (ss.jp1 -  ss.ij) / grid.dy
                                            : iosys.v[0] * (ss.ij  - ss.jm1) / grid.dy;
              rhs += dtTempAge * ((iosys.Sigma[0] / ice_rho) - UpEnthu - UpEnthv);
            }
            const PetscScalar nuw0 = (dtTempAge / fdz) * iosys.w[0];
            a0 = 1 - nuw0;
            a1 = nuw0;
          } else {
            // Dirichlet cond. for enthalpy at ice base
            rhs = iosys.Enth_s[0];
            a0  = 1.0;
            a1  = 0.0;
          }
          const PetscScalar
            alpha      = (iosys.Enth[0] < iosys.Enth_s[0] + warm_dE)
                            ? 1.0 - ((iosys.Enth[0] - iosys.Enth_s[0]) / warm_dE)
                            : 0.0;
          const PetscScalar C = ice_c * fdz / ice_k;
          rhs = (1.0 - alpha) * rhs + alpha * ( C * (hf_base + vRb(i,j)) );
          a0  = (1.0 - alpha) * a0  + alpha * 1.0,
          a1  = (1.0 - alpha) * a1  + alpha * (-1.0);

          if (!is_floating)    vbmr(i,j) *= 1.0 - alpha;

          ierr = iosys.setLevel0EqnThisColumn(a0,a1,rhs); CHKERRQ(ierr);
        }

        ierr = iosys.solveThisColumn(&Enthnew);
        if (ierr) reportColumnSolveError(ierr, iosys, "iceenthOnly", i, j);
        CHKERRQ(ierr);
        if (viewOneColumn && issounding(i,j)) {
          ierr = reportColumn(grid.com, iosys, "iceenthOnly", 
                              i, j, Enthnew, fMz); CHKERRQ(ierr);
        }
      }

      // thermodynamic basal melt rate causes water to be added to layer
      PetscScalar Hmeltnew = vHmelt(i,j);
      if (!is_floating) {
        Hmeltnew += vbmr(i,j) * dtTempAge;
      }

      // drain ice segments; has result that Enthnew[] is ice with at most
      //   omega_max liquid
      PetscScalar Hdrainedtotal = 0.0;
      for (PetscInt k=0; k < ks; k++) {
        if (EC->isLiquified(Enthnew[k],EC->getPressureFromDepth(vH(i,j) - fzlev[k]))) {
          liquifiedCount++;
        }

        // if there is liquid water already, thus temperate, consider whether there
        //   is enough to cause drainage;  UNACCOUNTED ENERGY LOSS IF E>E_l
        const PetscScalar p     = EC->getPressureFromDepth(vH(i,j) - fzlev[k]),
                          omega = EC->getWaterFractionLimited(Enthnew[k], p);
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
      if (!is_floating) {
        vbmr(i,j) += Hdrainedtotal / dtTempAge;
        Hmeltnew += Hdrainedtotal;
      }

      // Enthnew[] is finalized!:  apply bulge limiter and transfer column
      //   into Enthnew3; communication will occur later
      const PetscReal lowerEnthLimit = Enth_ks - bulgeEnthMax;
      for (PetscInt k=0; k < ks; k++) {
        if (Enthnew[k] < lowerEnthLimit) {
          *bulgeCount += 1;      // count the columns which have very large cold 
          Enthnew[k] = lowerEnthLimit;  // advection bulge ... and then actually
                                        // limit how low the enthalpy
        }
      }
      ierr = Enthnew3.setValColumnPL(i,j,Enthnew); CHKERRQ(ierr);

      // if no thermal layer then need to fill Tbnew[0] directly
      if (fMbz == 1) {
        if (is_floating) { // floating: get from PISMOceanModel
          Tbnew[0] = shelfbtemp(i,j);
        } else {                      // grounded: duplicate temp from ice
          ierr = EC->getAbsTemp(Enthnew[0],p_basal,Tbnew[0]); CHKERRQ(ierr);
        }
      }

      // Tbnew[] is finalized!:  transfer column into Tb3; no need for
      //    communication, even later
      ierr = Tb3.setValColumnPL(i,j,Tbnew); CHKERRQ(ierr);

      // finalize Hmelt value
      if (updateHmelt == PETSC_TRUE) {
        if (is_floating) {
          // if floating assume maximally saturated "till"
          // UNACCOUNTED MASS & ENERGY (LATENT) LOSS/GAIN (TO/FROM OCEAN)!!
          vHmelt(i,j) = hmelt_max;
        } else if (ks == 0) {
          vHmelt(i,j) = 0.0;  // no stored water on ice free land
        } else {
          // limit Hmelt to be in [0.0, hmelt_max]
          // UNACCOUNTED MASS & ENERGY (LATENT) LOSS (TO INFINITY AND BEYOND)!!
          vHmelt(i,j) = PetscMax(0.0, PetscMin(hmelt_max, Hmeltnew) );

          // refreeze case: if grounded base has become cold then put back ice at
          //   externally-set maximum rate; basal enthalpy not altered
          if ( (Enthnew[0] < iosys.Enth_s[0]) && (vHmelt(i,j) > 0.0) ) {
            if (vHmelt(i,j) > refreeze_rate * dtTempAge) {
              vbmr(i,j) -= refreeze_rate;
              vHmelt(i,j) -= refreeze_rate * dtTempAge;
            } else {
              // in this case we refreeze all available Hmelt
              vbmr(i,j) -= vHmelt(i,j) / dtTempAge;
              vHmelt(i,j) = 0.0;
            }
          }

        }
      }

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
  ierr = vRb.end_access(); CHKERRQ(ierr);
  ierr = vGhf.end_access(); CHKERRQ(ierr);
  ierr = vbmr.end_access(); CHKERRQ(ierr);

  ierr = Tb3.end_access(); CHKERRQ(ierr);
  ierr = u3.end_access(); CHKERRQ(ierr);
  ierr = v3.end_access(); CHKERRQ(ierr);
  ierr = w3.end_access(); CHKERRQ(ierr);
  ierr = Sigma3.end_access(); CHKERRQ(ierr);
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = Enthnew3.end_access(); CHKERRQ(ierr);

  delete [] Enthnew; delete [] Tbnew;  delete [] xcombined;

  *liquifiedVol = ((double) liquifiedCount) * fdz * grid.dx * grid.dy;
  return 0;
}


