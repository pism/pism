// Copyright (C) 2009-2012 Andreas Aschwanden and Ed Bueler and Constantine Khroulev
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
#include "varenthSystem.hh"
#include "DrainageCalculator.hh"
#include "Mask.hh"
#include "PISMStressBalance.hh"
#include "PISMSurface.hh"
#include "PISMOcean.hh"
#include "bedrockThermalUnit.hh"
#include "enthalpyConverter.hh"
#include "pism_options.hh"

//! \file iMenthalpy.cc Methods of IceModel which implement the enthalpy formulation of conservation of energy.


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
        const PetscScalar depth = vH(i,j) - grid.zlevels[k]; // FIXME issue #15
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
        const PetscScalar depth = vH(i,j) - grid.zlevels[k]; // FIXME issue #15
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
        const PetscScalar depth = vH(i,j) - grid.zlevels[k]; // FIXME issue #15
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
        const PetscScalar depth = vH(i,j) - grid.zlevels[k]; // FIXME issue #15
        CTSij[k] = EC->getCTS(Enthij[k], EC->getPressureFromDepth(depth));
      }
    }
  }
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = useForCTS.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  return 0;
}


//! Compute the CTS value of enthalpy in an ice column.
/*!
Return argument Enth_s[Mz] has the enthalpy value for the pressure-melting 
temperature at the corresponding z level.
 */
PetscErrorCode IceModel::getEnthalpyCTSColumn(PetscScalar p_air,
					      PetscScalar thk,
					      PetscInt ks,
					      PetscScalar **Enth_s) {

  for (PetscInt k = 0; k <= ks; k++) {
    const PetscScalar p = EC->getPressureFromDepth(thk - grid.zlevels_fine[k]); // FIXME issue #15
    (*Enth_s)[k] = EC->getEnthalpyCTS(p);
  }
  const PetscScalar Es_air = EC->getEnthalpyCTS(p_air);
  for (PetscInt k = ks+1; k < grid.Mz_fine; k++) {
    (*Enth_s)[k] = Es_air;
  }
  return 0;
}


//! Compute the lambda for BOMBPROOF.
/*!
See page \ref bombproofenth.
 */
PetscErrorCode IceModel::getlambdaColumn(PetscInt ks,
					 PetscScalar ice_rho_c,
                                         PetscScalar ice_k,
					 const PetscScalar *Enth,
					 const PetscScalar *Enth_s,
					 const PetscScalar *w,
					 PetscScalar *lambda) {

  *lambda = 1.0;  // start with centered implicit for more accuracy
  for (PetscInt k = 0; k <= ks; k++) {
    if (Enth[k] > Enth_s[k]) { // lambda = 0 if temperate ice present in column
      *lambda = 0.0;
    } else {
      const PetscScalar 
          denom = (PetscAbs(w[k]) + 0.000001/secpera) * ice_rho_c * grid.dz_fine;
      *lambda = PetscMin(*lambda, 2.0 * ice_k / denom);
    }
  }
  return 0;
}


//! Update ice enthalpy field based on conservation of energy.
/*!
This method is documented by the page \ref bombproofenth and by [\ref
AschwandenBuelerKhroulevBlatter].

This method updates IceModelVec3 vWork3d = vEnthnew, IceModelVec2S vbmr, and 
IceModelVec2S vbwat.  No communication of ghosts is done for any of these fields.

We use an instance of enthSystemCtx.

Regarding drainage, see [\ref AschwandenBuelerKhroulevBlatter] and references therein.
 */
PetscErrorCode IceModel::enthalpyAndDrainageStep(
                      PetscScalar* vertSacrCount, PetscScalar* liquifiedVol,
                      PetscScalar* bulgeCount) {
  PetscErrorCode  ierr;

  if (config.get_flag("do_cold_ice_methods")) {
    SETERRQ(grid.com, 1,
      "PISM ERROR:  enthalpyAndDrainageStep() called but do_cold_ice_methods==true\n");
  }

  const PetscReal dt_secs = dt_TempAge;

  // get fine grid levels in ice
  PetscInt    fMz = grid.Mz_fine;  
  PetscScalar fdz = grid.dz_fine;
  vector<double> &fzlev = grid.zlevels_fine;

  // essentially physical constants
  const PetscScalar
    p_air     = config.get("surface_pressure"),          // Pa
    ice_rho   = config.get("ice_density"),               // kg m-3
    L         = config.get("water_latent_heat_fusion");  // J kg-1

  // constants used in controlling numerical scheme
  const PetscScalar
    default_ice_k = config.get("ice_thermal_conductivity"),   // used in setting lambda
    default_ice_c = config.get("ice_specific_heat_capacity"), // used in setting lambda
    bulgeEnthMax  = config.get("enthalpy_cold_bulge_max");    // J kg-1

  // constants used in hydrology model
  const PetscScalar
    bwat_decay_rate = config.get("bwat_decay_rate"),   // m s-1
    bwat_max        = config.get("bwat_max");          // m

  DrainageCalculator dc(config);
  
  IceModelVec2S *Rb;
  IceModelVec3 *u3, *v3, *w3, *Sigma3;
  ierr = stress_balance->get_basal_frictional_heating(Rb); CHKERRQ(ierr);
  ierr = stress_balance->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr);
  ierr = stress_balance->get_volumetric_strain_heating(Sigma3); CHKERRQ(ierr); 

  PetscScalar *Enthnew;
  Enthnew = new PetscScalar[fMz];  // new enthalpy in column

  enthSystemCtx *esys;
  if (config.get_flag("use_temperature_dependent_thermal_conductivity") ||
      config.get_flag("use_linear_in_temperature_heat_capacity")) {
    esys  = new varenthSystemCtx(config, Enth3, fMz, "varenth", EC);
  } else {
    esys  = new enthSystemCtx(config, Enth3, fMz, "enth");
  }
  ierr = esys->initAllColumns(grid.dx, grid.dy, dt_secs, fdz); CHKERRQ(ierr);

  bool viewOneColumn;
  ierr = PISMOptionsIsSet("-view_sys", viewOneColumn); CHKERRQ(ierr);

  if (getVerbosityLevel() >= 4) {  // view: all column-independent constants correct?
    ierr = EC->viewConstants(NULL); CHKERRQ(ierr);
    ierr = esys->viewConstants(NULL, false); CHKERRQ(ierr);
  }

  // now get map-plane coupler fields: Dirichlet upper surface boundary and
  //    mass balance lower boundary under shelves
  if (surface != PETSC_NULL) {
    ierr = surface->ice_surface_temperature(artm);
    ierr = surface->ice_surface_liquid_water_fraction(liqfrac_surface); CHKERRQ(ierr);
    CHKERRQ(ierr);
  } else {
    SETERRQ(grid.com, 4,"PISM ERROR: surface == PETSC_NULL");
  }
  if (ocean != PETSC_NULL) {
    ierr = ocean->shelf_base_mass_flux(shelfbmassflux);
        CHKERRQ(ierr);
    ierr = ocean->shelf_base_temperature(shelfbtemp);
        CHKERRQ(ierr);
  } else {
    SETERRQ(grid.com, 5,"PISM ERROR: ocean == PETSC_NULL");
  }

  IceModelVec2S G0 = vWork2d[0];
  ierr = G0.set_attrs("internal","upward geothermal flux at z=0","W m-2", ""); CHKERRQ(ierr);
  ierr = G0.set_glaciological_units("mW m-2");
  if (btu) {
    ierr = btu->get_upward_geothermal_flux(G0); CHKERRQ(ierr);
  } else {
    SETERRQ(grid.com, 3,"PISM ERROR: PISMBedThermalUnit* btu == PETSC_NULL in enthalpyAndDrainageStep()");
  }

  ierr = artm.begin_access(); CHKERRQ(ierr);
  ierr = shelfbmassflux.begin_access(); CHKERRQ(ierr);
  ierr = shelfbtemp.begin_access(); CHKERRQ(ierr);

  // get other map-plane fields
  ierr = liqfrac_surface.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vbwat.begin_access(); CHKERRQ(ierr);
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

  MaskQuery mask(vMask);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      // for fine grid; this should *not* be replaced by call to grid.kBelowHeight()
      const PetscInt ks = static_cast<PetscInt>(floor(vH(i,j)/fdz));
#if (PISM_DEBUG==1)
      // check if ks is valid
      if ((ks < 0) || (ks >= grid.Mz_fine)) {
        PetscPrintf(grid.com,
                    "ERROR: ks = %d computed at i = %d, j = %d is invalid,"
                    " possibly because of invalid ice thickness.\n",
                    ks, i, j);
        SETERRQ(grid.com, 1, "invalid ks");
      }
#endif

      const bool ice_free_column = (ks == 0),
                 is_floating     = mask.ocean(i,j);

      // enthalpy and pressures at top of ice
      const PetscScalar p_ks = EC->getPressureFromDepth(vH(i,j) - fzlev[ks]); // FIXME issue #15
      PetscScalar Enth_ks;
      ierr = EC->getEnthPermissive(artm(i,j), liqfrac_surface(i,j), p_ks, Enth_ks); CHKERRQ(ierr);

      // deal completely with columns with no ice; enthalpy, vbwat, vbmr all need setting
      if (ice_free_column) {
        ierr = vWork3d.setColumn(i,j,Enth_ks); CHKERRQ(ierr);
        if (mask.floating_ice(i,j)) {
          // if floating then assume-maximally saturated till to avoid "shock"
          //   when grounding line advances
          vbwat(i,j) = bwat_max;
          vbmr(i,j) = shelfbmassflux(i,j);
        } else {
          // either truely no ice or grounded or both; either way zero-out subglacial fields
          vbwat(i,j) = 0.0;  // no stored water on ice free land
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

        ierr = Enth3.getValColumn(i,j,ks,esys->Enth); CHKERRQ(ierr);
        ierr = w3->getValColumn(i,j,ks,esys->w); CHKERRQ(ierr);

        ierr = getEnthalpyCTSColumn(p_air, vH(i,j), ks, &esys->Enth_s); CHKERRQ(ierr);

        PetscScalar lambda;
        ierr = getlambdaColumn(ks, ice_rho * default_ice_c, default_ice_k,
                               esys->Enth, esys->Enth_s, esys->w,
                               &lambda); CHKERRQ(ierr);
        if (lambda < 1.0)  *vertSacrCount += 1; // count columns with lambda < 1

        // if there is subglacial water, don't allow ice base enthalpy to be below
        // pressure-melting; that is, assume subglacial water is at the pressure-
        // melting temperature and enforce continuity of temperature
        if ((vbwat(i,j) > 0.0) && (esys->Enth[0] < esys->Enth_s[0])) { 
          esys->Enth[0] = esys->Enth_s[0];
        }

        const bool base_is_cold = (esys->Enth[0] < esys->Enth_s[0]);
        const PetscScalar p1 = EC->getPressureFromDepth(vH(i,j) - fdz); // FIXME issue #15
        const bool k1_istemperate = EC->isTemperate(esys->Enth[1], p1); // level  z = + \Delta z

        // can now determine melt, but only preliminarily because of drainage,
        //   from heat flux out of bedrock, heat flux into ice, and frictional heating
        if (is_floating) {
          vbmr(i,j) = shelfbmassflux(i,j);
        } else {
          if (base_is_cold) {
              vbmr(i,j) = 0.0;  // zero melt rate if cold base
          } else {
            const PetscScalar pbasal = EC->getPressureFromDepth(vH(i,j)); // FIXME issue #15
            PetscScalar hf_up;
            if (k1_istemperate) {
              const PetscScalar Tpmpbasal = EC->getMeltingTemp(pbasal);
              hf_up = - esys->k_from_T(Tpmpbasal) * (EC->getMeltingTemp(p1) - Tpmpbasal) / fdz;
            } else {
              PetscScalar Tbasal;
              ierr = EC->getAbsTemp(esys->Enth[0], pbasal, Tbasal); CHKERRQ(ierr);
              const PetscScalar Kbasal = esys->k_from_T(Tbasal) / EC->c_from_T(Tbasal);
              hf_up = - Kbasal * (esys->Enth[1] - esys->Enth[0]) / fdz;
            }

            // compute basal melt rate from flux balance; vbmr = - Mb / rho in
            //   efgis paper; after we compute it we make sure there is no
            //   refreeze if there is no available basal water
            vbmr(i,j) = ( (*Rb)(i,j) + G0(i,j) - hf_up ) / (ice_rho * L);

            if ((vbwat(i,j) <= 0) && (vbmr(i,j) < 0))
              vbmr(i,j) = 0.0;
          }
        }

        // now set-up for solve in ice; note esys->Enth[], esys->w[],
        //   esys->Enth_s[] are already filled
        ierr = esys->setIndicesAndClearThisColumn(i,j,ks); CHKERRQ(ierr);

        ierr = u3->getValColumn(i,j,ks,esys->u); CHKERRQ(ierr);
        ierr = v3->getValColumn(i,j,ks,esys->v); CHKERRQ(ierr);
        ierr = Sigma3->getValColumn(i,j,ks,esys->Sigma); CHKERRQ(ierr);

        ierr = esys->initThisColumn(isMarginal, lambda, vH(i, j)); CHKERRQ(ierr);
        ierr = esys->setBoundaryValuesThisColumn(Enth_ks); CHKERRQ(ierr);

        // determine lowest-level equation at bottom of ice; see decision chart
        //   in [\ref AschwandenBuelerKhroulevBlatter], and page documenting BOMBPROOF
        if (is_floating) {
          // floating base: Dirichlet application of known temperature from ocean
          //   coupler; assumes base of ice shelf has zero liquid fraction
          PetscScalar Enth0;
          ierr = EC->getEnthPermissive(shelfbtemp(i,j), 0.0, EC->getPressureFromDepth(vH(i,j)),
                                       Enth0); CHKERRQ(ierr);
          ierr = esys->setDirichletBasal(Enth0); CHKERRQ(ierr);
        } else if (base_is_cold) {
          // cold, grounded base (Neumann) case:  q . n = q_lith . n + F_b
          ierr = esys->setBasalHeatFlux(G0(i,j) + (*Rb)(i,j)); CHKERRQ(ierr);
        } else {
          // warm, grounded base case
          if (k1_istemperate) {
            // positive thickness of temperate ice; homogeneous Neumann case:  q . n = 0
            ierr = esys->setBasalHeatFlux(0.0); CHKERRQ(ierr);
          } else {
            // no thickness of temperate ice:  Dirichlet  H = H_s(pbasal)
            ierr = esys->setDirichletBasal(esys->Enth_s[0]); CHKERRQ(ierr);
          }
        }

        // solve the system
        PetscErrorCode pivoterr;
        ierr = esys->solveThisColumn(&Enthnew,pivoterr); CHKERRQ(ierr);
        if (pivoterr != 0) {
          ierr = PetscPrintf(PETSC_COMM_SELF,
            "\n\ntridiagonal solve of enthSystemCtx in enthalpyAndDrainageStep() FAILED at (%d,%d)\n"
                " with zero pivot position %d; viewing system to m-file ... \n",
            i, j, pivoterr); CHKERRQ(ierr);
          ierr = esys->reportColumnZeroPivotErrorMFile(pivoterr); CHKERRQ(ierr);
          SETERRQ(grid.com, 1,"PISM ERROR in enthalpyDrainageStep()\n");
        }
        if (viewOneColumn && issounding(i,j)) {
          ierr = PetscPrintf(PETSC_COMM_SELF,
            "\n\nin enthalpyAndDrainageStep(): viewing enthSystemCtx at (i,j)=(%d,%d) to m-file ... \n\n",
            i, j); CHKERRQ(ierr);
          ierr = esys->viewColumnInfoMFile(Enthnew, fMz); CHKERRQ(ierr);
        }

        // thermodynamic basal melt rate causes water to be added to layer
        PetscScalar bwatnew = vbwat(i,j);
        if (mask.grounded(i,j)) {
          bwatnew += vbmr(i,j) * dt_secs;
        }

        // drain ice segments by mechanism in [\ref AschwandenBuelerKhroulevBlatter],
        //   using DrainageCalculator dc
        PetscScalar Hdrainedtotal = 0.0;
        for (PetscInt k=0; k < ks; k++) {
          if (Enthnew[k] > esys->Enth_s[k]) { // avoid doing any more work if cold
            if (Enthnew[k] >= esys->Enth_s[k] + 0.5 * L) {
              liquifiedCount++; // count these rare events ...
              Enthnew[k] = esys->Enth_s[k] + 0.5 * L; //  but lose the energy
            }
            const PetscReal p = EC->getPressureFromDepth(vH(i,j) - fzlev[k]); // FIXME issue #15
            PetscReal omega;
            EC->getWaterFraction(Enthnew[k], p, omega);  // return code not checked
            if (omega > 0.01) {
              PetscReal fractiondrained = dc.get_drainage_rate(omega) * dt_secs; // pure number
              fractiondrained = PetscMin(fractiondrained, omega - 0.01); // only drain down to 0.01
              Hdrainedtotal += fractiondrained * fdz;  // always a positive contribution
              Enthnew[k] -= fractiondrained * L;
            }
          }
        }

        // in grounded case, add to both basal melt rate and bwat; if floating,
        // Hdrainedtotal is discarded because ocean determines basal melt rate
        if (mask.grounded(i,j)) {
          vbmr(i,j) += Hdrainedtotal / dt_secs;
          bwatnew += Hdrainedtotal;
        }

        // finalize Enthnew[]:  apply bulge limiter and transfer column
        //   into vWork3d; communication will occur later
        const PetscReal lowerEnthLimit = Enth_ks - bulgeEnthMax;
        for (PetscInt k=0; k < ks; k++) {
          if (Enthnew[k] < lowerEnthLimit) {
            *bulgeCount += 1;      // count the columns which have very large cold 
            Enthnew[k] = lowerEnthLimit;  // limit advection bulge ... enthalpy not too low
          }
        }
        ierr = vWork3d.setValColumnPL(i,j,Enthnew); CHKERRQ(ierr);

        // finalize bwat value
        bwatnew -= bwat_decay_rate * dt_secs;
        if (is_floating) {
          // if floating assume maximally saturated till to avoid "shock" if grounding line advances
          // UNACCOUNTED MASS & ENERGY (LATENT) LOSS/GAIN (TO/FROM OCEAN)!!
          vbwat(i,j) = bwat_max;
        } else {
          // limit bwat to be in [0.0, bwat_max]
          // UNACCOUNTED MASS & ENERGY (LATENT) LOSS (TO INFINITY AND BEYOND)!!
          vbwat(i,j) = PetscMax(0.0, PetscMin(bwat_max, bwatnew) );
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
  ierr = vbwat.end_access(); CHKERRQ(ierr);
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
  delete esys;

  *liquifiedVol = ((double) liquifiedCount) * fdz * grid.dx * grid.dy;
  return 0;
}

