// Copyright (C) 2009 Andreas Aschwanden and Ed Bueler and Constantine Khroulev
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

#include "iceEnthalpyModel.hh"
#include "enthColumnSystem.hh"
#include "enthalpyConverter.hh"


/*********** procedures for init ****************/

IceEnthalpyModel::IceEnthalpyModel(IceGrid &g) : IceModel(g) {

  doColdIceMethods = PETSC_FALSE;     // default is to actually use enthalpy for a polythermal model

  bmr_in_pore_pressure = PETSC_FALSE; // default to bwat-only model for pore water pressure
  bmr_enhance_scale = 0.10 / secpera; // default to say that basal melt rate must be on order of
                                      //   10 cm/a to start making major difference in pore pressure

  thk_affects_pore_pressure = PETSC_FALSE; // default is that basal water weakening of till
                                           //   (from increased pore water pressure) is un-affected
                                           //   by ice thickness, which is a surrogate for distance 
                                           //   to margin
  margin_pore_pressure_H_high  = 2000.0; // defaults for the mechanism in
  margin_pore_pressure_H_low   = 1000.0; //   getEffectivePressureOnTill() which reduces
  margin_pore_pressure_reduced = 0.97;  //   pore pressure near margin; H_high, H_low are thicknesses

}


PetscErrorCode IceEnthalpyModel::createVecs() {
  PetscErrorCode ierr;

  ierr = Enth3.create(grid, "enthalpy", true); CHKERRQ(ierr);
  // POSSIBLE standard name = land_ice_enthalpy
  ierr = Enth3.set_attrs(
     "model_state",
     "ice enthalpy (sensible plus latent heat, plus potential energy of pressure)",
     "J kg-1",
     ""); CHKERRQ(ierr);

  ierr = IceModel::createVecs(); CHKERRQ(ierr);

  // see IceModel::allocate_internal_objects(), which is where this should go
  ierr = EnthNew3.create(grid,"enthalpy_new",false); CHKERRQ(ierr); // global
  ierr = EnthNew3.set_attrs(
     "internal",
     "ice enthalpy; temporary space during timestep",
     "J kg-1",
     ""); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceEnthalpyModel::setFromOptions() {
  PetscErrorCode ierr;

  ierr = IceModel::setFromOptions(); CHKERRQ(ierr);

  // if set, use old IceModel::temperatureStep(), and set enthalpy as though
  //   ice is cold
  ierr = check_option("-cold", doColdIceMethods); CHKERRQ(ierr);

  // if set, add function of instantaneous basal melt rate to pore water pressure
  //   computed by the usual model, which is a function of stored till water thickness
  ierr = check_option("-bmr_enhance", bmr_in_pore_pressure); CHKERRQ(ierr);
  // bmr_enhance_scale is set by giving a basal melt rate in m a-1
  //   this amount of basal melt is then a scale for the basal melt rate
  //   weakening effect from -bmr_enhance;
  //   if -bmr_enhance_scale is set to some value, then bmr_in_pore_pressure
  //   is set to true;  default is "-bmr_enhance_scale 0.10";
  //   see getEffectivePressureOnTill()
  PetscScalar bmres;
  PetscTruth  bmres_set;
  ierr = PetscOptionsGetReal(PETSC_NULL, "-bmr_enhance_scale", &bmres, 
                             &bmres_set);  CHKERRQ(ierr);
  if (bmres_set == PETSC_TRUE) {
    bmr_in_pore_pressure = PETSC_TRUE;
    bmr_enhance_scale = bmres / secpera;
  }

  // if set, makes the thickness affect the pore_pressure; near margin there
  //   is a reduction in pore pressure, a conceptual drainage mechanism
  ierr = check_option("-thk_eff", thk_affects_pore_pressure); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL, "-margin_reduced",
                             &margin_pore_pressure_reduced, PETSC_NULL);  CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL, "-margin_H_high", 
                             &margin_pore_pressure_H_high, PETSC_NULL);  CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL, "-margin_H_low",
                             &margin_pore_pressure_H_low, PETSC_NULL);  CHKERRQ(ierr);

  // DEBUG:  report settings
  const int vlevel = 2;
  ierr = verbPrintf(vlevel, grid.com,
      "  IceEnthalpyModel::setFromOptions():\n"); CHKERRQ(ierr);
  ierr = verbPrintf(vlevel, grid.com,
      "    doColdIceMethods is %s\n", (doColdIceMethods == PETSC_TRUE) ? "TRUE" : "FALSE");
      CHKERRQ(ierr);
  ierr = verbPrintf(vlevel, grid.com,
      "    bmr_in_pore_pressure is %s\n", (bmr_in_pore_pressure == PETSC_TRUE) ? "TRUE" : "FALSE");
      CHKERRQ(ierr);
  if (bmr_in_pore_pressure == PETSC_TRUE) {
    ierr = verbPrintf(vlevel, grid.com,
      "      bmr_enhance_scale = %f m a-1\n", bmr_enhance_scale * secpera);
      CHKERRQ(ierr);
  }
  ierr = verbPrintf(vlevel, grid.com,
      "    thk_affects_pore_pressure is %s\n", (thk_affects_pore_pressure == PETSC_TRUE) ? "TRUE" : "FALSE");
      CHKERRQ(ierr);
  if (thk_affects_pore_pressure == PETSC_TRUE) {
    ierr = verbPrintf(vlevel, grid.com,
      "      margin_pore_pressure_reduced = %f\n", margin_pore_pressure_reduced);
      CHKERRQ(ierr);
    ierr = verbPrintf(vlevel, grid.com,
      "      margin_pore_pressure_H_high = %f m\n", margin_pore_pressure_H_high);
      CHKERRQ(ierr);
    ierr = verbPrintf(vlevel, grid.com,
      "      margin_pore_pressure_H_low = %f m\n", margin_pore_pressure_H_low);
      CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode IceEnthalpyModel::initFromFile(const char *fname) {
  PetscErrorCode  ierr;

  ierr = IceModel::initFromFile(fname); CHKERRQ(ierr);

  // kludge: add it to the dictionary now (as opposed to before the
  // initFromFile call above) so that penth does not stop if 'enthalpy' in not
  // found in the -i file:
  ierr = variables.add(Enth3); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
     "  entering IceEnthalpyModel::initFromFile() after base class version;\n"
     "  looking in '%s' for variable 'enthalpy' ... \n",fname);
     CHKERRQ(ierr);

  NCTool nc(&grid);
  ierr = nc.open_for_reading(fname); CHKERRQ(ierr);

  // Find the index of the last record in the file:
  int last_record;
  ierr = nc.get_dim_length("t", &last_record); CHKERRQ(ierr);
  last_record -= 1;

  bool enthalpy_exists=false;
  ierr = nc.find_variable("enthalpy", NULL, enthalpy_exists); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  if (enthalpy_exists) {
    ierr = Enth3.read(fname, last_record); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com,
      "  variable 'enthalpy' not found so setting it as cold ice, from temperature ...\n");
      CHKERRQ(ierr);
    ierr = setEnth3FromT3_ColdIce(); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode IceEnthalpyModel::bootstrapFromFile(const char *filename) {
  PetscErrorCode ierr;
  ierr = IceModel::bootstrapFromFile(filename); CHKERRQ(ierr);
  
  ierr = verbPrintf(2, grid.com, "continuing bootstrapping in IceEnthalpyModel ...\n"); CHKERRQ(ierr);
  ierr = verbPrintf(2, grid.com,
      "  ice enthalpy set from temperature, as cold ice (zero liquid fraction)\n");
      CHKERRQ(ierr);
  ierr = setEnth3FromT3_ColdIce(); CHKERRQ(ierr);
  ierr = verbPrintf(2, grid.com, "bootstrapping done (IceEnthalpyModel)\n"); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceEnthalpyModel::init_physics() {
  PetscErrorCode ierr;

  // let the base class create the ice and process its options:
  ierr = IceModel::init_physics(); CHKERRQ(ierr);

  if (grid.Mbz == 2) {
    SETERRQ(2,"IceEnthalpyModel does not allow grid.Mbz==2, because separation\n"
              "between ice/bedrock interface and application of geothermal flux\n"
              "is much easier.  Note grid.Mbz==1 and grid.Mbz>=3 are the allowed values.\n");
  }

  ierr = verbPrintf(2, grid.com,
      "  setting flow law to Glen-Paterson-Budd-Lliboutry-Duval type ...\n");
      CHKERRQ(ierr);
  if (ice != NULL)  delete ice;  // kill choice already made!
  iceFactory.setType(ICE_GPBLD); // new flowlaw which has dependence on enthalpy not temperature
  iceFactory.create(&ice);

  PolyThermalGPBLDIce *gpbldi = dynamic_cast<PolyThermalGPBLDIce*>(ice);
  if (gpbldi == NULL) {
    ThermoGlenIce *tgi = dynamic_cast<ThermoGlenIce*>(ice);
    if (tgi) {
      ierr = verbPrintf(2, grid.com,
        "  [flow law was actually set to ThermoGlenIce by IceEnthalpyModel ...]\n"); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2, grid.com,
        "  [WARNING: flow law unclear in IceEnthalpyModel ...]\n"); CHKERRQ(ierr);
    }
  }
  
  ierr = ice->printInfo(4);CHKERRQ(ierr);

  ierr = ice->setFromOptions();CHKERRQ(ierr);

  return 0;
}



/*********** procedures for read/write ****************/

PetscErrorCode IceEnthalpyModel::write_extra_fields(const char* filename) {
  PetscErrorCode ierr;

  ierr = Enth3.write(filename, NC_DOUBLE); CHKERRQ(ierr);//! Total code duplication with IceModel version, but checks flag doColdIceMethods and uses correct flow law.


  // also write omega = liquid water fraction
  //   we use EnthNew3 (global) as temporary, allocated space for this purpose
  ierr = verbPrintf(4, grid.com,
      "  writing liquid water fraction 'liquid_frac' from enthalpy ...\n"); CHKERRQ(ierr);
  ierr = setLiquidFracFromEnthalpy(EnthNew3); CHKERRQ(ierr);
  ierr = EnthNew3.write(filename, NC_FLOAT); CHKERRQ(ierr);

  // also write temp_pa = pressure-adjusted temp in Celcius
  //   again use EnthNew3 (global) as temporary, allocated space for this purpose
  ierr = verbPrintf(4, grid.com,
      "  writing pressure-adjusted ice temperature (deg C) 'temp_pa' ...\n"); CHKERRQ(ierr);
  ierr = setPATempFromEnthalpy(EnthNew3); CHKERRQ(ierr); // returns K
  ierr = EnthNew3.shift(- config.get("water_melting_temperature")); CHKERRQ(ierr); // make deg C
  ierr = EnthNew3.write(filename, NC_FLOAT); CHKERRQ(ierr);

  // write CTS position (unitless) if command line option -cts is given
  //   again use EnthNew3 (global) as temporary, allocated space for this purpose
  PetscTruth userWantsCTS;
  ierr = check_option("-cts", userWantsCTS); CHKERRQ(ierr);
  if (userWantsCTS == PETSC_TRUE) {
    ierr = verbPrintf(4, grid.com,
		      "  writing CTS position E/Es (-) 'cts' ...\n"); CHKERRQ(ierr);
    ierr = setCTSFromEnthalpy(EnthNew3); CHKERRQ(ierr); // returns K
    ierr = EnthNew3.write(filename, NC_FLOAT); CHKERRQ(ierr);
  }

  // reset attributes on EnthNew3, a temporary; probaly not needed
  ierr = EnthNew3.set_name("enthalpy_new"); CHKERRQ(ierr);
  ierr = EnthNew3.set_attrs(
     "internal",
     "ice enthalpy; temporary space during timestep",
     "J kg-1",
     ""); CHKERRQ(ierr);

  // also write bfrict = friction heating, in units mW m-2 like geothermal heat
  ierr = verbPrintf(4, grid.com,
      "  writing basal frictional heating 'bfrict' ...\n"); CHKERRQ(ierr);
  ierr = vRb.write(filename, NC_FLOAT); CHKERRQ(ierr);

  // also write bmelt = ice basal melt rate in m a-1
  ierr = verbPrintf(4, grid.com,
      "  writing ice basal melt rate 'bmelt' ...\n"); CHKERRQ(ierr);
  ierr = vbasalMeltRate.write(filename, NC_FLOAT); CHKERRQ(ierr);

  return 0;
}

/*********** setting fields ****************/

//! Compute Enth3 from temperature T3 by assuming the ice has zero liquid fraction.
/*!
Forces temperature to at most the pressure melting point, before computing
the enthalpy for that temperature and zero liquid fraction.
 */
PetscErrorCode IceEnthalpyModel::setEnth3FromT3_ColdIce() {
  PetscErrorCode ierr;

  EnthalpyConverter EC(config);
  
  PetscScalar **H;
  ierr = T3.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = vH.get_array(H); CHKERRQ(ierr);

  PetscScalar *Tij, *Enthij; // columns of these values
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = T3.getInternalColumn(i,j,&Tij); CHKERRQ(ierr);
      ierr = Enth3.getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
      for (PetscInt k=0; k<grid.Mz; ++k) {
        const PetscScalar depth = H[i][j] - grid.zlevels[k];
        // because of how getPressureFromDepth() works, the
        //   energy content in the air is set to the value ice would have if it a chunk of it
        //   occupied the air; the atmosphere actually has much lower energy content;
        //   done this way for regularity (i.e. dEnth/dz computations)
        ierr = EC.getEnthPermissive(Tij[k],0.0,EC.getPressureFromDepth(depth), Enthij[k] );
           CHKERRQ(ierr);
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


//! Compute the ice temperature corresponding to Enth3, and put in Tnew3; use just after Enth3 is determined.
/*!
Does not communicate.  Ghosts will be invalid, but "T3.endGhostCommTransfer(Tnew3)" in
IceModel::temperatureAgeStep() will have desired effect.
 */
PetscErrorCode IceEnthalpyModel::setTnew3FromEnth3() {
  PetscErrorCode ierr;

  EnthalpyConverter EC(config);

  PetscScalar **thickness;
  PetscScalar *Tij, *Enthij; // columns of these values
  ierr = Tnew3.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = vH.get_array(thickness); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = Tnew3.getInternalColumn(i,j,&Tij); CHKERRQ(ierr);
      ierr = Enth3.getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
      for (PetscInt k=0; k<grid.Mz; ++k) {
        const PetscScalar depth = thickness[i][j] - grid.zlevels[k];
        ierr = EC.getAbsTemp(Enthij[k],EC.getPressureFromDepth(depth), Tij[k]); CHKERRQ(ierr);
      }
    }
  }
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = Tnew3.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  return 0;
}


//! Compute the liquid fraction corresponding to Enth3, and put in a global IceModelVec3 provided by user.
PetscErrorCode IceEnthalpyModel::setLiquidFracFromEnthalpy(IceModelVec3 &useForLiquidFrac) {
  PetscErrorCode ierr;

  ierr = useForLiquidFrac.set_name("liqfrac"); CHKERRQ(ierr);
  ierr = useForLiquidFrac.set_attrs(
     "diagnostic",
     "liquid water fraction in ice (between 0 and 1)",
     "",
     ""); CHKERRQ(ierr);

  EnthalpyConverter EC(config);

  PetscScalar **thickness;
  PetscScalar *omegaij, *Enthij; // columns of these values
  ierr = useForLiquidFrac.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = vH.get_array(thickness); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = useForLiquidFrac.getInternalColumn(i,j,&omegaij); CHKERRQ(ierr);
      ierr = Enth3.getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
      for (PetscInt k=0; k<grid.Mz; ++k) {
        const PetscScalar depth = thickness[i][j] - grid.zlevels[k];
        ierr = EC.getWaterFraction(Enthij[k],EC.getPressureFromDepth(depth), omegaij[k]);
          CHKERRQ(ierr);
      }
    }
  }
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = useForLiquidFrac.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  // communication not done; we allow global IceModelVec3s as useForLiquidFrac

  return 0;
}


//! Compute the pressure-adjusted temperature corresponding to Enth3, and put in a global IceModelVec3 provided by user.
PetscErrorCode IceEnthalpyModel::setPATempFromEnthalpy(IceModelVec3 &useForPATemp) {
  PetscErrorCode ierr;

  ierr = useForPATemp.set_name("temp_pa"); CHKERRQ(ierr);
  ierr = useForPATemp.set_attrs(
     "diagnostic",
     "pressure-adjusted ice temperature (degrees above pressure-melting point)",
     "deg_C",
     ""); CHKERRQ(ierr);

  EnthalpyConverter EC(config);

  PetscScalar **thickness;
  PetscScalar *Tpaij, *Enthij; // columns of these values
  ierr = useForPATemp.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = vH.get_array(thickness); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = useForPATemp.getInternalColumn(i,j,&Tpaij); CHKERRQ(ierr);
      ierr = Enth3.getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
      for (PetscInt k=0; k<grid.Mz; ++k) {
        const PetscScalar depth = thickness[i][j] - grid.zlevels[k];
        ierr = EC.getPATemp(Enthij[k],EC.getPressureFromDepth(depth), Tpaij[k]);
          CHKERRQ(ierr);
      }
    }
  }
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = useForPATemp.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  // communication not done; we allow global IceModelVec3s as useForPATemp
  return 0;
}



//! Compute CTS = E/E_s from Enth3 and Es, and put in a global IceModelVec3 provided by user.
/*!
The actual cold-temperate transition surface (CTS) is the level set CTS = E/E_s(p) = 1.
 */
PetscErrorCode IceEnthalpyModel::setCTSFromEnthalpy(IceModelVec3 &useForCTS) {
  PetscErrorCode ierr;

  ierr = useForCTS.set_name("cts"); CHKERRQ(ierr);
  ierr = useForCTS.set_attrs(
     "diagnostic",
     "cts = E/E_s(p), so cold-temperate transition surface is at cts = 1",
     "",
     ""); CHKERRQ(ierr);

  EnthalpyConverter EC(config);

  PetscScalar **thickness;
  PetscScalar *CTSij, *Enthij; // columns of these values
  ierr = useForCTS.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = vH.get_array(thickness); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = useForCTS.getInternalColumn(i,j,&CTSij); CHKERRQ(ierr);
      ierr = Enth3.getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
      for (PetscInt k=0; k<grid.Mz; ++k) {
        const PetscScalar depth = thickness[i][j] - grid.zlevels[k];
        CTSij[k] = EC.getCTS(Enthij[k], EC.getPressureFromDepth(depth));
      }
    }
  }
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = useForCTS.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  return 0;
}


/*********** modified reporting to standard out ******************/

PetscErrorCode IceEnthalpyModel::energyAgeStats(
                    PetscScalar ivol, PetscScalar iarea, bool /* useHomoTemp */, 
                    PetscScalar &gmeltfrac, PetscScalar &gtemp0, PetscScalar &gorigfrac) {
  PetscErrorCode  ierr;
  PetscScalar     **H, **Enthbase, *tau;
  PetscScalar     meltarea, temp0, origvol;
  
  EnthalpyConverter EC(config);

  // put basal ice enthalpy in vWork2d[0]
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.getHorSlice(vWork2d[0], 0.0); CHKERRQ(ierr);  // z=0 slice
  ierr = Enth3.end_access(); CHKERRQ(ierr);

  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = vWork2d[0].get_array(Enthbase); CHKERRQ(ierr);
  ierr = tau3.begin_access(); CHKERRQ(ierr);

  const PetscScalar   a = grid.dx * grid.dy * 1e-3 * 1e-3; // area unit (km^2)
  const PetscScalar   currtime = grid.year * secpera;

  meltarea = 0.0; temp0 = 0.0; origvol = 0.0;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (H[i][j] > 0) {

        // accumulate area of base which is at melt point
        if (EC.isTemperate(Enthbase[i][j], EC.getPressureFromDepth(H[i][j]) ))  
          meltarea += a;

        // accumulate volume of ice which is original
        ierr = tau3.getInternalColumn(i,j,&tau); CHKERRQ(ierr);
        const PetscInt  ks = grid.kBelowHeight(H[i][j]);
        for (PetscInt k=1; k<=ks; k++) {
          // ice is original if it is at least one year older than current time
          if (0.5*(tau[k-1]+tau[k]) > currtime + secpera)
            origvol += a * 1.0e-3 * (grid.zlevels[k] - grid.zlevels[k-1]);
        }

      }

      // if you happen to be at center, record absolute basal temp there
      if (i == (grid.Mx - 1)/2 && j == (grid.My - 1)/2) {
        ierr = EC.getAbsTemp(Enthbase[i][j],EC.getPressureFromDepth(H[i][j]), temp0);
          CHKERRQ(ierr);
      }
    }
  }
  
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vWork2d[0].end_access(); CHKERRQ(ierr);
  ierr = tau3.end_access(); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&meltarea, &gmeltfrac, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&origvol,  &gorigfrac, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&temp0,    &gtemp0,    grid.com); CHKERRQ(ierr);

  // normalize fractions correctly
  if (ivol > 0.0)    gorigfrac = gorigfrac / ivol;
  else gorigfrac = 0.0;
  if (iarea > 0.0)   gmeltfrac = gmeltfrac / iarea;
  else gmeltfrac = 0.0;

  return 0;
}


/*********** timestep routines ****************/

PetscErrorCode IceEnthalpyModel::temperatureStep(
     PetscScalar* vertSacrCount, PetscScalar* bulgeCount) {
  PetscErrorCode ierr;
  if (doColdIceMethods==PETSC_TRUE) {
    ierr = verbPrintf(4,grid.com,
      "    [IceEnthalpyModel::temperatureStep(): ENTHALPY IS OFF. CALLING IceModel::temperatureStep()]\n");
    CHKERRQ(ierr);

    ierr = IceModel::temperatureStep(vertSacrCount,bulgeCount);  CHKERRQ(ierr);

    // start & complete communication: THIS IS REDUNDANT WITH temperatureAgeStep(),
    //    BUT NEEDED TO GET UPDATED ENTHALPY AT END OF TIME STEP
    ierr = T3.beginGhostCommTransfer(Tnew3); CHKERRQ(ierr);
    ierr = T3.endGhostCommTransfer(Tnew3); CHKERRQ(ierr);

    ierr = setEnth3FromT3_ColdIce();  CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(4,grid.com,
      "    [IceEnthalpyModel::temperatureStep(): ENTHALPY IS ON. CALLING IceEnthalpyModel::enthalpyAndDrainageStep()]\n");
      CHKERRQ(ierr);
    // new enthalpy values go in EnthNew3; also updates (and communicates) Hmelt
    // enthalpyStep() is in iceEnthalpyModel.cc
    ierr = enthalpyAndDrainageStep(vertSacrCount,bulgeCount);  CHKERRQ(ierr);

    // start & complete communication
    ierr = Enth3.beginGhostCommTransfer(EnthNew3); CHKERRQ(ierr);
    ierr = Enth3.endGhostCommTransfer(EnthNew3); CHKERRQ(ierr);

    ierr = setTnew3FromEnth3();  CHKERRQ(ierr);  // temperatureAgeStep() ASSUMES Tnew3 valid
  }
  return 0;
}


PetscErrorCode IceEnthalpyModel::enthalpyAndDrainageStep(PetscScalar* vertSacrCount, PetscScalar* bulgeCount) {
  PetscErrorCode  ierr;

  if (doColdIceMethods==PETSC_TRUE) {
    SETERRQ(1,"\n\n    IceEnthalpyModel::enthalpyAndDrainageStep() called but doColdIceMethods==true ... ending\n");
  }

  // set up fine grid in ice and bedrock
  PetscInt    fMz, fMbz;
  PetscScalar fdz, *fzlev, fdzb, *fzblev;

  ierr = grid.get_fine_vertical_grid(fMz, fMbz, fdz, fdzb, fzlev, fzblev); CHKERRQ(ierr);

  ierr = verbPrintf(4,grid.com,
    "\n  [entering enthalpyAndDrainageStep(); fMz = %d, fdz = %5.3f, fMbz = %d, fdzb = %5.3f]",
    fMz, fdz, fMbz, fdzb); CHKERRQ(ierr);

  PetscTruth viewOneColumn;
  ierr = check_option("-view_sys", viewOneColumn); CHKERRQ(ierr);

  PetscTruth viewOneRedoColumn;
  ierr = check_option("-view_redo_sys", viewOneRedoColumn); CHKERRQ(ierr);

  EnthalpyConverter EC(config);
  if (getVerbosityLevel() >= 4) { ierr = EC.viewConstants(NULL); CHKERRQ(ierr); }

  enthSystemCtx system(fMz,fMbz);
  system.dx              = grid.dx;
  system.dy              = grid.dy;
  system.dtTemp          = dtTempAge; // same time step for temp and age, currently
  system.dzEQ            = fdz;
  system.dzbEQ           = fdzb;
  system.ice_rho         = config.get("ice_density"); // ice->rho;
  system.ice_c           = config.get("ice_specific_heat_capacity"); // ice->c_p;
  system.ice_k           = config.get("ice_thermal_conductivity"); // ice->k;
  system.ice_nu          = config.get("enthalpy_temperate_diffusivity"); // diffusion in temperate ice
  system.bed_thermal_rho = config.get("bedrock_thermal_density");
  system.bed_thermal_c   = config.get("bedrock_thermal_specific_heat_capacity");
  system.bed_thermal_k   = config.get("bedrock_thermal_conductivity");

  // space for solution of system; length = fMz + fMbz - 1
  const PetscInt k0 = fMbz - 1;
  PetscScalar *x;
  x = new PetscScalar[fMz + k0];

  // constants (needed after solution of system, in insertion phase, or for boundary values)
  const PetscScalar rho_c_I   = system.ice_rho * system.ice_c,
                    p_air     = config.get("surface_pressure"),          // Pa
                    L         = config.get("water_latent_heat_fusion"),  // J kg-1
                    omega_max = config.get("liquid_water_fraction_max"); // pure
  
  // this is bulge limit constant in J kg-1; is max amount by which ice
  //   enthalpy can be lower than surface temperature (as an enthalpy);
  //   value is enthalpy change equivalent to change in cold ice temp by 15 K
  const PetscScalar bulgeMaxTemp = 15.0,
                    bulgeMaxEnth = system.ice_c * bulgeMaxTemp;

  // pointers to values in current column
  PetscScalar *Enthnew, *Tb, *Tbnew;

  system.u     = new PetscScalar[fMz];
  system.v     = new PetscScalar[fMz];
  system.w     = new PetscScalar[fMz];
  system.Sigma = new PetscScalar[fMz];
  system.Enth  = new PetscScalar[fMz];
  system.Enth_s= new PetscScalar[fMz];
  system.Enth_b= new PetscScalar[fMbz];

  Enthnew      = new PetscScalar[fMz];
  Tb           = new PetscScalar[fMbz];
  Tbnew        = new PetscScalar[fMbz];

  // system needs access to Enth3 for planeStar()
  system.Enth3 = &Enth3;

  // checks that all needed constants and pointers got set:
  ierr = system.initAllColumns(); CHKERRQ(ierr);

  PetscScalar *xredo;
  xredo = new PetscScalar[fMbz];
  bedrockOnlySystemCtx bedredosystem(fMbz);
  bedredosystem.dtTemp          = dtTempAge; // same time step for temp and age, currently
  bedredosystem.dzbEQ           = fdzb;
  bedredosystem.bed_thermal_rho = config.get("bedrock_thermal_density");
  bedredosystem.bed_thermal_c   = config.get("bedrock_thermal_specific_heat_capacity");
  bedredosystem.bed_thermal_k   = config.get("bedrock_thermal_conductivity");
  bedredosystem.T_b             = new PetscScalar[fMbz];
  // checks that all needed constants and pointers got set:
  ierr = bedredosystem.initAllColumns(); CHKERRQ(ierr);

  // now get map-plane coupler fields
  IceModelVec2 *pccTs, *pccsbt, *pccsbmf;
  if (atmosPCC != PETSC_NULL) {
    // call sets pccTs to point to IceModelVec2 with current surface temps
    ierr = atmosPCC->updateSurfTempAndProvide(
              grid.year, dtTempAge / secpera, pccTs);
              CHKERRQ(ierr);
  } else {
    SETERRQ(3,"PISM ERROR: atmosPCC == PETSC_NULL");
  }
  if (oceanPCC != PETSC_NULL) {
    ierr = oceanPCC->updateShelfBaseTempAndProvide(
              grid.year, dt / secpera, pccsbt);
              CHKERRQ(ierr);
    ierr = oceanPCC->updateShelfBaseMassFluxAndProvide(
              grid.year, dt / secpera, pccsbmf);
              CHKERRQ(ierr);
  } else {
    SETERRQ(4,"PISM ERROR: oceanPCC == PETSC_NULL");
  }
  PetscScalar  **Ts, **Tshelfbase, **bmr_float;
  ierr = pccTs->get_array(Ts);  CHKERRQ(ierr);
  ierr = pccsbt->get_array(Tshelfbase);  CHKERRQ(ierr);
  ierr = pccsbmf->get_array(bmr_float);  CHKERRQ(ierr);

  // get other map-plane fields
  PetscScalar  **H, **Ghf, **mask, **Hmelt, **Rb, **basalMeltRate;
  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = vHmelt.get_array(Hmelt); CHKERRQ(ierr);
  ierr = vbasalMeltRate.get_array(basalMeltRate); CHKERRQ(ierr);
  ierr = vMask.get_array(mask); CHKERRQ(ierr);
  ierr = vRb.get_array(Rb); CHKERRQ(ierr);
  ierr = vGhf.get_array(Ghf); CHKERRQ(ierr);

  // these are accessed a column at a time
  ierr = u3.begin_access(); CHKERRQ(ierr);
  ierr = v3.begin_access(); CHKERRQ(ierr);
  ierr = w3.begin_access(); CHKERRQ(ierr);
  ierr = Sigma3.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = EnthNew3.begin_access(); CHKERRQ(ierr);
  ierr = Tb3.begin_access(); CHKERRQ(ierr);

  PetscScalar max_hmelt = config.get("max_hmelt");

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      // for fine grid; this should *not* be replaced by call to grid.kBelowHeight()
      const PetscInt  ks = static_cast<PetscInt>(floor(H[i][j]/fdz));

      // enthalpy and pressures at boundaries
      const PetscScalar p_basal = EC.getPressureFromDepth(H[i][j]),
                        p_ks   = EC.getPressureFromDepth(H[i][j] - fzlev[ks]);
      PetscScalar Enth_air, Enth_ks, Enth_shelfbase;
      ierr = EC.getEnthPermissive(Ts[i][j], 0.0, p_air, Enth_air ); CHKERRQ(ierr);
      // in theory we could have a water fraction at k=ks level, but for
      //   now there is no case where we have that:
      ierr = EC.getEnthPermissive(Ts[i][j], 0.0, p_ks, Enth_ks ); CHKERRQ(ierr);
      // at underside of ice shelf, set enthalpy to that of max liquid water
      //   temperate ice; probably does not make much difference anyway because of
      //   no upward/downward liquid water transport:
      ierr = EC.getEnthPermissive(Tshelfbase[i][j], omega_max, p_basal,
                              Enth_shelfbase); CHKERRQ(ierr);

      if (k0+ks>0) { // if there are enough points in bedrock&ice to bother ...
        ierr = system.setIndicesThisColumn(i,j,ks); CHKERRQ(ierr);
        ierr = Tb3.getValColumnPL(i,j,fMbz,fzblev,Tb); CHKERRQ(ierr);

        if (grid.ice_vertical_spacing == EQUAL) {
          ierr = u3.getValColumnPL(i,j,fMz,fzlev,system.u); CHKERRQ(ierr);
          ierr = v3.getValColumnPL(i,j,fMz,fzlev,system.v); CHKERRQ(ierr);
          ierr = w3.getValColumnPL(i,j,fMz,fzlev,system.w); CHKERRQ(ierr);
          ierr = Sigma3.getValColumnPL(i,j,fMz,fzlev,system.Sigma); CHKERRQ(ierr);
          ierr = Enth3.getValColumnPL(i,j,fMz,fzlev,system.Enth); CHKERRQ(ierr);
        } else {
          // slower, but right for not-equal spaced
          ierr = u3.getValColumnQUAD(i,j,fMz,fzlev,system.u); CHKERRQ(ierr);
          ierr = v3.getValColumnQUAD(i,j,fMz,fzlev,system.v); CHKERRQ(ierr);
          ierr = w3.getValColumnQUAD(i,j,fMz,fzlev,system.w); CHKERRQ(ierr);
          ierr = Sigma3.getValColumnQUAD(i,j,fMz,fzlev,system.Sigma); CHKERRQ(ierr);
          ierr = Enth3.getValColumnQUAD(i,j,fMz,fzlev,system.Enth); CHKERRQ(ierr);
        }

        // go through column and find appropriate lambda for BOMBPROOF;
        //   at the same time fill system.Enth_s[];
        PetscScalar lambda = 1.0;  // start with centered implicit for more accuracy
        system.Enth_s[0] = EC.getEnthalpyCTS( p_basal );
        for (PetscInt k = 1; k <= ks; k++) {
          system.Enth_s[k] = EC.getEnthalpyCTS( EC.getPressureFromDepth(H[i][j] - fzlev[k]) );
          if (system.Enth[k] > system.Enth_s[k]) {
            // if there is a liquid water fraction we will switch to upwind;
            //    conductivity goes to zero
            lambda = 0.0;
          } else {
            const PetscScalar denom = (PetscAbs(system.w[k]) + 0.000001/secpera)
                                       * rho_c_I * fdz;
            lambda = PetscMin(lambda, 2.0 * system.ice_k / denom);
          }
        }
        if (lambda < 1.0)  *vertSacrCount += 1; // count columns with lambda < 1
        for (PetscInt k = ks+1; k < fMz; k++) {
          system.Enth_s[k] = EC.getEnthalpyCTS(p_air);
        }

        // go though bedrock and set system enthalpy from temperature
        for (PetscInt k=0; k < fMbz; k++) {
          ierr = EC.getEnthBedrock(Tb[k], system.Enth_b[k]); CHKERRQ(ierr);
        }

        // if isMarginal then only do vertical conduction for ice;
        //   will ignore advection and strain heating if isMarginal
        const bool isMarginal = checkThinNeigh(H[i+1][j],H[i+1][j+1],H[i][j+1],H[i-1][j+1],
                                               H[i-1][j],H[i-1][j-1],H[i][j-1],H[i+1][j-1]);

        ierr = system.setSchemeParamsThisColumn(mask[i][j], isMarginal, lambda); CHKERRQ(ierr);

        // set boundary values for tridiagonal system
        ierr = system.setSurfaceBoundaryValuesThisColumn(Enth_ks); CHKERRQ(ierr);
        ierr = system.setBasalBoundaryValuesThisColumn(
                 Ghf[i][j],Enth_shelfbase,Rb[i][j]); CHKERRQ(ierr);

        // solve the system for this column: x will contain new enthalpy in ice and in bedrock
        //   that is, x[k] is always an enthalpy
        ierr = system.solveThisColumn(&x); // no CHKERRQ(ierr) immediately because:
        if (ierr > 0) {
          SETERRQ3(2,
            "Tridiagonal solve failed at (%d,%d) with zero pivot position %d.\n", i, j, ierr);
        } else { CHKERRQ(ierr); }

        // diagnostic/debug
        if (viewOneColumn) {
          if ((i==id) && (j==jd)) {
            ierr = verbPrintf(1,grid.com,
              "\nin IceEnthalpyModel::enthalpyDrainageStep();\n"
                "   fMz = %d, fdz = %5.3f, fMbz = %d, fdzb = %5.3f, k0 = %d, ks = %d\n\n",
              fMz, fdz, fMbz, fdzb, k0, ks); CHKERRQ(ierr);
            ierr = verbPrintf(1,grid.com,
              "viewing system and solution at (i,j)=(%d,%d):\n", i, j); CHKERRQ(ierr);
            ierr = system.viewConstants(NULL); CHKERRQ(ierr);
            ierr = system.viewSystem(NULL,"system"); CHKERRQ(ierr);
            ierr = system.viewColumnValues(NULL, x, fMz+k0, "solution x"); CHKERRQ(ierr);
          }
        }

      }

      // top ice level; no possibility of drainage (a separate modeling issue!)
      Enthnew[ks] = x[k0 + ks];

      // now that enthalpy is known, check for and correct any extreme advection bulges
      bool bulgeLimiterTripped = false;
      for (PetscInt k=0; k < ks; k++) {
        Enthnew[k] = x[k0 + k];
        if (Enthnew[k] < Enthnew[ks] - bulgeMaxEnth) {
          Enthnew[k] = Enthnew[ks] - bulgeMaxEnth;
          bulgeCount++;
          bulgeLimiterTripped = true;
        }
      }

      // prepare for melting/refreezing
      PetscScalar Hmeltnew = Hmelt[i][j];

      // drain ice segments
      for (PetscInt k=1; k < ks; k++) {
        // modifies last two arguments, generally:
        ierr = drainageToBaseModelEnth(EC, L, omega_max, H[i][j], fzlev[k], fdz,
                                       Enthnew[k], Hmeltnew); CHKERRQ(ierr);
      }

      // drain ice/rock interface (or base of ice shelf) segment
      if (ks > 0) {
        if (PismModMask(mask[i][j]) == MASK_FLOATING) {
          // if the ice is floating then mass and energy balance
          //   is the responsibility of the PISMOceanCoupler
          Enthnew[0] = Enth_shelfbase;
          Hmeltnew = 0.0;
        } else {
          // use the drainage model for the basal ice segment to determine both
          //   Hmeltnew and Enthnew[0]; modifies last two arguments, generally:
          ierr = drainageToBaseModelEnth(EC, L, omega_max, H[i][j], 0.0, fdz,
                                         Enthnew[0], Hmeltnew); CHKERRQ(ierr);
        }
      } else {
        Hmeltnew = 0.0; // no stored water if no ice present
        Enthnew[0] = Enth_air;
      }

      // set enthalpy to energy content of surface ice; this is a regularity
      //   issue not an atmosphere model!
      for (PetscInt k=ks+1; k<fMz; k++) {
        Enthnew[k] = Enth_air;
      }

      // transfer column into EnthNew3; communication later
      ierr = EnthNew3.setValColumnPL(i,j,fMz,fzlev,Enthnew); CHKERRQ(ierr);

      // recover temperature in bedrock from enthalpy solution
      for (PetscInt k=0; k < k0; k++) {
        Tbnew[k] = EC.getAbsTempBedrock(x[k]);
      }

      // bottom of ice is top of bedrock when grounded, so
      //   bedrock value T(z=0) should match ice value T(z=0);
      //   when floating just match ocean temp provided by PISMOceanCoupler
      if (PismModMask(mask[i][j]) == MASK_FLOATING) { // top of bedrock sees ocean
          Tbnew[k0] = Tshelfbase[i][j];
      } else {
        if (ks > 0) { // grounded ice present
          // get ice temperature at z=0; enforces continuity of temperature
          ierr = EC.getAbsTemp(Enthnew[0], p_basal, Tbnew[k0]); CHKERRQ(ierr);
        } else {      // no significant ice; top of bedrock sees atmosphere
          Tbnew[k0] = Ts[i][j];
        }
      }

      // in temperate base case, or if bulge limiter was tripped, 
      //    redo bedrock temperature solution
      // FIXME:  and put energy back into melting
      if (fMbz > 1) {
        if (    (bulgeLimiterTripped)
             || (    (PismModMask(mask[i][j]) != MASK_FLOATING)
                  && (EC.isTemperate(Enthnew[0], p_basal )) ) ) {
          for (PetscInt k=0; k < fMbz; k++) {
            bedredosystem.T_b[k] = Tb[k];
          }
          PetscScalar T_basal_new;
          ierr = EC.getAbsTemp(Enthnew[0], p_basal, T_basal_new); CHKERRQ(ierr);
          ierr = bedredosystem.setTopBoundaryValueThisColumn(T_basal_new); CHKERRQ(ierr);
          ierr = bedredosystem.setBasalBoundaryValueThisColumn(Ghf[i][j]); CHKERRQ(ierr);
          // solve the system
          ierr = bedredosystem.solveThisColumn(&xredo); // no CHKERRQ(ierr) immediately because:
          if (ierr > 0) {
            SETERRQ3(2,
              "bedredosystem tridiagonal solve failed at (%d,%d)\n"
              "   with zero pivot position %d.\n", i, j, ierr);
          } else { CHKERRQ(ierr); }
          // insert
          for (PetscInt k=0; k < fMbz; k++) {
            Tbnew[k] = xredo[k];
          }
          // FIXME:  need to subtract Tbnew - xredo and put into melting
  
          // diagnostic/debug
          if (viewOneRedoColumn) {
            if ((i==id) && (j==jd)) {
              ierr = verbPrintf(1,PETSC_COMM_SELF,
                "\n\nviewing bedredosystem and solution at (i,j)=(%d,%d):\n", i, j);
                CHKERRQ(ierr);
              ierr = bedredosystem.viewConstants(NULL); CHKERRQ(ierr);
              ierr = bedredosystem.viewSystem(NULL,"system"); CHKERRQ(ierr);
              ierr = bedredosystem.viewColumnValues(NULL, xredo, fMbz, "solution xredo");
                CHKERRQ(ierr);
            }
          }
        }  
      }  // if (fMbz > 1)

      // transfer column into Tb3; no need for communication, even later
      ierr = Tb3.setValColumnPL(i,j,fMbz,fzblev,Tbnew); CHKERRQ(ierr);

      // basalMeltRate[][] is rate of mass loss from bottom of ice
      if (PismModMask(mask[i][j]) == MASK_FLOATING) {
        // rate of mass loss at bottom of ice shelf;  can be negative (marine freeze-on)
        basalMeltRate[i][j] = bmr_float[i][j]; // set by PISMOceanCoupler
      } else {
        // rate of change of Hmelt[][];  can be negative (till water freeze-on)
        basalMeltRate[i][j] = (Hmeltnew - Hmelt[i][j]) / dtTempAge;
      }

      // finalize Hmelt value
      if (PismModMask(mask[i][j]) == MASK_FLOATING) {
        // if floating assume maximally saturated till
        // UNACCOUNTED MASS & ENERGY (LATENT) LOSS/GAIN (TO/FROM OCEAN)!!
        Hmelt[i][j] = max_hmelt;
      } else {
        // limit Hmelt to be in [0.0, max_hmelt]
        // UNACCOUNTED MASS & ENERGY (LATENT) LOSS (TO INFINITY AND BEYOND)!!
        Hmelt[i][j] = PetscMax(0.0, PetscMin(max_hmelt, Hmeltnew) );
      }

    }
  }

  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vHmelt.end_access(); CHKERRQ(ierr);
  ierr = vRb.end_access(); CHKERRQ(ierr);
  ierr = vGhf.end_access(); CHKERRQ(ierr);
  ierr = vbasalMeltRate.end_access(); CHKERRQ(ierr);

  ierr = pccTs->end_access(); CHKERRQ(ierr);
  ierr = pccsbt->end_access();  CHKERRQ(ierr);
  ierr = pccsbmf->end_access();  CHKERRQ(ierr);

  ierr = Tb3.end_access(); CHKERRQ(ierr);
  ierr = u3.end_access(); CHKERRQ(ierr);
  ierr = v3.end_access(); CHKERRQ(ierr);
  ierr = w3.end_access(); CHKERRQ(ierr);
  ierr = Sigma3.end_access(); CHKERRQ(ierr);
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = EnthNew3.end_access(); CHKERRQ(ierr);


  delete [] system.u;     delete [] system.v;     delete [] system.w;
  delete [] system.Sigma; delete [] system.Enth;  delete [] system.Enth_s; 
  delete [] system.Enth_b;

  delete [] bedredosystem.T_b;
  delete [] xredo;

  delete [] x;
  delete [] Tb;     delete [] Tbnew;   delete [] Enthnew;
  delete [] fzlev;  delete [] fzblev;

  return 0;
}


//! Move some of the liquid water fraction in a column segment [z,z+dz] to the base according to heuristics.
/*!
Generally this procedure alters the values in the last two arguments,
'enthalpy' and 'Hmelt'.  There is conservation of energy here because the 
enthalpy lost/gained by the ice becomes latent heat lost/gained by the 
layer of basal water.

Heuristic: Once liquid water fraction exceeds a cap, all of it goes to the base.
Follows \ref Greve97Greenland and references therein.

Interesting note:  If the basal ice is cold and there is available water
(Hmelt > 0.0) then ice will freeze on, causing a negative basal melt rate which can 
enter into the mass continuity equation, and we bring the lowest ice layer (basal ice)
up to temperate.
 */
PetscErrorCode IceEnthalpyModel::drainageToBaseModelEnth(EnthalpyConverter &EC,
                PetscScalar L, PetscScalar omega_max,
                PetscScalar thickness, PetscScalar z, PetscScalar dz,
                PetscScalar &enthalpy, PetscScalar &Hmelt) {
  PetscErrorCode ierr;

  if (allowAboveMelting == PETSC_TRUE) {
    SETERRQ(1,"IceEnthalpyModel::drainageToBaseModelEnth() called but allowAboveMelting==TRUE");
  }

  // no change to either enthalpy or basal layer thickness in this case
  if (updateHmelt == PETSC_FALSE)  return 0;

  const PetscScalar p     = EC.getPressureFromDepth(thickness - z);

  // if there is liquid water already, thus temperate, consider whether there
  //   is enough to cause drainage
  const PetscScalar omega = EC.getWaterFractionLimited(enthalpy, p);
  if (omega > 1.0e-6) {
    const PetscScalar abovecap = omega - omega_max;
    if (abovecap > 0.0) {
      // alternative, but different if E > E_l: enthalpy -= abovecap * L;
      //   following form is safer, but gives UNACCOUNTED ENERGY LOSS IF E>E_l
      ierr = EC.getEnthAtWaterFraction(omega_max, p, enthalpy); CHKERRQ(ierr);
      Hmelt    += abovecap * dz;   // ice-equivalent water thickness change
    }
    return 0; // done with temperate case
  }
  
  // if cold and in the basal layer, consider whether there is available 
  //   water to freeze on
  // FIXME: think about ocean case!?
  if ((z >= -1.0e-6) && (z <= 1.0e-6)) {
    // only consider freeze-on if column segment is at base of ice;
    //   E_s = getEnthalpyCTS(config, p)
    const PetscScalar dEnth_to_reach_temperate = EC.getEnthalpyCTS(p) - enthalpy;
    if (dEnth_to_reach_temperate > 0.0) {
      // if below E_s, then freeze on, and bring up enthalpy to E_s if
      //   enough water is available; first quantity is also
      //   ((rho Hmelt dx dy) * L) / (rho dx dy dz)
      const PetscScalar
         dEnth_available = (Hmelt / dz) * L,
         dEnth_added     = PetscMin(dEnth_available, dEnth_to_reach_temperate);
      enthalpy += dEnth_added;
      Hmelt    -= (dEnth_added * dz) / L;  // will cause negative basal melt rate;
                                           //   enters into mass continuity
    }
  }

  return 0;
}


//! Compute effective pressure on till using thickness of stored till water and basal melt rate.
/*!
The effective pressure on the till is 
    \f[   N = \rho g H - p_w   \f]
where \f$\rho g H\f$ is the ice over-burden pressure (in the shallow approximation),
and  \f$p_w\f$ is the modeled pore water pressure.

This procedure provides a simple model of pore water pressure \f$p_w\f$.  It is
a function of the thickness of the basal stored water plus (optionally)
of the basal melt rate.

Input \c bwat is thickness of basal water.  Input \c bmr is the basal melt rate.
Because both \c bwat and \c bmr are zero at points where base of ice is frozen,
the resulting effective pressure on the till, from this routine, equals
the overburden pressure.

Option \c -plastic_pwfrac controls parameter till_pw_fraction.

Option \c -bmr_enhance turns on the basal melt rate dependency in pore pressure.
Option \c -bmr_enhance_scale sets the value.  Default is equivalent to
\c -bmr_enhance_scale 0.10, for 10 cm/a as a significant enough level of
basal melt rate to cause a big weakening effect.

Need \f$0 \le\f$ \c bwat \f$\le\f$ \c max_hmelt before calling this.  There is
no error checking.

Compare the porewater pressure computed by formula (4) in [\ref Ritzetal2001],
where the pressure is a function of sea level and bed elevation.
Also, the method using "elevation of the bed at the grounding line"
as in Lingle&Brown 1987 is not implementable because that elevation is
at an unknowable location.  We are not doing a flow line model!
 */
PetscScalar IceEnthalpyModel::getEffectivePressureOnTill(PetscScalar thk, PetscScalar bwat,
							 PetscScalar bmr, PetscScalar frac,
							 PetscScalar max_hmelt) const {

  const PetscScalar
    p_overburden = ice->rho * standard_gravity * thk;

  // base model for pore water pressure;  note  0 <= p_pw <= frac * p_overburden
  //   because  0 <= bwat <= max_hmelt
  PetscScalar  p_pw = frac * (bwat / max_hmelt) * p_overburden;

  if (bmr_in_pore_pressure) {
    // more weakening from instantaneous basal melt rate; additional
    //   pore water pressure means reduction of effective pressure and of tau_c
    //   note  (additional) <= (1.0 - frac) * p_overburden so  0 <= p_pw <= p_overburden
    p_pw += ( 1.0 - exp( - PetscMax(0.0,bmr) / bmr_enhance_scale ) ) * (1.0 - frac) * p_overburden;
  }

  if (thk_affects_pore_pressure) {
    // ice thickness is surrogate for distance to margin; near margin the till
    //   is presumably better drained so we reduce the pore pressure
    const PetscScalar re    = margin_pore_pressure_reduced,
                      Hhigh = margin_pore_pressure_H_high,
                      Hlow  = margin_pore_pressure_H_low;
    if (thk < Hhigh) {
      if (thk <= Hlow) {
        p_pw *= re;
      } else { // case Hlow < thk < Hhigh;
               //   use linear to connect (Hlow, reduced * p_pw) to (Hhigh, 1.0 * p_w)
        p_pw *= re + (1.0 - re) * (thk - Hlow) / (Hhigh - Hlow);
      }
    }
  }

  return p_overburden - p_pw;
}


//! Update the till yield stress for the pseudo-plastic SSA model.
/*!
Updates based on stored till water and basal melt rate.  We implement
formula (2.4) in [\ref SchoofStream],
    \f[   \tau_c = \mu (\rho g H - p_w), \f]
where \f$\tau_c\f$ is the till yield stress, \f$\rho g H\f$ is the ice over-burden
pressure (in the shallow approximation), \f$p_w\f$ is the modeled
pore water pressure, and \f$\mu\f$ is a strength coefficient for the mineral till
(at least, it is independent of \f$p_w\f$).  The difference
    \f[   N = \rho g H - p_w   \f]
is the effective pressure on the till.
 
We modify Schoof's formula by (possibly) adding a small till cohesion \f$c_0\f$
and by expressing the coefficient as the tangent of a till friction angle

\f$\varphi\f$:
    \f[   \tau_c = c_0 + (\tan \varphi) N. \f]
See [\ref Paterson] table 8.1) regarding values of \f$c_0\f$.
Option  \c -plastic_c0 controls it.

The main issue is the model for pore water pressure \f$p_w\f$ when
computing \f$N\f$.  See getEffectivePressureOnTill().

See [\ref BBssasliding] for a discussion of a complete model using these tools.

Note that IceModel::updateSurfaceElevationAndMask() also
checks whether do_plastic_till is true and if so it sets all mask points to
DRAGGING.

FIXME: Should be renamed "updateYieldStressUsingBasalWater()"?
 */
PetscErrorCode IceEnthalpyModel::updateYieldStressFromHmelt() {
  PetscErrorCode  ierr;

  bool do_plastic_till = config.get_flag("do_plastic_till");
  // only makes sense when do_plastic_till == TRUE
  if (do_plastic_till == PETSC_FALSE) {
    SETERRQ(1,"do_plastic_till == PETSC_FALSE but updateYieldStressFromHmelt() called");
  }

  if (holdTillYieldStress == PETSC_FALSE) { // usual case: use Hmelt to determine tauc
    PetscScalar **mask, **tauc, **H, **Hmelt, **bmr, **tillphi; 

    PetscScalar till_pw_fraction = config.get("till_pw_fraction"),
      till_c_0 = config.get("till_c_0") * 1e3, // convert from kPa to Pa
      till_mu = tan((pi/180.0)*config.get("default_till_phi")),
      max_hmelt = config.get("max_hmelt");

    ierr =          vMask.get_array(mask);    CHKERRQ(ierr);
    ierr =          vtauc.get_array(tauc);    CHKERRQ(ierr);
    ierr =             vH.get_array(H);       CHKERRQ(ierr);
    ierr =         vHmelt.get_array(Hmelt);   CHKERRQ(ierr);
    ierr = vbasalMeltRate.get_array(bmr);     CHKERRQ(ierr);
    ierr =       vtillphi.get_array(tillphi); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        if (PismModMask(mask[i][j]) == MASK_FLOATING) {
          tauc[i][j] = 0.0;  
        } else if (H[i][j] == 0.0) {
          tauc[i][j] = 1000.0e3;  // large yield stress of 1000 kPa = 10 bar if no ice
        } else { // grounded and there is some ice
          const PetscScalar
            N = getEffectivePressureOnTill(H[i][j], Hmelt[i][j], bmr[i][j], till_pw_fraction,
					   max_hmelt);
          if (useConstantTillPhi == PETSC_TRUE) {
            tauc[i][j] = till_c_0 + N * till_mu;
          } else {
            tauc[i][j] = till_c_0 + N * tan((pi/180.0) * tillphi[i][j]);
          }
        }
      }
    }
    ierr =          vMask.end_access(); CHKERRQ(ierr);
    ierr =          vtauc.end_access(); CHKERRQ(ierr);
    ierr =             vH.end_access(); CHKERRQ(ierr);
    ierr =       vtillphi.end_access(); CHKERRQ(ierr);
    ierr = vbasalMeltRate.end_access(); CHKERRQ(ierr);
    ierr =         vHmelt.end_access(); CHKERRQ(ierr);
  }

  return 0;
}


