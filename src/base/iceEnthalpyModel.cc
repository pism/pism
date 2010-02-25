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

#include "iceEnthalpyModel.hh"
#include "bedrockOnlySystem.hh"
#include "iceenthOnlySystem.hh"

#include "enthColumnSystem.hh"

#include "enthalpyConverter.hh"


/*********** procedures for init ****************/

IceEnthalpyModel::IceEnthalpyModel(IceGrid &g, NCConfigVariable &conf,
     NCConfigVariable &conf_overrides) : IceModel(g, conf, conf_overrides) {

  doColdIceMethods = PETSC_FALSE; // default is to USE enthalpy
}


PetscErrorCode IceEnthalpyModel::setFromOptions() {
  PetscErrorCode ierr;

  ierr = IceModel::setFromOptions(); CHKERRQ(ierr);

  // if set, use old IceModel::temperatureStep(), and set enthalpy as though
  //   ice is cold
  ierr = check_option("-cold", doColdIceMethods); CHKERRQ(ierr);

  // DEBUG:  report settings
  const int vlevel = 2;
  ierr = verbPrintf(vlevel, grid.com,
      "  IceEnthalpyModel::setFromOptions():\n"); CHKERRQ(ierr);
  ierr = verbPrintf(vlevel, grid.com,
      "    doColdIceMethods is %s\n", 
      (doColdIceMethods == PETSC_TRUE) ? "TRUE" : "FALSE");
      CHKERRQ(ierr);
  ierr = verbPrintf(vlevel, grid.com,
      "    config:bmr_enhance_basal_water_pressure is %s\n",
      config.get_flag("bmr_enhance_basal_water_pressure") ? "TRUE" : "FALSE");
      CHKERRQ(ierr);
  ierr = verbPrintf(vlevel, grid.com,
      "    config:thk_eff_basal_water_pressure is %s\n",
      config.get_flag("thk_eff_basal_water_pressure") ? "TRUE" : "FALSE");
      CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceEnthalpyModel::createVecs() {
  PetscErrorCode ierr;

  ierr = IceModel::createVecs(); CHKERRQ(ierr);

  ierr = Enth3.create(grid, "enthalpy", true); CHKERRQ(ierr);
  // POSSIBLE standard name = land_ice_enthalpy
  ierr = Enth3.set_attrs(
     "model_state",
     "ice enthalpy (includes sensible heat, latent heat, pressure)",
     "J kg-1",
     ""); CHKERRQ(ierr);
  ierr = variables.add(Enth3); CHKERRQ(ierr);

  // see IceModel::allocate_internal_objects(), which is where this should go
  ierr = EnthNew3.create(grid,"enthalpy_new",false); CHKERRQ(ierr); // global
  ierr = EnthNew3.set_attrs(
     "internal",
     "ice enthalpy; temporary space during timestep",
     "J kg-1",
     ""); CHKERRQ(ierr);
  // work vectors with "internal" are not in PISMVars IceModel::variables

  return 0;
}


PetscErrorCode IceEnthalpyModel::init_physics() {
  PetscErrorCode ierr;

  // let the base class create the ice and process its options:
  ierr = IceModel::init_physics(); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
      "  setting flow law to Glen-Paterson-Budd-Lliboutry-Duval type ...\n");
      CHKERRQ(ierr);
  if (ice != NULL)  delete ice;  // kill choice already made
  iceFactory.setType(ICE_GPBLD); // new flowlaw which has dependence on enthalpy
                                 //   not temperature
  iceFactory.create(&ice);
  PolyThermalGPBLDIce *gpbldi = dynamic_cast<PolyThermalGPBLDIce*>(ice);
  if (gpbldi == NULL) {
    ThermoGlenIce *tgi = dynamic_cast<ThermoGlenIce*>(ice);
    if (tgi) {
      ierr = verbPrintf(2, grid.com,
        "    [flow law was actually set to ThermoGlenIce by IceEnthalpyModel]\n");
        CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(1, grid.com,
        "PISM WARNING: flow law unclear in IceEnthalpyModel ...\n"); CHKERRQ(ierr);
    }
  }
  ierr = ice->setFromOptions();CHKERRQ(ierr);
  ierr = ice->printInfo(4);CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceEnthalpyModel::initFromFile(const char *filename) {
  PetscErrorCode ierr;

  // option to override default behavior: if -init_from_temp or
  //   -init_from_temp_and_liqfrac are set then we *MODIFY* the input file
  bool initfromT, initfromTandOm;
  ierr = check_option("-init_from_temp", initfromT); CHKERRQ(ierr);
  ierr = check_option("-init_from_temp_and_liqfrac", initfromTandOm); CHKERRQ(ierr);
  if ((!initfromT) && (!initfromTandOm)) {
    // if neither option is set, proceed with normal initialization
    ierr = IceModel::initFromFile(filename); CHKERRQ(ierr);
    return 0;
  }

  // -init_from_temp was set, calculate enthalpy from temperature (i.e. cold-ice only)
  if (initfromT) {
    ierr = verbPrintf(2, grid.com,
      "  option -init_from_temp seen ... IceEnthalpyModel doing SPECIAL ACTIONS:\n"
      "      reading ice temperature and thickness from %s ...\n",
      filename); CHKERRQ(ierr);  
    NCTool nc(&grid);
    ierr = nc.open_for_reading(filename); CHKERRQ(ierr);
    int last_record;
    ierr = nc.get_dim_length("t", &last_record); CHKERRQ(ierr);
    last_record -= 1;
    // computation of Enth3 will need thickness
    ierr = vH.read(filename,last_record); CHKERRQ(ierr);
    ierr = T3.read(filename,last_record); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);
  
    ierr = verbPrintf(2, grid.com,
      "      computing enthalpy from ice temperature and thickness ...\n");
    CHKERRQ(ierr);  
    ierr = setEnth3FromT3_ColdIce(); CHKERRQ(ierr);

    ierr = verbPrintf(2, grid.com,
      "      MODIFYING file '%s' by writing enthalpy ...\n",filename);
    CHKERRQ(ierr);  
    ierr = nc.open_for_writing(filename, true, false); CHKERRQ(ierr);
    ierr = Enth3.write(filename,NC_DOUBLE); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);

  } else if (initfromTandOm) {
    // -init_from_and_liqfrac was set
    IceModelVec3 Liqfrac3;
    ierr = Liqfrac3.create(grid, "liqfrac", false); CHKERRQ(ierr);

    ierr = verbPrintf(2, grid.com,
      "  option -init_from_temp_and_liqfrac seen ... IceEnthalpyModel doing SPECIAL ACTIONS:\n"
      "      reading ice temperature, liquid water fraction and thickness from %s ...\n",
      filename); CHKERRQ(ierr);  
    NCTool nc(&grid);
    ierr = nc.open_for_reading(filename); CHKERRQ(ierr);
    int last_record;
    ierr = nc.get_dim_length("t", &last_record); CHKERRQ(ierr);
    last_record -= 1;
    // computation of Enth3 will need thickness, temperature and liquid water fraction
    ierr = vH.read(filename,last_record); CHKERRQ(ierr);
    ierr = T3.read(filename,last_record); CHKERRQ(ierr);
    ierr = Liqfrac3.read(filename,last_record); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);

    ierr = verbPrintf(2, grid.com,
      "      computing enthalpy from ice temperature and thickness ...\n");
    CHKERRQ(ierr);  
    ierr = setEnth3FromT3AndLiqfrac3(Liqfrac3); CHKERRQ(ierr);

    ierr = verbPrintf(2, grid.com,
      "      MODIFYING file '%s' by writing enthalpy ...\n",filename);
    CHKERRQ(ierr);  
    ierr = nc.open_for_writing(filename, true, false); CHKERRQ(ierr);
    ierr = Enth3.write(filename,NC_DOUBLE); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);
  }

  // now init; this will actually re-read enthalpy because of
  //   "variables.add(Enth3)" above
  ierr = IceModel::initFromFile(filename); CHKERRQ(ierr);
  
  return 0;
}


PetscErrorCode IceEnthalpyModel::bootstrapFromFile(const char *filename) {
  PetscErrorCode ierr;

  ierr = IceModel::bootstrapFromFile(filename); CHKERRQ(ierr);
  
  ierr = verbPrintf(2, grid.com, 
    "continuing bootstrapping in IceEnthalpyModel ...\n"); CHKERRQ(ierr);
  ierr = verbPrintf(2, grid.com,
    "  ice enthalpy set from temperature, as cold ice (zero liquid fraction)\n");
    CHKERRQ(ierr);
  ierr = setEnth3FromT3_ColdIce(); CHKERRQ(ierr);
  ierr = verbPrintf(2, grid.com, 
    "bootstrapping done (IceEnthalpyModel)\n"); CHKERRQ(ierr);

  return 0;
}


/*********** procedures for read/write ****************/

PetscErrorCode IceEnthalpyModel::write_extra_fields(const char* filename) {
  PetscErrorCode ierr;

  ierr = IceModel::write_extra_fields(filename); CHKERRQ(ierr);
  
  ierr = Enth3.write(filename, NC_DOUBLE); CHKERRQ(ierr);

  // also write omega = liquid water fraction
  //   we use EnthNew3 (global) as temporary, allocated space for this purpose
  ierr = verbPrintf(4, grid.com,
      "  writing liquid water fraction 'liquid_frac' from enthalpy ...\n"); CHKERRQ(ierr);
  ierr = setLiquidFracFromEnthalpy(EnthNew3); CHKERRQ(ierr);
  ierr = EnthNew3.write(filename, NC_FLOAT); CHKERRQ(ierr);

  // also write temp_pa = pressure-adjusted temp in Celcius
  //   again use EnthNew3 (global) as temporary, allocated space
  ierr = verbPrintf(4, grid.com,
      "  writing pressure-adjusted ice temperature (deg C) 'temp_pa' ...\n");
      CHKERRQ(ierr);
  ierr = setPATempFromEnthalpy(EnthNew3); CHKERRQ(ierr); // sets to Kelvin
  ierr = EnthNew3.shift(- config.get("water_melting_temperature"));// make deg C
      CHKERRQ(ierr);
  ierr = EnthNew3.write(filename, NC_FLOAT); CHKERRQ(ierr);

  // write CTS position (unitless) if command line option -cts is given
  //   again use EnthNew3 (global) as temporary, allocated space
  bool userWantsCTS;
  ierr = check_option("-cts", userWantsCTS); CHKERRQ(ierr);
  if (userWantsCTS) {
    ierr = verbPrintf(4, grid.com,
      "  writing CTS scalar field 'cts'  (= E/Es) ...\n"); CHKERRQ(ierr);
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
        // because of how getPressureFromDepth() works, the energy content in
        //   the air is set to the value ice would have if it a chunk of it
        //   occupied the air; the atmosphere actually has much lower energy
        //   content; done this way for regularity (i.e. dEnth/dz computations)
        ierr = EC.getEnthPermissive(Tij[k],0.0,EC.getPressureFromDepth(depth),
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
PetscErrorCode IceEnthalpyModel::setEnth3FromT3AndLiqfrac3(
                                          IceModelVec3 &Liqfrac3) {
  PetscErrorCode ierr;
  EnthalpyConverter EC(config);
  
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
        ierr = EC.getEnthPermissive(Tij[k],Liqfracij[k],
                      EC.getPressureFromDepth(depth), Enthij[k]); CHKERRQ(ierr);
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


//! Compute the ice temperature corresponding to Enth3, and put in Tnew3; use just after Enth3 is determined.
/*!
Does not communicate.  Ghosts will be invalid, but 
"T3.endGhostCommTransfer(Tnew3)" in IceModel::temperatureAgeStep() will have
the desired effect.
 */
PetscErrorCode IceEnthalpyModel::setTnew3FromEnth3() {
  PetscErrorCode ierr;
  EnthalpyConverter EC(config);

  PetscScalar *Tij, *Enthij; // columns of these values
  ierr = Tnew3.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = Tnew3.getInternalColumn(i,j,&Tij); CHKERRQ(ierr);
      ierr = Enth3.getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
      for (PetscInt k=0; k<grid.Mz; ++k) {
        const PetscScalar depth = vH(i,j) - grid.zlevels[k];
        ierr = EC.getAbsTemp(Enthij[k],EC.getPressureFromDepth(depth), Tij[k]); 
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
  ierr = Tnew3.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  return 0;
}


//! Compute the liquid fraction corresponding to Enth3, and put in a global IceModelVec3 provided by user.
PetscErrorCode IceEnthalpyModel::setLiquidFracFromEnthalpy(IceModelVec3 &useForLiquidFrac) {
  PetscErrorCode ierr;
  EnthalpyConverter EC(config);

  ierr = useForLiquidFrac.set_name("liqfrac"); CHKERRQ(ierr);
  ierr = useForLiquidFrac.set_attrs(
     "diagnostic",
     "liquid water fraction in ice (between 0 and 1)",
     "",
     ""); CHKERRQ(ierr);

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
        ierr = EC.getWaterFraction(Enthij[k],EC.getPressureFromDepth(depth),
                                   omegaij[k]); CHKERRQ(ierr);
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
  EnthalpyConverter EC(config);

  ierr = useForPATemp.set_name("temp_pa"); CHKERRQ(ierr);
  ierr = useForPATemp.set_attrs(
     "diagnostic",
     "pressure-adjusted ice temperature (degrees above pressure-melting point)",
     "deg_C",
     ""); CHKERRQ(ierr);

  PetscScalar *Tpaij, *Enthij; // columns of these values
  ierr = useForPATemp.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = useForPATemp.getInternalColumn(i,j,&Tpaij); CHKERRQ(ierr);
      ierr = Enth3.getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
      for (PetscInt k=0; k<grid.Mz; ++k) {
        const PetscScalar depth = vH(i,j) - grid.zlevels[k];
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
  EnthalpyConverter EC(config);

  ierr = useForCTS.set_name("cts"); CHKERRQ(ierr);
  ierr = useForCTS.set_attrs(
     "diagnostic",
     "cts = E/E_s(p), so cold-temperate transition surface is at cts = 1",
     "",
     ""); CHKERRQ(ierr);

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
        CTSij[k] = EC.getCTS(Enthij[k], EC.getPressureFromDepth(depth));
      }
    }
  }
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = useForCTS.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceEnthalpyModel::energyStats(PetscScalar iarea, bool /*useHomoTemp*/, 
                                             PetscScalar &gmeltfrac, PetscScalar &gtemp0) {
  PetscErrorCode  ierr;
  EnthalpyConverter EC(config);

  ierr = vH.begin_access(); CHKERRQ(ierr);
  // put basal ice enthalpy in vWork2d[0]
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.getHorSlice(vWork2d[0], 0.0); CHKERRQ(ierr);  // z=0 slice
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  PetscScalar **Enthbase;
  ierr = vWork2d[0].get_array(Enthbase); CHKERRQ(ierr);

  const PetscScalar a = grid.dx * grid.dy * 1e-3 * 1e-3; // area unit (km^2)
  PetscScalar  meltarea = 0.0, temp0 = 0.0;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vH(i,j) > 0) {
        // accumulate area of base which is at melt point
        if (EC.isTemperate(Enthbase[i][j], EC.getPressureFromDepth(vH(i,j)) ))  
          meltarea += a;
      }
      // if you happen to be at center, record absolute basal temp there
      if (i == (grid.Mx - 1)/2 && j == (grid.My - 1)/2) {
        ierr = EC.getAbsTemp(Enthbase[i][j],EC.getPressureFromDepth(vH(i,j)), temp0);
          CHKERRQ(ierr);
      }
    }
  }
  
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vWork2d[0].end_access(); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&meltarea, &gmeltfrac, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&temp0,    &gtemp0,    grid.com); CHKERRQ(ierr);

  // normalize fraction correctly
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
      "    [IceEnthalpyModel::temperatureStep(): ENTHALPY IS OFF.\n"
      "         CALLING IceModel::temperatureStep()]\n");
    CHKERRQ(ierr);

    ierr = IceModel::temperatureStep(vertSacrCount,bulgeCount);  CHKERRQ(ierr);

    // start & complete communication: THIS IS REDUNDANT WITH temperatureAgeStep(),
    //    BUT NEEDED TO GET UPDATED ENTHALPY AT END OF TIME STEP
    ierr = T3.beginGhostCommTransfer(Tnew3); CHKERRQ(ierr);
    ierr = T3.endGhostCommTransfer(Tnew3); CHKERRQ(ierr);

    ierr = setEnth3FromT3_ColdIce();  CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(4,grid.com,
      "    [IceEnthalpyModel::temperatureStep(): ENTHALPY IS ON.\n"
      "         CALLING IceEnthalpyModel::enthalpyAndDrainageStep()]\n");
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


//! Update enthalpy field based on conservation of energy in ice and bedrock.
PetscErrorCode IceEnthalpyModel::enthalpyAndDrainageStep(
                      PetscScalar* vertSacrCount, PetscScalar* /*bulgeCount*/) {
  PetscErrorCode  ierr;
  EnthalpyConverter EC(config);

  if (doColdIceMethods==PETSC_TRUE) {
    SETERRQ(1,
       "\n\nPISM ERROR:  IceEnthalpyModel::enthalpyAndDrainageStep() called but\n"
           "             doColdIceMethods==true ... ending\n");
  }

  PetscInt    fMz, fMbz;  // number of fine grid levels in ice and bedrock, resp
  PetscScalar fdz, *fzlev, fdzb, *fzblev;
  ierr = grid.get_fine_vertical_grid(fMz, fMbz, fdz, fdzb, fzlev, fzblev); CHKERRQ(ierr);

  bedrockOnlySystemCtx bosys(config, fMbz);
  ierr = bosys.initAllColumns(dtTempAge, fdzb); CHKERRQ(ierr);
  ierr = bosys.viewConstants(NULL, false); CHKERRQ(ierr);

  iceenthOnlySystemCtx ieosys(config, Enth3, fMz);
  ierr = ieosys.initAllColumns(grid.dx, grid.dy, dtTempAge, fdz); CHKERRQ(ierr);
  ierr = ieosys.viewConstants(NULL, false); CHKERRQ(ierr);

  enthSystemCtx system(fMz,fMbz);
  system.dx      = grid.dx;
  system.dy      = grid.dy;
  system.dtTemp  = dtTempAge; // same time step for temp and age, currently
  system.dzEQ    = fdz;
  system.dzbEQ   = fdzb;
  system.ice_rho = config.get("ice_density");
  system.ice_c   = config.get("ice_specific_heat_capacity");
  system.ice_k   = config.get("ice_thermal_conductivity");
  system.ice_nu  = config.get("enthalpy_temperate_diffusivity"); // diffusion in temperate ice
  system.bed_rho = config.get("bedrock_thermal_density");
  system.bed_c   = config.get("bedrock_thermal_specific_heat_capacity");
  system.bed_k   = config.get("bedrock_thermal_conductivity");

  // space for solution of system
  PetscScalar *x;
  x = new PetscScalar[fMbz + 1 + fMz];

  const PetscScalar
    p_air     = config.get("surface_pressure"),          // Pa
    L         = config.get("water_latent_heat_fusion"),  // J kg-1
    omega_max = config.get("liquid_water_fraction_max"), // pure
    max_hmelt = config.get("max_hmelt");

  // pointers to values in current column
  PetscScalar *Enthnew, *Tbnew;

  system.u     = new PetscScalar[fMz];
  system.v     = new PetscScalar[fMz];
  system.w     = new PetscScalar[fMz];
  system.Sigma = new PetscScalar[fMz];
  system.Enth_s= new PetscScalar[fMz];  // enthalpy of pressure-melting-point

  system.Enth  = new PetscScalar[fMz];
  system.Tb    = new PetscScalar[fMbz];

  Enthnew      = new PetscScalar[fMz];
  Tbnew        = new PetscScalar[fMbz];

  // system needs access to Enth3 for planeStar()
  system.Enth3 = &Enth3;

  // checks that all needed constants and pointers got set:
  ierr = system.initAllColumns(); CHKERRQ(ierr);

  bool viewOneColumn;
  ierr = check_option("-view_sys", viewOneColumn); CHKERRQ(ierr);

  if (getVerbosityLevel() >= 4) {  // view: all constants correct at this point?
    ierr = EC.viewConstants(NULL); CHKERRQ(ierr);
    ierr = system.viewConstants(NULL,false); CHKERRQ(ierr);
  }

  // now get map-plane coupler fields: Dirichlet upper surface boundary and
  //    mass balance lower boundary
  if (surface != PETSC_NULL) {
    ierr = surface->ice_surface_temperature(grid.year, dtTempAge / secpera, artm); CHKERRQ(ierr);
  } else {
    SETERRQ(4,"PISM ERROR: surface == PETSC_NULL");
  }

  if (ocean != PETSC_NULL) {
    ierr = ocean->shelf_base_mass_flux(grid.year, dt / secpera, shelfbmassflux); CHKERRQ(ierr);
  } else {
    SETERRQ(4,"PISM ERROR: ocean == PETSC_NULL");
  }

  ierr = artm.begin_access(); CHKERRQ(ierr);
  ierr = shelfbmassflux.begin_access(); CHKERRQ(ierr);

  // get other map-plane fields
  PetscScalar  **H, **bmr;
  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = vHmelt.begin_access(); CHKERRQ(ierr);
  ierr = vbasalMeltRate.get_array(bmr); CHKERRQ(ierr);
  ierr = vRb.begin_access(); CHKERRQ(ierr);
  ierr = vGhf.begin_access(); CHKERRQ(ierr);

  ierr = vMask.begin_access(); CHKERRQ(ierr);

  // these are accessed a column at a time
  ierr = u3.begin_access(); CHKERRQ(ierr);
  ierr = v3.begin_access(); CHKERRQ(ierr);
  ierr = w3.begin_access(); CHKERRQ(ierr);
  ierr = Sigma3.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = EnthNew3.begin_access(); CHKERRQ(ierr);
  ierr = Tb3.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      // for fine grid; this should *not* be replaced by call to grid.kBelowHeight()
      const PetscInt  ks = static_cast<PetscInt>(floor(H[i][j]/fdz));

      // enthalpy and pressures at boundaries of ice
      const PetscScalar p_basal = EC.getPressureFromDepth(H[i][j]),
                        p_ks    = EC.getPressureFromDepth(H[i][j] - fzlev[ks]);
      PetscScalar Enth_air, Enth_ks;
      ierr = EC.getEnthPermissive(artm(i,j), 0.0, p_air, Enth_air); CHKERRQ(ierr);
      // in theory we could have a water fraction at k=ks level, but for
      //   now there is no case where we have that:
      ierr = EC.getEnthPermissive(artm(i,j), 0.0, p_ks,  Enth_ks); CHKERRQ(ierr);

      ierr = system.setIndicesAndClearThisColumn(i,j,ks); CHKERRQ(ierr);

      ierr = Tb3.getValColumnPL(i,j,fMbz,fzblev,bosys.Tb); CHKERRQ(ierr);
      ierr = Tb3.getValColumnPL(i,j,fMbz,fzblev,system.Tb); CHKERRQ(ierr);

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

      // at each level in the ice, compute the enthalpy that the CTS would have
      //   there (i.e. the enthalpy for the pressure-melting temperature there)
      system.Enth_s[0] = EC.getEnthalpyCTS(p_basal);
      for (PetscInt k = 1; k <= ks; k++) {
        system.Enth_s[k] = EC.getEnthalpyCTS(EC.getPressureFromDepth(H[i][j]-fzlev[k]));
      }
      for (PetscInt k = ks+1; k < fMz; k++) {
        system.Enth_s[k] = EC.getEnthalpyCTS(p_air);
      }

      ierr = bosys.setBoundaryValuesThisColumn(EC.getMeltingTemp(p_basal), vGhf(i,j)); CHKERRQ(ierr); // FIXME

      // go through column and find appropriate lambda for BOMBPROOF
      PetscScalar lambda = 1.0;  // start with centered implicit for more accuracy
      for (PetscInt k = 0; k <= ks; k++) {
        const PetscScalar 
          kconduct = (system.Enth[k] > system.Enth_s[k]) ? 0.0 : system.ice_k,
          denom = (PetscAbs(system.w[k]) + 0.000001/secpera)
                                     * system.ice_rho * system.ice_c * fdz;
        lambda = PetscMin(lambda, 2.0 * kconduct / denom); // so lambda = 0 in temperate ice
      }
      if (lambda < 1.0)  *vertSacrCount += 1; // count columns with lambda < 1

      // if isMarginal then only do vertical conduction for ice;
      //   will ignore advection and strain heating if isMarginal
      const bool isMarginal = checkThinNeigh(
                                 H[i+1][j],H[i+1][j+1],H[i][j+1],H[i-1][j+1],
                                 H[i-1][j],H[i-1][j-1],H[i][j-1],H[i+1][j-1]);
      ierr = system.setSchemeParamsThisColumn(
               vMask.is_floating(i,j), isMarginal, lambda); CHKERRQ(ierr);

      const PetscScalar
	hf = vMask.is_floating(i,j) ? system.ice_rho * L * shelfbmassflux(i,j) 
                                      : vGhf(i,j);
      ierr = system.setBoundaryValuesThisColumn(Enth_ks, hf, vRb(i,j)); CHKERRQ(ierr);

      // diagnostic/debug
      if (viewOneColumn && issounding(i,j)) {
        ierr = verbPrintf(1,grid.com,
          "just before solving system at (i,j)=(%d,%d):\n", i, j); CHKERRQ(ierr);
        ierr = system.viewConstants(NULL,true); CHKERRQ(ierr);
      }

      // solve the system for this column: x will contain new enthalpy in ice
      //   and new temperature in bedrock
      ierr = system.solveThisColumn(&x); // no CHKERRQ(ierr) immediately because:
      if (ierr > 0) {
        char fname[PETSC_MAX_PATH_LEN];
        snprintf(fname, PETSC_MAX_PATH_LEN, "enthsys_i%d_j%d_zeropivot%d.m",
                 i,j,ierr);
        PetscPrintf(grid.com,
          "\n\ntridiagonal solve in enthalpyAndDrainageStep() failed at (%d,%d)\n"
          "   with zero pivot position %d\n"
          "   viewing system to file %s ... \n",
          i, j, ierr, fname);
        PetscViewer viewer;
        ierr = PetscViewerCreate(grid.com, &viewer);CHKERRQ(ierr);
        ierr = PetscViewerSetType(viewer, PETSC_VIEWER_ASCII);CHKERRQ(ierr);
        ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
        ierr = PetscViewerFileSetName(viewer, fname);CHKERRQ(ierr);
        ierr = system.viewSystem(viewer,"system"); CHKERRQ(ierr);
        ierr = PetscViewerDestroy(viewer); CHKERRQ(ierr);
        PetscPrintf(grid.com, "\n   ENDING ...\n");
        PetscEnd();
      } else { CHKERRQ(ierr); }

      // diagnostic/debug
      if (viewOneColumn && issounding(i,j)) {
        char fname[PETSC_MAX_PATH_LEN];
        snprintf(fname, PETSC_MAX_PATH_LEN, "enthsyssoln_i%d_j%d.m",i,j);
        ierr = verbPrintf(1,grid.com,
            "viewing system and solution at (i,j)=(%d,%d) to file %s\n\n",
            i, j, fname); CHKERRQ(ierr);
        PetscViewer    viewer;
        ierr = PetscViewerCreate(grid.com, &viewer);CHKERRQ(ierr);
        ierr = PetscViewerSetType(viewer, PETSC_VIEWER_ASCII);CHKERRQ(ierr);
        ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
        ierr = PetscViewerFileSetName(viewer, fname);CHKERRQ(ierr);
        ierr = system.viewSystem(viewer,"system"); CHKERRQ(ierr);
        ierr = system.viewColumnValues(viewer, x, fMbz+fMz, "solution x");
            CHKERRQ(ierr);
        ierr = PetscViewerDestroy(viewer); CHKERRQ(ierr);
      }

      // start from bottom and work up column, using solution vector x[]
      //   to determine bedrock temperature, basal melt rate, and ice enthalpy

      // temperature in bedrock now complete
      for (PetscInt k=0; k < fMbz; k++) {
        Tbnew[k] = x[k];
      }
      // transfer column into Tb3; no need for communication, even later
      ierr = Tb3.setValColumnPL(i,j,fMbz,fzblev,Tbnew); CHKERRQ(ierr);

      // determine melt rate if any
      PetscScalar Hmeltnew = vHmelt(i,j);  // prepare for melting/refreezing
      if (ks > 0) {
        // three cases: see tables in doc/doxy/html/bombproofenth.html
        if (vMask.is_floating(i,j)) {
          bmr[i][j] = shelfbmassflux(i,j); 
          // and ignor x[fMbz] which contains only known ocean melt rate
          // ... and no contribute to stored Hmelt, which will be set to zero
        } else if (system.Enth[0] >= system.Enth_s[0]) {
          // grounded melting case; add to stored melt water
          bmr[i][j] = x[fMbz] / (system.ice_rho * L);
          Hmeltnew += dtTempAge * bmr[i][j];
        } else {
          // ignor x[Mbz] which contains only heat flux into ice
          bmr[i][j] = 0.0;
        }
      } else {
        bmr[i][j] = 0.0;
      }

      // enthalpy is *temporarily* known; drainage will change
      for (PetscInt k=0; k < fMz; k++) {
        Enthnew[k] = x[fMbz+1 + k];
      }

      // drain ice segments; alters Enthnew[]; adds to basal melt rate too
      PetscScalar Hdrainedtotal = 0.0;
      for (PetscInt k=0; k < ks; k++) {
        // modifies last two arguments, generally:
        PetscScalar dHdrained;
        ierr = drainageToBaseModelEnth(EC, omega_max, H[i][j], fzlev[k], fdz,
                                       Enthnew[k], dHdrained); CHKERRQ(ierr);
        Hdrainedtotal += dHdrained;
      }
      bmr[i][j] += Hdrainedtotal / dtTempAge;
      Hmeltnew += Hdrainedtotal;

      // transfer column into EnthNew3; communication later
      ierr = EnthNew3.setValColumnPL(i,j,fMz,fzlev,Enthnew); CHKERRQ(ierr);

      if (updateHmelt == PETSC_TRUE) {
        // finalize Hmelt value
        if (vMask.is_floating(i,j)) {
          // if floating assume maximally saturated "till"
          // UNACCOUNTED MASS & ENERGY (LATENT) LOSS/GAIN (TO/FROM OCEAN)!!
          vHmelt(i,j) = max_hmelt;
        } if (ks == 0) {
          vHmelt(i,j) = 0.0;  // no stored water on ice free land
        } else {
          // limit Hmelt to be in [0.0, max_hmelt]
          // UNACCOUNTED MASS & ENERGY (LATENT) LOSS (TO INFINITY AND BEYOND)!!
          vHmelt(i,j) = PetscMax(0.0, PetscMin(max_hmelt, Hmeltnew) );
        }
      }

    }
  }

  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vHmelt.end_access(); CHKERRQ(ierr);
  ierr = vRb.end_access(); CHKERRQ(ierr);
  ierr = vGhf.end_access(); CHKERRQ(ierr);
  ierr = vbasalMeltRate.end_access(); CHKERRQ(ierr);

  ierr = artm.end_access(); CHKERRQ(ierr);
  ierr = shelfbmassflux.end_access(); CHKERRQ(ierr);

  ierr = Tb3.end_access(); CHKERRQ(ierr);
  ierr = u3.end_access(); CHKERRQ(ierr);
  ierr = v3.end_access(); CHKERRQ(ierr);
  ierr = w3.end_access(); CHKERRQ(ierr);
  ierr = Sigma3.end_access(); CHKERRQ(ierr);
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = EnthNew3.end_access(); CHKERRQ(ierr);

  delete [] system.u;     delete [] system.v;     delete [] system.w;
  delete [] system.Sigma; delete [] system.Enth;  delete [] system.Enth_s; 
  delete [] system.Tb;

  delete [] x;
  delete [] Tbnew;   delete [] Enthnew;
  delete [] fzlev;  delete [] fzblev;

  return 0;
}


//! Move some of the liquid water fraction in a column segment [z,z+dz] to the base.
/*!
Generally this procedure alters the values in the last two arguments,
'enthalpy' and 'Hdrained'.  The latter argument is the ice-equivalent water
thickness which is moved to the bed by drainage.

Heuristic: Once liquid water fraction exceeds a cap, all of it goes to the base.
Follows \ref Greve97Greenland and references therein.
 */
PetscErrorCode IceEnthalpyModel::drainageToBaseModelEnth(
                EnthalpyConverter &EC,
                PetscScalar omega_max, PetscScalar thickness,
                PetscScalar z, PetscScalar dz,
                PetscScalar &enthalpy, PetscScalar &Hdrained) {
  PetscErrorCode ierr;

  if (allowAboveMelting == PETSC_TRUE) {
    SETERRQ(1,"IceEnthalpyModel::drainageToBaseModelEnth() called\n"
              "   BUT allowAboveMelting==TRUE");
  }

  // if there is liquid water already, thus temperate, consider whether there
  //   is enough to cause drainage;  UNACCOUNTED ENERGY LOSS IF E>E_l
  const PetscScalar p     = EC.getPressureFromDepth(thickness - z),
                    omega = EC.getWaterFractionLimited(enthalpy, p);
  if (omega > omega_max) {
    // drain water:
    Hdrained = (omega - omega_max) * dz;
    // update enthalpy because omega == omega_max now:
    ierr = EC.getEnthAtWaterFraction(omega_max, p, enthalpy); CHKERRQ(ierr);
  }
  return 0;
}

