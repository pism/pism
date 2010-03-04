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
#include "enthalpyConverter.hh"
#include "bedrockOnlySystem.hh"
#include "iceenthOnlySystem.hh"
#include "combinedSystem.hh"


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


//!  If this gets called then we need to extend the IceModelVec3s owned by IceEnthalpyModel.
PetscErrorCode IceEnthalpyModel::check_maximum_thickness_hook(const int old_Mz) {
  PetscErrorCode  ierr;

  // We use surface temperatures to extend Enth3 and Enthnew3. We get them from the
  // PISMSurfaceModel.
  if (surface != PETSC_NULL) {
    ierr = surface->ice_surface_temperature(grid.year, 0.0, artm); CHKERRQ(ierr);
  } else {
    SETERRQ(1,"PISM ERROR: surface == PETSC_NULL");
  }

  // vWork2d[0] will have the enthalpy of the air, very close to the value of
  //   Enth_ks in enthalpyAndDrainageStep() below
  EnthalpyConverter EC(config);
  ierr = vWork2d[0].begin_access(); CHKERRQ(ierr);
  ierr = artm.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = EC.getEnthPermissive(
         artm(i,j),0.0,EC.getPressureFromDepth(0.0),vWork2d[0](i,j));
         CHKERRQ(ierr);
    }
  }
  ierr = vWork2d[0].end_access(); CHKERRQ(ierr);
  ierr = artm.end_access(); CHKERRQ(ierr);

  // Model state 3D vectors:
  ierr = Enth3.extend_vertically(old_Mz, vWork2d[0]); CHKERRQ(ierr);

  // Work 3D vectors:
  ierr = EnthNew3.extend_vertically(old_Mz, vWork2d[0]); CHKERRQ(ierr);

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
    PetscScalar myLiquifiedVol = 0.0;
    ierr = enthalpyAndDrainageStep(vertSacrCount,&myLiquifiedVol);  CHKERRQ(ierr);

    PetscScalar gLiquifiedVol;
    ierr = PetscGlobalSum(&myLiquifiedVol, &gLiquifiedVol, grid.com); CHKERRQ(ierr);
    if (gLiquifiedVol > 0.0) {
      ierr = verbPrintf(1,grid.com,
        "\n IceEnthalpyModel WARNING: fully-liquified cells detected: volume liquified = %.3f km^3\n",
        gLiquifiedVol / 1.0e9); CHKERRQ(ierr);
    }
    

    // start & complete communication
    ierr = Enth3.beginGhostCommTransfer(EnthNew3); CHKERRQ(ierr);
    ierr = Enth3.endGhostCommTransfer(EnthNew3); CHKERRQ(ierr);

    ierr = setTnew3FromEnth3();  CHKERRQ(ierr);  // temperatureAgeStep() ASSUMES Tnew3 valid
  }
  return 0;
}


//! Compute the CTS value of enthalpy in an ice column, and the lambda for BOMBPROOF.
/*!
Return argument Enth_s[Mz] has the enthalpy value for the pressure-melting 
temperature at the corresponding z level.
 */
PetscErrorCode getEnthalpyCTSColumn(
      const NCConfigVariable &config, const EnthalpyConverter &EC, 
      const PetscInt Mz, const PetscScalar dzEQ, const PetscScalar *zlev,
      const PetscScalar thk, const PetscInt ks,
      const PetscScalar *Enth, const PetscScalar *w,
      PetscScalar *lambda, PetscScalar **Enth_s) {

  *lambda = 1.0;  // start with centered implicit for more accuracy

  const PetscScalar
      ice_rho_c = config.get("ice_density") * config.get("ice_specific_heat_capacity"),
      ice_k     = config.get("ice_thermal_conductivity");
  for (PetscInt k = 0; k <= ks; k++) {
    (*Enth_s)[k] = EC.getEnthalpyCTS(EC.getPressureFromDepth(thk - zlev[k]));

    if (Enth[k] > (*Enth_s)[k]) { // lambda = 0 if temperate ice present in column
      *lambda = 0.0;
    } else {
      const PetscScalar 
          denom = (PetscAbs(w[k]) + 0.000001/secpera) * ice_rho_c * dzEQ;
      *lambda = PetscMin(*lambda, 2.0 * ice_k / denom);
    }
  }

  const PetscScalar p_air = config.get("surface_pressure");
  for (PetscInt k = ks+1; k < Mz; k++) {
    (*Enth_s)[k] = EC.getEnthalpyCTS(p_air);
  }

  return 0;
}


PetscErrorCode reportColumnSolveError(
    const PetscErrorCode solve_ierr, MPI_Comm com, columnSystemCtx &sys, 
    const char *prefix, const PetscInt i, const PetscInt j) {

  char fname[PETSC_MAX_PATH_LEN];
  snprintf(fname, PETSC_MAX_PATH_LEN, "%s_i%d_j%d_zeropivot%d.m",
           prefix,i,j,solve_ierr);
  PetscErrorCode ierr;
  ierr = PetscPrintf(com,
    "\n\ntridiagonal solve in enthalpyAndDrainageStep(), for %sSystemCtx,\n"
        "   failed at (%d,%d) with zero pivot position %d\n"
        "   viewing system to file %s ... \n",
        prefix, i, j, solve_ierr, fname); CHKERRQ(ierr);
  PetscViewer viewer;
  ierr = PetscViewerCreate(com, &viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer, PETSC_VIEWER_ASCII);CHKERRQ(ierr);
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, fname);CHKERRQ(ierr);

  ierr = sys.viewSystem(viewer,"system"); CHKERRQ(ierr);

  ierr = PetscViewerDestroy(viewer); CHKERRQ(ierr);
  ierr = PetscPrintf(com, "\n   ENDING ...\n"); CHKERRQ(ierr);
  PetscEnd();
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
        "   viewing system to file %s ... \n",
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

This method modifies IceModelVec3 EnthNew3, IceModelVec3Bedrock Tb3,
IceModelVec2 vBasalMeltRate, and IceModelVec2 vHmelt.  No communication of
ghosts is done for any of these fields.
 */
PetscErrorCode IceEnthalpyModel::enthalpyAndDrainageStep(
                      PetscScalar* vertSacrCount, PetscScalar* liquifiedVol) {
  PetscErrorCode  ierr;

  if (doColdIceMethods==PETSC_TRUE) {
    SETERRQ(1,
       "\n\nPISM ERROR:  IceEnthalpyModel::enthalpyAndDrainageStep() called but\n"
           "             doColdIceMethods==true ... ending\n");
  }

  PetscInt    fMz, fMbz;  // number of fine grid levels in ice and bedrock, resp
  PetscScalar fdz, *fzlev, fdzb, *fzblev;
  ierr = grid.get_fine_vertical_grid(fMz, fMbz, fdz, fdzb, fzlev, fzblev); CHKERRQ(ierr);

  if (fMbz == 2) {
    SETERRQ(1,"method in IceEnthalpyModel does not make sense when fMbz == 2;\n"
              "   fMbz==1 and fMbz>2 are allowed\n");
  }

  EnthalpyConverter EC(config);

  const PetscScalar
    p_air     = config.get("surface_pressure"),
    ice_rho   = config.get("ice_density"),
    ice_c     = config.get("ice_specific_heat_capacity"),
    ice_k     = config.get("ice_thermal_conductivity"),
    L         = config.get("water_latent_heat_fusion"),  // J kg-1
    omega_max = config.get("liquid_water_fraction_max"), // pure
    max_hmelt = config.get("max_hmelt");                 // m

  PetscScalar *Enthnew, *Tbnew;
  Enthnew = new PetscScalar[fMz];  // new enthalpy in column
  Tbnew   = new PetscScalar[fMbz]; // new bedrock temperature in column

  combinedSystemCtx    cbsys(config, Enth3, fMz, fMbz);
  ierr = cbsys.initAllColumns(grid.dx, grid.dy, dtTempAge, fdz, fdzb); CHKERRQ(ierr);
  // space for solution when ice and bedrock are combined in one system
  PetscScalar *xcombined;
  xcombined = new PetscScalar[fMbz + fMz - 1];

  bedrockOnlySystemCtx bosys(config, fMbz);
  ierr = bosys.initAllColumns(dtTempAge, fdzb); CHKERRQ(ierr);

  iceenthOnlySystemCtx iosys(config, Enth3, fMz);
  ierr = iosys.initAllColumns(grid.dx, grid.dy, dtTempAge, fdz); CHKERRQ(ierr);

  bool viewOneColumn;
  ierr = check_option("-view_sys", viewOneColumn); CHKERRQ(ierr);

  // FIXME: verbosity failure?: option "-verbose 4" does not generate true here?
  if (getVerbosityLevel() >= 4) {  // view: all column-independent constants correct?
    ierr = EC.viewConstants(NULL); CHKERRQ(ierr);
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
  ierr = vbasalMeltRate.begin_access(); CHKERRQ(ierr);
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

  PetscInt liquifiedCount = 0;

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      // for fine grid; this should *not* be replaced by call to grid.kBelowHeight()
      const PetscInt ks = static_cast<PetscInt>(floor(vH(i,j)/fdz));
      // ignor advection and strain heating in ice if isMarginal
      const bool isMarginal = checkThinNeigh(
                                 vH(i+1,j),vH(i+1,j+1),vH(i,j+1),vH(i-1,j+1),
                                 vH(i-1,j),vH(i-1,j-1),vH(i,j-1),vH(i+1,j-1)  );

      // enthalpy and pressures at boundaries of ice
      const PetscScalar p_basal = EC.getPressureFromDepth(vH(i,j)),
                        p_ks    = EC.getPressureFromDepth(vH(i,j) - fzlev[ks]);
      PetscScalar Enth_air, Enth_ks;
      ierr = EC.getEnthPermissive(artm(i,j), 0.0, p_air, Enth_air); CHKERRQ(ierr);
      // in theory we could have a water fraction at k=ks level, but for
      //   now there is no case where we have that:
      ierr = EC.getEnthPermissive(artm(i,j), 0.0, p_ks,  Enth_ks); CHKERRQ(ierr);

      ierr = Enth3.getValColumn(i,j,fMz,fzlev,iosys.Enth); CHKERRQ(ierr);
      ierr = w3.getValColumn(i,j,fMz,fzlev,iosys.w); CHKERRQ(ierr);

      PetscScalar lambda;
      ierr = getEnthalpyCTSColumn(config, EC, fMz, fdz, fzlev,
                                  vH(i,j), ks, iosys.Enth, iosys.w,
                                  &lambda, &iosys.Enth_s); CHKERRQ(ierr);

      if (lambda < 1.0)  *vertSacrCount += 1; // count columns with lambda < 1

      // major decision: is cold base and grounded and has bedrock layer?:
      if ( (iosys.Enth[0] < iosys.Enth_s[0]) && (fMbz > 1) && (!vMask.is_floating(i,j)) ) {

        // ***** COLD BASE, GROUNDED CASE WITH BEDROCK *****
        ierr = cbsys.setIndicesAndClearThisColumn(i,j,ks); CHKERRQ(ierr);

        ierr = copyColumn(iosys.Enth,cbsys.Enth,fMz); CHKERRQ(ierr);
        ierr = copyColumn(iosys.Enth_s,cbsys.Enth_s,fMz); CHKERRQ(ierr);
        ierr = u3.getValColumn(i,j,fMz,fzlev,cbsys.u); CHKERRQ(ierr);
        ierr = v3.getValColumn(i,j,fMz,fzlev,cbsys.v); CHKERRQ(ierr);
        ierr = copyColumn(iosys.w,cbsys.w,fMz); CHKERRQ(ierr);
        ierr = Sigma3.getValColumn(i,j,fMz,fzlev,cbsys.Sigma); CHKERRQ(ierr);
        ierr = Tb3.getValColumn(i,j,fMbz,fzblev,cbsys.Tb); CHKERRQ(ierr);

        ierr = cbsys.setSchemeParamsThisColumn(isMarginal, lambda); CHKERRQ(ierr);
        ierr = cbsys.setBoundaryValuesThisColumn(Enth_ks, vGhf(i,j), vRb(i,j)); CHKERRQ(ierr);

        ierr = cbsys.solveThisColumn(&xcombined);
        if (ierr) {
          reportColumnSolveError(ierr, grid.com, cbsys, "combined", i, j);
        }
        if (viewOneColumn && issounding(i,j)) {
          ierr = reportColumn(grid.com, cbsys, "combined", i, j, xcombined, fMbz+fMz-1);
            CHKERRQ(ierr);
        }
        // break result x[] of combined system between Enthnew[fMz] and Tbnew[fMbz]
        for (PetscInt k = 0; k < fMbz-1; k++) {
          Tbnew[k] = xcombined[k];
        }
        // at this point we need a temperature from ice that could in extreme
        //   situations *be fully melted*; thus we catch the return code and
        //   and count this phenomenon
        ierr = EC.getAbsTemp(xcombined[fMbz-1], p_basal, Tbnew[fMbz-1]);
        if (ierr==1) { // return code of 1 means block of ice melted completely
          liquifiedCount++;
        } else CHKERRQ(ierr);
        for (PetscInt k = 0; k < fMz; k++) {
          Enthnew[k] = xcombined[k + fMbz-1];
        }

        vbasalMeltRate(i,j) = 0.0;  // zero melt rate if cold base

      } else {

        // ***** OTHER CASES *****
        // ***** BEDROCK ONLY SOLVE *****
        PetscScalar hf_base;
        if (fMbz > 1) { // deal with bedrock layer first, if present
          // case of temperate bed and a bedrock layer
          ierr = bosys.setIndicesAndClearThisColumn(i,j,-1); CHKERRQ(ierr);  

          ierr = Tb3.getValColumn(i,j,fMbz,fzblev,bosys.Tb);
                   CHKERRQ(ierr);

          const PetscScalar Tbtop = (vMask.is_floating(i,j)) ? shelfbtemp(i,j)
                                                             : EC.getMeltingTemp(p_basal);
          ierr = bosys.setBoundaryValuesThisColumn(Tbtop, vGhf(i,j)); CHKERRQ(ierr);

          ierr = bosys.solveThisColumn(&Tbnew);
          if (ierr) {
            reportColumnSolveError(ierr, grid.com, bosys, "bedrockOnly", i, j);
          }
          if (viewOneColumn && issounding(i,j)) {
            ierr = reportColumn(grid.com, bosys, "bedrockOnly", i, j, Tbnew, fMbz);
                CHKERRQ(ierr);
          }

          hf_base = bosys.extractHeatFluxFromSoln(Tbnew);
        } else {
          hf_base = vGhf(i,j);
        }

        // can determine melt now from heat flux out of base, etc.
        if (vMask.is_floating(i,j)) {
          vbasalMeltRate(i,j) = shelfbmassflux(i,j);
        } else {
          if (iosys.Enth[0] < iosys.Enth_s[0]) {
            vbasalMeltRate(i,j) = 0.0;  // zero melt rate if cold base
          } else {
            vbasalMeltRate(i,j) = ( hf_base + vRb(i,j) ) / (ice_rho * L);
          }
        }

        // ***** ICE ONLY SOLVE *****
        // now set-up for solve in ice; note iosys.Enth[], iosys.w[], iosys.Enth_s[]
        //   are already filled
        ierr = iosys.setIndicesAndClearThisColumn(i,j,ks); CHKERRQ(ierr);

        ierr = u3.getValColumn(i,j,fMz,fzlev,iosys.u); CHKERRQ(ierr);
        ierr = v3.getValColumn(i,j,fMz,fzlev,iosys.v); CHKERRQ(ierr);
        ierr = Sigma3.getValColumn(i,j,fMz,fzlev,iosys.Sigma); CHKERRQ(ierr);

        ierr = iosys.setSchemeParamsThisColumn(isMarginal, lambda); CHKERRQ(ierr);
        ierr = iosys.setBoundaryValuesThisColumn(Enth_ks); CHKERRQ(ierr);
        if (iosys.Enth[0] < iosys.Enth_s[0]) {
          // cold base case: ice base equation says heat flux is known
          const PetscScalar C = ice_c * fdz / ice_k;
          ierr = iosys.setLevel0EqnThisColumn(
                   1.0,-1.0,C * (hf_base + vRb(i,j))); CHKERRQ(ierr);
        } else {
          // determine lowest-level equation by vert vel at bottom of ice
          if (iosys.w[0] < 0.0) {
            // outflow "boundary condition"
            const PetscScalar nuw0 = (dtTempAge / fdz) * iosys.w[0];
            ierr = iosys.setLevel0EqnThisColumn(
                     1 - nuw0, nuw0, iosys.Enth[0]); CHKERRQ(ierr);            
          } else {
            // Dirichlet cond. for enthalpy at ice base
            ierr = iosys.setLevel0EqnThisColumn(
                     1.0,0.0,iosys.Enth_s[0]); CHKERRQ(ierr);
          }
        }

        ierr = iosys.solveThisColumn(&Enthnew);
        if (ierr) {
          reportColumnSolveError(ierr, grid.com, iosys, "iceenthOnly", i, j);
        }
        if (viewOneColumn && issounding(i,j)) {
          ierr = reportColumn(grid.com, iosys, "iceenthOnly", i, j, Enthnew, fMz);
              CHKERRQ(ierr);
        }
      }

      // basal melt rate causes water to be added to layer
      PetscScalar Hmeltnew = vHmelt(i,j);
      if (!vMask.is_floating(i,j)) {
        Hmeltnew += vbasalMeltRate(i,j) * dtTempAge;
      }

      // drain ice segments; alters Enthnew[]; adds to both basal melt rate and Hmelt;
      //    has side-effect that Enthnew[] is ice with at most omega_max liquid
      //    fraction
      PetscScalar Hdrainedtotal = 0.0;
      for (PetscInt k=0; k < ks; k++) {
        PetscScalar dHdrained = 0.0;
        if (EC.isLiquified(Enthnew[k],EC.getPressureFromDepth(vH(i,j) - fzlev[k]))) {
          liquifiedCount++;
        }
        // modifies last two arguments, generally:
        ierr = drainageToBaseModelEnth(EC, omega_max, vH(i,j), fzlev[k], fdz,
                                       Enthnew[k], dHdrained); CHKERRQ(ierr);
        Hdrainedtotal += dHdrained;  // always a positive contribution
      }
      if (!vMask.is_floating(i,j)) {
        vbasalMeltRate(i,j) += Hdrainedtotal / dtTempAge;
        Hmeltnew += Hdrainedtotal;
      }

      // transfer column into EnthNew3; communication later
      ierr = EnthNew3.setValColumnPL(i,j,fMz,fzlev,Enthnew); CHKERRQ(ierr);

      // if no thermal layer then need to fill Tb directly
      if (fMbz == 1) {
        if (vMask.is_floating(i,j)) { // floating: get from PISMOceanModel
          Tbnew[0] = shelfbtemp(i,j);
        } else {                      // grounded: duplicate temp from ice
          ierr = EC.getAbsTemp(Enthnew[0],EC.getPressureFromDepth(vH(i,j)), Tbnew[0]);
                    CHKERRQ(ierr);
        }
      }

      // transfer column into Tb3; no need for communication, even later
      ierr = Tb3.setValColumnPL(i,j,fMbz,fzblev,Tbnew); CHKERRQ(ierr);

      // finalize Hmelt value
      if (updateHmelt == PETSC_TRUE) {
        if (vMask.is_floating(i,j)) {
          // if floating assume maximally saturated "till"
          // UNACCOUNTED MASS & ENERGY (LATENT) LOSS/GAIN (TO/FROM OCEAN)!!
          vHmelt(i,j) = max_hmelt;
        } else if (ks == 0) {
          vHmelt(i,j) = 0.0;  // no stored water on ice free land
        } else {
          // limit Hmelt to be in [0.0, max_hmelt]
          // UNACCOUNTED MASS & ENERGY (LATENT) LOSS (TO INFINITY AND BEYOND)!!
          vHmelt(i,j) = PetscMax(0.0, PetscMin(max_hmelt, Hmeltnew) );
        }
      }

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
  ierr = vbasalMeltRate.end_access(); CHKERRQ(ierr);

  ierr = Tb3.end_access(); CHKERRQ(ierr);
  ierr = u3.end_access(); CHKERRQ(ierr);
  ierr = v3.end_access(); CHKERRQ(ierr);
  ierr = w3.end_access(); CHKERRQ(ierr);
  ierr = Sigma3.end_access(); CHKERRQ(ierr);
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = EnthNew3.end_access(); CHKERRQ(ierr);

  delete [] Enthnew; delete [] Tbnew;  delete [] xcombined;
  delete [] fzlev;   delete [] fzblev;

  *liquifiedVol = ((double) liquifiedCount) * fdz * grid.dx * grid.dy;
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
  } else {
    Hdrained = 0.0;
  }
  return 0;
}

