// Copyright (C) 2009 Andreas Aschwanden and Ed Bueler
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


/*!
This constructor just sets flow law factor for nonzero water content, from
\ref AschwandenBlatter2009 and \ref LliboutryDuval1985.
 */
PolyThermalGPBLDIce::PolyThermalGPBLDIce(MPI_Comm c,const char pre[]) : ThermoGlenIce(c,pre) {
  EC = NULL;
  T_0 = 273.15;               // default overridden through config interface
  water_frac_coeff = 184.0;   // default overridden through config interface
}


PetscErrorCode PolyThermalGPBLDIce::setFromConfig(NCConfigVariable *config) {
  if (config == NULL) {
    PetscPrintf(PETSC_COMM_WORLD,
      "\n\n\n  PolyThermalGPBLDIce ERROR in setFromConfig():\n"
       "    config == NULL ... \n\n");
    PetscEnd();
  }
  EC = new EnthalpyConverter(config);
  T_0  = config->get("water_melting_temperature");    // K
  water_frac_coeff = config->get("gpbld_water_frac_coeff");                
  return 0;
}


PolyThermalGPBLDIce::~PolyThermalGPBLDIce() {
  delete EC;
}


PetscErrorCode PolyThermalGPBLDIce::setFromOptions() {
  PetscErrorCode ierr;

  ierr = ThermoGlenIce::setFromOptions(); CHKERRQ(ierr);
  
  ierr = PetscOptionsBegin(comm,prefix,"PolyThermalGPBLDIce options",NULL);CHKERRQ(ierr);
  {
    ierr = PetscOptionsReal("-ice_gpbld_water_frac_coeff",
      "coefficient of softness factor in temperate ice, as function of liquid water fraction (no units)",
      "",water_frac_coeff,&water_frac_coeff,NULL);CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PolyThermalGPBLDIce::view(PetscViewer viewer) const {
  PetscErrorCode ierr;

  ierr = ThermoGlenIce::view(viewer); CHKERRQ(ierr);
  
  PetscTruth iascii;
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(comm,&viewer);CHKERRQ(ierr);
  }
  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii);CHKERRQ(ierr);
  if (iascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"<\nderived PolyThermalGPBLDIce object (%s)\n",prefix);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  T_0             =%10.3f (K)\n",T_0);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  water_frac_coeff=%10.1f\n",water_frac_coeff);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,">\n",water_frac_coeff);CHKERRQ(ierr);
  } else {
    SETERRQ(1,"No binary viewer for this object\n");
  }
  return 0;
}


//! The softness factor in the Glen-Paterson-Budd-Lliboutry-Duval flow law.  For constitutive law form.
/*!
This is a modification of Glen-Paterson-Budd ice, which is ThermoGlenIce.  In particular, if
\f$A()\f$ is the softness factor for ThermoGlenIce, if \f$E\f$ is the enthalpy, and \f$p\f$ is
the pressure then the softness we compute is
   \f[A = A(T_{pa}(E,p))(1+184\omega).\f]
The pressure-melting temperature \f$T_{pa}(E,p)\f$ is computed by getPATemp().
 */
PetscScalar PolyThermalGPBLDIce::softnessParameterFromEnth(
                PetscScalar enthalpy, PetscScalar pressure) const {
  if (EC == NULL) {
    PetscPrintf(PETSC_COMM_WORLD,"EC is NULL in PolyThermalGPBLDIce::flowFromEnth()... ending\n");
    PetscEnd();
  }
  PetscScalar E_s, E_l;
  EC->getEnthalpyInterval(pressure, E_s, E_l);
  if (enthalpy <= E_s) {       // cold ice
    return softnessParameter( EC->getPATemp(enthalpy,pressure) ); // uses ThermoGlenIce formula
  } else if (enthalpy < E_l) { // temperate ice
    const PetscScalar omega = EC->getWaterFraction(enthalpy,pressure);
    // next line implements eqn (23) in \ref AschwandenBlatter2009
    return softnessParameter(T_0) * (1.0 + water_frac_coeff * omega);  // uses ThermoGlenIce formula
  } else { // liquid water not allowed
    PetscPrintf(PETSC_COMM_WORLD,
      "\n\n\n  PISM ERROR in PolyThermalGlenPBLDIce::flow(): liquid water not allowed; ending ... \n\n");
    PetscEnd();
    return 0.0;
  }
}


//! The hardness factor in the Paterson-Budd-Lliboutry-Duval flow law.  For viscosity form.
PetscScalar PolyThermalGPBLDIce::hardnessParameterFromEnth(
                PetscScalar enthalpy, PetscScalar pressure) const {
  return pow(softnessParameterFromEnth(enthalpy,pressure), -1.0/n);
}


//! Glen-Paterson-Budd-Lliboutry-Duval flow law itself.
PetscScalar PolyThermalGPBLDIce::flowFromEnth(
                PetscScalar stress, PetscScalar enthalpy, PetscScalar pressure, PetscScalar /* gs */) const {
  return softnessParameterFromEnth(enthalpy,pressure) * pow(stress,n-1);
}


PetscScalar PolyThermalGPBLDIce::effectiveViscosityColumnFromEnth(
                PetscScalar thickness,  PetscInt kbelowH, const PetscScalar *zlevels,
                PetscScalar u_x,  PetscScalar u_y, PetscScalar v_x,  PetscScalar v_y,
                const PetscScalar *enthalpy1, const PetscScalar *enthalpy2) const {
  if (EC == NULL) {
    PetscPrintf(PETSC_COMM_WORLD,
       "EC is NULL in PolyThermalGPBLDIce::effectiveViscosityColumnFromEnth()... ending\n");
    PetscEnd();
  }

  // result is \nu_e H, i.e. viscosity times thickness; B is really hardness times thickness
  // integrates the hardness parameter using the trapezoid rule.
  PetscScalar B = 0;
  if (kbelowH > 0) {
    PetscScalar dz = zlevels[1] - zlevels[0];
    B += 0.5 * dz * hardnessParameterFromEnth( 0.5 * (enthalpy1[0] + enthalpy2[0]),
                                               EC->getPressureFromDepth(thickness) );
    for (PetscInt m=1; m < kbelowH; m++) {
      const PetscScalar dzNEXT = zlevels[m+1] - zlevels[m],
                        depth  = thickness - 0.5 * (zlevels[m+1] + zlevels[m]);
      B += 0.5 * (dz + dzNEXT) * hardnessParameterFromEnth( 0.5 * (enthalpy1[m] + enthalpy2[m]),
                                                            EC->getPressureFromDepth(depth) );
      dz = dzNEXT;
    }
    // use last dz from for loop
    const PetscScalar depth  = 0.5 * (thickness - zlevels[kbelowH]);
    B += 0.5 * dz * hardnessParameterFromEnth( 0.5 * (enthalpy1[kbelowH] + enthalpy2[kbelowH]),
                                               EC->getPressureFromDepth(depth) );
  }
  const PetscScalar alpha = secondInvariant(u_x, u_y, v_x, v_y);
  return 0.5 * B * pow(schoofReg + alpha, (1-n)/(2*n));
}


#define ICE_GPBLD      "gpbld"
//! Create new kind of ice, PolyThermalGPBLDIce.  For IceFactory registration.
static PetscErrorCode create_gpbld(MPI_Comm comm,const char pre[],IceType **i) {
  *i = new (PolyThermalGPBLDIce)(comm,pre);  return 0;
}


/*********** procedures for init ****************/

IceEnthalpyModel::IceEnthalpyModel(IceGrid &g) : IceModel(g) {
  doColdIceMethods = false;
}


PetscErrorCode IceEnthalpyModel::createVecs() {
  PetscErrorCode ierr;

  ierr = Enth3.create(grid, "enthalpy", true); CHKERRQ(ierr);
  // POSSIBLE standard name = land_ice_enthalpy
  ierr = Enth3.set_attrs(
     "model_state",
     "ice enthalpy (sensible heat plus latent heat of liquid fraction)",
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


PetscErrorCode IceEnthalpyModel::init_physics() {
  PetscErrorCode ierr;

  // let the base class create the ice and process its options:
  ierr = IceModel::init_physics(); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
      "  setting flow law to Glen-Paterson-Budd-Lliboutry-Duval type ...\n");
      CHKERRQ(ierr);
  ierr = iceFactory.registerType(ICE_GPBLD, &create_gpbld);
  if (ierr != 0) {
    PetscPrintf(grid.com,
       "FAILURE OF iceFactory.registerType() ... return value %d ... ending ....\n",ierr);
    PetscEnd();
  }
  CHKERRQ(ierr);
  if (ice != NULL)  delete ice;  // kill choice already made!
  iceFactory.setType(ICE_GPBLD); // new flowlaw which has dependence on enthalpy not temperature
  iceFactory.create(&ice);

  PolyThermalGPBLDIce *gpbldi = dynamic_cast<PolyThermalGPBLDIce*>(ice);
  if (gpbldi) {
    gpbldi->setFromConfig(&config);
  } else {
    ThermoGlenIce *tgi = dynamic_cast<ThermoGlenIce*>(ice);
    if (tgi) {
      ierr = verbPrintf(2, grid.com,
        "  [flow law was actually set to ThermoGlenIce by IceEnthalpyModel ...]\n"); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2, grid.com,
        "  [WARNING: flow law unclear in IceEnthalpyModel ...]\n"); CHKERRQ(ierr);
    }
  }
  
  ierr = ice->printInfo(4);CHKERRQ(ierr); // DEBUG

  ierr = ice->setFromOptions();CHKERRQ(ierr);

  return 0;
}



/*********** procedures for read/write ****************/

PetscErrorCode IceEnthalpyModel::write_extra_fields(const char filename[]) {
  PetscErrorCode ierr;

  if (doColdIceMethods) { // in this case, just update Enth3 to reflect
                                  // temperature in ice at final time
    ierr = verbPrintf(2, grid.com,
      "  using temperature to set enthalpy for writing (as cold ice) ...\n");
      CHKERRQ(ierr);
    ierr = setEnth3FromT3_ColdIce(); CHKERRQ(ierr);
  }
  ierr = Enth3.write(filename, NC_DOUBLE); CHKERRQ(ierr);

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

  // reset attributes on EnthNew3, a temporary; probaly not needed
  ierr = EnthNew3.set_name("enthalpy_new"); CHKERRQ(ierr);
  ierr = EnthNew3.set_attrs(
     "internal",
     "ice enthalpy; temporary space during timestep",
     "J kg-1",
     ""); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceEnthalpyModel::initFromFile(const char *fname) {
  PetscErrorCode  ierr;

  ierr = IceModel::initFromFile(fname); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
     "entering IceEnthalpyModel::initFromFile() after base class version;\n"
     "  looking in '%s' for variable 'enthalpy' ... \n",fname);
     CHKERRQ(ierr);

  NCTool nc(&grid);
  ierr = nc.open_for_reading(fname); CHKERRQ(ierr);

/* if we were to require "enthalpy" to be present then the code would be simpler:
  ierr = Enth3.read(fname, last_record); CHKERRQ(ierr);
*/

  grid_info g;
  ierr = nc.get_grid_info(g); CHKERRQ(ierr);
  bool enthExists=false;
  ierr = nc.find_variable("enthalpy", NULL, enthExists); CHKERRQ(ierr);

  if (enthExists) {
    // act like we are regridding the variable
    double *zlevs = NULL, *zblevs = NULL; // NULLs correspond to 2D-only regridding
    if ((g.z_len != 0) && (g.zb_len != 0)) {
      ierr = nc.get_vertical_dims(zlevs, zblevs); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(1, grid.com,
         "PISM ERROR: -i file does not look right; at least one of 'z' and 'zb' is absent in '%s'.\n",
         fname); CHKERRQ(ierr);
      PetscEnd();
    }
    ierr = nc.close(); CHKERRQ(ierr);
    LocalInterpCtx lic(g, zlevs, zblevs, grid);
    ierr = Enth3.regrid(fname, lic, true); CHKERRQ(ierr);  // at this point, it is critical
  } else {
    ierr = verbPrintf(2, grid.com,
      "  variable 'enthalpy' not found so setting it as cold ice, from temperature ...\n");
      CHKERRQ(ierr);
    ierr = setEnth3FromT3_ColdIce(); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode IceEnthalpyModel::regrid() {
  PetscErrorCode ierr;
  
  ierr = IceModel::regrid(); CHKERRQ(ierr);
  
  PetscTruth regridVarsSet, regrid_from_set;
  char filename[PETSC_MAX_PATH_LEN], regridVars[PETSC_MAX_PATH_LEN];
  NCTool nc(&grid);

  // Get the regridding file name:
  ierr = PetscOptionsGetString(PETSC_NULL, "-regrid_from", filename, PETSC_MAX_PATH_LEN,
                               &regrid_from_set); CHKERRQ(ierr);

  // Return if no regridding is requested:
  if (!regrid_from_set) return 0;

  ierr = PetscOptionsGetString(PETSC_NULL, "-regrid_vars", regridVars,
                               PETSC_MAX_PATH_LEN, &regridVarsSet); CHKERRQ(ierr);
  if (regridVarsSet == PETSC_FALSE) {
    strcpy(regridVars, "");
  }
  if (!strchr(regridVars, 'y'))   return 0;

  // create "local interpolation context" from dimensions, limits, and lengths
  //   extracted from regridFile, and from information about the part of the
  //   grid owned by this processor

  ierr = nc.open_for_reading(filename);
  
  grid_info g;
  // Note that after this call g.z_len and g.zb_len are zero if the
  // corresponding dimension does not exist.
  ierr = nc.get_grid_info(g); CHKERRQ(ierr);

  double *zlevs = NULL, *zblevs = NULL; // NULLs correspond to 2D-only regridding
  if ((g.z_len != 0) && (g.zb_len != 0)) {
    ierr = nc.get_vertical_dims(zlevs, zblevs); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com,
		      "PISM WARNING: at least one of 'z' and 'zb' is absent in '%s'.\n"
		      "              3D regridding is disabled.\n",
		      filename);
    CHKERRQ(ierr);
  }
  ierr = nc.close(); CHKERRQ(ierr);

  { // explicit scoping means destructor will be called for lic
    LocalInterpCtx lic(g, zlevs, zblevs, grid);

    if (lic.regrid_2d_only) {
      ierr = verbPrintf(2, grid.com,
           "IceEnthalpyModel WARNING: unable to regrid enthalpy ('y') from NetCDF file '%s'!\n",
           filename); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2,grid.com, 
         "IceEnthalpyModel regridding enthalpy ('y') from NetCDF file '%s'\n", 
         filename); CHKERRQ(ierr);
      ierr = Enth3.regrid(filename, lic, true); CHKERRQ(ierr);
    }
  }

  // Note that deleting a NULL pointer is safe.
  delete [] zlevs;  delete [] zblevs;
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

  EnthalpyConverter EC(&config);
  
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
        Enthij[k] = EC.getEnthPermissive(Tij[k],0.0, EC.getPressureFromDepth(depth) );
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

  EnthalpyConverter EC(&config);

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
        Tij[k] = EC.getAbsTemp(Enthij[k], EC.getPressureFromDepth(depth));
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

  ierr = useForLiquidFrac.set_name("liquid_frac"); CHKERRQ(ierr);
  ierr = useForLiquidFrac.set_attrs(
     "diagnostic",
     "liquid water fraction in ice (between 0 and 1)",
     "",
     ""); CHKERRQ(ierr);

  EnthalpyConverter EC(&config);

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
        omegaij[k] = EC.getWaterFraction(Enthij[k], EC.getPressureFromDepth(depth));
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
     "pressure-adjusted ice temperature",
     "deg_C",
     ""); CHKERRQ(ierr);

  EnthalpyConverter EC(&config);

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
        Tpaij[k] = EC.getPATemp(Enthij[k], EC.getPressureFromDepth(depth));
      }
    }
  }
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = useForPATemp.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  // communication not done; we allow global IceModelVec3s as useForPATemp
  return 0;
}


/*********** velocity routines in which new flow law gets used ****************/

//! Total code duplication with IceModel version, but checks flag doColdIceMethods and uses correct flow law.
PetscErrorCode IceEnthalpyModel::velocitySIAStaggered() {
  PetscErrorCode  ierr;

  PetscScalar *delta, *I, *J, *K, *Sigma;
  delta = new PetscScalar[grid.Mz];
  I = new PetscScalar[grid.Mz];
  J = new PetscScalar[grid.Mz];
  K = new PetscScalar[grid.Mz];
  Sigma = new PetscScalar[grid.Mz];

  PetscScalar **h_x[2], **h_y[2], **H, **uvbar[2];

  PetscScalar *Tij, *Toffset, *ageij, *ageoffset;

  const bool usetau3 = (IceTypeUsesGrainSize(ice) && (realAgeForGrainSize == PETSC_TRUE));

  const PetscTruth usesGrainSize = IceTypeUsesGrainSize(ice);

  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = vWork2d[0].get_array(h_x[0]); CHKERRQ(ierr);
  ierr = vWork2d[1].get_array(h_x[1]); CHKERRQ(ierr);
  ierr = vWork2d[2].get_array(h_y[0]); CHKERRQ(ierr);
  ierr = vWork2d[3].get_array(h_y[1]); CHKERRQ(ierr);
  ierr = vuvbar[0].get_array(uvbar[0]); CHKERRQ(ierr);
  ierr = vuvbar[1].get_array(uvbar[1]); CHKERRQ(ierr);

  ierr = T3.begin_access(); CHKERRQ(ierr);
  if (usetau3) {
    ierr = tau3.begin_access(); CHKERRQ(ierr);
  }
  ierr = w3.begin_access(); CHKERRQ(ierr);
  ierr = Istag3[0].begin_access(); CHKERRQ(ierr);
  ierr = Istag3[1].begin_access(); CHKERRQ(ierr);
  ierr = Sigmastag3[0].begin_access(); CHKERRQ(ierr);
  ierr = Sigmastag3[1].begin_access(); CHKERRQ(ierr);

  PetscScalar *Enthij, *Enthoffset;
  PolyThermalGPBLDIce *gpbldi = NULL;
  if (!doColdIceMethods) {
    gpbldi = dynamic_cast<PolyThermalGPBLDIce*>(ice);
    if (!gpbldi) {
      PetscPrintf(grid.com,
        "doColdIceMethods==false in IceEnthalpyMethod::velocitySIAStaggered()\n"
        "   but not using PolyThermalGPBLDIce ... ending ....\n");
      PetscEnd();
    }
    ierr = Enth3.begin_access(); CHKERRQ(ierr);
  }

  // staggered grid computation of: I, J, Sigma
  for (PetscInt o=0; o<2; o++) {
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
        // staggered point: o=0 is i+1/2, o=1 is j+1/2,
        //   (i,j) and (i+oi,j+oj) are reg grid neighbors of staggered pt:
        const PetscInt     oi = 1-o, oj=o;
        const PetscScalar  slope = (o==0) ? h_x[o][i][j] : h_y[o][i][j];
        const PetscScalar  thickness = 0.5 * (H[i][j] + H[i+oi][j+oj]);

        if (thickness > 0) {
          ierr = T3.getInternalColumn(i,j,&Tij); CHKERRQ(ierr);
          ierr = T3.getInternalColumn(i+oi,j+oj,&Toffset); CHKERRQ(ierr);
          if (usetau3) {
            ierr = tau3.getInternalColumn(i,j,&ageij); CHKERRQ(ierr);
            ierr = tau3.getInternalColumn(i+oi,j+oj,&ageoffset); CHKERRQ(ierr);
          }

          if (!doColdIceMethods) {
            ierr = Enth3.getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
            ierr = Enth3.getInternalColumn(i+oi,j+oj,&Enthoffset); CHKERRQ(ierr);
          }

          // does validity check for thickness:
          const PetscInt      ks = grid.kBelowHeight(thickness);
          const PetscScalar   alpha =
                  sqrt(PetscSqr(h_x[o][i][j]) + PetscSqr(h_y[o][i][j]));

          I[0] = 0;   J[0] = 0;   K[0] = 0;
          for (PetscInt k=0; k<=ks; ++k) {
            const PetscScalar   pressure = ice->rho * earth_grav * (thickness - grid.zlevels[k]);
            PetscScalar flow,grainsize = constantGrainSize;
            if (usesGrainSize && realAgeForGrainSize) {
              grainsize = grainSizeVostok(0.5 * (ageij[k] + ageoffset[k]));
            }
            // If the flow law does not use grain size, it will just ignore it, no harm there
            if (doColdIceMethods) {
              flow = ice->flow(alpha * pressure, 0.5 * (Tij[k] + Toffset[k]), pressure, grainsize);
            } else {
              flow = gpbldi->flowFromEnth(alpha * pressure, 0.5 * (Enthij[k] + Enthoffset[k]), 
                                          pressure, grainsize);
            }

            delta[k] = 2.0 * pressure * enhancementFactor * flow;

            // for Sigma, ignore mask value and assume SHEET; will be overwritten
            // by correctSigma() in iMssa.cc
            Sigma[k] = delta[k] * PetscSqr(alpha) * pressure;

            if (k>0) { // trapezoid rule for I[k] and K[k]
              const PetscScalar dz = grid.zlevels[k] - grid.zlevels[k-1];
              I[k] = I[k-1] + 0.5 * dz * (delta[k-1] + delta[k]);
              K[k] = K[k-1] + 0.5 * dz * (grid.zlevels[k-1] * delta[k-1]
                                          + grid.zlevels[k] * delta[k]);
              J[k] = grid.zlevels[k] * I[k] - K[k];
            }
          }
          for (PetscInt k=ks+1; k<grid.Mz; ++k) { // above the ice
            Sigma[k] = 0.0;
            I[k] = I[ks];
            J[k] = grid.zlevels[k] * I[ks];
          }

          // diffusivity for deformational flow (vs basal diffusivity, incorporated in ub,vb)
          const PetscScalar  Dfoffset = J[ks] + (thickness - grid.zlevels[ks]) * I[ks];

          // vertically-averaged SIA-only velocity, sans sliding;
          //   note uvbar[0][i][j] is  u  at right staggered point (i+1/2,j)
          //   but uvbar[1][i][j] is  v  at up staggered point (i,j+1/2)
          uvbar[o][i][j] = - Dfoffset * slope / thickness;

          ierr = Istag3[o].setValColumnPL(i,j,grid.Mz,grid.zlevels,I); CHKERRQ(ierr);
          ierr = Sigmastag3[o].setValColumnPL(i,j,grid.Mz,grid.zlevels,Sigma); CHKERRQ(ierr);
        } else {  // zero thickness case
          uvbar[o][i][j] = 0;
          ierr = Istag3[o].setColumn(i,j,0.0); CHKERRQ(ierr);
          ierr = Sigmastag3[o].setColumn(i,j,0.0); CHKERRQ(ierr);
        }
      } // o
    } // j
  } // i

  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vuvbar[0].end_access(); CHKERRQ(ierr);
  ierr = vuvbar[1].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[0].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[1].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[2].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[3].end_access(); CHKERRQ(ierr);

  ierr = T3.end_access(); CHKERRQ(ierr);
  if (usetau3) {
    ierr = tau3.end_access(); CHKERRQ(ierr);
  }
  ierr = w3.end_access(); CHKERRQ(ierr);
  ierr = Sigmastag3[0].end_access(); CHKERRQ(ierr);
  ierr = Sigmastag3[1].end_access(); CHKERRQ(ierr);
  ierr = Istag3[0].end_access(); CHKERRQ(ierr);
  ierr = Istag3[1].end_access(); CHKERRQ(ierr);

  if (!doColdIceMethods) {
    ierr = Enth3.end_access(); CHKERRQ(ierr);
  }

  delete [] delta;   delete [] I;   delete [] J;   delete [] K;   delete [] Sigma;

  return 0;
}


PetscErrorCode IceEnthalpyModel::computeEffectiveViscosity(IceModelVec2 vNuH[2], PetscReal epsilon) {
  PetscErrorCode ierr;

  if (leaveNuHAloneSSA == PETSC_TRUE) {
    return 0;
  }

  //CHECK_NOT_SSA_EXTERNAL(ssa);
  if (ssa) {SETERRQ(1,"This should not be called when the external SSA solver is active");}

  if (useConstantNuHForSSA == PETSC_TRUE) {
    // Intended only for debugging, this treats the entire domain as though it was the strength extension
    // (i.e. strength does not even depend on thickness)
    PetscReal nuH = ssaStrengthExtend.notional_strength();
    ierr = vNuH[0].set(nuH); CHKERRQ(ierr);
    ierr = vNuH[1].set(nuH); CHKERRQ(ierr);
    return 0;
  }

  // We need to compute integrated effective viscosity (\bar\nu * H).
  // It is locally determined by the strain rates and temperature field.
  PetscScalar *Tij, *Toffset, **H, **nuH[2], **u, **v;
  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = T3.begin_access(); CHKERRQ(ierr);
  ierr = vNuH[0].get_array(nuH[0]); CHKERRQ(ierr);
  ierr = vNuH[1].get_array(nuH[1]); CHKERRQ(ierr);

  ierr = vubarSSA.get_array(u); CHKERRQ(ierr);
  ierr = vvbarSSA.get_array(v); CHKERRQ(ierr);

  PetscScalar *Enthij, *Enthoffset;
  PolyThermalGPBLDIce *gpbldi = NULL;
  if (!doColdIceMethods) {
    gpbldi = dynamic_cast<PolyThermalGPBLDIce*>(ice);
    if (!gpbldi) {
      PetscPrintf(grid.com,
        "doColdIceMethods==false in IceEnthalpyMethod::computeEffectiveViscosity()\n"
        "   but not using PolyThermalGPBLDIce ... ending ....\n");
      PetscEnd();
    }
    ierr = Enth3.begin_access(); CHKERRQ(ierr);
  }

  for (PetscInt o=0; o<2; ++o) {
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        if (H[i][j] < ssaStrengthExtend.min_thickness_for_extension()) {
          // Extends strength of SSA (i.e. nuH coeff) into the ice free region.  Does not add or subtract ice mass.
          nuH[o][i][j] = ssaStrengthExtend.notional_strength();
        } else {
          const PetscInt      oi = 1-o, oj=o;
          const PetscScalar   dx = grid.dx,
                              dy = grid.dy;
          PetscScalar u_x, u_y, v_x, v_y;
          // Check the offset to determine how to differentiate velocity
          if (o == 0) {
            u_x = (u[i+1][j] - u[i][j]) / dx;
            u_y = (u[i][j+1] + u[i+1][j+1] - u[i][j-1] - u[i+1][j-1]) / (4*dy);
            v_x = (v[i+1][j] - v[i][j]) / dx;
            v_y = (v[i][j+1] + v[i+1][j+1] - v[i][j-1] - v[i+1][j-1]) / (4*dy);
          } else {
            u_x = (u[i+1][j] + u[i+1][j+1] - u[i-1][j] - u[i-1][j+1]) / (4*dx);
            u_y = (u[i][j+1] - u[i][j]) / dy;
            v_x = (v[i+1][j] + v[i+1][j+1] - v[i-1][j] - v[i-1][j+1]) / (4*dx);
            v_y = (v[i][j+1] - v[i][j]) / dy;
          }
          const PetscScalar myH = 0.5 * (H[i][j] + H[i+oi][j+oj]);

          if (doColdIceMethods) {
            ierr = T3.getInternalColumn(i,j,&Tij); CHKERRQ(ierr);
            ierr = T3.getInternalColumn(i+oi,j+oj,&Toffset); CHKERRQ(ierr);
            nuH[o][i][j] = ice->effectiveViscosityColumn(
                                myH, grid.kBelowHeight(myH), grid.zlevels,
                                u_x, u_y, v_x, v_y, Tij, Toffset);
          } else {
            ierr = Enth3.getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
            ierr = Enth3.getInternalColumn(i+oi,j+oj,&Enthoffset); CHKERRQ(ierr);
            nuH[o][i][j] = gpbldi->effectiveViscosityColumnFromEnth(
                                myH, grid.kBelowHeight(myH), grid.zlevels,
                                u_x, u_y, v_x, v_y, Enthij, Enthoffset);
          }

          if (! finite(nuH[o][i][j]) || false) {
            ierr = PetscPrintf(grid.com, "nuH[%d][%d][%d] = %e\n", o, i, j, nuH[o][i][j]);
              CHKERRQ(ierr);
            ierr = PetscPrintf(grid.com, "  u_x, u_y, v_x, v_y = %e, %e, %e, %e\n",
                               u_x, u_y, v_x, v_y);
              CHKERRQ(ierr);
          }

          // We ensure that nuH is bounded below by a positive constant.
          nuH[o][i][j] += epsilon;
        }
      }
    }
  }
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = T3.end_access(); CHKERRQ(ierr);
  ierr = vNuH[0].end_access(); CHKERRQ(ierr);
  ierr = vNuH[1].end_access(); CHKERRQ(ierr);
  ierr = vubarSSA.end_access(); CHKERRQ(ierr);
  ierr = vvbarSSA.end_access(); CHKERRQ(ierr);

  if (!doColdIceMethods) {
    ierr = Enth3.end_access(); CHKERRQ(ierr);
  }

  // Some communication
  ierr = vNuH[0].beginGhostComm(); CHKERRQ(ierr);
  ierr = vNuH[0].endGhostComm(); CHKERRQ(ierr);
  ierr = vNuH[1].beginGhostComm(); CHKERRQ(ierr);
  ierr = vNuH[1].endGhostComm(); CHKERRQ(ierr);
  return 0;
}


/*********** timestep routines ****************/

PetscErrorCode IceEnthalpyModel::temperatureStep(
     PetscScalar* vertSacrCount, PetscScalar* bulgeCount) {
  PetscErrorCode ierr;
  if (doColdIceMethods) {
    ierr = verbPrintf(4,grid.com,
      "    [IceEnthalpyModel::temperatureStep(): ENTHALPY IS OFF. CALLING IceModel::temperatureStep()]\n");
      CHKERRQ(ierr);
    ierr = IceModel::temperatureStep(vertSacrCount,bulgeCount);  CHKERRQ(ierr);
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

    ierr = setTnew3FromEnth3();  CHKERRQ(ierr);  // INEFFICIENT, BUT temperatureAgeStep() ASSUMES Tnew3 valid
  }
  return 0;
}


PetscErrorCode IceEnthalpyModel::enthalpyAndDrainageStep(PetscScalar* vertSacrCount, PetscScalar* bulgeCount) {
  PetscErrorCode  ierr;

  if (doColdIceMethods) {
    ierr = PetscPrintf(grid.com,
      "\n\n    IceEnthalpyModel::enthalpyAndDrainageStep() called but doColdIceMethods==true ... ending\n");
      CHKERRQ(ierr);
    PetscEnd();
  }

  // set up fine grid in ice and bedrock
  PetscInt    fMz, fMbz;
  PetscScalar fdz, *fzlev, fdzb, *fzblev;
  ierr = grid.getFineEqualVertCounts(fMz,fMbz); CHKERRQ(ierr);
  fzlev = new PetscScalar[fMz];
  fzblev = new PetscScalar[fMbz];
  ierr = grid.getFineEqualVertLevs(fMz,fMbz,fdz,fdzb,fzlev,fzblev); CHKERRQ(ierr);

  ierr = verbPrintf(4,grid.com,
    "\n  [entering enthalpyAndDrainageStep(); fMz = %d, fdz = %5.3f, fMbz = %d, fdzb = %5.3f]",
    fMz, fdz, fMbz, fdzb); CHKERRQ(ierr);

  PetscTruth viewOneColumn;
  ierr = check_option("-view_sys", viewOneColumn); CHKERRQ(ierr);

  EnthalpyConverter EC(&config);
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
  system.bed_thermal_rho = config.get("bedrock_thermal_density"); // bed_thermal.rho;
  system.bed_thermal_c   = config.get("bedrock_thermal_specific_heat_capacity"); // bed_thermal.c_p;
  system.bed_thermal_k   = config.get("bedrock_thermal_conductivity"); // bed_thermal.k;

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
  
  // now get map-plane coupler fields
  IceModelVec2 *pccTs, *pccsbt, *pccsbmf;
  if (atmosPCC != PETSC_NULL) {
    // call sets pccTs to point to IceModelVec2 with current surface temps
    ierr = atmosPCC->updateSurfTempAndProvide(
              grid.year, dtTempAge / secpera, &info_coupler, pccTs);
              CHKERRQ(ierr);
  } else {
    SETERRQ(1,"PISM ERROR: atmosPCC == PETSC_NULL");
  }
  if (oceanPCC != PETSC_NULL) {
    ierr = oceanPCC->updateShelfBaseTempAndProvide(
              grid.year, dt / secpera, &info_coupler, pccsbt);
              CHKERRQ(ierr);
    ierr = oceanPCC->updateShelfBaseMassFluxAndProvide(
              grid.year, dt / secpera, &info_coupler, pccsbmf);
              CHKERRQ(ierr);
  } else {
    SETERRQ(1,"PISM ERROR: oceanPCC == PETSC_NULL");
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

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      // for fine grid; this should *not* be replaced by call to grid.kBelowHeight()
      const PetscInt  ks = static_cast<PetscInt>(floor(H[i][j]/fdz));

      // enthalpy at boundaries
      const PetscScalar
          Enth_air       = EC.getEnthPermissive(Ts[i][j], 0.0, p_air ),
          // in theory we could have a water fraction at k=ks level, but for
          //   now there is no case where we have that:
          Enth_ks        = EC.getEnthPermissive(Ts[i][j], 0.0,
                                                EC.getPressureFromDepth(H[i][j] - fzlev[ks]) ),
          // at underside of ice shelf, set enthalpy to that of max liquid water
          //   temperate ice; probably does not make much difference anyway because of
          //   no upward/downward liquid water transport:
          Enth_shelfbase = EC.getEnthPermissive(Tshelfbase[i][j],
                                                config.get("liquid_water_fraction_max"), 
                                                EC.getPressureFromDepth(H[i][j]) );

      if (k0+ks>0) { // if there are enough points in bedrock&ice to bother ...
        ierr = system.setIndicesThisColumn(i,j,ks); CHKERRQ(ierr);
        ierr = Tb3.getValColumn(i,j,fMbz,fzblev,Tb); CHKERRQ(ierr);

        if (grid.vertical_spacing == EQUAL) {
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
        system.Enth_s[0] = EC.getEnthalpyCTS( EC.getPressureFromDepth(H[i][j]) );
        for (PetscInt k = 1; k <= ks; k++) {
          system.Enth_s[k] = EC.getEnthalpyCTS( EC.getPressureFromDepth(H[i][j] - fzlev[k]) );
          if (system.Enth[k] > system.Enth_s[k]) {
            // if there is a liquid water fraction we will switch to upwind; conductivity goes to zero
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
          system.Enth_b[k] = EC.getEnthBedrock(system.Enth[0], Tb[k0], Tb[k]);
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

#if 0
#define iSHOW 25
#define jSHOW 25

if ((i==iSHOW) && (j==jSHOW)) {
  ierr = verbPrintf(1,grid.com,
     "\n\n k0 = %d,  Tb[k0] = %12.5f, Tb[0] = %12.5f, Enth[0] = %12.5f, Enth_b[0] = %12.5f\n",
     k0,Tb[k0],Tb[0],system.Enth[0],system.Enth_b[0]); CHKERRQ(ierr);
  ierr = verbPrintf(1,grid.com,"\n\n SHOWING COLUMN (i,j)=(%d,%d) values:\n",i,j); CHKERRQ(ierr);
  ierr = system.viewColumnValues(NULL,fzblev, fMbz, "z values (elevation) in bedrock"); CHKERRQ(ierr);
  ierr = system.viewColumnValues(NULL,system.Enth, fMz, "ice enthalpy Enth[] before solve"); CHKERRQ(ierr);
  ierr = system.viewColumnValues(NULL,Tb, fMbz, "bedrock temperature Tb[] before solve"); CHKERRQ(ierr);
  ierr = system.viewColumnValues(NULL,system.Enth_b, fMbz, "bedrock enthalpy Enth_b[] before solve"); CHKERRQ(ierr);
}
#endif

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
            ierr = system.viewSystem(NULL,"system"); CHKERRQ(ierr);
            ierr = system.viewColumnValues(NULL, x, fMz+k0, "solution x"); CHKERRQ(ierr);
          }
        }


      }

      // prepare for melting/refreezing
      PetscScalar Hmeltnew = Hmelt[i][j];

      // top ice level; no possibility of drainage (a separate modeling issue!)
      Enthnew[ks] = x[k0 + ks];

      // insert solution for generic ice segments, including draining to base
      for (PetscInt k=ks-1; k >= 1; k--) {
        Enthnew[k] = x[k0 + k];
        // modifies last two arguments, generally:
        ierr = drainageToBaseModelEnth(EC, L, omega_max, H[i][j], fzlev[k], fdz,
                                       Enthnew[k], Hmeltnew); CHKERRQ(ierr);
      }

      // insert solution for ice/rock interface (or base of ice shelf) segment
      if (ks > 0) {
        Enthnew[0] = x[k0 + 0];
        if (PismModMask(mask[i][j]) != MASK_FLOATING) {
          // we only use the drainage model for the basal ice segment if 
          //   the ice is grounded; if the ice is floating then mass and energy balance
          //   is the responsibility of the PISMOceanCoupler
          // FIXME: must this be called to avoid unreasonable Enth values anyway?
          // modifies last two arguments, generally:
          ierr = drainageToBaseModelEnth(EC, L, omega_max, H[i][j], 0.0, fdz,
                                         Enthnew[0], Hmeltnew); CHKERRQ(ierr);
        }
      } else {
        Hmeltnew = 0.0; // no stored water if no ice present
        // Enthnew[0] = Enthnew[ks] already set
      }

      // bottom of ice is top of bedrock when grounded, so
      //   T(z=0) at top of bedrock should match enthalpy at z=0;
      //   when floating just match ocean temp provided by PISMOceanCoupler
      if (PismModMask(mask[i][j]) == MASK_FLOATING) { // top of bedrock sees ocean
          Tbnew[k0] = Tshelfbase[i][j];
      } else {
        if (ks > 0) { // grounded ice present
          Tbnew[k0] = EC.getAbsTemp(Enthnew[0], EC.getPressureFromDepth(H[i][j]) );
        } else {      // no significant ice; top of bedrock sees atmosphere
          Tbnew[k0] = Ts[i][j];
        }
      }

      // insert generic bedrock segments solution, which refers to z=0 segment
      for (PetscInt k=k0-1; k >= 0; k--) {
        //Tbnew[k] = EC.getAbsTempBedrock(Enthnew[0], Tbnew[k0], x[k]);
        Tbnew[k] = EC.getAbsTempBedrock(system.Enth[0], Tb[k0], x[k]);  // use same pt in E-T space as before
      }

      // transfer column into Tb3; neighboring columns will not reference so no need for communication
      ierr = Tb3.setValColumn(i,j,fMbz,fzblev,Tbnew); CHKERRQ(ierr);

      // now that enthalpy is known in top layer, check for (and correct) any extreme advection bulges
      for (PetscInt k=0; k < ks; k++) {
        if (Enthnew[k] < Enthnew[ks] - bulgeMaxEnth) {
          Enthnew[k] = Enthnew[ks] - bulgeMaxEnth;
          bulgeCount++;
        }
      }

      // set enthalpy to energy content of surface ice; this is a regularity issue not an atmosphere model!
      for (PetscInt k=ks+1; k<fMz; k++) {
        Enthnew[k] = Enth_air;
      }

      // transfer column into EnthNew3; communication later
      ierr = EnthNew3.setValColumnPL(i,j,fMz,fzlev,Enthnew); CHKERRQ(ierr);

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
        // eliminate basal lubrication water if floating
        // UNACCOUNTED MASS & ENERGY (LATENT) LOSS (TO OCEAN)!!
        Hmelt[i][j] = 0.0;
      } else {
        // limit Hmelt by default max
        // UNACCOUNTED MASS & ENERGY (LATENT) LOSS (TO INFINITY AND BEYOND)!!
        Hmelt[i][j] = PetscMin(Hmelt_max, Hmeltnew);
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
                PetscScalar thickness, PetscScalar z, const PetscScalar dz,
                PetscScalar &enthalpy, PetscScalar &Hmelt) {

  if (allowAboveMelting == PETSC_TRUE) {
    SETERRQ(1,"drainageToBaseModelEnth() called but allowAboveMelting==TRUE");
  }

  // no change to either enthalpy or basal layer thickness in this case
  if (updateHmelt == PETSC_FALSE)  return 0;

  const PetscScalar p     = EC.getPressureFromDepth(thickness - z),
                    omega = EC.getWaterFraction(enthalpy, p);

  // if there is liquid water already, thus temperate, consider whether there
  //   is enough to cause drainage
  if (omega > 1.0e-6) {
    const PetscScalar abovecap = omega - omega_max;
    if (abovecap > 0.0) {
      enthalpy -= abovecap * L;
      Hmelt    += abovecap * dz;   // ice-equivalent water thickness change
    }
    return 0; // done with temperate case
  }
  
  // if cold and in the basal layer, consider whether there is available 
  //   water to freeze on
  // FIXME: think about ocean case!?
  if ((z >= -1.0e-6) && (z <= 1.0e-6)) {
    // only consider freeze-on if column segment is at base of ice; E_s = getEnthalpyCTS(config, p)
    const PetscScalar dEnth_to_reach_temperate = EC.getEnthalpyCTS(p) - enthalpy;
    if (dEnth_to_reach_temperate > 0.0) {
      // if below E_s, then freeze on, and bring up enthalpy to E_s if enough water is available
      const PetscScalar dEnth_available = (Hmelt / dz) * L, // = ((rho Hmelt dx dy) * L) / (rho dx dy dz)
                        dEnth_added     = PetscMin(dEnth_available, dEnth_to_reach_temperate);
      enthalpy += dEnth_added;
      Hmelt    -= (dEnth_added * dz) / L;  // will cause negative basal melt rate,
                                           //   which can enter into mass continuity
    }
  }

  return 0;
}

