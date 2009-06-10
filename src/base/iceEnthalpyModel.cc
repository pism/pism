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

#define DEBUGVERB 2


static PetscScalar getPressureFromDepth(NCConfigVariable &config, PetscScalar depth) {
  const PetscScalar p_air = config.get("surface_pressure"); // Pa
  if (depth <= 0.0) { // at or above surface of ice
    return p_air;
  } else {
    const PetscScalar g     = config.get("earth_gravity"),
                      rho_i = config.get("ice_density");
    return p_air + rho_i * g * depth;
  }
}


//! Get enthalpy at phase transition endpoints, and pressure melting temperature, from pressure p.
/*! From \ref AschwandenBlatter2009,
     \f[ T_m(p) = T_0 - \beta p, \f]
     \f[ H_l(p) = c_w T_m(p), \f]
     \f[ H_s(p) = -L + H_l(p). \f]
Also returns H_s.
 */ 
static PetscScalar getEnthalpyInterval(NCConfigVariable &config, PetscScalar p, 
                                       PetscScalar &T_m, PetscScalar &H_l, PetscScalar &H_s) {
  const PetscScalar T_0  = config.get("water_melting_temperature"),    // K
                    beta = config.get("beta_CC"),                      // K Pa-1
                    c_w  = config.get("water_specific_heat_capacity"), // J kg-1 K-1
                    L    = config.get("water_latent_heat_fusion");     // J kg-1
  T_m = T_0 - beta * p;
  H_l = c_w * T_m;
  H_s = - L + H_l;
  return H_s;
}


//! Get absolute ice temperature (K) from enthalpy H and pressure p.
/*! From \ref AschwandenBlatter2009,
     \f[ T=T(H,p) = \begin{cases} 
                       c_i^{-1} (H-H_s(p)) + T_m(p), & H < H_s(p), \\
                       T_m(p), &                       H_s(p) \le H < H_l(p).
                    \end{cases} \f]

We do not allow liquid water (i.e. water fraction \f$\omega=1.0\f$) so we fail if
\f$H \ge H_l(p)\f$.
  
See getEnthalpyInterval() for calculation of \f$T_m(p)\f$, \f$H_l(p)\f$, and \f$H_s(p)\f$. 
 */
static PetscScalar getAbsTemp(NCConfigVariable &config, PetscScalar H, PetscScalar p) {
  PetscScalar T_m, H_l, H_s;
  getEnthalpyInterval(config, p, T_m, H_l, H_s);
  // implement T part of eqn (12) in AB2009, but bonk if liquid water
  if (H < H_s) {
    const PetscScalar c_i = config.get("ice_specific_heat_capacity");   // J kg-1 K-1
    return ((H - H_s) / c_i) + T_m;
  } else if (H < H_l) { // two cases in (12)
    return T_m;
  } else {
    PetscPrintf(PETSC_COMM_WORLD,
      "\n\n\n  PISM ERROR in getAbsTemp():\n"
            "    enthalpy equals or exceeds that of liquid water; ending ... \n\n");
    PetscEnd();
    return T_m;
  }
}


//! Get pressure-adjusted ice temperature (K) from enthalpy H and pressure p.
/*! See getAbsTemp(), and compare:
     \f[ T_pa = T_pa(H,p) = \begin{cases} 
                       c_i^{-1} (H-H_s(p)) + T_0, & H < H_s(p), \\
                       T_0,                       & H_s(p) \le H < H_l(p).
                    \end{cases} \f]
 */
static PetscScalar getPATemp(NCConfigVariable &config, PetscScalar H, PetscScalar p) {
  PetscScalar T_m, H_l, H_s;
  getEnthalpyInterval(config, p, T_m, H_l, H_s);
  const PetscScalar T_0  = config.get("water_melting_temperature");     // K
  if (H < H_s) {
    const PetscScalar c_i = config.get("ice_specific_heat_capacity");   // J kg-1 K-1
    return ((H - H_s) / c_i) + T_0;
  } else if (H < H_l) { // two cases in (12)
    return T_0;
  } else {
    PetscPrintf(PETSC_COMM_WORLD,
      "\n\n\n  PISM ERROR in getPATemp():\n"
            "    enthalpy equals or exceeds that of liquid water; ending ... \n\n");
    PetscEnd();
    return T_0;
  }
}


//! Get liquid water fraction from enthalpy H and pressure p.
/*! From \ref AschwandenBlatter2009,
   \f[ \omega=\omega(H,p) = \begin{cases}
                               0.0,            & H \le H_s(p), \\
                               (H-H_s(p)) / L, & H_s(p) < H < H_l(p).
                            \end{cases} \f]

We do not allow liquid water (i.e. water fraction \f$\omega=1.0\f$) so we fail if
\f$H \ge H_l(p)\f$.
  
See getEnthalpyInterval() for calculation of \f$T_m(p)\f$, \f$H_l(p)\f$, and \f$H_s(p)\f$. 
 */
static PetscScalar getWaterFraction(NCConfigVariable &config, PetscScalar H, PetscScalar p) {
  PetscScalar T_m, H_l, H_s;
  getEnthalpyInterval(config, p, T_m, H_l, H_s);
  // implement omega part of eqn (12) in AB2009, but bonk if liquid water
  if (H <= H_s) { // two cases in (12)
    return 0.0;
  } else if (H < H_l) {
    const PetscScalar L = config.get("water_latent_heat_fusion");     // J kg-1
    return (H - H_s) / L;
  } else {
    PetscPrintf(PETSC_COMM_WORLD,
      "\n\n\n  PISM ERROR in getWaterFraction():\n"
            "    enthalpy equals or exceeds that of liquid water; ending ... \n\n");
    PetscEnd();
    return 1.0;
  }
}


//! Compute enthalpy from absolute temperature, liquid water fraction, and pressure p.
static PetscScalar getEnth(NCConfigVariable &config, PetscScalar T, PetscScalar omega, PetscScalar p) {
  if ((omega < 0.0) || (1.0 < omega)) {
    PetscPrintf(PETSC_COMM_WORLD,
      "\n\n\n  PISM ERROR in getEnth(): water fraction omega not in range [0,1]; ending ... \n\n");
    PetscEnd();
  }
  const PetscScalar T_0 = config.get("water_melting_temperature");    // K
  if (T > T_0 + 0.000001) {
    PetscPrintf(PETSC_COMM_WORLD,
      "\n\n\n  PISM ERROR in getEnth(): T exceeds T_0 so we have liquid water; ending ... \n\n");
    PetscEnd();
  }
  PetscScalar T_m, H_l, H_s;
  getEnthalpyInterval(config, p, T_m, H_l, H_s);
  const PetscScalar c_w = config.get("water_specific_heat_capacity"), // J kg-1 K-1
                    c_i = config.get("ice_specific_heat_capacity");   // J kg-1 K-1
  const PetscScalar c = (1.0 - omega) * c_i + omega * c_w;
  return H_s + c * (T - T_0);
}


/*! 
This constructor just sets flow law factor for nonzero water content, from 
\ref AschwandenBlatter2009 and \ref LliboutryDuval1985.
 */
PolyThermalGPBLDIce::PolyThermalGPBLDIce(MPI_Comm c,const char pre[]) : ThermoGlenIce(c,pre) {
  config = NULL;
  water_frac_coeff = 184.0;   // FIXME:  should also come through config interface
}


PetscErrorCode PolyThermalGPBLDIce::view(PetscViewer viewer) const {
  PetscErrorCode ierr;
  PetscTruth iascii;
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(comm,&viewer);CHKERRQ(ierr);
  }
  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii);CHKERRQ(ierr);
  if (iascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"<\nPolyThermalGPBLDIce object (%s)\n",prefix);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  v_schoof=%4f m/a L_schoof=%4f km\n",
                                  schoofVel*secpera,schoofLen/1e3);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  water_frac_coeff=%4f\n>\n",water_frac_coeff);CHKERRQ(ierr);
  } else {
    SETERRQ(1,"No binary viewer for this object\n");
  }
  return 0;
}


//! The factor in the Paterson-Budd-Lliboutry-Duval flow law.  For constitutive law form.
PetscScalar PolyThermalGPBLDIce::softnessParameterFromEnth(PetscScalar enthalpy, PetscScalar pressure) const {
  if (config == NULL) { 
    PetscPrintf(PETSC_COMM_WORLD,"config ptr is NULL in PolyThermalGPBLDIce::flowFromEnth()... ending\n");
    PetscEnd();
  }
  PetscScalar T_m, H_l, H_s;
  getEnthalpyInterval(*config, pressure, T_m, H_l, H_s);
  if (enthalpy <= H_s) {       // cold ice
    return softnessParameter( getPATemp(*config,enthalpy,pressure) ); // uses ThermoGlenIce formula
  } else if (enthalpy < H_l) { // temperate ice
    const PetscScalar T_0  = config->get("water_melting_temperature"),    // K
                      omega = getWaterFraction(*config,enthalpy,pressure);
    // next line implements eqn (23) in \ref AschwandenBlatter2009
    return softnessParameter(T_0) * (1.0 + water_frac_coeff * omega);  // uses ThermoGlenIce formula 
  } else { // liquid water not allowed
    PetscPrintf(PETSC_COMM_WORLD,
      "\n\n\n  PISM ERROR in PolyThermalGlenPBLDIce::flow(): liquid water not allowed; ending ... \n\n");
    PetscEnd();
    return 0.0;
  }
}


//! The factor in the Paterson-Budd-Lliboutry-Duval flow law.  For viscosity form.
PetscScalar PolyThermalGPBLDIce::hardnessParameterFromEnth(PetscScalar enthalpy, PetscScalar pressure) const {
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
  if (config == NULL) { 
    PetscPrintf(PETSC_COMM_WORLD,
       "config ptr is NULL in PolyThermalGPBLDIce::effectiveViscosityColumnFromEnth()... ending\n");
    PetscEnd();
  }

  // DESPITE NAME, does *not* return effective viscosity.
  // The result is \nu_e H, i.e. viscosity times thickness.
  // B is really hardness times thickness.

  // Integrate the hardness parameter using the trapezoid rule.
  PetscScalar B = 0;
  if (kbelowH > 0) {
    PetscScalar dz = zlevels[1] - zlevels[0];
    B += 0.5 * dz * hardnessParameterFromEnth( 0.5 * (enthalpy1[0] + enthalpy2[0]),
                                               getPressureFromDepth(*config, thickness) );
    for (PetscInt m=1; m < kbelowH; m++) {
      const PetscScalar dzNEXT = zlevels[m+1] - zlevels[m],
                        depth  = thickness - 0.5 * (zlevels[m+1] + zlevels[m]);
      B += 0.5 * (dz + dzNEXT) * hardnessParameterFromEnth( 0.5 * (enthalpy1[m] + enthalpy2[m]),
                                                            getPressureFromDepth(*config, depth) );
      dz = dzNEXT;
    }
    // use last dz from for loop
    const PetscScalar depth  = 0.5 * (thickness - zlevels[kbelowH]);
    B += 0.5 * dz * hardnessParameterFromEnth( 0.5 * (enthalpy1[kbelowH] + enthalpy2[kbelowH]),
                                               getPressureFromDepth(*config, depth) );
  }
  const PetscScalar alpha = secondInvariant(u_x, u_y, v_x, v_y);
  return 0.5 * B * pow(schoofReg + alpha, (1-n)/(2*n));
}



/*********** for registering new kind of ice with IceFactory ****************/

#define ICE_GPBLD      "gpbld"    

static PetscErrorCode create_gpbld(MPI_Comm comm,const char pre[],IceType **i) {
  // *i = new (PolyThermalGPBLDIce)(comm,pre);  return 0;
  *i = new (ThermoGlenIce)(comm,pre);  return 0; // FIXME:  for now this is just old flow law
}


/*********** procedures for init ****************/

IceEnthalpyModel::IceEnthalpyModel(IceGrid &g) : IceModel(g) {
  doColdIceTemperatureStep = true;   // for start, default to no actual enthalpy computation;
                                     // just read and write additional enthalpy field to and from file
}


PetscErrorCode IceEnthalpyModel::createVecs() {
  PetscErrorCode ierr;

  ierr = Enth3.create(grid, "enthalpy", true); CHKERRQ(ierr);
  // PROPOSED standard name = land_ice_enthalpy
  ierr = Enth3.set_attrs(
     "model_state",
     "ice enthalpy (sensible heat plus latent heat of liquid fraction)",
     "J kg-1", 
     ""); CHKERRQ(ierr);

  ierr = IceModel::createVecs(); CHKERRQ(ierr);
  
  // see IceModel::allocate_internal_objects(), which is where this should go
  ierr = EnthNew3.create(grid,"enthalpy_new",false); CHKERRQ(ierr);
  ierr = EnthNew3.set_attrs(
     "internal",
     "ice enthalpy; temporary during update",
     "J kg-1",
     ""); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceEnthalpyModel::init_physics() {
  PetscErrorCode ierr;

  // let the base class create the ice and process its options:
  ierr = IceModel::init_physics(); CHKERRQ(ierr);

  ierr = iceFactory.registerType(ICE_GPBLD, &create_gpbld);
  if (ierr != 0) {
    PetscPrintf(grid.com,
       "FAILURE OF iceFactory.registerType() ... return value %d ... ending ....\n",ierr);
    PetscEnd();
  }
  CHKERRQ(ierr);

  if (ice != NULL)  delete ice;  // kill choice already made!

  iceFactory.setType(ICE_GPBLD);    // new flowlaw which has dependence on enthalpy not temperature
  iceFactory.create(&ice);

#if 0
  PolyThermalGPBLDIce *gpbldi = dynamic_cast<PolyThermalGPBLDIce*>(ice);
  if (gpbldi) {
    gpbldi->config = &config;
  } else {
    PetscPrintf(grid.com,
       "not using PolyThermalGPBLDIce ... ending ....\n");
    PetscEnd();
  }
  ierr = ice->printInfo(1);CHKERRQ(ierr); // DEBUG
#endif

  ierr = ice->setFromOptions();CHKERRQ(ierr);

  return 0;
}



/*********** procedures for read/write ****************/

PetscErrorCode IceEnthalpyModel::write_extra_fields(const char filename[]) {
  PetscErrorCode ierr;
  
  if (doColdIceTemperatureStep) { // in this case, just update Enth3 to reflect 
                                  // temperature in ice at final time
    ierr = verbPrintf(DEBUGVERB, grid.com, 
      "  using temperature to set enthalpy for writing (as cold ice) ...\n");
      CHKERRQ(ierr);
    ierr = setEnth3FromT3_ColdIce(); CHKERRQ(ierr);
  }
  ierr = Enth3.write(filename, NC_DOUBLE); CHKERRQ(ierr);

  // also write omega = liquid water fraction; we use EnthNew3 as temporary space,
  //   for this purpose
  ierr = verbPrintf(DEBUGVERB, grid.com, 
      "  writing liquid water fraction 'liquid_frac' from enthalpy ...\n"); CHKERRQ(ierr);
  ierr = EnthNew3.set_name("liquid_frac"); CHKERRQ(ierr);
  ierr = EnthNew3.set_attrs(
     "diagnostic",
     "liquid water fraction in ice; 0 <= omega <= 1",
     "", 
     ""); CHKERRQ(ierr);
  PetscScalar **thickness;
  PetscScalar *omegaij, *Enthij; // columns of these values
  ierr = EnthNew3.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = vH.get_array(thickness); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = EnthNew3.getInternalColumn(i,j,&omegaij); CHKERRQ(ierr);
      ierr = Enth3.getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
      for (PetscInt k=0; k<grid.Mz; ++k) {
        const PetscScalar depth = thickness[i][j] - grid.zlevels[k];
        omegaij[k] = getWaterFraction(config, Enthij[k], getPressureFromDepth(config, depth));
      }
    }
  }
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = EnthNew3.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = EnthNew3.write(filename, NC_DOUBLE); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceEnthalpyModel::initFromFile(const char *fname) {
  PetscErrorCode  ierr;

  ierr = IceModel::initFromFile(fname); CHKERRQ(ierr);

  ierr = verbPrintf(DEBUGVERB, grid.com, 
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
    ierr = verbPrintf(DEBUGVERB, grid.com, 
      "  variable 'enthalpy' not found so setting it as cold ice, from temperature ...\n");
      CHKERRQ(ierr);
    ierr = setEnth3FromT3_ColdIce(); CHKERRQ(ierr);
  }
    
  return 0;
}



/*********** setting fields ****************/

PetscErrorCode IceEnthalpyModel::setEnth3FromT3_ColdIce() {
  PetscErrorCode ierr;

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
        if (depth > 0.0) { // in ice
          Enthij[k] = getEnth(config,Tij[k],0.0,getPressureFromDepth(config,depth));
        } else {
          Enthij[k] = 0.0;  // set enthalpy in air to zero
        }
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


/*********** timestep routines ****************/

PetscErrorCode IceEnthalpyModel::temperatureStep(
     PetscScalar* vertSacrCount, PetscScalar* bulgeCount) {
  PetscErrorCode ierr;
  if (doColdIceTemperatureStep) {
    ierr = verbPrintf(DEBUGVERB,grid.com,
      "     IceEnthalpyModel::temperatureStep(): CALLING IceModel::temperatureStep()\n"); CHKERRQ(ierr);
    ierr = IceModel::temperatureStep(vertSacrCount,bulgeCount);  CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(DEBUGVERB,grid.com,
      "     IceEnthalpyModel::temperatureStep(): CALLING IceEnthalpyModel::enthalpyStep()\n"); CHKERRQ(ierr);
    // new enthalpy values go in EnthNew3; also updates (and communicates) Hmelt
    ierr = enthalpyStep(vertSacrCount,bulgeCount);  CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode IceEnthalpyModel::enthalpyStep(PetscScalar* vertSacrCount, PetscScalar* bulgeCount) {

  //FIXME: introduce serious code duplication from IceModel::temperatureStep()
  verbPrintf(1,grid.com,
    "\n\n    IceEnthalpyModel::enthalpyStep(): NOT IMPLEMENTED ... ending\n");
  PetscEnd();
  return 0;
}


PetscErrorCode IceEnthalpyModel::temperatureAgeStep() {
  PetscErrorCode  ierr;

  ierr = verbPrintf(DEBUGVERB,grid.com,
    "\n  [IceEnthalpyModel::temperatureAgeStep():  ENTERING; DOING IceModel::temperatureAgeStep() FIRST\n");
    CHKERRQ(ierr);
  
  ierr = IceModel::temperatureAgeStep(); CHKERRQ(ierr);

  if (doColdIceTemperatureStep) {
    ierr = verbPrintf(DEBUGVERB,grid.com,
      "   IceEnthalpyModel::temperatureAgeStep(): ENTHALPY IS OFF.  DONE.]\n"); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(DEBUGVERB,grid.com,
      "   IceEnthalpyModel::temperatureAgeStep(): ENTHALPY IS ON.  COMMUNICATING ENTHALPY]\n"); CHKERRQ(ierr);

    // start & complete communication
    ierr = Enth3.beginGhostCommTransfer(EnthNew3); CHKERRQ(ierr);
    ierr = Enth3.endGhostCommTransfer(EnthNew3); CHKERRQ(ierr);
  }
  return 0;
}

