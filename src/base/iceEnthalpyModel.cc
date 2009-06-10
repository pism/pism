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


// FIXME:  NEED AN ICETYPE WHICH DEPENDS ON WATER AND COMBINES PATERSON-BUDD WITH EQN
//         (23) IN AB09; IGNORES EQN (22) IN AB09S


/*********** procedures for init ****************/

IceEnthalpyModel::IceEnthalpyModel(IceGrid &g) : IceModel(g) {
  doColdIceTemperatureStep = true;   // for start, default to no actual enthalpy computation;
                                     // just read and write additional enthalpy field to and from file
}


PetscErrorCode IceEnthalpyModel::createVecs() {
  PetscErrorCode ierr;

  ierr = Enth3.create(grid, "enthalpy", true); CHKERRQ(ierr);
  // PROPOSED standard name = land_ice_enthalpy
  ierr = Enth3.set_attrs("model_state",
                         "ice enthalpy (sensible heat per mass plus latent heat content of liquid fraction)",
		         "J kg-1", 
		         ""); CHKERRQ(ierr);

  ierr = IceModel::createVecs(); CHKERRQ(ierr);
  
  // see IceModel::allocate_internal_objects()
  ierr = EnthNew3.create(grid,"enthalpy_new",false); CHKERRQ(ierr);
  ierr = EnthNew3.set_attrs("internal",
                            "ice enthalpy; temporary during update",
                            "J kg-1",
                            ""); CHKERRQ(ierr);

  return 0;
}


/*********** procedures for read/write ****************/

PetscErrorCode IceEnthalpyModel::write_extra_fields(const char filename[]) {
  PetscErrorCode ierr;
  ierr = Enth3.write(filename, NC_DOUBLE); CHKERRQ(ierr);
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
    ierr = setEnth3FromTemp_ColdIce(); CHKERRQ(ierr);
  }
    
  return 0;
}


/*********** scalar functions ****************/

PetscScalar IceEnthalpyModel::getPressureFromDepth(PetscScalar depth) {
  const PetscScalar p_air = config.get("surface_pressure"); // Pa
  if (depth <= 0.0) { // at or above surface of ice
    return p_air;
  } else {
    const PetscScalar g     = config.get("earth_gravity"),
                      rho_i = config.get("ice_density");
    return p_air + rho_i * g * depth;
  }
}


PetscScalar IceEnthalpyModel::get_H_s(PetscScalar p, 
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


PetscScalar IceEnthalpyModel::getAbsTemp(PetscScalar H, PetscScalar p) {
  PetscScalar T_m, H_l, H_s;
  get_H_s(p, T_m, H_l, H_s);
  // implement T part of eqn (12) in AB2009, but bonk if liquid water
  if (H < H_s) {
    const PetscScalar c_i = config.get("ice_specific_heat_capacity");   // J kg-1 K-1
    return ((H - H_s) / c_i) + T_m;
  } else if (H < H_l) { // two cases in (12)
    return T_m;
  } else {
    PetscPrintf(grid.com,
      "\n\n\n  PISM ERROR in getAbsTemp():\n"
            "    enthalpy equals or exceeds that of liquid water; ending ... \n\n");
    PetscEnd();
    return T_m;
  }
}


PetscScalar IceEnthalpyModel::getWaterFraction(PetscScalar H, PetscScalar p) {
  PetscScalar T_m, H_l, H_s;
  get_H_s(p, T_m, H_l, H_s);
  // implement omega part of eqn (12) in AB2009, but bonk if liquid water
  if (H <= H_s) { // two cases in (12)
    return 0.0;
  } else if (H < H_l) {
    const PetscScalar L = config.get("water_latent_heat_fusion");     // J kg-1
    return (H - H_s) / L;
  } else {
    PetscPrintf(grid.com,
      "\n\n\n  PISM ERROR in getWaterFraction():\n"
            "    enthalpy equals or exceeds that of liquid water; ending ... \n\n");
    PetscEnd();
    return 1.0;
  }
}


PetscScalar IceEnthalpyModel::getEnth(PetscScalar T, PetscScalar omega, PetscScalar p) {
  if ((omega < 0.0) || (1.0 < omega)) {
    PetscPrintf(grid.com,
      "\n\n\n  PISM ERROR in getEnth(): water fraction omega not in range [0,1]; ending ... \n\n");
    PetscEnd();
  }
  const PetscScalar T_0 = config.get("water_melting_temperature");    // K
  if (T > T_0 + 0.000001) {
    PetscPrintf(grid.com,
      "\n\n\n  PISM ERROR in getEnth(): T exceeds T_0 so we have liquid water; ending ... \n\n");
    PetscEnd();
  }
  PetscScalar T_m, H_l, H_s;
  get_H_s(p, T_m, H_l, H_s);
  const PetscScalar c_w = config.get("water_specific_heat_capacity"), // J kg-1 K-1
                    c_i = config.get("ice_specific_heat_capacity");   // J kg-1 K-1
  const PetscScalar c = (1.0 - omega) * c_i + omega * c_w;
  return H_s + c * (T - T_0);
}


PetscScalar IceEnthalpyModel::getEnth_pa(PetscScalar T_pa, PetscScalar omega, PetscScalar p) {
  if ((omega < 0.0) || (1.0 < omega)) {
    PetscPrintf(grid.com,
      "\n\n\n  PISM ERROR in getEnth(): water fraction omega not in range [0,1]; ending ... \n\n");
    PetscEnd();
  }
  PetscScalar T_m, H_l, H_s;
  get_H_s(p, T_m, H_l, H_s);
  if (T_pa > T_m + 0.000001) {
    PetscPrintf(grid.com,
      "\n\n\n  PISM ERROR in getEnth_pa(): T_pa exceeds T_m so we have liquid water; ending ... \n\n");
    PetscEnd();
  }
  const PetscScalar c_w = config.get("water_specific_heat_capacity"), // J kg-1 K-1
                    c_i = config.get("ice_specific_heat_capacity");   // J kg-1 K-1
  const PetscScalar c = (1.0 - omega) * c_i + omega * c_w;
  return H_s + c * (T_pa - T_m);
}


/*********** setting fields ****************/

PetscErrorCode IceEnthalpyModel::setEnth3FromTemp_ColdIce() {
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
        const PetscScalar z = grid.zlevels[k];
        PetscScalar depth = (z < H[i][j]) ? H[i][j] - z : 0.0;
        Enthij[k] = getEnth(Tij[k],0.0,getPressureFromDepth(depth));
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

  //FIXME: introduce serious code duplication from IceModel::temperatureStep(), but actually computes
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

