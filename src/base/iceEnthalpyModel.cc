// Copyright (C) 2009 Andreas Aschwandend and Ed Bueler
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

IceEnthalpyModel::IceEnthalpyModel(IceGrid &g) : IceModel(g) {
}


PetscErrorCode IceEnthalpyModel::createVecs() {
  PetscErrorCode ierr;

  ierr = Enth3.create(grid, "enthalpy", true); CHKERRQ(ierr);
  // PROPOSED standard name = land_ice_enthalpy
  ierr = Enth3.set_attrs("model_state",
                         "ice enthalpy (sensible and latent heat content per unit volume)",
		         "J m-3", 
		         ""); CHKERRQ(ierr);
  ierr = Enth3.set_attr("valid_min", 0.0); CHKERRQ(ierr);

  ierr = IceModel::createVecs(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceEnthalpyModel::write_extra_fields(const char filename[]) {
  PetscErrorCode ierr;
  ierr = Enth3.write(filename, NC_DOUBLE); CHKERRQ(ierr);
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
    ierr = setEnth3FromTemp_ColdIce(); CHKERRQ(ierr);
  }
    
  return 0;
}


PetscErrorCode IceEnthalpyModel::setEnth3FromTemp_ColdIce() {
  PetscErrorCode ierr;

  PetscScalar **surfelev, *Tij, *Enthij;
  ierr = vH.get_array(surfelev); CHKERRQ(ierr);
  ierr = T3.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = T3.getInternalColumn(i,j,&Tij); CHKERRQ(ierr);
      ierr = Enth3.getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
      for (PetscInt k=0; k<grid.Mz; ++k) {
        Enthij[k] = getAbsTemp(Tij[k],surfelev[i][j] - grid.zlevels[k]);
      }
    }
  }
  
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = T3.end_access(); CHKERRQ(ierr);
  ierr = vh.end_access(); CHKERRQ(ierr);

  ierr = Enth3.beginGhostComm(); CHKERRQ(ierr);
  ierr = Enth3.endGhostComm(); CHKERRQ(ierr);
  return 0;
}


PetscScalar IceEnthalpyModel::getAbsTemp(PetscScalar enth, PetscScalar depth) {
  const PetscScalar Tm     = ice->meltingTemp - ice->beta_CC_grad * depth,
                    c_w    = 2.008e3, // should be for pure water
                    Hl_atp = c_w * (Tm - ice->meltingTemp), // generally negative; in Celcius
                    Hs     = - ice->latentHeat + Hl_atp;
  return (1.0/ice->c_p) * (enth - Hs) + Tm;  // eqn 5.12 in AB09
}

