// Copyright (C) 2008--2011 Ed Bueler and Constantine Khroulev
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

#include "iceCalvBCModel.hh"


IceCalvBCModel::IceCalvBCModel(IceGrid &g, NCConfigVariable &conf, NCConfigVariable &conf_overrides, int mytest)
  : IceExactSSAModel(g, conf, conf_overrides, mytest) {
}


IceCalvBCModel::~IceCalvBCModel() {}

PetscErrorCode IceCalvBCModel::createVecs() {
  PetscErrorCode ierr;

  ierr = IceExactSSAModel::createVecs(); CHKERRQ(ierr);
  
  // additionally create this stuff
  ierr = vsmoothCFmask.create(grid, "smoothCFmask", true); CHKERRQ(ierr); // stencil = BOX
  ierr = vsmoothCFmask.set_attrs(  // no pism_intent; no standard_name
           "NONE", "smoothed version of (thk>0) mask, for calving front normal direction", 
           "", ""); CHKERRQ(ierr);  
  ierr = vnCF[0].create(grid, "nCF_x", true); CHKERRQ(ierr);
  ierr = vnCF[0].set_attrs(  // no pism_intent; no standard_name
           "NONE", "x component of calving front normal direction", 
           "m-1", ""); CHKERRQ(ierr);  
  ierr = vnCF[1].create(grid, "nCF_y", true); CHKERRQ(ierr);
  ierr = vnCF[1].set_attrs(  // no pism_intent; no standard_name
           "NONE", "y component of calving front normal direction", 
           "m-1", ""); CHKERRQ(ierr);  
  return 0;
}

PetscErrorCode IceCalvBCModel::compute_nCF() {
  PetscErrorCode ierr;

  PetscScalar     **cfmask, **ncf[2], **H;

  // smoothCFmask is smoothed version of (H>100) mask
  ierr = vsmoothCFmask.get_array(cfmask); CHKERRQ(ierr);
  ierr = vH.get_array(H); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (H[i][j] <= 0.001)       cfmask[i][j] = 0.0;
      else if  (H[i][j] >= 100.0) cfmask[i][j] = 1.0;
      else                        cfmask[i][j] = H[i][j] / 100.0;
    }
  }
  ierr = vsmoothCFmask.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  // smooth by iterated averaging
  IceModelVec2S vnCFOLD = vnCF[0];
  PetscScalar **cfOLD;
  const PetscInt smoothStages = 3;
  for (PetscInt k=0; k<smoothStages; ++k) {
    // communicate/copy to get valid ghost, into "OLD" copy
    ierr = vsmoothCFmask.beginGhostComm(vnCFOLD); CHKERRQ(ierr);
    ierr = vsmoothCFmask.endGhostComm(vnCFOLD); CHKERRQ(ierr);
    // now smooth
    ierr = vsmoothCFmask.get_array(cfmask); CHKERRQ(ierr);
    ierr = vnCFOLD.get_array(cfOLD); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        cfmask[i][j] = 
           0.075 * cfOLD[i-1][j+1] + 0.15 * cfOLD[i][j+1] + 0.075 * cfOLD[i+1][j+1] 
         + 0.15  * cfOLD[i-1][j]   + 0.1  * cfOLD[i][j]   + 0.15  * cfOLD[i+1][j] 
         + 0.075 * cfOLD[i-1][j-1] + 0.15 * cfOLD[i][j-1] + 0.075 * cfOLD[i+1][j-1]; 
      }
    }
    ierr = vsmoothCFmask.end_access(); CHKERRQ(ierr);
    ierr = vnCFOLD.end_access(); CHKERRQ(ierr);
  }

  // communicate one more time for ghosts for gradient calculation
  ierr = vsmoothCFmask.beginGhostComm(); CHKERRQ(ierr);
  ierr = vsmoothCFmask.endGhostComm(); CHKERRQ(ierr);

  // calving front normal direction is negative of gradient of smoothCFmask:
  ierr = vsmoothCFmask.get_array(cfmask); CHKERRQ(ierr);
  ierr = vnCF[0].get_array(ncf[0]); CHKERRQ(ierr);
  ierr = vnCF[1].get_array(ncf[1]); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ncf[0][i][j] = - (cfmask[i+1][j] - cfmask[i-1][j]) / (2.0 * grid.dx);
      ncf[1][i][j] = - (cfmask[i][j+1] - cfmask[i][j-1]) / (2.0 * grid.dy);
    }
  } 
  ierr = vsmoothCFmask.end_access(); CHKERRQ(ierr);
  ierr = vnCF[0].end_access(); CHKERRQ(ierr);
  ierr = vnCF[1].end_access(); CHKERRQ(ierr);

  // it is reasonable to scale to improve condition number; what scaling?
//  ierr = vsmoothCFmask.scale(1.0); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceCalvBCModel::writeCFfields(const char* default_filename) {
  PetscErrorCode ierr;
  
  char filename[PETSC_MAX_PATH_LEN];
  PetscBool o_set;
  ierr = PetscOptionsGetString(PETSC_NULL, "-o", filename, PETSC_MAX_PATH_LEN, &o_set); CHKERRQ(ierr);

  // Use the default if the output file name was not given:
  if (!o_set)
    strncpy(filename, default_filename, PETSC_MAX_PATH_LEN);

  ierr = vsmoothCFmask.write(filename, NC_DOUBLE); CHKERRQ(ierr);
  ierr = vnCF[0].write(filename, NC_DOUBLE); CHKERRQ(ierr);
  ierr = vnCF[1].write(filename, NC_DOUBLE); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceCalvBCModel::assembleSSAMatrix(const bool includeBasalShear,
                                 IceModelVec2S vNuH[2], Mat A) {
  const PetscInt  Mx=grid.Mx, My=grid.My, M=2*My;
  const PetscScalar   dx=grid.dx, dy=grid.dy;
  const PetscScalar   one = 1.0;
  PetscErrorCode  ierr;
  PetscScalar     **nuH[2], **tauc, **cfmask;
  PISMVector2 **uv;

  // clear it out
  ierr = MatZeroEntries(A); CHKERRQ(ierr);

  PetscReal beta_shelves_drag_too = config.get("beta_shelves_drag_too");

  /* matrix assembly loop */
  ierr = vMask.begin_access();  CHKERRQ(ierr);
  ierr = vtauc.get_array(tauc); CHKERRQ(ierr);
  ierr = vel_ssa.get_array(uv); CHKERRQ(ierr);
  ierr = vNuH[0].get_array(nuH[0]); CHKERRQ(ierr);
  ierr = vNuH[1].get_array(nuH[1]); CHKERRQ(ierr);
  ierr = vsmoothCFmask.get_array(cfmask); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscInt J = 2*j;
      const PetscInt rowU = i*M + J;
      const PetscInt rowV = i*M + J+1;
      if (vMask.value(i,j) == MASK_SHEET) {
        // set diagonal entry to one; RHS entry will be known (e.g. SIA) velocity;
        //   this is where Dirichlet boundary value to SSA is set
        ierr = MatSetValues(A, 1, &rowU, 1, &rowU, &one, INSERT_VALUES); CHKERRQ(ierr);
        ierr = MatSetValues(A, 1, &rowV, 1, &rowV, &one, INSERT_VALUES); CHKERRQ(ierr);
      } else {
        const PetscInt im = (i + Mx - 1) % Mx,
                       ip = (i + 1) % Mx,
                       Jm = 2 * ((j + My - 1) % My),
                       Jp = 2 * ((j + 1) % My);
        const PetscScalar dx2 = dx*dx, 
                          d4 = dx*dy*4,
                          dy2 = dy*dy;
        /* Provide shorthand for the following staggered coefficients
        *      c11
        *  c00     c01
        *      c10
        * Note that the positive i (x) direction is right and the positive j (y)
        * direction is up.  Thickness H is incorporated here because 
        * nuH[][][] is actually viscosity times thickness. */
        const PetscScalar c00 = nuH[0][i-1][j],
                          c01 = nuH[0][i][j],
                          c10 = nuH[1][i][j-1],
                          c11 = nuH[1][i][j];
        const PetscInt stencilSize = 13;
        /* The locations of the stencil points for the U equation */
        const PetscInt colU[] = {
          /*       */ i*M+Jp,                 // U neighbors: x-dx,y+dy   x,y+dy   x+dx,y+dy
          im*M+Jp+1,  i*M+Jp+1,   ip*M+Jp+1,  // V neighbors:     "        "           "
          im*M+J,     i*M+J,      ip*M+J,     // U neighbors: x-dx,y      x,y      x+dx,y
          im*M+J+1,               ip*M+J+1,   // V neighbors:     "        "           "
          /*       */ i*M+Jm,                 // U neighbors: x-dx,y-dy   x,y-dy   x+dx,y-dy
          im*M+Jm+1,  i*M+Jm+1,   ip*M+Jm+1}; // V neighbors:     "        "           "
        /* The locations of the stencil points for the V equation */
        const PetscInt colV[] = {
          im*M+Jp,        i*M+Jp,     ip*M+Jp,// < same scheme as for colU >
          /*           */ i*M+Jp+1,
          im*M+J,                     ip*M+J,
          im*M+J+1,       i*M+J+1,    ip*M+J+1,
          im*M+Jm,        i*M+Jm,     ip*M+Jm,
          /*           */ i*M+Jm+1 };
 
        if ( (vMask.value(i,j) == MASK_FLOATING)
             && (cfmask[i][j] > 0.01) && (cfmask[i][j] < 0.99) ) {
          // FIXME: put in the correct vals
          PetscScalar valU[] = {
            /*               */ -c11/dy2,
            (2*c00+c11)/d4,     2*(c00-c01)/d4,                 -(2*c01+c11)/d4,
            -4*c00/dx2,         4*(c01+c00)/dx2+(c11+c10)/dy2,  -4*c01/dx2,
            (c11-c10)/d4,                                       (c10-c11)/d4,
            /*               */ -c10/dy2,
            -(2*c00+c10)/d4,    2*(c01-c00)/d4,                 (2*c01+c10)/d4 };
  
          // values for generic PDE discretization point; "v" eqn
          PetscScalar valV[] = {
            (2*c11+c00)/d4,     (c00-c01)/d4,                   -(2*c11+c01)/d4,
            /*               */ -4*c11/dy2,
            2*(c11-c10)/d4,                                     2*(c10-c11)/d4,
            -c00/dx2,           4*(c11+c10)/dy2+(c01+c00)/dx2,  -c01/dx2,
            -(2*c10+c00)/d4,    (c01-c00)/d4,                   (2*c10+c01)/d4,
            /*               */ -4*c10/dy2 };

          // now actually set the values in the matrix 
          ierr = MatSetValues(A, 1, &rowU, stencilSize, colU, valU, INSERT_VALUES); CHKERRQ(ierr);
          ierr = MatSetValues(A, 1, &rowV, stencilSize, colV, valV, INSERT_VALUES); CHKERRQ(ierr);
        } else {
          // values for generic PDE discretization point; "u" eqn
          PetscScalar valU[] = {
            /*               */ -c11/dy2,
            (2*c00+c11)/d4,     2*(c00-c01)/d4,                 -(2*c01+c11)/d4,
            -4*c00/dx2,         4*(c01+c00)/dx2+(c11+c10)/dy2,  -4*c01/dx2,
            (c11-c10)/d4,                                       (c10-c11)/d4,
            /*               */ -c10/dy2,
            -(2*c00+c10)/d4,    2*(c01-c00)/d4,                 (2*c01+c10)/d4 };
  
          // values for generic PDE discretization point; "v" eqn
          PetscScalar valV[] = {
            (2*c11+c00)/d4,     (c00-c01)/d4,                   -(2*c11+c01)/d4,
            /*               */ -4*c11/dy2,
            2*(c11-c10)/d4,                                     2*(c10-c11)/d4,
            -c00/dx2,           4*(c11+c10)/dy2+(c01+c00)/dx2,  -c01/dx2,
            -(2*c10+c00)/d4,    (c01-c00)/d4,                   (2*c10+c01)/d4,
            /*               */ -4*c10/dy2 };

          /* Dragging ice experiences friction at the bed determined by the
           *    basalDrag[x|y]() methods.  These may be a plastic, pseudo-plastic,
           *    or linear friction law according to basal->drag(), which gets called
           *    by basalDragx(),basalDragy().  */
          if ((includeBasalShear) && (vMask.value(i,j) == MASK_DRAGGING_SHEET)) {
            // Dragging is done implicitly (i.e. on left side of SSA eqns for u,v).
            valU[5] += basalDragx(tauc, uv, i, j);
            valV[7] += basalDragy(tauc, uv, i, j);
          }

          // make shelf drag a little bit if desired
          if ((shelvesDragToo == PETSC_TRUE) && (vMask.value(i,j) == MASK_FLOATING)) {
            valU[5] += beta_shelves_drag_too;
            valV[7] += beta_shelves_drag_too;
          }

          // now actually set the values in the matrix 
          ierr = MatSetValues(A, 1, &rowU, stencilSize, colU, valU, INSERT_VALUES); CHKERRQ(ierr);
          ierr = MatSetValues(A, 1, &rowV, stencilSize, colV, valV, INSERT_VALUES); CHKERRQ(ierr);

        }
      }
    }
  }
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vsmoothCFmask.end_access(); CHKERRQ(ierr);

  ierr = vel_ssa.end_access(); CHKERRQ(ierr);
  ierr = vtauc.end_access(); CHKERRQ(ierr);

  ierr = vNuH[0].end_access(); CHKERRQ(ierr);
  ierr = vNuH[1].end_access(); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  return 0;
}



PetscErrorCode IceCalvBCModel::assembleSSARhs(Vec rhs) {
  PetscErrorCode  ierr;
  const PetscInt  M=2*grid.My;

  PetscScalar  **h, **H, **taudx, **taudy,
               **cfmask, **ncf[2];

  double ocean_rho = config.get("sea_water_density"),
    ice_rho = config.get("ice_density");

  // first, find and store normal direction to calving front;
  //   used in forming both rhs and A; also get access to it here
  ierr = compute_nCF(); CHKERRQ(ierr);
  ierr = vsmoothCFmask.get_array(cfmask); CHKERRQ(ierr);
  ierr = vnCF[0].get_array(ncf[0]); CHKERRQ(ierr);
  ierr = vnCF[1].get_array(ncf[1]); CHKERRQ(ierr);
  const PetscScalar Gamma = 0.5 * (1.0 - ice_rho / ocean_rho) * ice_rho * standard_gravity;

  // clear right-hand side
  ierr = VecSet(rhs, 0.0); CHKERRQ(ierr);

  // get driving stress components
  ierr = computeDrivingStress(vWork2d[0],vWork2d[1]); CHKERRQ(ierr); // in iMutil.cc
  ierr = vWork2d[0].get_array(taudx); CHKERRQ(ierr);
  ierr = vWork2d[1].get_array(taudy); CHKERRQ(ierr);

  /* rhs (= right-hand side) assembly loop */
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vh.get_array(h); CHKERRQ(ierr);
  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = uvbar.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscInt J = 2*j,
                     rowU = i*M + J,
                     rowV = i*M + J+1;
      if (vMask.value(i,j) == MASK_SHEET) {
        // non-SSA case: set the velocity; in some cases this is Dirichlet cond
        //   for SSA region
        ierr = VecSetValue(rhs, rowU, 0.5*(uvbar(i-1,j,0) + uvbar(i,j,0)),
                           INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValue(rhs, rowV, 0.5*(uvbar(i,j-1,1) + uvbar(i,j,1)),
                           INSERT_VALUES); CHKERRQ(ierr);
      } else if ( (vMask.value(i,j) == MASK_FLOATING)
                  && (cfmask[i][j] > 0.01) && (cfmask[i][j] < 0.99) ) {
        // calving front case: use info in vnCF[2]
        const PetscScalar GammaHsqr = Gamma * H[i][j] * H[i][j];
        ierr = VecSetValue(rhs, rowU, GammaHsqr * ncf[0][i][j],
                           INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValue(rhs, rowV, GammaHsqr * ncf[1][i][j],
                           INSERT_VALUES); CHKERRQ(ierr);
      } else {
        // general SSA case: get driving stress for right hand side
	ierr = VecSetValue(rhs, rowU, taudx[i][j], INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(rhs, rowV, taudy[i][j], INSERT_VALUES); CHKERRQ(ierr);          
      }
    }
  }
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vh.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = uvbar.end_access(); CHKERRQ(ierr);

  ierr = vWork2d[0].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[1].end_access(); CHKERRQ(ierr);

  ierr = vsmoothCFmask.end_access(); CHKERRQ(ierr);
  ierr = vnCF[0].end_access(); CHKERRQ(ierr);
  ierr = vnCF[1].end_access(); CHKERRQ(ierr);

  // now actually build Vec for rhs
  ierr = VecAssemblyBegin(rhs); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(rhs); CHKERRQ(ierr);
  return 0;
}

