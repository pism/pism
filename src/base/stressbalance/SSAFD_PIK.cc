// Copyright (C) 2011 Torsten Albrecht and Constantine Khroulev
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

#include "SSAFD.hh"

SSA *SSAFD_PIKFactory(IceGrid &g, IceBasalResistancePlasticLaw &b,
                IceFlowLaw &i, EnthalpyConverter &ec,
                const NCConfigVariable &c)
{
  return new SSAFD_PIK(g,b,i,ec,c);
}

/*! */
PetscErrorCode SSAFD_PIK::init(PISMVars &vars) {
  PetscErrorCode ierr;
  ierr = SSAFD::init(vars); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
                    "  [... including PIK CFBC implementation]\n"); CHKERRQ(ierr);

  return 0;
}

/*! */
PetscErrorCode SSAFD_PIK::assemble_matrix(bool include_basal_shear, Mat A) {
  PetscErrorCode ierr;

  ierr = verbPrintf(3,grid.com, "SSAFD_PIK:assemble_matrix is called\n"); CHKERRQ(ierr);

  // put the new matrix assembly here

  const PetscScalar   dx=grid.dx, dy=grid.dy;
  IceModelVec2V vel = velocity;         // a shortcut

  ierr = MatZeroEntries(A); CHKERRQ(ierr);

  /* matrix assembly loop */

  ierr = nuH.begin_access(); CHKERRQ(ierr);
  ierr = tauc->begin_access(); CHKERRQ(ierr);
  ierr = mask->begin_access(); CHKERRQ(ierr);
  ierr = vel.begin_access(); CHKERRQ(ierr);

  if (vel_bc && bc_locations) {
    ierr = bc_locations->begin_access(); CHKERRQ(ierr);
  }

  //IceModelVec2S &thk = *thickness;
  ierr = thickness->begin_access(); CHKERRQ(ierr);
  //ierr = bed->begin_access(); CHKERRQ(ierr);

  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {

      PetscScalar H_ij = (*thickness)(i,j),
        H_e = (*thickness)(i + 1,j),
        H_w = (*thickness)(i - 1,j),
        H_n = (*thickness)(i,j + 1),
        H_s = (*thickness)(i,j - 1);

      bool onIcefreeOcean = H_ij <= 1.0;

      //defined so far via ice thickness
      bool atBoundary = H_ij > 100.0 && (H_e <= 1.0 || H_w <= 1.0 || H_s <= 1.0 || H_n <= 1.0);

      if (vel_bc && bc_locations && bc_locations->as_int(i,j) == 1) {
        // set diagonal entry to one; RHS entry will be known (e.g. SIA) velocity;
        //   this is where boundary value to SSA is set
        ierr = set_diagonal_matrix_entry(A, i, j, scaling); CHKERRQ(ierr);
        continue;
      }

      if (onIcefreeOcean) { // vanish ice velocities on the ice free ocean
        ierr = set_diagonal_matrix_entry(A, i, j, scaling); CHKERRQ(ierr);
        continue;
      }

      if (atBoundary) {
        const PetscScalar dx2 = dx*dx, d4 = dx*dy*4, dy2 = dy*dy;

        const PetscScalar c_w = nuH(i - 1, j, 0);
        const PetscScalar c_e = nuH(i, j, 0);
        const PetscScalar c_s = nuH(i, j - 1, 1);
        const PetscScalar c_n = nuH(i, j, 1);

        //in i-direction
        PetscInt aMn = 1, aPn = 1, aMM = 1, aPP = 1, aMs = 1, aPs = 1;

        //in j-direction
        PetscInt bPw = 1, bPP = 1, bPe = 1, bMw = 1, bMM = 1, bMe = 1;

        //direct adjacent neighbors
        if (H_e <= 1.0) aPP = 0;
        if (H_w <= 1.0) aMM = 0;
        if (H_n <= 1.0) bPP = 0;
        if (H_s <= 1.0) bMM = 0;

        //neighbors in the corners
        PetscScalar H_ne = (*thickness)(i + 1,j + 1),
          H_se = (*thickness)(i + 1,j - 1),
          H_nw = (*thickness)(i - 1,j + 1),
          H_sw = (*thickness)(i - 1,j - 1);

        //for the each single boundary to decide, which derivative to drop
        if (H_n <= 1.0 || H_ne <= 1.0) aPn = 0;
        if (H_e <= 1.0 || H_ne <= 1.0) bPe = 0;
        if (H_e <= 1.0 || H_se <= 1.0) bMe = 0;
        if (H_s <= 1.0 || H_se <= 1.0) aPs = 0;
        if (H_s <= 1.0 || H_sw <= 1.0) aMs = 0;
        if (H_w <= 1.0 || H_sw <= 1.0) bMw = 0;
        if (H_w <= 1.0 || H_nw <= 1.0) bPw = 0;
        if (H_n <= 1.0 || H_nw <= 1.0) aMn = 0;

        const PetscInt sten = 14;//one more than default
        MatStencil  row, col[sten];

        // This complicated stencil is the result of adding factors (0 if boundary or 1 if not) for every single derivatice of the SSA equations across one of the 4 direct or 8 neighboring grid cell boundaries.
        // If this ice grid cell was surrounded only by other ice grid cells, we would get the standard stencil of the SSA.
        // Derivatives across the 4 direct neighbors are replaced by hydrostatic pressure on the rhs.

        PetscScalar valU[] = {

          /*                                           */ -bPP*c_n/dy2,
          (2*bPw*aMM*c_w+aMn*bPP*c_n)/d4,                 (2*bPP*(c_w*aMM-c_e*aPP)+c_n*bPP*(aPn-aMn))/d4,                -(2*bPe*aPP*c_e+aPn*bPP*c_n)/d4,
          -4*aMM*c_w/dx2,                                  4*(aPP*c_e+aMM*c_w)/dx2+(bPP*c_n+bMM*c_s)/dy2,                 -4*aPP*c_e/dx2,
          (aMM*(bPP*c_n-bMM*c_s)+2*c_w*aMM*(bMw-bPw))/d4, (2*(c_e*aPP-c_w*aMM)*(bPP-bMM)+(c_n*bPP-c_s*bMM)*(aPP-aMM))/d4, (aPP*(c_s*bMM-c_n*bPP)+2*c_e*aPP*(bPe-bMe))/d4,
          /*                                           */ -bMM*c_s/dy2,
          -(2*bMw*aMM*c_w+aMs*bMM*c_s)/d4,                (2*bMM*(aPP*c_e-c_w*aMM)+c_s*bMM*(aMs-aPs))/d4,                 (2*bMe*aPP*c_e+aPs*bMM*c_s)/d4};


        PetscScalar valV[] = {

          (2*aMn*bPP*c_n+bPw*aMM*c_w)/d4,                 (bPP*(c_w*aMM-aPP*c_e)+2*c_n*bPP*(aPn-aMn))/d4,                -(2*aPn*bPP*c_n+bPe*aPP*c_e)/d4,
          /*                                           */ -4*bPP*c_n/dy2,
          (2*aMM*(bPP*c_n-c_s*bMM)+c_w*aMM*(bMw-bPw))/d4, (2*(bPP*c_n-c_s*bMM)*(aPP-aMM)+(aPP*c_e-c_w*aMM)*(bPP-bMM))/d4, (2*aPP*(c_s*bMM-c_n*bPP)+c_e*aPP*(bPe-bMe))/d4,
          -aMM*c_w/dx2,                                    4*(bPP*c_n+bMM*c_s)/dy2+(aPP*c_e+aMM*c_w)/dx2,                 -aPP*c_e/dx2,
          -(2*aMs*bMM*c_s+bMw*aMM*c_w)/d4,                (bMM*(aPP*c_e-c_w*aMM)+2*c_s*bMM*(aMs-aPs))/d4,                 (2*aPs*bMM*c_s+bMe*aPP*c_e)/d4,
          /*                                           */ -4*bMM*c_s/dy2 };



    	if (include_basal_shear && (mask->as_int(i,j) == MASK_GROUNDED)) {
          // Dragging is done implicitly (i.e. on left side of SSA eqns for u,v).
          valU[5] += basal.drag((*tauc)(i,j), vel(i,j).u, vel(i,j).v);
          valV[8] += basal.drag((*tauc)(i,j), vel(i,j).u, vel(i,j).v);// without bc it would be the 7th point of the stencil
    	}


        // build "u" equation: NOTE TRANSPOSE
     	row.j = i; row.i = j; row.c = 0;
     	const PetscInt UI[] = {
          /*       */ i,
          i-1,        i,          i+1,
          i-1,        i,          i+1,
          i-1,        i,          i+1,
          /*       */ i,
          i-1,        i,          i+1};
     	const PetscInt UJ[] = {
          /*       */ j+1,
          j+1,        j+1,        j+1,
          j,          j,          j,
          j,          j,          j,
          /*       */ j-1,
          j-1,        j-1,        j-1};
     	const PetscInt UC[] = {
          /*       */ 0,
          1,          1,          1,
          0,          0,          0,
          1,          1,          1,
          /*       */ 0,
          1,          1,          1};
     	for (PetscInt m=0; m<sten; m++) {
          col[m].j = UI[m]; col[m].i = UJ[m], col[m].c = UC[m];
     	}
     	ierr = MatSetValuesStencil(A,1,&row,sten,col,valU,INSERT_VALUES); CHKERRQ(ierr);

        // build "v" equation: NOTE TRANSPOSE
        row.j = i; row.i = j; row.c = 1;
        const PetscInt VI[] = {
          i-1,        i,          i+1,
          /*       */ i,
          i-1,        i,          i+1,
          i-1,        i,          i+1,
          i-1,        i,          i+1,
          /*       */ i};
        const PetscInt VJ[] = {
          j+1,        j+1,        j+1,
          /*       */ j+1,
          j,          j,          j,
          j,          j,          j,
          j-1,        j-1,        j-1,
          /*       */ j-1};
        const PetscInt VC[] = {
          0,          0,          0,
          /*       */ 1,
          0,          0,          0,
          1,          1,          1,
          0,          0,          0,
          /*       */ 1};
        for (PetscInt m=0; m<sten; m++) {
          col[m].j = VI[m]; col[m].i = VJ[m], col[m].c = VC[m];
        }
        ierr = MatSetValuesStencil(A,1,&row,sten,col,valV,INSERT_VALUES); CHKERRQ(ierr);

	//the rest is just a copy of SSAFD
      } else {
        const PetscScalar dx2 = dx*dx, d4 = dx*dy*4, dy2 = dy*dy;
        /* Provide shorthand for the following staggered coefficients  nu H:
         *      c_n
         *  c_w     c_e
         *      c_s
         * Note that the positive i (x) direction is right and the positive j (y)
         * direction is up. */
        const PetscScalar c_w = nuH(i-1,j,0);
        const PetscScalar c_e = nuH(i,j,0);
        const PetscScalar c_s = nuH(i,j-1,1);
        const PetscScalar c_n = nuH(i,j,1);

        const PetscInt sten = 13;
        MatStencil  row, col[sten];

        /* start with the values at the points */
        PetscScalar valU[] = {
          /*               */ -c_n/dy2,
          (2*c_w+c_n)/d4,     2*(c_w-c_e)/d4,                 -(2*c_e+c_n)/d4,
          -4*c_w/dx2,         4*(c_e+c_w)/dx2+(c_n+c_s)/dy2,  -4*c_e/dx2,
          (c_n-c_s)/d4,                                       (c_s-c_n)/d4,
          /*               */ -c_s/dy2,
          -(2*c_w+c_s)/d4,    2*(c_e-c_w)/d4,                 (2*c_e+c_s)/d4 };
        PetscScalar valV[] = {
          (2*c_n+c_w)/d4,     (c_w-c_e)/d4,                   -(2*c_n+c_e)/d4,
          /*               */ -4*c_n/dy2,
          2*(c_n-c_s)/d4,                                     2*(c_s-c_n)/d4,
          -c_w/dx2,           4*(c_n+c_s)/dy2+(c_e+c_w)/dx2,  -c_e/dx2,
          -(2*c_s+c_w)/d4,    (c_e-c_w)/d4,                   (2*c_s+c_e)/d4,
          /*               */ -4*c_s/dy2 };

        /* Dragging ice experiences friction at the bed determined by the
         *    basalDrag[x|y]() methods.  These may be a plastic, pseudo-plastic,
         *    or linear friction law according to basal->drag(), which gets called
         *    by basalDragx(),basalDragy().  */
        if (include_basal_shear && (mask->as_int(i,j) == MASK_GROUNDED)) {
          // Dragging is done implicitly (i.e. on left side of SSA eqns for u,v).
          valU[5] += basal.drag((*tauc)(i,j), vel(i,j).u, vel(i,j).v);
          valV[7] += basal.drag((*tauc)(i,j), vel(i,j).u, vel(i,j).v);
        }

        // build "u" equation: NOTE TRANSPOSE
        row.j = i; row.i = j; row.c = 0;
        const PetscInt UI[] = {
          /*       */ i,
          i-1,        i,          i+1,
          i-1,        i,          i+1,
          i-1,                    i+1,
          /*       */ i,
          i-1,        i,          i+1};
        const PetscInt UJ[] = {
          /*       */ j+1,
          j+1,        j+1,        j+1,
          j,          j,          j,
          j,                      j,
          /*       */ j-1,
          j-1,        j-1,        j-1};
        const PetscInt UC[] = {
          /*       */ 0,
          1,          1,          1,
          0,          0,          0,
          1,                      1,
          /*       */ 0,
          1,          1,          1};
        for (PetscInt m=0; m<sten; m++) {
          col[m].j = UI[m]; col[m].i = UJ[m], col[m].c = UC[m];
        }
        ierr = MatSetValuesStencil(A,1,&row,sten,col,valU,INSERT_VALUES); CHKERRQ(ierr);

        // build "v" equation: NOTE TRANSPOSE
        row.j = i; row.i = j; row.c = 1;
        const PetscInt VI[] = {
          i-1,        i,          i+1,
          /*       */ i,
          i-1,                    i+1,
          i-1,        i,          i+1,
          i-1,        i,          i+1,
          /*       */ i};
        const PetscInt VJ[] = {
          j+1,        j+1,        j+1,
          /*       */ j+1,
          j,                      j,
          j,          j,          j,
          j-1,        j-1,        j-1,
          /*       */ j-1};
        const PetscInt VC[] = {
          0,          0,          0,
          /*       */ 1,
          0,                      0,
          1,          1,          1,
          0,          0,          0,
          /*       */ 1};
        for (PetscInt m=0; m<sten; m++) {
          col[m].j = VI[m]; col[m].i = VJ[m], col[m].c = VC[m];
        }
        ierr = MatSetValuesStencil(A,1,&row,sten,col,valV,INSERT_VALUES); CHKERRQ(ierr);
      }
    }
  }

  if (vel_bc && bc_locations) {
    ierr = bc_locations->end_access(); CHKERRQ(ierr);
  }

  ierr = vel.end_access(); CHKERRQ(ierr);
  ierr = mask->end_access(); CHKERRQ(ierr);
  ierr = tauc->end_access(); CHKERRQ(ierr);
  ierr = nuH.end_access(); CHKERRQ(ierr);

  ierr = thickness->end_access(); CHKERRQ(ierr);
  //ierr = bed->end_access(); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  return 0;
}


/*! */
PetscErrorCode SSAFD_PIK::assemble_rhs(Vec rhs) {
  PetscErrorCode ierr;

  ierr = verbPrintf(3, grid.com, "SSAFD_PIK:assemble_rhs is called\n"); CHKERRQ(ierr);

  const double dx = grid.dx, dy = grid.dy;
  PISMVector2 **rhs_uv;

  ierr = VecSet(rhs, 0.0); CHKERRQ(ierr);

  // get driving stress components
  ierr = compute_driving_stress(taud); CHKERRQ(ierr);

  ierr = taud.begin_access(); CHKERRQ(ierr);
  ierr = DAVecGetArray(SSADA, rhs, &rhs_uv); CHKERRQ(ierr);

  if (vel_bc && bc_locations) {
    ierr = vel_bc->begin_access(); CHKERRQ(ierr);
    ierr = bc_locations->begin_access(); CHKERRQ(ierr);
  }

  ierr = thickness->begin_access(); CHKERRQ(ierr);
  ierr = bed->begin_access(); CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {
      double H_ij = (*thickness)(i, j),
        H_e = (*thickness)(i + 1, j),
        H_w = (*thickness)(i - 1, j),
        H_n = (*thickness)(i, j + 1),
        H_s = (*thickness)(i, j - 1);

      bool ice_free = (H_ij <= 1.0);
      bool boundary = H_ij > 100.0 &&
        (H_e <= 1.0 || H_w <= 1.0 || H_s <= 1.0 || H_n <= 1.0);

      if (vel_bc && (bc_locations->as_int(i, j) == 1)) {
        rhs_uv[i][j].u = scaling * (*vel_bc)(i, j).u;
        rhs_uv[i][j].v = scaling * (*vel_bc)(i, j).v;
        continue;
      }

      if (ice_free) {
        rhs_uv[i][j].u = 0.0;
        rhs_uv[i][j].v = 0.0;
        continue;
      }

      if (boundary) {
        PetscInt aMM = 1, aPP = 1, bMM = 1, bPP = 1;
        //direct adjacent neighbors
        if (H_e <= 1.0) aPP = 0;
        if (H_w <= 1.0) aMM = 0;
        if (H_n <= 1.0) bPP = 0;
        if (H_s <= 1.0) bMM = 0;

        const double standard_gravity = config.get("standard_gravity");
        double ocean_rho = config.get("sea_water_density");

        const double ice_pressure = ice.rho * standard_gravity * H_ij;
        double ocean_pressure;

        const double h_grounded = (*bed)(i,j) + H_ij,
          h_floating = sea_level + (1.0 - ice.rho / ocean_rho) * H_ij,
          H_ij2 = H_ij*H_ij;

        double h_ij = 0.0, tdx = 0.0, tdy = 0.0;

        if (H_ij > 0.0 && (*bed)(i,j) < (sea_level - (ice.rho / ocean_rho) * H_ij)) {
          //calving front boundary condition for floating shelf
          ocean_pressure = 0.5 * ice.rho * standard_gravity * (1 - (ice.rho / ocean_rho))*H_ij2;
          // this is not really the ocean_pressure, but the difference between
          // ocean_pressure and isotrop.normal stresses (=pressure) from within
          // the ice
          h_ij = h_floating;
        } else {
          h_ij = h_grounded;
          if( (*bed)(i,j) >= sea_level){//boundary condition for cliff --> zero stress = ocean_pressure
            ocean_pressure = 0.5 * ice.rho * standard_gravity * H_ij2;
            // this is not 'zero' because the isotrop.normal stresses
            // (=pressure) from within the ice figures on RHS
          }else{//boundary condition for marine terminating glacier
            ocean_pressure = 0.5 * ice.rho * standard_gravity *
              (H_ij2 - (ocean_rho / ice.rho)*(sea_level - (*bed)(i,j))*(sea_level - (*bed)(i,j)));
          }
        }


        //here we take the direct gradient at the boundary (not centered)
        if (aPP == 0 && aMM == 1) tdx = ice_pressure*h_ij / dx;
        else if (aMM == 0 && aPP == 1) tdx = -ice_pressure*h_ij / dx;
        else if (aPP == 0 && aMM == 0) tdx = 0; //in case of some kind of ice nose, or ice bridge

        if (bPP == 0 && bMM == 1) tdy = ice_pressure*h_ij / dy;
        else if (bMM == 0 && bPP == 1) tdy = -ice_pressure*h_ij / dy;
        else if (bPP == 0 && bMM == 0) tdy = 0;

        rhs_uv[i][j].u = tdx - (aMM - aPP)*ocean_pressure / dx;
        rhs_uv[i][j].v = tdy - (bMM - bPP)*ocean_pressure / dy;

        continue;
      }

      // usual case: use already computed driving stress
      rhs_uv[i][j].u = taud(i,j).u;
      rhs_uv[i][j].v = taud(i,j).v;
    }
  }
  ierr = thickness->end_access(); CHKERRQ(ierr);
  ierr = bed->end_access(); CHKERRQ(ierr);

  if (vel_bc) {
    ierr = bc_locations->end_access(); CHKERRQ(ierr);
    ierr = vel_bc->end_access(); CHKERRQ(ierr);
  }

  ierr = taud.end_access(); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(SSADA,rhs,&rhs_uv); CHKERRQ(ierr);

  return 0;
}

