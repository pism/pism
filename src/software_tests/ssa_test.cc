// Copyright (C) 2010 Ed Bueler
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

// TO DO:  make std out from velocitySSA() immediate by using PetscPrintf
// (and then run on greenland example again)

/*

Example uses, from a PISM directory:

$ pismv -test I -My 25 -Mx 3 -o testI.nc
$ ssa_test -i testI.nc

$ pismv -test J -Mx 40 -My 40 -o testJ.nc
$ ssa_test -i testJ.nc

$ cd example/pst/
$ ./pst.sh 8 >> out.pst &  # takes a while; produces P1.nc; can be shortened
$ ssa_test -i P1.nc


*/

static char help[] =
  "\nSSA_TEST\n"
  "  Testing program for SSA, time-independent calculations separate from\n"
  "  IceModel.  Also may be used in a PISM software (regression) test.\n\n";

#include <cmath>
#include <cstdio>
#include <string>
#include <petscksp.h>
#include "../base/pism_const.hh"
#include "../base/grid.hh"
#include "../base/iceModelVec.hh"
#include "../base/flowlaw_factory.hh" // IceFlowLawFactory, IceFlowLaw
#include "../base/materials.hh" // SSAStrengthExtension, IceBasalResistancePlasticLaw
#include "../base/nc_util.hh"
#include "../base/PISMIO.hh"
#include "../base/NCVariable.hh"



// LARGE BLOCKS OF COMMENTS IN pism-dev/src/base/iMssa.cc ARE REMOVED FROM THIS
// STAND-ALONE, WHICH HAS LOTS OF CODE DUPLICATION ...


PetscErrorCode trivialMove(IceGrid &grid, DA &SSADA, Vec &SSAX, 
                           IceModelVec2V &vel_ssa) {
  PetscErrorCode  ierr;
  PISMVector2 **Xuv;
  ierr = vel_ssa.begin_access(); CHKERRQ(ierr);
  ierr = DAVecGetArray(SSADA,SSAX,&Xuv); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      vel_ssa(i,j).u = Xuv[i][j].u;
      vel_ssa(i,j).v = Xuv[i][j].v;
    }
  }
  ierr = DAVecRestoreArray(SSADA,SSAX,&Xuv); CHKERRQ(ierr);
  ierr = vel_ssa.end_access(); CHKERRQ(ierr);
  return 0;
}



PetscErrorCode computeEffectiveViscosity(
    IceGrid &grid, NCConfigVariable &config, IceFlowLaw *ice,
    SSAStrengthExtension &ssaStrengthExtend, IceModelVec2S &vH, 
    IceModelVec2V &vel_ssa, PetscReal hardness, PetscReal epsilon,
    IceModelVec2S vNuH[2]) {
  PetscErrorCode ierr;

  bool use_constant_nuh_for_ssa = config.get_flag("use_constant_nuh_for_ssa");
  if (use_constant_nuh_for_ssa) {
    // Intended only for debugging, this treats the entire domain as though
    // it were the strength extension (i.e. strength does not depend on thickness)
    PetscReal nuH = ssaStrengthExtend.notional_strength();
    ierr = vNuH[0].set(nuH); CHKERRQ(ierr);
    ierr = vNuH[1].set(nuH); CHKERRQ(ierr);
    return 0;
  }

  // We need to compute integrated effective viscosity (\bar\nu * H).
  // It is locally determined by the strain rates and temperature field.
  PetscScalar **nuH[2];
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vNuH[0].get_array(nuH[0]); CHKERRQ(ierr);
  ierr = vNuH[1].get_array(nuH[1]); CHKERRQ(ierr);

  PISMVector2 **uv;
  ierr = vel_ssa.get_array(uv); CHKERRQ(ierr);

  const PetscScalar   dx = grid.dx, dy = grid.dy;

  for (PetscInt o=0; o<2; ++o) {
    const PetscInt oi = 1 - o, oj=o;
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        const PetscScalar H = 0.5 * (vH(i,j) + vH(i+oi,j+oj));
        if (H < ssaStrengthExtend.min_thickness_for_extension()) {
          // Extends strength of SSA (i.e. nuH coeff) into the ice free region.
          //  Does not add or subtract ice mass.
          nuH[o][i][j] = ssaStrengthExtend.notional_strength();
        } else {
          PetscScalar u_x, u_y, v_x, v_y;
          // Check the offset to determine how to differentiate velocity
          if (o == 0) {
            u_x = (uv[i+1][j].u - uv[i][j].u) / dx;
            u_y = (uv[i][j+1].u + uv[i+1][j+1].u - uv[i][j-1].u - uv[i+1][j-1].u) / (4*dy);
            v_x = (uv[i+1][j].v - uv[i][j].v) / dx;
            v_y = (uv[i][j+1].v + uv[i+1][j+1].v - uv[i][j-1].v - uv[i+1][j-1].v) / (4*dy);
          } else {
            u_x = (uv[i+1][j].u + uv[i+1][j+1].u - uv[i-1][j].u - uv[i-1][j+1].u) / (4*dx);
            u_y = (uv[i][j+1].u - uv[i][j].u) / dy;
            v_x = (uv[i+1][j].v + uv[i+1][j+1].v - uv[i-1][j].v - uv[i-1][j+1].v) / (4*dx);
            v_y = (uv[i][j+1].v - uv[i][j].v) / dy;
          }

          nuH[o][i][j] = H * ice->effectiveViscosity(hardness, u_x, u_y, v_x, v_y);

          if (! finite(nuH[o][i][j]) || false) {
            ierr = PetscPrintf(grid.com, "nuH[%d][%d][%d] = %e\n", o, i, j, nuH[o][i][j]);
              CHKERRQ(ierr); 
            ierr = PetscPrintf(grid.com, "  u_x, u_y, v_x, v_y = %e, %e, %e, %e\n", 
                               u_x, u_y, v_x, v_y);
              CHKERRQ(ierr);
          }
          
          // We ensure that nuH is bounded below by a positive constant.
          nuH[o][i][j] += epsilon;
        } // end of if (vH(i,j) < ssaStrengthExtend.min_thickness_for_extension()) { ... } else {
      } // j
    } // i
  } // o
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vNuH[0].end_access(); CHKERRQ(ierr);
  ierr = vNuH[1].end_access(); CHKERRQ(ierr);

  ierr = vel_ssa.end_access(); CHKERRQ(ierr);

  // Some communication
  ierr = vNuH[0].beginGhostComm(); CHKERRQ(ierr);
  ierr = vNuH[0].endGhostComm(); CHKERRQ(ierr);
  ierr = vNuH[1].beginGhostComm(); CHKERRQ(ierr);
  ierr = vNuH[1].endGhostComm(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode testConvergenceOfNu(IceGrid &grid, 
                  IceModelVec2S vNuH[2], IceModelVec2S vNuHOld[2],
                  PetscReal *norm, PetscReal *normChange) {
  PetscErrorCode  ierr;
  PetscReal nuNorm[2], nuChange[2];
  const PetscScalar area = grid.dx * grid.dy;
#define MY_NORM     NORM_1

  // Test for change in nu
  ierr = vNuHOld[0].add(-1, vNuH[0]); CHKERRQ(ierr);
  ierr = vNuHOld[1].add(-1, vNuH[1]); CHKERRQ(ierr);

  ierr = vNuHOld[0].norm(MY_NORM, nuChange[0]); CHKERRQ(ierr);
  nuChange[0] *= area;
  ierr = vNuHOld[1].norm(MY_NORM, nuChange[1]); CHKERRQ(ierr);
  nuChange[1] *= area;

  *normChange = sqrt(PetscSqr(nuChange[0]) + PetscSqr(nuChange[1]));

  ierr = vNuH[0].norm(MY_NORM, nuNorm[0]); CHKERRQ(ierr);
  nuNorm[0] *= area;
  ierr = vNuH[1].norm(MY_NORM, nuNorm[1]); CHKERRQ(ierr);
  nuNorm[1] *= area;

  *norm = sqrt(PetscSqr(nuNorm[0]) + PetscSqr(nuNorm[1]));
  return 0;
}


PetscErrorCode assembleSSAMatrix(
    IceGrid &grid, NCConfigVariable &config, IceBasalResistancePlasticLaw &basal,
    IceModelVec2S &vtauc, IceModelVec2Mask &vMask, IceModelVec2V &vel_ssa,
    bool includeBasalShear, IceModelVec2S vNuH[2], Mat &A) {
  PetscErrorCode  ierr;

  const PetscScalar   dx=grid.dx, dy=grid.dy;
  // next constant not too sensitive, but must match value in assembleSSARhs():
  const PetscScalar   scaling = 1.0e9;  // comparable to typical beta for an ice stream
  PetscScalar     **nuH[2], **tauc;
  PISMVector2     **uvssa;

  ierr = MatZeroEntries(A); CHKERRQ(ierr);

  /* matrix assembly loop */
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vtauc.get_array(tauc); CHKERRQ(ierr);
  ierr = vel_ssa.get_array(uvssa); CHKERRQ(ierr);
  ierr = vNuH[0].get_array(nuH[0]); CHKERRQ(ierr);
  ierr = vNuH[1].get_array(nuH[1]); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PismMask mask_value = vMask.value(i,j);
      if (mask_value == MASK_SHEET) {
        // set diagonal entry to one; RHS entry will be known (e.g. SIA) velocity;
        //   this is where boundary value to SSA is set
        MatStencil  row, col;
        row.j = i; row.i = j; row.c = 0;
        col.j = i; col.i = j; col.c = 0;
        ierr = MatSetValuesStencil(A,1,&row,1,&col,&scaling,INSERT_VALUES); CHKERRQ(ierr);
        row.c = 1;
        col.c = 1;
        ierr = MatSetValuesStencil(A,1,&row,1,&col,&scaling,INSERT_VALUES); CHKERRQ(ierr);
      } else {
        const PetscScalar dx2 = dx*dx, d4 = dx*dy*4, dy2 = dy*dy;
        /* Provide shorthand for the following staggered coefficients  nu H:
        *      c11
        *  c00     c01
        *      c10
        * Note that the positive i (x) direction is right and the positive j (y)
        * direction is up. */
        const PetscScalar c00 = nuH[0][i-1][j];
        const PetscScalar c01 = nuH[0][i][j];
        const PetscScalar c10 = nuH[1][i][j-1];
        const PetscScalar c11 = nuH[1][i][j];

        const PetscInt sten = 13;
        MatStencil  row, col[sten];

        /* start with the values at the points */
        PetscScalar valU[] = {
          /*               */ -c11/dy2,
          (2*c00+c11)/d4,     2*(c00-c01)/d4,                 -(2*c01+c11)/d4,
          -4*c00/dx2,         4*(c01+c00)/dx2+(c11+c10)/dy2,  -4*c01/dx2,
          (c11-c10)/d4,                                       (c10-c11)/d4,
          /*               */ -c10/dy2,
          -(2*c00+c10)/d4,    2*(c01-c00)/d4,                 (2*c01+c10)/d4 };
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
        if ((includeBasalShear) && (mask_value == MASK_DRAGGING_SHEET)) {
          // Dragging is done implicitly (i.e. on left side of SSA eqns for u,v).
          const PetscScalar
            mydrag = basal.drag(tauc[i][j], uvssa[i][j].u, uvssa[i][j].v);
          valU[5] += mydrag;
          valV[7] += mydrag;
          //valU[5] += basalDragx(tauc, uvssa, i, j);
          //valV[7] += basalDragy(tauc, uvssa, i, j);
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
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vel_ssa.end_access(); CHKERRQ(ierr);
  ierr = vtauc.end_access(); CHKERRQ(ierr);
  ierr = vNuH[0].end_access(); CHKERRQ(ierr);
  ierr = vNuH[1].end_access(); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode computeDrivingStress(
    IceGrid &grid, NCConfigVariable &config, PetscReal n, PetscReal rho,
    IceModelVec2S &vh, IceModelVec2S &vH, 
    IceModelVec2S &vbed, IceModelVec2Mask &vMask,
    IceModelVec2S &vtaudx, IceModelVec2S &vtaudy) {
  PetscErrorCode ierr;

  const PetscScalar etapow  = (2.0 * n + 2.0)/n,  // = 8/3 if n = 3
                    invpow  = 1.0 / etapow,  // = 3/8
                    dinvpow = (- n - 2.0) / (2.0 * n + 2.0); // = -5/8
  const PetscScalar minThickEtaTransform = 5.0; // m
  const PetscScalar dx=grid.dx, dy=grid.dy;

  bool compute_surf_grad_inward_ssa = config.get_flag("compute_surf_grad_inward_ssa");
  bool use_eta = (config.get_string("surface_gradient_method") == "eta");
  double standard_gravity = config.get("standard_gravity");

  ierr =    vh.begin_access();    CHKERRQ(ierr);
  ierr =    vH.begin_access();  CHKERRQ(ierr);
  ierr =  vbed.begin_access();  CHKERRQ(ierr);
  ierr = vMask.begin_access();  CHKERRQ(ierr);
  ierr = vtaudx.begin_access(); CHKERRQ(ierr);
  ierr = vtaudy.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar pressure = rho * standard_gravity * vH(i,j);
      if (pressure <= 0.0) {
        vtaudx(i,j) = 0.0;
        vtaudy(i,j) = 0.0;
      } else {
        PetscScalar h_x = 0.0, h_y = 0.0;
        // FIXME: we need to handle grid periodicity correctly.
        if (vMask.is_grounded(i,j) && (use_eta == true)) {
	        // in grounded case, differentiate eta = H^{8/3} by chain rule
          if (vH(i,j) > 0.0) {
            const PetscScalar myH = (vH(i,j) < minThickEtaTransform) ?
	                                  minThickEtaTransform : vH(i,j);
            const PetscScalar eta = pow(myH, etapow), factor = invpow * pow(eta, dinvpow);
            h_x = factor * (pow(vH(i+1,j),etapow) - pow(vH(i-1,j),etapow)) / (2*dx);
            h_y = factor * (pow(vH(i,j+1),etapow) - pow(vH(i,j-1),etapow)) / (2*dy);
          }
          // now add bed slope to get actual h_x,h_y
          // FIXME: there is no reason to assume user's bed is periodized
          h_x += vbed.diff_x(i,j);
          h_y += vbed.diff_y(i,j);
        } else {  // floating or eta transformation is not used
          if (compute_surf_grad_inward_ssa) {
            h_x = vh.diff_x_p(i,j);
            h_y = vh.diff_y_p(i,j);
          } else {
            h_x = vh.diff_x(i,j);
            h_y = vh.diff_y(i,j);
          }
        }

        vtaudx(i,j) = - pressure * h_x;
        vtaudy(i,j) = - pressure * h_y;
      }
    }
  }
  ierr =   vbed.end_access(); CHKERRQ(ierr);
  ierr =     vh.end_access(); CHKERRQ(ierr);
  ierr =     vH.end_access(); CHKERRQ(ierr);
  ierr =  vMask.end_access(); CHKERRQ(ierr);
  ierr = vtaudx.end_access(); CHKERRQ(ierr);
  ierr = vtaudy.end_access(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode assembleSSARhs(
    IceGrid &grid, NCConfigVariable &config, PetscReal n, PetscReal rho,
    DA &SSADA, IceModelVec2S *vWork2d, IceModelVec2Stag &uvbar,
    IceModelVec2S &vh, IceModelVec2S &vH,
    IceModelVec2S &vbed, IceModelVec2Mask &vMask,
    Vec &rhs) {
  PetscErrorCode  ierr;

  // next constant not too sensitive, but must match value in assembleSSAMatrix():
  const PetscScalar   scaling = 1.0e9;  // comparable to typical beta for an ice stream;

  ierr = VecSet(rhs, 0.0); CHKERRQ(ierr);

  // get driving stress components
  ierr = computeDrivingStress(grid, config, n, rho,
                              vh, vH, vbed, vMask, 
                              vWork2d[0], vWork2d[1]); CHKERRQ(ierr);

  PetscScalar     **taudx, **taudy;
  PISMVector2     **rhs_uv;
  ierr = vWork2d[0].get_array(taudx); CHKERRQ(ierr);
  ierr = vWork2d[1].get_array(taudy); CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = uvbar.begin_access(); CHKERRQ(ierr);
  ierr = DAVecGetArray(SSADA,rhs,&rhs_uv); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vMask.value(i,j) == MASK_SHEET) {
        rhs_uv[i][j].u = scaling * 0.5*(uvbar(i-1,j,0) + uvbar(i,j,0));
        rhs_uv[i][j].v = scaling * 0.5*(uvbar(i,j-1,1) + uvbar(i,j,1));
      } else {	// usual case: use already computed driving stress
        rhs_uv[i][j].u = taudx[i][j];
        rhs_uv[i][j].v = taudy[i][j];
      }
    }
  }
  ierr = DAVecRestoreArray(SSADA,rhs,&rhs_uv); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = uvbar.end_access(); CHKERRQ(ierr);
  ierr = vWork2d[0].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[1].end_access(); CHKERRQ(ierr);

  ierr = VecAssemblyBegin(rhs); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(rhs); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode velocitySSA(
    IceGrid &grid, NCConfigVariable &config, IceFlowLaw *ice,
    SSAStrengthExtension &ssaStrengthExtend, IceBasalResistancePlasticLaw &basal,
    IceModelVec2S &vtauc, 
    IceModelVec2S &vh, IceModelVec2S &vH,
    IceModelVec2S &vbed, IceModelVec2Mask &vMask,
    IceModelVec2Stag &uvbar, IceModelVec2S *vWork2d,
    IceModelVec2S vNuH[2], IceModelVec2V &vel_ssa_old,
    DA &SSADA, Vec &SSAX, Vec &SSARHS, Mat &A, KSP &SSAKSP, 
    IceModelVec2V &vel_ssa, PetscInt *numiter) {

  PetscErrorCode ierr;

  PetscPrintf(grid.com, "SSA stress balance computation:\n");

  IceModelVec2S vNuHOld[2] = {vWork2d[2], vWork2d[3]};
  PetscReal   norm, normChange;
  PetscInt    its;
  KSPConvergedReason  reason;

  const PetscReal hardness = 1.9e8;  // Pa s^{1/3}; see p. 49 of MacAyeal et al 1996
  
  PetscReal ssaRelativeTolerance = config.get("ssa_relative_convergence"),
            epsilon              = config.get("epsilon_ssa");

  PetscInt ssaMaxIterations = static_cast<PetscInt>(config.get("max_iterations_ssa"));
  
  ierr = vel_ssa.copy_to(vel_ssa_old); CHKERRQ(ierr);

  // computation of RHS only needs to be done once; does not depend on solution;
  //   but matrix changes under nonlinear iteration (loop over k below)
  ierr = assembleSSARhs(grid, config, ice->exponent(), ice->rho, SSADA, vWork2d, uvbar,
                        vh, vH, vbed, vMask, SSARHS); CHKERRQ(ierr);

  for (PetscInt l=0; ; ++l) { // iterate with increasing regularization parameter
    ierr = computeEffectiveViscosity(
                grid, config, ice, ssaStrengthExtend,
                vH, vel_ssa, hardness, epsilon, vNuH); CHKERRQ(ierr);
    //ierr = update_nu_viewers(vNuH); CHKERRQ(ierr);
    // iterate on effective viscosity: "outer nonlinear iteration":
    for (PetscInt k=0; k<ssaMaxIterations; ++k) { 
      PetscPrintf(grid.com, "  %d,%2d:", l, k);
    
      // in preparation of measuring change of effective viscosity:
      ierr = vNuH[0].copy_to(vNuHOld[0]); CHKERRQ(ierr);
      ierr = vNuH[1].copy_to(vNuHOld[1]); CHKERRQ(ierr);

      // assemble (or re-assemble) matrix, which depends on updated viscosity
      ierr = assembleSSAMatrix(grid, config, basal, vtauc, vMask, vel_ssa, 
                               true, vNuH, A); CHKERRQ(ierr);
      PetscPrintf(grid.com, "A:");

      // call PETSc to solve linear system by iterative method; "inner linear iteration"
      ierr = KSPSetOperators(SSAKSP, A, A, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
      ierr = KSPSolve(SSAKSP, SSARHS, SSAX); CHKERRQ(ierr); // SOLVE

      // report to standard out about iteration
      ierr = KSPGetConvergedReason(SSAKSP, &reason); CHKERRQ(ierr);
      if (reason < 0) {
        ierr = PetscPrintf(grid.com, 
            "\n\n\nPISM ERROR:  KSPSolve() reports 'diverged'; reason = %d = '%s';\n"
                  "  see PETSc man page for KSPGetConvergedReason();   ENDING ...\n\n",
            reason,KSPConvergedReasons[reason]); CHKERRQ(ierr);
        PetscEnd();
      }
      ierr = KSPGetIterationNumber(SSAKSP, &its); CHKERRQ(ierr);
      PetscPrintf(grid.com, "S:%d,%d: ", its, reason);

      // Communicate so that we have stencil width for evaluation of effective
      //   viscosity on next "outer" iteration (and geometry etc. if done):
      ierr = trivialMove(grid, SSADA, SSAX, vel_ssa); CHKERRQ(ierr);
      ierr = vel_ssa.beginGhostComm(); CHKERRQ(ierr);
      ierr = vel_ssa.endGhostComm(); CHKERRQ(ierr);

      // update viscosity and check for viscosity convergence
      ierr = computeEffectiveViscosity(
                grid, config, ice, ssaStrengthExtend,
                vH, vel_ssa, hardness, epsilon, vNuH); CHKERRQ(ierr);
      ierr = testConvergenceOfNu(grid, vNuH, vNuHOld, &norm, &normChange); CHKERRQ(ierr);
      PetscPrintf(grid.com, "|nu|_2, |Delta nu|_2/|nu|_2 = %10.3e %10.3e\n", 
                         norm, normChange/norm);

      *numiter = k + 1;
      if (norm == 0 || normChange / norm < ssaRelativeTolerance) goto done;

    } // end of the "outer loop" (index: k)

    if (epsilon > 0.0) {
       // this has no units; epsilon goes up by this ratio when previous value failed
       const PetscScalar DEFAULT_EPSILON_MULTIPLIER_SSA = 4.0;
       ierr = PetscPrintf(grid.com,
			 "WARNING: Effective viscosity not converged after %d iterations\n"
			 "\twith epsilon=%8.2e. Retrying with epsilon * %8.2e.\n",
			 ssaMaxIterations, epsilon, DEFAULT_EPSILON_MULTIPLIER_SSA);
       CHKERRQ(ierr);

       ierr = vel_ssa.copy_from(vel_ssa_old); CHKERRQ(ierr);
       epsilon *= DEFAULT_EPSILON_MULTIPLIER_SSA;
    } else {
       SETERRQ1(1, 
         "Effective viscosity not converged after %d iterations; epsilon=0.0.\n"
         "  Stopping.                \n", 
         ssaMaxIterations);
    }

  } // end of the "outer outer loop" (index: l)

  done:

  PetscPrintf(grid.com, "... =%5d outer iterations\n", *numiter);
  return 0;
}


int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm    com;  // won't be used except for rank,size
  PetscMPIInt rank, size;

  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);

  com = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);
  
  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  {  
    NCConfigVariable config, overrides;
    ierr = init_config(com, rank, config, overrides); CHKERRQ(ierr);

    PetscTruth i_set, usage_set, help_set;
    ierr = PetscOptionsHasName(PETSC_NULL, "-usage", &usage_set); CHKERRQ(ierr);
    ierr = PetscOptionsHasName(PETSC_NULL, "-help", &help_set); CHKERRQ(ierr);
    if ((usage_set==PETSC_TRUE) || (help_set==PETSC_TRUE)) {
      PetscPrintf(com,
        "\nusage of SSA_TEST:\n"
          "  1) create a PISM output file foo.nc with variables thk, usurf, bed, tauc\n"
          "  2) do 'ssa_test -i foo.nc' or\n"
          "  3) or do 'mpiexec -n 2 ssa_test -i foo.nc -display :0'\n\n");
    }
    
    // get input file name and open
    char filename[PETSC_MAX_PATH_LEN];
    ierr = PetscOptionsGetString(PETSC_NULL, "-i", filename, 
			       PETSC_MAX_PATH_LEN, &i_set); CHKERRQ(ierr);
    if (i_set==PETSC_FALSE) {
      PetscPrintf(com,
        "\nSSA_TEST ERROR:  requires PISM-written NetCDF file as input ... ending!\n\n");
      PetscEnd();
    }

    IceGrid grid(com, rank, size, config);
    ierr = PetscPrintf(grid.com, 
        "SSA_TEST\n  initializing from NetCDF file '%s'...\n", filename); CHKERRQ(ierr);
    PISMIO pio(&grid);
    ierr = pio.get_grid(filename); CHKERRQ(ierr); // fails if filename not present
    grid.start_year = grid.year;

    grid.compute_nprocs();
    grid.compute_ownership_ranges();
    ierr = grid.compute_horizontal_spacing(); CHKERRQ(ierr);
    ierr = grid.createDA(); CHKERRQ(ierr);

    ierr = setVerbosityLevel(5); CHKERRQ(ierr);
    ierr = grid.printInfo(1); CHKERRQ(ierr);
    //ierr = grid.printVertLevels(1); CHKERRQ(ierr); 

    // SSA solve vars; note pieces of the SSA Velocity routine are defined in iMssa.cc
    KSP SSAKSP;
    Mat SSAStiffnessMatrix;
    Vec SSAX, SSARHS;  // Global vectors for solution of the linear system
    DA  SSADA;

    // allocate SSA tools
    // mimic IceGrid::createDA() with TRANSPOSE :
    PetscInt dof=2, stencilwidth=1;
    ierr = DACreate2d(grid.com, DA_XYPERIODIC, DA_STENCIL_BOX,
                    grid.My, grid.Mx,
                    grid.Ny, grid.Nx,
                    dof, stencilwidth,
                    grid.procs_y, grid.procs_x,
                    &SSADA); CHKERRQ(ierr);

    ierr = DACreateGlobalVector(SSADA, &SSAX); CHKERRQ(ierr);
    ierr = VecDuplicate(SSAX, &SSARHS); CHKERRQ(ierr);

    ierr = DAGetMatrix(SSADA, MATMPIAIJ, &SSAStiffnessMatrix); CHKERRQ(ierr);
    ierr = MatSetFromOptions(SSAStiffnessMatrix);CHKERRQ(ierr);

    ierr = KSPCreate(grid.com, &SSAKSP); CHKERRQ(ierr);
    // the default PC type somehow is ILU, which now fails (?) while block jacobi
    //   seems to work; runtime options can override (see test J in vfnow.py)
    PC pc;
    ierr = KSPGetPC(SSAKSP,&pc); CHKERRQ(ierr);
    ierr = PCSetType(pc,PCBJACOBI); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(SSAKSP); CHKERRQ(ierr);

    IceFlowLaw *ice = NULL;
    IceFlowLawFactory ice_factory(com, NULL, config);
    string ice_type = ICE_GPBLD;
    ice_factory.setType(ICE_GPBLD); // set the default type
    ierr = ice_factory.setFromOptions(); CHKERRQ(ierr);
    ice_factory.create(&ice);
    //ierr = ice->printInfo(1); CHKERRQ(ierr);

    IceBasalResistancePlasticLaw basal(
           config.get("plastic_regularization") / secpera, 
           config.get_flag("do_pseudo_plastic_till"),
           config.get("pseudo_plastic_q"),
           config.get("pseudo_plastic_uthreshold") / secpera);
    //ierr = basal.printInfo(1,grid.com); CHKERRQ(ierr);

    SSAStrengthExtension  ssaStrengthExtend;

    IceModelVec2S vh, vH, vbed, vtauc;
    IceModelVec2Mask vMask;
    IceModelVec2V vel_ssa, vel_ssa_old;
    IceModelVec2Stag uvbar;
    const PetscInt WIDE_STENCIL = 2;

    // ice upper surface elevation
    ierr = vh.create(grid, "usurf", true, WIDE_STENCIL); CHKERRQ(ierr);
    ierr = vh.set_attrs("diagnostic", "ice upper surface elevation",
          "m", "surface_altitude"); CHKERRQ(ierr);
    // land ice thickness
    ierr = vH.create(grid, "thk", true, WIDE_STENCIL); CHKERRQ(ierr);
    ierr = vH.set_attrs("model_state", "land ice thickness",
          "m", "land_ice_thickness"); CHKERRQ(ierr);
    ierr = vH.set_attr("valid_min", 0.0); CHKERRQ(ierr);
    // bedrock surface elevation
    ierr = vbed.create(grid, "topg", true, WIDE_STENCIL); CHKERRQ(ierr);
    ierr = vbed.set_attrs("model_state", "bedrock surface elevation",
          "m", "bedrock_altitude"); CHKERRQ(ierr);
    // yield stress for basal till (plastic or pseudo-plastic model)
    ierr = vtauc.create(grid, "tauc", false); CHKERRQ(ierr);
    ierr = vtauc.set_attrs("diagnostic", 
          "yield stress for basal till (plastic or pseudo-plastic model)",
          "Pa", ""); CHKERRQ(ierr);

    // grounded_dragging_floating integer mask
    ierr = vMask.create(grid, "mask", true, WIDE_STENCIL); CHKERRQ(ierr);
    ierr = vMask.set_attrs("model_state", "grounded_dragging_floating integer mask",
			 "", ""); CHKERRQ(ierr);
    vector<double> mask_values(6);
    mask_values[0] = MASK_ICE_FREE_BEDROCK;
    mask_values[1] = MASK_SHEET;
    mask_values[2] = MASK_DRAGGING_SHEET;
    mask_values[3] = MASK_FLOATING;
    mask_values[4] = MASK_ICE_FREE_OCEAN;
    mask_values[5] = MASK_OCEAN_AT_TIME_0;
    ierr = vMask.set_attr("flag_values", mask_values); CHKERRQ(ierr);
    ierr = vMask.set_attr("flag_meanings",
			"ice_free_bedrock sheet dragging_sheet floating ice_free_ocean ocean_at_time_zero");
		  CHKERRQ(ierr);
    vMask.output_data_type = NC_BYTE;

    ierr = vel_ssa.create(grid, "bar_ssa", true, WIDE_STENCIL); // components are ubar_ssa and vbar_ssa
    ierr = vel_ssa.set_attrs(
       "internal_restart", "SSA model ice velocity in the X direction",
       "m s-1", "", 0); CHKERRQ(ierr);
    ierr = vel_ssa.set_attrs(
       "internal_restart", "SSA model ice velocity in the Y direction",
       "m s-1", "", 1); CHKERRQ(ierr);
    ierr = vel_ssa.set_glaciological_units("m year-1"); CHKERRQ(ierr);
    ierr = vel_ssa.set(0.0); CHKERRQ(ierr);

    ierr = vel_ssa_old.create(grid, "bar_ssa_old", true, WIDE_STENCIL); CHKERRQ(ierr);
    ierr = vel_ssa_old.set_attrs(
      "internal", "latest SSA velocities for rapid re-solve of SSA equations", "","");
      CHKERRQ(ierr);

    ierr = uvbar.create(grid, "uvbar", true); CHKERRQ(ierr);
    ierr = uvbar.set_attrs("internal", 
			 "vertically averaged ice velocity, on staggered grid offset in X direction,"
			 " from SIA, in the X direction",
			 "m s-1", "", 0); CHKERRQ(ierr);
    ierr = uvbar.set_attrs("internal", 
			 "vertically averaged ice velocity, on staggered grid offset in Y direction,"
			 " from SIA, in the Y direction",
			 "m s-1", "", 1); CHKERRQ(ierr);

    // 2d work vectors
    const int nWork2d = 4;
    IceModelVec2S vWork2d[nWork2d];
    for (int j = 0; j < nWork2d; j++) {
      char namestr[30], longnamestr[70];
      snprintf(namestr, sizeof(namestr),
               "work_vector_%d", j);
      snprintf(longnamestr, sizeof(longnamestr), 
               "temporary working space 2d vector %d", j);
      ierr = vWork2d[j].create(grid, namestr, true, WIDE_STENCIL); CHKERRQ(ierr);
      ierr = vWork2d[j].set_attrs(
         "internal", longnamestr, "","");
      CHKERRQ(ierr);
    }
  
    // read fields
    NCTool nc(grid.com, grid.rank);
    ierr = nc.open_for_reading(filename); CHKERRQ(ierr);
    int last_record;
    ierr = nc.get_dim_length("t", &last_record); CHKERRQ(ierr);
    last_record -= 1;
    ierr = vh.read(filename, last_record); CHKERRQ(ierr);
    ierr = vH.read(filename, last_record); CHKERRQ(ierr);
    ierr = vbed.read(filename, last_record); CHKERRQ(ierr);
    ierr = vMask.read(filename, last_record); CHKERRQ(ierr);
    ierr = vtauc.read(filename, last_record); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);

    bool show = true;
    const PetscInt  window = 400;
    if (show) {
      //ierr = vh.view(window);  CHKERRQ(ierr);
      ierr = vH.view(window);  CHKERRQ(ierr);
      //ierr = vbed.view(window);  CHKERRQ(ierr);
      ierr = vtauc.view(window);  CHKERRQ(ierr);
      //ierr = vMask.view(window);  CHKERRQ(ierr);
      PetscPrintf(grid.com,"  before SSA: showing fields in X windows for 5 seconds ...\n");
      ierr = PetscSleep(5); CHKERRQ(ierr);
    }

    IceModelVec2S vNuDefault[2] = {vWork2d[0], vWork2d[1]};
    PetscInt numiter = -999;
    ierr = velocitySSA(grid, config, ice, ssaStrengthExtend, basal,
                       vtauc, vh, vH, vbed, vMask, uvbar, vWork2d,
                       vNuDefault, vel_ssa_old,
                       SSADA, SSAX, SSARHS, SSAStiffnessMatrix, SSAKSP,
                       vel_ssa, &numiter); CHKERRQ(ierr); 

    show = true;
    if (show) {
      ierr = vel_ssa.view(window);  CHKERRQ(ierr);
      PetscPrintf(grid.com,"[after SSA: showing components of velocity solution in X windows for 5 seconds ...]\n");
      ierr = PetscSleep(5); CHKERRQ(ierr);
    }

    if (ice != NULL)  delete ice;
    ierr = KSPDestroy(SSAKSP); CHKERRQ(ierr);
    ierr = MatDestroy(SSAStiffnessMatrix); CHKERRQ(ierr);
    ierr = VecDestroy(SSAX); CHKERRQ(ierr);
    ierr = VecDestroy(SSARHS); CHKERRQ(ierr);
    ierr = DADestroy(SSADA);CHKERRQ(ierr);

  }
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
