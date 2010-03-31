// Copyright (C) 2004--2010 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <cmath>
#include <cstring>
#include <petscda.h>
#include "iceModel.hh"

//! Each step of SSA uses previously saved values to start iteration; zero them here to start.
PetscErrorCode IceModel::initSSA() {
  PetscErrorCode ierr;
  if (!have_ssa_velocities) {
    ierr = vel_ssa.set(0.0); CHKERRQ(ierr);
  }
  return 0;
}


//! Compute the product of the effective viscosity \f$\nu\f$ and ice thickness \f$H\f$ for the SSA model.
/*! 
In PISM the product \f$\nu H\f$ can be
  - constant, or
  - can be computed with a constant ice hardness \f$\bar B\f$ (temperature-independent)
    but with dependence of the viscosity on the strain rates, or 
  - it can depend on the strain rates and have a vertically-averaged ice hardness.

The flow law in ice stream and ice shelf regions must, for now, be a temperature-dependent Glen
law.  This is the only flow law we know how to convert to ``viscosity form''.  (More general 
forms like Goldsby-Kohlstedt are not yet inverted.)  The viscosity form is
   \f[ \nu(T^*,D) = \frac{1}{2} B(T^*) D^{(1/n)-1}\, D_{ij} \f]
where 
   \f[  D_{ij} = \frac{1}{2} \left(\frac{\partial U_i}{\partial x_j} +
                                   \frac{\partial U_j}{\partial x_i}\right) \f]
is the strain rate tensor and \f$B\f$ is an ice hardness related to 
the ice softness \f$A(T^*)\f$ by
   \f[ B(T^*)=A(T^*)^{-1/n}  \f]
in the case of a temperature dependent Glen-type law.  (Here \f$T^*\f$ is the pressure-adjusted temperature.)

The effective viscosity is then
   \f[ \nu = \frac{\bar B}{2} \left[\left(\frac{\partial u}{\partial x}\right)^2 + 
                               \left(\frac{\partial v}{\partial y}\right)^2 + 
                               \frac{\partial u}{\partial x} \frac{\partial v}{\partial y} + 
                               \frac{1}{4} \left(\frac{\partial u}{\partial y}
                                                 + \frac{\partial v}{\partial x}\right)^2
                               \right]^{(1-n)/(2n)}                                                \f]
where in the temperature-dependent case
   \f[ \bar B = \frac{1}{H}\,\int_b^h B(T^*)\,dz\f]
This integral is approximately computed by the trapezoid rule.

In fact the integral is regularized as described in [\ref SchoofStream].
The regularization constant \f$\epsilon\f$ is an argument to this procedure.

Also we put \f$\bar\nu H = \f$\c constantNuHForSSA anywhere the ice is thinner
than \c min_thickness_SSA.  The geometry is not changed, but this has the effect 
of producing a shelf extension in ice free ocean, which affects the driving stress
and the force balance at the calving front.
 */
PetscErrorCode IceModel::computeEffectiveViscosity(IceModelVec2S vNuH[2], PetscReal epsilon) {
  PetscErrorCode ierr;

  if (leaveNuHAloneSSA == PETSC_TRUE) {
    return 0;
  }

  bool use_constant_nuh_for_ssa = config.get_flag("use_constant_nuh_for_ssa");
  if (use_constant_nuh_for_ssa) {
    // Intended only for debugging, this treats the entire domain as though it was the strength extension
    // (i.e. strength does not even depend on thickness)
    PetscReal nuH = ssaStrengthExtend.notional_strength();
    ierr = vNuH[0].set(nuH); CHKERRQ(ierr);
    ierr = vNuH[1].set(nuH); CHKERRQ(ierr);
    return 0;
  }

  // We need to compute integrated effective viscosity (\bar\nu * H).
  // It is locally determined by the strain rates and temperature field.
  PetscScalar *Tij = NULL, *Toffset = NULL, **H, **nuH[2];
  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = vNuH[0].get_array(nuH[0]); CHKERRQ(ierr);
  ierr = vNuH[1].get_array(nuH[1]); CHKERRQ(ierr);

  PISMVector2 **uv;
  ierr = vel_ssa.get_array(uv); CHKERRQ(ierr);

  PetscScalar *Enthij, *Enthoffset;
  PolyThermalGPBLDIce *gpbldi = NULL;
  if (config.get_flag("do_cold_ice_methods") == false) {
    gpbldi = dynamic_cast<PolyThermalGPBLDIce*>(ice);
    if (!gpbldi) {
      PetscPrintf(grid.com,
        "do_cold_ice_methods == false in IceModel::computeEffectiveViscosity()\n"
        "   but not using PolyThermalGPBLDIce ... ending ....\n");
      PetscEnd();
    }
  } else {
    Tij = new PetscScalar[grid.Mz];
    Toffset = new PetscScalar[grid.Mz];
  }

  ierr = Enth3.begin_access(); CHKERRQ(ierr);

  for (PetscInt o=0; o<2; ++o) {
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        if (H[i][j] < ssaStrengthExtend.min_thickness_for_extension()) {
          // Extends strength of SSA (i.e. nuH coeff) into the ice free region.
          //  Does not add or subtract ice mass.
          nuH[o][i][j] = ssaStrengthExtend.notional_strength();
        } else {
          const PetscInt      oi = 1-o, oj=o;
          const PetscScalar   dx = grid.dx, 
                              dy = grid.dy;
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
          const PetscScalar myH = 0.5 * (H[i][j] + H[i+oi][j+oj]);

	  ierr = Enth3.getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
	  ierr = Enth3.getInternalColumn(i+oi,j+oj,&Enthoffset); CHKERRQ(ierr);

          if (config.get_flag("do_cold_ice_methods") == true) {
	    for (int k = 0; k < grid.Mz; ++k) {
	      ierr = EC->getAbsTemp(Enthij[k],
				    EC->getPressureFromDepth(H[i][j]-grid.zlevels[k]),
				    Tij[k]); CHKERRQ(ierr);
	      ierr = EC->getAbsTemp(Enthoffset[k],
				    EC->getPressureFromDepth(H[i][j]-grid.zlevels[k]),
				    Toffset[k]); CHKERRQ(ierr);
	    }

            nuH[o][i][j] = ice->effectiveViscosityColumn(
                                myH, grid.kBelowHeight(myH), grid.zlevels,
                                u_x, u_y, v_x, v_y, Tij, Toffset);
          } else {
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
        } // end of if (H[i][j] < ssaStrengthExtend.min_thickness_for_extension()) { ... } else {
      } // j
    } // i
  } // o
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vNuH[0].end_access(); CHKERRQ(ierr);
  ierr = vNuH[1].end_access(); CHKERRQ(ierr);

  ierr = vel_ssa.end_access(); CHKERRQ(ierr);

  ierr = Enth3.end_access(); CHKERRQ(ierr);

  delete[] Tij;
  delete[] Toffset;

  // Some communication
  ierr = vNuH[0].beginGhostComm(); CHKERRQ(ierr);
  ierr = vNuH[0].endGhostComm(); CHKERRQ(ierr);
  ierr = vNuH[1].beginGhostComm(); CHKERRQ(ierr);
  ierr = vNuH[1].endGhostComm(); CHKERRQ(ierr);
  return 0;
}



/*!
Compares saved to current product of vertically-averaged viscosity times height,
the quantity denoted \f$\bar\nu H\f$.  This comparison is used to determine
if the outer iteration is converged.

Verification and PST experiments
suggest that an \f$L^1\f$ criterion for convergence is best.  For verification
there seems to be little difference, presumably because the solutions are smooth
and the norms are roughly equivalent on a subspace of smooth functions.  For PST,
the \f$L^1\f$ criterion gives faster runs with essentially the same results.
Presumably that is because rapid (temporal and spatial) variation in 
\f$\bar\nu H\f$ occurs at margins, occupying very few horizontal grid cells.
For the significant (e.g.~in terms of flux) parts of the flow, it is o.k. to ignore
a bit of bad behavior at these few places, and \f$L^1\f$ ignores it more than
\f$L^2\f$ (much less \f$L^\infty\f$, which might not work at all).
 */
PetscErrorCode IceModel::testConvergenceOfNu(IceModelVec2S vNuH[2], IceModelVec2S vNuHOld[2],
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


//! Assemble the matrix (left-hand side) for the numerical approximation of the SSA equations.
/*! 
COMMENT FIXME:  Finish plastic till paper first.  Then fix these comments to match.

The SSA equations are in their clearest form
    \f[ - \frac{\partial T_{ij}}{\partial x_j} + \tau_{(b)i} = f_i \f]
where \f$i,j\f$ range over \f$x,y\f$, \f$T_{ij}\f$ is a depth-integrated viscous stress tensor 
(i.e. equation (2.6) in [\ref SchoofStream] , and following [\ref Morland], 
and \f$\tau_{(b)i}\f$ are the components of the basal shear stress.  
Also \f$f_i\f$ is the driving shear stress \f$f_i = - \rho g H \frac{\partial h}{\partial x_i}\f$.  
These equations determine velocity in a more-or-less elliptic equation manner.  Here \f$H\f$ 
is the ice thickness and \f$h\f$ is the elevation of the surface of the ice.

More specifically, the SSA equations are
\latexonly
\def\ddt#1{\ensuremath{\frac{\partial #1}{\partial t}}}
\def\ddx#1{\ensuremath{\frac{\partial #1}{\partial x}}}
\def\ddy#1{\ensuremath{\frac{\partial #1}{\partial y}}}
\begin{align*}
  - 2 \ddx{}\left[\nu H \left(2 \ddx{u} + \ddy{v}\right)\right]
        - \ddy{}\left[\nu H \left(\ddy{u} + \ddx{v}\right)\right]
        + \tau_{(b)x}  &=  - \rho g H \ddx{h}, \\
    - \ddx{}\left[\nu H \left(\ddy{u} + \ddx{v}\right)\right]
      - 2 \ddy{}\left[\nu H \left(\ddx{u} + 2 \ddy{v}\right)\right]
        + \tau_{(b)y}  &=  - \rho g H \ddy{h}, 
\end{align*}
\endlatexonly
 where \f$u\f$ is the \f$x\f$-component of the velocity and \f$v\f$ is the \f$y\f$-component 
of the velocity.  Note \f$\nu\f$ is the vertically-averaged effective viscosity of the ice.  

For ice shelves \f$\tau_{(b)i} = 0\f$ [\ref MacAyealetal].  For ice streams with a basal 
till modelled as a plastic material, \f$\tau_{(b)i} = \tau_c u_i/|\mathbf{u}|\f$ where 
\f$\mathbf{u} = (u,v)\f$, \f$|\mathbf{u}| = \left(u^2 + v^2\right)^{1/2}\f$, and \f$\tau_c\f$ 
is the yield stress of the till [\ref SchoofStream].  More generally,
ice streams can be modeled with a pseudo-plastic basal till.  This includes assuming the
basal till is a linearly-viscous material, \f$\tau_{(b)i} = \beta u_i\f$ where \f$\beta\f$ 
is the basal drag (friction) parameter [\ref HulbeMacAyeal].
See IceBasalResistancePlasticLaw::drag().

Note that the basal shear stress appears on the \em left side of the above system.  
We believe this is crucial, because of its effect on the spectrum of the linear 
approximations of each stage.  The effect on spectrum is clearest in the linearly-viscous
till case (i.e. [\ref HulbeMacAyeal]) but there seems to be an analogous effect in the 
plastic till case [\ref SchoofStream].

This method assembles the matrix for the left side of the SSA equations.  The numerical method 
is finite difference.  In particular [FIXME: explain f.d. approxs, esp. mixed derivatives]
 */
PetscErrorCode IceModel::assembleSSAMatrix(bool includeBasalShear, IceModelVec2S vNuH[2], Mat A) {
  const PetscInt  Mx=grid.Mx, My=grid.My, M=2*My;
  const PetscScalar   dx=grid.dx, dy=grid.dy;
  // next constant not too sensitive, but must match value in assembleSSARhs():
  const PetscScalar   scaling = 1.0e9;  // comparable to typical beta for an ice stream
  PetscErrorCode  ierr;
  PetscScalar     **nuH[2], **tauc;
  PISMVector2     **uvssa;

  ierr = MatZeroEntries(A); CHKERRQ(ierr);

  PetscReal beta_shelves_drag_too = config.get("beta_shelves_drag_too");

  /* matrix assembly loop */
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vtauc.get_array(tauc); CHKERRQ(ierr);

  ierr = vel_ssa.get_array(uvssa); CHKERRQ(ierr);

  ierr = vNuH[0].get_array(nuH[0]); CHKERRQ(ierr);
  ierr = vNuH[1].get_array(nuH[1]); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscInt J = 2*j;
      const PetscInt rowU = i*M + J;
      const PetscInt rowV = i*M + J+1;
      const PismMask mask_value = vMask.value(i,j);
      if (mask_value == MASK_SHEET) {
        // set diagonal entry to one; RHS entry will be known (e.g. SIA) velocity;
        //   this is where boundary value to SSA is set
        ierr = MatSetValues(A, 1, &rowU, 1, &rowU, &scaling, INSERT_VALUES); CHKERRQ(ierr);
        ierr = MatSetValues(A, 1, &rowV, 1, &rowV, &scaling, INSERT_VALUES); CHKERRQ(ierr);
      } else {
        const PetscInt im = (i + Mx - 1) % Mx, ip = (i + 1) % Mx;
        const PetscInt Jm = 2 * ((j + My - 1) % My), Jp = 2 * ((j + 1) % My);
        const PetscScalar dx2 = dx*dx, d4 = dx*dy*4, dy2 = dy*dy;
        /* Provide shorthand for the following staggered coefficients
        *      c11
        *  c00     c01
        *      c10
        * Note that the positive i (x) direction is right and the positive j (y)
        * direction is up. */

        // thickness H is incorporated here because nuH[][][] is actually 
        //   viscosity times thickness
        const PetscScalar c00 = nuH[0][i-1][j];
        const PetscScalar c01 = nuH[0][i][j];
        const PetscScalar c10 = nuH[1][i][j-1];
        const PetscScalar c11 = nuH[1][i][j];

        const PetscInt stencilSize = 13;
        /* The locations of the stencil points for the U equation */
        const PetscInt colU[] = {
          /*       */ i*M+Jp,
          im*M+Jp+1,  i*M+Jp+1,   ip*M+Jp+1,
          im*M+J,     i*M+J,      ip*M+J,
          im*M+J+1,               ip*M+J+1,
          /*       */ i*M+Jm,
          im*M+Jm+1,  i*M+Jm+1,   ip*M+Jm+1};
        /* The values at those points */
        PetscScalar valU[] = {
          /*               */ -c11/dy2,
          (2*c00+c11)/d4,     2*(c00-c01)/d4,                 -(2*c01+c11)/d4,
          -4*c00/dx2,         4*(c01+c00)/dx2+(c11+c10)/dy2,  -4*c01/dx2,
          (c11-c10)/d4,                                       (c10-c11)/d4,
          /*               */ -c10/dy2,
          -(2*c00+c10)/d4,    2*(c01-c00)/d4,                 (2*c01+c10)/d4 };

        /* The locations of the stencil points for the V equation */
        const PetscInt colV[] = {
          im*M+Jp,        i*M+Jp,     ip*M+Jp,
          /*           */ i*M+Jp+1,
          im*M+J,                     ip*M+J,
          im*M+J+1,       i*M+J+1,    ip*M+J+1,
          im*M+Jm,        i*M+Jm,     ip*M+Jm,
          /*           */ i*M+Jm+1 };
        /* The values at those points */
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
          valU[5] += basalDragx(tauc, uvssa, i, j);
          valV[7] += basalDragy(tauc, uvssa, i, j);
        }

        // make shelf drag a little bit if desired
        if ((shelvesDragToo == PETSC_TRUE) && (mask_value == MASK_FLOATING)) {
          //ierr = verbPrintf(1,grid.com,"... SHELF IS DRAGGING ..."); CHKERRQ(ierr);
          valU[5] += beta_shelves_drag_too;
          valV[7] += beta_shelves_drag_too;
        }

        ierr = MatSetValues(A, 1, &rowU, stencilSize, colU, valU, INSERT_VALUES); CHKERRQ(ierr);
        ierr = MatSetValues(A, 1, &rowV, stencilSize, colV, valV, INSERT_VALUES); CHKERRQ(ierr);
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


//! Computes the right-hand side ("rhs") of the linear problem for the SSA equations.
/*! 
The right side of the SSA equations is just the driving stress term
   \f[ - \rho g H \nabla h \f]
because, in particular, the basal stress is put on the left side of the system.  
(See comment for assembleSSAMatrix().)  

The surface slope is computed by centered difference onto the regular grid,
which may use periodic ghosting.  (Optionally the surface slope can be computed 
by only differencing into the grid for points at the edge; see test I.)

Note that the grid points with mask value MASK_SHEET correspond to the trivial 
equations
   \f[ \bar u_{ij} = \frac{uvbar[0][i-1][j] + uvbar[0][i][j]}{2}, \f]
and similarly for \f$\bar v_{ij}\f$.  That is, the vertically-averaged horizontal
velocity is already known for these points because it was computed (on the staggered
grid) using the SIA.

Regarding the first argument, if \c surfGradInward == \c true then we differentiate
the surface \f$h(x,y)\f$ inward from edge of grid at the edge of the grid.  This allows 
certain to make sense on a period grid.  This applies to Test I and Test J, in particular.
 */
PetscErrorCode IceModel::assembleSSARhs(Vec rhs) {
  const PetscInt  M=2*grid.My;

  // next constant not too sensitive, but must match value in assembleSSAMatrix():
  const PetscScalar   scaling = 1.0e9;  // comparable to typical beta for an ice stream;

  PetscErrorCode  ierr;
  PetscScalar     **uvbar[2], **taudx, **taudy;

  ierr = VecSet(rhs, 0.0); CHKERRQ(ierr);

  // get driving stress components
  ierr = computeDrivingStress(vWork2d[0],vWork2d[1]); CHKERRQ(ierr); // in iMgeometry.cc
  ierr = vWork2d[0].get_array(taudx); CHKERRQ(ierr);
  ierr = vWork2d[1].get_array(taudy); CHKERRQ(ierr);

  /* rhs (= right-hand side) assembly loop */
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vuvbar[0].get_array(uvbar[0]); CHKERRQ(ierr);
  ierr = vuvbar[1].get_array(uvbar[1]); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscInt J = 2*j;
      const PetscInt rowU = i*M + J;
      const PetscInt rowV = i*M + J+1;
      if (vMask.value(i,j) == MASK_SHEET) {
        ierr = VecSetValue(rhs, rowU, scaling * 0.5*(uvbar[0][i-1][j] + uvbar[0][i][j]),
                           INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValue(rhs, rowV, scaling * 0.5*(uvbar[1][i][j-1] + uvbar[1][i][j]),
                           INSERT_VALUES); CHKERRQ(ierr);
      } else {
	// usual case: use already computed driving stress
	ierr = VecSetValue(rhs, rowU, taudx[i][j], INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(rhs, rowV, taudy[i][j], INSERT_VALUES); CHKERRQ(ierr);          
      }
    }
  }
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vuvbar[0].end_access(); CHKERRQ(ierr);
  ierr = vuvbar[1].end_access(); CHKERRQ(ierr);

  ierr = vWork2d[0].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[1].end_access(); CHKERRQ(ierr);

  ierr = VecAssemblyBegin(rhs); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(rhs); CHKERRQ(ierr);
  return 0;
}


//! Move the solution to the SSA system into a \c Vec which is \c DA distributed.
/*!
Since the parallel layout of the vector \c x in the \c KSP representation of the 
linear SSA system does not in general have anything to do with the \c DA-based
vectors, for the rest of \c IceModel, we must scatter the entire vector to all 
processors.  In particular, this procedure runs once for each outer iteration
in order to put the velocities where needed to compute effective viscosity.
 */
PetscErrorCode IceModel::moveVelocityToDAVectors(Vec x) {
  const PetscInt  M = 2 * grid.My;
  PetscErrorCode  ierr;
  PetscScalar     *uv;
  Vec             xLoc = SSAXLocal;
  PISMVector2     **uv_da;

  ierr = VecScatterBegin(SSAScatterGlobalToLocal, x, xLoc, 
           INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterEnd(SSAScatterGlobalToLocal, x, xLoc, 
           INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = VecGetArray(xLoc, &uv); CHKERRQ(ierr);
  ierr = vel_ssa.get_array(uv_da); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      uv_da[i][j].u = uv[i*M + 2*j];
      uv_da[i][j].v = uv[i*M + 2*j+1];
    }
  }
  ierr = VecRestoreArray(xLoc, &uv); CHKERRQ(ierr);
  ierr = vel_ssa.end_access(); CHKERRQ(ierr);

  // Communicate so that we have stencil width for evaluation of effective viscosity (and geometry)
  ierr = vel_ssa.beginGhostComm(); CHKERRQ(ierr);
  ierr = vel_ssa.endGhostComm(); CHKERRQ(ierr);
  return 0;
}


//! Compute the vertically-averaged horizontal velocity from the shallow shelf approximation (SSA).
/*!
This is the main procedure implementing the SSA.  

The outer loop (over index \c k) is the nonlinear iteration.  In this loop the effective 
viscosity is computed by computeEffectiveViscosity() and then the linear system is 
set up and solved.

This procedure creates a PETSC \c KSP, it calls assembleSSAMatrix() and assembleSSARhs() 
to store the linear system in the \c KSP, and then calls the PETSc procedure KSPSolve()
to solve the linear system.  

Solving the linear system is also a loop, an iteration, but it occurs
inside KSPSolve().  This inner loop is controlled by PETSc but the user can set the option 
<tt>-ksp_rtol</tt>, in particular, and that tolerance controls when the inner loop terminates.

Note that <tt>-ksp_type</tt> can be used to choose the \c KSP.  This will set 
which type of linear iterative method is used.  The KSP is important because 
          omfv = 
the eigenvalues of the linearized SSA are not well understood, but they 
determine the convergence of this linear, inner, iteration.  The default KSP 
is GMRES(30).

Note that <tt>-ksp_pc</tt> will set which preconditioner to use; the default 
is ILU.  A well-chosen preconditioner can put the eigenvalues in the right
place so that the KSP can converge quickly.  The preconditioner is also
important because it will behave differently on different numbers of
processors.  If the user wants the results of SSA calculations to be 
independent of the number of processors, then <tt>-ksp_pc none</tt> should
be used.

If you want to test different KSP methods, it may be helpful to see how many 
iterations were necessary.  Use <tt>-ksp_monitor</tt> or <tt>-d k</tt>.
Initial testing implies that CGS takes roughly half the iterations of 
GMRES(30), but is not significantly faster because the iterations are each 
roughly twice as slow.  Furthermore, ILU and BJACOBI seem roughly equivalent
as preconditioners.

The outer loop terminates when the effective viscosity is no longer changing 
much, according to the tolerance set by the option <tt>-ssa_rtol</tt>.  (The 
outer loop also terminates when a maximum number of iterations is exceeded.)

In truth there is an outer outer loop (over index \c l).  This one manages an
attempt to over-regularize the effective viscosity so that the nonlinear
iteration (the "outer" loop over \c k) has a chance to converge.

  // We need to save the velocity from the last time step since we may have to
  // restart the iteration with larger values of epsilon.

 */
PetscErrorCode IceModel::velocitySSA(PetscInt *numiter) {
  PetscErrorCode ierr;
  bool dosnes;
  IceModelVec2S vNuDefault[2] = {vWork2d[0], vWork2d[1]}; // already allocated space

  ierr = PISMOptionsIsSet("-ssa_bueler", dosnes); CHKERRQ(ierr);
  if (dosnes == PETSC_TRUE) {
    ierr = velocitySSA_SNES(vNuDefault, numiter); CHKERRQ(ierr);
  } else {
    ierr = velocitySSA(vNuDefault, numiter); CHKERRQ(ierr);
  }

  return 0;
}


//! Call this one directly if control over allocation of vNuH[2] is needed (e.g. test J).
/*!
Generally use velocitySSA(PetscInt*) unless you have a vNuH[2] already stored away.
 */
PetscErrorCode IceModel::velocitySSA(IceModelVec2S vNuH[2], PetscInt *numiter) {
  PetscErrorCode ierr;
  KSP ksp = SSAKSP;
  Mat A = SSAStiffnessMatrix;
  Vec x = SSAX, rhs = SSARHS; // solve  A x = rhs
  IceModelVec2S vNuHOld[2] = {vWork2d[2], vWork2d[3]};
  PetscReal   norm, normChange;
  PetscInt    its;
  KSPConvergedReason  reason;

  stdout_ssa = "";
  
  PetscReal ssaRelativeTolerance = config.get("ssa_relative_convergence"),
            epsilon              = config.get("epsilon_ssa");

  PetscInt ssaMaxIterations = static_cast<PetscInt>(config.get("max_iterations_ssa"));
  
  ierr = vel_ssa.copy_to(vel_ssa_old); CHKERRQ(ierr);

  // computation of RHS only needs to be done once; does not depend on solution;
  //   but matrix changes under nonlinear iteration (loop over k below)
  ierr = assembleSSARhs(rhs); CHKERRQ(ierr);

  for (PetscInt l=0; ; ++l) { // iterate with increasing regularization parameter
    ierr = computeEffectiveViscosity(vNuH, epsilon); CHKERRQ(ierr);
    ierr = update_nu_viewers(vNuH, vNuHOld, true); CHKERRQ(ierr);
    // iterate on effective viscosity: "outer nonlinear iteration":
    for (PetscInt k=0; k<ssaMaxIterations; ++k) { 
    
      // in preparation of measuring change of effective viscosity:
      ierr = vNuH[0].copy_to(vNuHOld[0]); CHKERRQ(ierr);
      ierr = vNuH[1].copy_to(vNuHOld[1]); CHKERRQ(ierr);

      // reform matrix, which depends on updated viscosity
      if (getVerbosityLevel() > 2) {
        char tempstr[50] = "";  snprintf(tempstr,50, "  %d,%2d:", l, k);
        stdout_ssa += tempstr;
      }
      ierr = assembleSSAMatrix(true, vNuH, A); CHKERRQ(ierr);
      if (getVerbosityLevel() > 2)
        stdout_ssa += "A:";

      // call PETSc to solve linear system by iterative method; "inner linear iteration"
      ierr = KSPSetOperators(ksp, A, A, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
      ierr = KSPSolve(ksp, rhs, x); CHKERRQ(ierr); // SOLVE
      ierr = KSPGetConvergedReason(ksp, &reason); CHKERRQ(ierr);
      if (reason < 0) {
        ierr = verbPrintf(1,grid.com, 
            "\n\n\nPISM ERROR:  KSPSolve() reports 'diverged'; reason = %d = '%s';\n"
                  "  see PETSc man page for KSPGetConvergedReason();   ENDING ...\n\n",
            reason,KSPConvergedReasons[reason]); CHKERRQ(ierr);
        PetscEnd();
      }
      ierr = KSPGetIterationNumber(ksp, &its); CHKERRQ(ierr);
      if (getVerbosityLevel() > 2) {
        char tempstr[50] = "";  snprintf(tempstr,50, "S:%d,%d: ", its, reason);
        stdout_ssa += tempstr;
      }

      // finish iteration and report to standard out
      ierr = moveVelocityToDAVectors(x); CHKERRQ(ierr);
      ierr = computeEffectiveViscosity(vNuH, epsilon); CHKERRQ(ierr);
      ierr = testConvergenceOfNu(vNuH, vNuHOld, &norm, &normChange); CHKERRQ(ierr);
      if (getVerbosityLevel() > 2) {
        // verbose output includes norm and change for each outer iteration
        char tempstr[100] = "";
        snprintf(tempstr,100, "|nu|_2, |Delta nu|_2/|nu|_2 = %10.3e %10.3e\n", 
                         norm, normChange/norm);
        stdout_ssa += tempstr;
      }

      ierr = update_nu_viewers(vNuH, vNuHOld, true); CHKERRQ(ierr);
      *numiter = k + 1;
      if (norm == 0 || normChange / norm < ssaRelativeTolerance) goto done;

    } // end of the "outer loop" (index: k)

    if (epsilon > 0.0) {
       // this has no units; epsilon goes up by this ratio when previous value failed
       const PetscScalar DEFAULT_EPSILON_MULTIPLIER_SSA = 4.0;
       ierr = verbPrintf(1,grid.com,
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

  if (getVerbosityLevel() > 2) {
    char tempstr[50] = "";
    snprintf(tempstr,50, "... =%5d outer iterations", *numiter);
    stdout_ssa += tempstr;
  } else if (getVerbosityLevel() == 2) {
    // at default verbosity, just record last normchange and iterations
    char tempstr[50] = "";
    snprintf(tempstr,50, "%5d outer iterations", *numiter);
    stdout_ssa += tempstr;
  }
  if (getVerbosityLevel() >= 2)
    stdout_ssa = "  SSA: " + stdout_ssa;
  if (ssaSystemToASCIIMatlab == PETSC_TRUE) {
    ierr = writeSSAsystemMatlab(vNuH); CHKERRQ(ierr);
  }

  return 0;
}


//! At all SSA points, update the velocity field.
/*!
Once the vertically-averaged velocity field is computed by the SSA, this 
procedure updates the three-dimensional horizontal velocities \f$u\f$ and
\f$v\f$.  Note that \f$w\f$ gets updated later by 
vertVelocityFromIncompressibility().  The three-dimensional velocity field
is needed, for example, so that the temperature equation can include advection.
Basal velocities also get updated.

Here is where the flag do_superpose controlled by option <tt>-super</tt> applies.
If do_superpose is true then the just-computed velocity \f$v\f$ from the SSA is
combined, in convex combination, to the stored velocity \f$u\f$ from the SIA 
computation:
   \f[U = f(|v|)\, u + \left(1-f(|v|)\right)\, v.\f]
Here
   \f[ f(|v|) = 1 - (2/\pi) \arctan(10^{-4} |v|^2) \f]
is a function which decreases smoothly from 1 for \f$|v| = 0\f$ to 0 as
\f$|v|\f$ becomes significantly larger than 100 m/a.
 */
PetscErrorCode IceModel::broadcastSSAVelocity(bool updateVelocityAtDepth) {

  PetscErrorCode ierr;
  PetscScalar **ubar, **vbar, **uvbar[2];
  PetscScalar *u, *v;
  
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vubar.get_array(ubar); CHKERRQ(ierr);
  ierr = vvbar.get_array(vbar); CHKERRQ(ierr);

  PISMVector2 **uvssa, **bvel;
  ierr = vel_ssa.get_array(uvssa); CHKERRQ(ierr);
  ierr = vel_basal.get_array(bvel); CHKERRQ(ierr);

  ierr = u3.begin_access(); CHKERRQ(ierr);
  ierr = v3.begin_access(); CHKERRQ(ierr);
  ierr = vuvbar[0].get_array(uvbar[0]); CHKERRQ(ierr);
  ierr = vuvbar[1].get_array(uvbar[1]); CHKERRQ(ierr);

  const PetscScalar inC_fofv = 1.0e-4 * PetscSqr(secpera),
                    outC_fofv = 2.0 / pi;

  bool do_superpose = config.get_flag("do_superpose");

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vMask.value(i,j) != MASK_SHEET) {
        // combine velocities if desired (and not floating)
        const bool addVels = ( do_superpose && (vMask.value(i,j) == MASK_DRAGGING_SHEET) );
        PetscScalar fv = 0.0, omfv = 1.0;  // case of formulas below where ssa
                                           // speed is infinity; i.e. when !addVels
                                           // we just pass through the SSA velocity
        if (addVels) {
          const PetscScalar c2 = PetscSqr(uvssa[i][j].u) + PetscSqr(uvssa[i][j].v);
          omfv = outC_fofv * atan(inC_fofv * c2);
          fv = 1.0 - omfv;
        }

        // update 3D velocity; u,v were from SIA
        if (updateVelocityAtDepth) {
          ierr = u3.getInternalColumn(i,j,&u); CHKERRQ(ierr); // returns pointer
          ierr = v3.getInternalColumn(i,j,&v); CHKERRQ(ierr);
          for (PetscInt k=0; k<grid.Mz; ++k) {
            u[k] = (addVels) ? fv * u[k] + omfv * uvssa[i][j].u : uvssa[i][j].u;
            v[k] = (addVels) ? fv * v[k] + omfv * uvssa[i][j].v : uvssa[i][j].v;
          }
        }

        // update basal velocity; ub,vb were from SIA
        bvel[i][j].u = (addVels) ? fv * bvel[i][j].u + omfv * uvssa[i][j].u : uvssa[i][j].u;
        bvel[i][j].v = (addVels) ? fv * bvel[i][j].v + omfv * uvssa[i][j].v : uvssa[i][j].v;
        
        // also update ubar,vbar by adding SIA contribution, interpolated from 
        //   staggered grid
        const PetscScalar ubarSIA = 0.5*(uvbar[0][i-1][j] + uvbar[0][i][j]),
                          vbarSIA = 0.5*(uvbar[1][i][j-1] + uvbar[1][i][j]);
        ubar[i][j] = (addVels) ? fv * ubarSIA + omfv * uvssa[i][j].u : uvssa[i][j].u;
        vbar[i][j] = (addVels) ? fv * vbarSIA + omfv * uvssa[i][j].v : uvssa[i][j].v;

      }
    }
  }

  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vubar.end_access(); CHKERRQ(ierr);
  ierr = vvbar.end_access(); CHKERRQ(ierr);
  ierr = vel_ssa.end_access(); CHKERRQ(ierr);
  ierr = vel_basal.end_access(); CHKERRQ(ierr);
  ierr = u3.end_access(); CHKERRQ(ierr);
  ierr = v3.end_access(); CHKERRQ(ierr);
  ierr = vuvbar[0].end_access(); CHKERRQ(ierr);
  ierr = vuvbar[1].end_access(); CHKERRQ(ierr);

  return 0;
}


//! At SSA points, correct the previously-computed basal frictional heating.
/*!
Ice shelves have zero basal friction heating.
 */
PetscErrorCode IceModel::correctBasalFrictionalHeating() {
  PetscErrorCode  ierr;
  PetscScalar **Rb, **tauc;

  bool use_ssa_velocity = config.get_flag("use_ssa_velocity");

  PISMVector2 **bvel;
  ierr = vel_basal.get_array(bvel); CHKERRQ(ierr);
  ierr = vRb.get_array(Rb); CHKERRQ(ierr);
  ierr = vtauc.get_array(tauc); CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vMask.is_floating(i,j)) {
        Rb[i][j] = 0.0;
      }
      if ((vMask.value(i,j) == MASK_DRAGGING_SHEET) && use_ssa_velocity) {
        // note basalDrag[x|y]() produces a coefficient, not a stress;
        //   uses *updated* ub,vb if do_superpose == TRUE
        const PetscScalar 
	  basal_stress_x = - basalDragx(tauc, bvel, i, j) * bvel[i][j].u,
	  basal_stress_y = - basalDragy(tauc, bvel, i, j) * bvel[i][j].v;
	Rb[i][j] = - basal_stress_x * bvel[i][j].u - basal_stress_y * bvel[i][j].v;
      }
      // otherwise leave SIA-computed value alone
    }
  }

  ierr = vel_basal.end_access(); CHKERRQ(ierr);
  ierr = vtauc.end_access(); CHKERRQ(ierr);
  ierr = vRb.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);

  return 0;
}


//! At SSA points, correct the previously-computed volume strain heating (dissipation heating).
/*!
Documented in \ref BBssasliding.
 */
PetscErrorCode IceModel::correctSigma() {
  PetscErrorCode  ierr;
  PetscScalar **H;
  PetscScalar *Sigma, *E;

  double enhancement_factor = config.get("enhancement_factor");

  bool do_superpose = config.get_flag("do_superpose");

  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  
  PISMVector2 **uvssa;
  ierr = vel_ssa.get_array(uvssa); CHKERRQ(ierr);

  ierr = Sigma3.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);

  const PetscScalar dx = grid.dx, 
                    dy = grid.dy;
  // next constant is the form of regularization used by C. Schoof 2006 "A variational
  // approach to ice streams" J Fluid Mech 556 pp 227--251
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vMask.value(i,j) != MASK_SHEET) {
        // note ubar_ssa and vbar_ssa in vel_ssa *are* communicated for differencing by last
        //   call to moveVelocityToDAVectors()
        // apply glaciological-superposition-to-low-order if desired (and not floating)
        bool addVels = ( do_superpose && (vMask.value(i,j) == MASK_DRAGGING_SHEET) );
        PetscScalar fv = 0.0, omfv = 1.0;  // case of formulas below where ssa
                                           // speed is infinity; i.e. when !addVels
                                           // we just pass through the SSA velocity
        if (addVels) {
          const PetscScalar c2peryear = PetscSqr(uvssa[i][j].u*secpera)
                                           + PetscSqr(uvssa[i][j].v*secpera);
          omfv = (2.0/pi) * atan(1.0e-4 * c2peryear);
          fv = 1.0 - omfv;
        }
        const PetscScalar 
                u_x   = (uvssa[i+1][j].u - uvssa[i-1][j].u)/(2*dx),
                u_y   = (uvssa[i][j+1].u - uvssa[i][j-1].u)/(2*dy),
                v_x   = (uvssa[i+1][j].v - uvssa[i-1][j].v)/(2*dx),
                v_y   = (uvssa[i][j+1].v - uvssa[i][j-1].v)/(2*dy),
                D2ssa = PetscSqr(u_x) + PetscSqr(v_y) + u_x * v_y
                          + PetscSqr(0.5*(u_y + v_x));
        // get valid pointers to column of Sigma, T values
        ierr = Sigma3.getInternalColumn(i,j,&Sigma); CHKERRQ(ierr);
        ierr = Enth3.getInternalColumn(i,j,&E); CHKERRQ(ierr);
        const PetscInt ks = grid.kBelowHeight(H[i][j]);
        for (PetscInt k=0; k<ks; ++k) {
          // use hydrostatic pressure; presumably this is not quite right in context
          //   of shelves and streams; here we hard-wire the Glen law
	  PetscScalar T;
	  ierr = EC->getAbsTemp(E[k], EC->getPressureFromDepth(H[i][j]-grid.zlevels[k]), T); CHKERRQ(ierr);
          const PetscScalar
            n_glen  = ice->exponent(),
            Sig_pow = (1.0 + n_glen) / (2.0 * n_glen),
            Tstar   = T + ice->beta_CC_grad * (H[i][j] - grid.zlevels[k]),
            // Use pressure-adjusted temperature and account for the enhancement factor.
            //   Note, enhancement factor is not used in SSA anyway.
            //   Should we get rid of it completely?  If not, what is most consistent here?
            BofT    = ice->hardnessParameter(Tstar) * pow(enhancement_factor,-1/n_glen);
          if (addVels) {
            const PetscScalar D2sia = pow(Sigma[k] / (2 * BofT), 1.0 / Sig_pow);
            Sigma[k] = 2.0 * BofT * pow(fv*fv*D2sia + omfv*omfv*D2ssa, Sig_pow);
          } else { // floating (or grounded SSA sans super)
            Sigma[k] = 2.0 * BofT * pow(D2ssa, Sig_pow);
          }
        }
        for (PetscInt k=ks+1; k<grid.Mz; ++k) {
          Sigma[k] = 0.0;
        }
      }
      // otherwise leave SIA-computed value alone
    }
  }

  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);

  ierr = vel_ssa.end_access(); CHKERRQ(ierr);

  ierr = Sigma3.end_access(); CHKERRQ(ierr);
  ierr = Enth3.end_access(); CHKERRQ(ierr);

  return 0;
}

