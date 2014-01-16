/* Copyright (C) 2013, 2014 Jed Brown and the PISM Authors

   This file is part of PISM.

   PISM is free software; you can redistribute it and/or modify it under the
   terms of the GNU General Public License as published by the Free Software
   Foundation; either version 3 of the License, or (at your option) any later
   version.

   PISM is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
   FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
   details.

   You should have received a copy of the GNU General Public License
   along with PISM; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <assert.h>
#include <petscsnes.h>
#include <petscmat.h>
#include <petscdmda.h>
#include <petscsys.h>

#if defined __SSE2__
#  include <emmintrin.h>
#endif

/* The SSE2 kernels are only for PetscScalar=double on architectures that support it */
#define USE_SSE2_KERNELS (!defined NO_SSE2 && !defined PETSC_USE_COMPLEX && !defined PETSC_USE_REAL_SINGLE && defined __SSE2__)

#include "FE3DTools.h"

#include "Blatter_implementation.h"

/*! \brief Compute the square of a number */
static PetscScalar Sqr(PetscScalar a) {return a*a;}

/*! \file Blatter_implementation.c

  This file contains the implementation of the \f$Q_1\f$ 3D FEM solver for
  the Blatter-Pattyn stress balance system.

  Define the Blatter effective strain rate tensor \f$M\f$
  \f[
  M =
  \left(
  \begin{array}{lll}
  4 u_x + 2 v_y & u_y + v_x & u_z\\
  u_y + v_x & 2 u_x + 4 v_y & v_z
  \end{array}
  \right).
  \f]

  Then the Blatter system of equations reads
  \f[
  \nabla\cdot (\eta M) + \rho g \nabla s = 0,
  \f]
  or
  \f{align*}{
  \left(\eta (4 u_x + 2 v_y)\right)_x + \left(\eta (u_y + v_x)\right)_y + \left(\eta u_z\right)_z + \rho g s_x &= 0,\\
  \left(\eta (2 u_x + 4 v_y)\right)_y + \left(\eta (u_y + v_x)\right)_x + \left(\eta v_z\right)_z + \rho g s_y &= 0,
  \f}
  where \f$\eta\f$ is the effective viscosity
  \f[
  \eta = \frac B 2 \left(\gamma + \epsilon^2\right)^{\displaystyle\frac{1-n}{2n}}
  \f]
  and \f$B\f$ and \f$\epsilon\f$ are ice hardness and the regularizing parameter, respectively.

  Also, \f$\gamma\f$ is the second invariant of the strain rate tensor, defined by
  \f[
  \gamma(u,v) = u_{x}^2 + v_{y}^2 + v_{y} \cdot u_{x} + \frac14{(u_{y}+v_{x})^2} + \frac14{u_{z}^2} + \frac14{v_{z}^2}
  \f]
  See compute_nonlinearity() for the computation of \f$\gamma\f$ and \f$\eta\f$.
*/

/*!
  There is one compile-time option:
  `NO_SSE2`: If the host supports SSE2, we use integration code that
  has been vectorized with SSE2 intrinsics, unless this macro is
  defined. The intrinsics speed up integration by about 30% on my
  architecture (P8700, gcc-4.5 snapshot) (JB).

*/

typedef PetscErrorCode (*DMDASNESJacobianLocal)(DMDALocalInfo*, void*, Mat, Mat, MatStructure*, void*);
typedef PetscErrorCode (*DMDASNESFunctionLocal)(DMDALocalInfo*, void*, void*, void*);

static void compute_surface_gradient(PetscReal dchi[4][4][2], const PrmNode parameters[],
                                     PetscReal dx, PetscReal dy, PetscReal ds[4][2]);

static void compute_nodal_z_coordinates(const PrmNode parameters[], PetscInt k, PetscInt zm, PetscReal zn[]);

static PetscErrorCode BlatterQ1_restriction_hook(DM fine,
                                                 Mat mrestrict, Vec rscale, Mat inject,
                                                 DM coarse, void *ctx);


/*! \brief Set up the DM and allocate storage for model parameters on the current grid level.
 */
static PetscErrorCode BlatterQ1_setup_level(BlatterQ1Ctx *ctx, DM dm)
{
  PetscErrorCode ierr;
  PetscInt refinelevel, coarsenlevel, level, Mx, My, Mz, mx, my, stencil_width;
  DMDAStencilType stencil_type;
  const PetscInt *lx, *ly;

  /* Storage (DMDA and Vec) for 2D parameters */
  DM da2prm;
  Vec parameters;
  /* Storage (DMDA and Vec) for 3D ice hardness */
  DM da3;
  Vec hardness;

  MPI_Comm comm;
  ierr = PetscObjectGetComm((PetscObject)dm, &comm); CHKERRQ(ierr);

  PetscFunctionBegin;
  ierr = DMDAGetInfo(dm, NULL,	/* dimensions */
                     &Mz, &My, &Mx, /* grid size */
                     0, &my, &mx,   /* number of processors in each direction */
                     NULL,          /* number of degrees of freedom */
                     &stencil_width,
                     NULL, NULL, NULL, /* types of ghost nodes at the boundary */
                     &stencil_type); CHKERRQ(ierr);

  ierr = DMDAGetOwnershipRanges(dm, NULL, &ly, &lx); CHKERRQ(ierr);

  ierr = DMGetRefineLevel(dm, &refinelevel); CHKERRQ(ierr);
  ierr = DMGetCoarsenLevel(dm, &coarsenlevel); CHKERRQ(ierr);
  level = refinelevel - coarsenlevel;

  ierr = DMDACreate2d(comm,
                      DMDA_BOUNDARY_PERIODIC,
                      DMDA_BOUNDARY_PERIODIC,
                      stencil_type,
                      My, Mx, my, mx, /* grid size */
                      sizeof(PrmNode)/sizeof(PetscScalar), /* number of parameters per map-plane location */
                      stencil_width,
                      ly, /* number of nodes per processor*/
                      lx, /* ditto */
                      &da2prm); CHKERRQ(ierr);

  {
    ierr = PetscPrintf(comm,
                       "Level %D domain size (m) %8.2g x %8.2g, num elements %3d x %3d x %3d (%8d), size (m) %g x %g\n",
                       level, ctx->Lx, ctx->Ly, Mx, My, Mz, Mx*My*Mz, ctx->Lx / Mx, ctx->Ly / My); CHKERRQ(ierr);
  }

  ierr = DMCreateLocalVector(da2prm, &parameters); CHKERRQ(ierr);

  ierr = PetscObjectCompose((PetscObject)dm, "DMDA_2D", (PetscObject)da2prm); CHKERRQ(ierr);
  ierr = PetscObjectCompose((PetscObject)dm, "DMDA_2D_Vec", (PetscObject)parameters); CHKERRQ(ierr);

  ierr = DMDestroy(&da2prm); CHKERRQ(ierr);
  ierr = VecDestroy(&parameters); CHKERRQ(ierr);


  ierr = DMDACreate3d(comm,
                      DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_PERIODIC, DMDA_BOUNDARY_PERIODIC,
                      DMDA_STENCIL_BOX,
                      Mz, My, Mx,
                      1, my, mx, /* number of processors in z, y, x directions. (Always *one* in the z-direction.) */
                      1,         /* number of degrees of freedom per node (one) */
                      1,         /* stencil width */
                      NULL, ly, lx, /* number of nodes per processor */
                      &da3); CHKERRQ(ierr);

  ierr = DMCreateLocalVector(da3, &hardness); CHKERRQ(ierr);

  ierr = PetscObjectCompose((PetscObject)dm, "DMDA_3D", (PetscObject)da3); CHKERRQ(ierr);
  ierr = PetscObjectCompose((PetscObject)dm, "DMDA_3D_Vec", (PetscObject)hardness); CHKERRQ(ierr);

  ierr = DMDestroy(&da3); CHKERRQ(ierr);
  ierr = VecDestroy(&hardness); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


/*! \brief Create the restriction matrix.
 *
 * The result of this call is attached to `dm_fine` under `mat_name`.
 *
 * \param[in] dm_fine DM corresponding to the fine grid
 * \param[in] dm_coarse DM corresponding to the coarse grid
 * \param[in] dm_name name of the DM ("DMDA_2D" or "DMDA_3D")
 * \param[in] mat_name name to use when attaching the interpolation matrix to `dm_fine`
 */
static PetscErrorCode BlatterQ1_create_restriction(DM dm_fine, DM dm_coarse,
                                                     const char dm_name[],
                                                     const char mat_name[]) {
  PetscErrorCode ierr;
  DM da2_fine, da2_coarse;
  Mat mat;

  /* 1. get the DM for parameters from the fine grid DM */
  ierr = PetscObjectQuery((PetscObject)dm_fine, dm_name,
                          (PetscObject*)&da2_fine); CHKERRQ(ierr);
  if (!da2_fine) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
                          "No %s composed with given DMDA", dm_name);

  /* 2. get the DM for parameters from the coarse grid DM */
  ierr = PetscObjectQuery((PetscObject)dm_coarse, dm_name,
                          (PetscObject*)&da2_coarse); CHKERRQ(ierr);
  if (!da2_coarse) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
                            "No %s composed with given DMDA", dm_name);

  /* call DMCreateInterpolation */
  ierr = DMCreateInterpolation(da2_coarse, da2_fine,
                               &mat, PETSC_NULL); CHKERRQ(ierr);

  /* attach to the fine grid DM */
  ierr = PetscObjectCompose((PetscObject)dm_fine, mat_name,
                            (PetscObject)mat); CHKERRQ(ierr);
  ierr = MatDestroy(&mat); CHKERRQ(ierr);

  return 0;
}

/*! \brief Grid coarsening hook.
 *
 * This hook is called *once* when SNES sets up the next coarse level.
 *
 * This hook does three things:
 * - Set up the DM for the newly created coarse level.
 * - Set up the matrix type on the coarsest level to allow using
 *   direct solvers for the coarse problem.
 * - Set up the interpolation matrix that will be used by the
 *   restriction hook to set model parameters on the new coarse level.
 *
 * See BlatterQ1_restriction_hook().
 */
static PetscErrorCode BlatterQ1_coarsening_hook(DM dm_fine, DM dm_coarse, void *ctx)
{
  PetscErrorCode ierr;
  BlatterQ1Ctx *blatter_ctx = (BlatterQ1Ctx*)ctx;
  PetscInt rlevel, clevel;

  PetscFunctionBegin;
  ierr = BlatterQ1_setup_level(blatter_ctx, dm_coarse); CHKERRQ(ierr);

  ierr = DMGetRefineLevel(dm_coarse, &rlevel); CHKERRQ(ierr);
  ierr = DMGetCoarsenLevel(dm_coarse, &clevel); CHKERRQ(ierr);
  if (rlevel-clevel == 0) {
    ierr = DMSetMatType(dm_coarse, MATAIJ); CHKERRQ(ierr);
  }

  ierr = DMCoarsenHookAdd(dm_coarse, BlatterQ1_coarsening_hook, BlatterQ1_restriction_hook,
                          ctx); CHKERRQ(ierr);

  /* create 2D interpolation */
  ierr = BlatterQ1_create_restriction(dm_fine, dm_coarse,
                                      "DMDA_2D", "DMDA_2D_Restriction"); CHKERRQ(ierr);

  /* create 3D interpolation */
  ierr = BlatterQ1_create_restriction(dm_fine, dm_coarse,
                                      "DMDA_3D", "DMDA_3D_Restriction"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*! \brief Restrict model parameters from the `fine` grid onto the `coarse` grid.
 *
 * This function uses the restriction matrix created by BlatterQ1_coarsening_hook().
 */
static PetscErrorCode BlatterQ1_restrict(DM fine, DM coarse,
                                         const char dm_name[],
                                         const char mat_name[],
                                         const char vec_name[])
{
  PetscErrorCode ierr;
  Vec X_fine, X_fine_global, X_coarse, X_coarse_global;
  DM da2_fine, da2_coarse;
  Mat mat;

  PetscFunctionBegin;

  /* get the restriction matrix from the fine grid DM */
  ierr = PetscObjectQuery((PetscObject)fine, mat_name,
                          (PetscObject*)&mat); CHKERRQ(ierr);
  if (!mat) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
                    "No DMDA_Restriction composed with given DMDA");

  /* get the 2D DMDA from the fine grid DM */
  ierr = PetscObjectQuery((PetscObject)fine, dm_name,
                          (PetscObject*)&da2_fine); CHKERRQ(ierr);
  if (!da2_fine) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
                         "No DMDA_Vec composed with given DMDA");

  /* get the storage vector from the fine grid DM */
  ierr = PetscObjectQuery((PetscObject)fine, vec_name,
                          (PetscObject*)&X_fine); CHKERRQ(ierr);
  if (!X_fine) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
                       "No DMDA_Vec composed with given DMDA");

  /* get the 2D DMDA from the coarse grid DM */
  ierr = PetscObjectQuery((PetscObject)coarse, dm_name,
                          (PetscObject*)&da2_coarse); CHKERRQ(ierr);
  if (!da2_coarse) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
                           "No DMDA_Vec composed with given DMDA");

  /* get the storage vector from the coarse grid DM */
  ierr = PetscObjectQuery((PetscObject)coarse, vec_name,
                          (PetscObject*)&X_coarse); CHKERRQ(ierr);
  if (!X_coarse) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
                         "No DMDA_Vec composed with given DMDA");

  ierr = DMGetGlobalVector(da2_fine, &X_fine_global); CHKERRQ(ierr);
  ierr = DMGetGlobalVector(da2_coarse, &X_coarse_global); CHKERRQ(ierr);

  ierr = DMLocalToGlobalBegin(da2_fine, X_fine, INSERT_VALUES, X_fine_global); CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(da2_fine, X_fine, INSERT_VALUES, X_fine_global); CHKERRQ(ierr);

  ierr = MatRestrict(mat, X_fine_global, X_coarse_global); CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(da2_coarse, X_coarse_global, INSERT_VALUES, X_coarse); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da2_coarse, X_coarse_global, INSERT_VALUES, X_coarse); CHKERRQ(ierr);

  ierr = DMRestoreGlobalVector(da2_fine, &X_fine_global); CHKERRQ(ierr);
  ierr = DMRestoreGlobalVector(da2_coarse, &X_coarse_global); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


/*! \brief Restrict model parameters to the grid corresponding to the next coarse level.

  This function computes the product
  \f[ X_{\text{coarse}} = R X_{\text{fine}}, \f]

  extracting necessary information from DM objects managing fine and coarse grids.

  This is called once per SNESSolve().
 */
static PetscErrorCode BlatterQ1_restriction_hook(DM fine,
                                                 Mat mrestrict, Vec rscale, Mat inject,
                                                 DM coarse, void *ctx)
{
  PetscErrorCode ierr;

  /* Get rid of "unused argument" warnings: */
  (void) mrestrict;
  (void) rscale;
  (void) inject;
  (void) ctx;

  ierr = BlatterQ1_restrict(fine, coarse,
                            "DMDA_2D", "DMDA_2D_Restriction", "DMDA_2D_Vec"); CHKERRQ(ierr);

  ierr = BlatterQ1_restrict(fine, coarse,
                            "DMDA_3D", "DMDA_3D_Restriction", "DMDA_3D_Vec"); CHKERRQ(ierr);

  return 0;
}

/*! \brief Get the pointer to the 2D array storing models parameters.

  The input argument da is the 3D DM for the current level, but we
  need the 2D DM managing 2D parameters, so we need to extract it
  first.

  \param[in] da the DM managed by the SNES object
  \param[out] prm pointer to the array
*/
PetscErrorCode BlatterQ1_begin_2D_parameter_access(DM da, PrmNode ***prm)
{
  PetscErrorCode ierr;
  DM             da2prm;
  Vec            X;
  PetscFunctionBegin;

  ierr = PetscObjectQuery((PetscObject)da, "DMDA_2D",
                          (PetscObject*)&da2prm); CHKERRQ(ierr);
  if (!da2prm) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
                       "No DMDA_2D composed with given DMDA");

  ierr = PetscObjectQuery((PetscObject)da, "DMDA_2D_Vec",
                          (PetscObject*)&X); CHKERRQ(ierr);
  if (!X) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
                  "No DMDA_2D_Vec composed with given DMDA");

  ierr = DMDAVecGetArray(da2prm, X, prm); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*! \brief Restore the 2D array of model parameters.

  See BlatterQ1_begin_2D_parameter_access() for details.
*/
PetscErrorCode BlatterQ1_end_2D_parameter_access(DM da, PrmNode ***prm)
{
  PetscErrorCode ierr;
  DM             da2prm;
  Vec            X;
  PetscFunctionBegin;

  ierr = PetscObjectQuery((PetscObject)da, "DMDA_2D",
                          (PetscObject*)&da2prm); CHKERRQ(ierr);
  if (!da2prm) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
                       "No DMDA_2D composed with given DMDA");

  ierr = PetscObjectQuery((PetscObject)da, "DMDA_2D_Vec",
                          (PetscObject*)&X); CHKERRQ(ierr);
  if (!X) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
                  "No DMDA_2D_Vec composed with given DMDA");

  ierr = DMDAVecRestoreArray(da2prm, X, prm); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*! \brief Get the 3D array of ice hardness.
 *
 * This is similar to BlatterQ1_begin_2D_parameter_access(), but for ice
 * hardness.
 */
PetscErrorCode BlatterQ1_begin_hardness_access(DM da, PetscScalar ****hardness)
{
  PetscErrorCode ierr;
  DM             da3;
  Vec            X;
  PetscFunctionBegin;

  ierr = PetscObjectQuery((PetscObject)da, "DMDA3",
                          (PetscObject*)&da3); CHKERRQ(ierr);
  if (!da3) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
                       "No DMDA3 composed with given DMDA");

  ierr = PetscObjectQuery((PetscObject)da, "DMDA3_Vec",
                          (PetscObject*)&X); CHKERRQ(ierr);
  if (!X) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
                  "No DMDA3_Vec composed with given DMDA");

  ierr = DMDAVecGetArray(da3, X, hardness); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*! \brief Restore the 3D array of ice hardness.

  See BlatterQ1_begin_hardness_access() for details.
*/
PetscErrorCode BlatterQ1_end_hardness_access(DM da, PetscScalar ****hardness)
{
  PetscErrorCode ierr;
  DM             da3;
  Vec            X;
  PetscFunctionBegin;

  ierr = PetscObjectQuery((PetscObject)da, "DMDA3",
                          (PetscObject*)&da3); CHKERRQ(ierr);
  if (!da3) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
                       "No DMDA3 composed with given DMDA");

  ierr = PetscObjectQuery((PetscObject)da, "DMDA3_Vec",
                          (PetscObject*)&X); CHKERRQ(ierr);
  if (!X) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
                  "No DMDA3_Vec composed with given DMDA");

  ierr = DMDAVecRestoreArray(da3, X, hardness); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


/*! \brief Set the initial guess. */
/*!
 * FIXME: we need a callback similar to drag() and viscosity().
 */
static PetscErrorCode BlatterQ1_initial_guess(SNES snes, Vec X, void* ctx)
{
  PetscErrorCode ierr;

  /* Get rid of "unused argument" warnings: */
  (void) snes;
  (void) ctx;

  PetscFunctionBegin;

  ierr = VecSet(X, 0.0); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*! \brief Compute the viscosity non-linearity and related quantities.
 *
 * This function is called once for each quadrature point. It uses
 * nodal values of \f$u\f$, \f$v\f$ and the basis expansion of \f$u\f$
 * and \f$v\f$ to compute \f$u\f$, \f$v\f$ at the quadrature point
 *
 * \f[ u(q_j) = \sum_{i=1}^8\phi_i(q_j)\cdot u_j \f]
 *
 * \f[ v(q_j) = \sum_{i=1}^8\phi_i(q_j)\cdot v_j \f]
 *
 * Here \f$q_j\f$ is a quadrature point, \f$u_j\f$ is a nodal value of \f$u\f$.
 *
 * In addition to this, it uses the map from the reference element to a
 * physical element to compute partial derivatives of u,v with respect to x,y.
 *
 * \c gamma is the second invariant \f$ \gamma \f$.
 *
 * \param[in]  ctx Solver's "application context". Provides `ctx->nonlinear.viscosity()`.
 * \param[in]  velocity nodal values of horizontal velocity
 * \param[in]  phi  values of global basis functions \f$\phi\f$ (at the current quadrature point).
 * \param[in]  dphi values of derivatives of global basis functions \f$\phi\f$ (with respect to \f$x,y,z\f$)
 * \param[out] u,v  components of the horizontal velocity (at the current quadrature point)
 * \param[out] du,dv partial derivatives of `u` and `v`
 * \param[out] eta  effective viscosity \f$\eta\f$ (at the current quadrature point)
 * \param[out] deta derivative of eta with respect to gamma \f$\frac{\partial \eta}{\partial \gamma}\f$
 */
static void compute_nonlinearity(BlatterQ1Ctx *ctx,
                                 const Node velocity[restrict],
                                 const PetscReal phi[restrict],
                                 PetscReal dphi[restrict][3],
                                 PetscScalar *restrict u,
                                 PetscScalar *restrict v,
                                 PetscScalar du[restrict],
                                 PetscScalar dv[restrict],
                                 PetscReal *eta,
                                 PetscReal *deta)
{
  PetscInt l, ll;
  PetscScalar second_invariant;

  du[0] = du[1] = du[2] = 0;
  dv[0] = dv[1] = dv[2] = 0;
  *u = 0;
  *v = 0;
  for (l=0; l<8; l++) {
    *u += phi[l] * velocity[l].u;
    *v += phi[l] * velocity[l].v;
    for (ll=0; ll<3; ll++) {
      du[ll] += dphi[l][ll] * velocity[l].u;
      dv[ll] += dphi[l][ll] * velocity[l].v;
    }
  }
  second_invariant = Sqr(du[0]) + Sqr(dv[1]) + du[0]*dv[1] + 0.25*Sqr(du[1]+dv[0]) + 0.25*Sqr(du[2]) + 0.25*Sqr(dv[2]);
  {				/* FIXME: this needs to be an argument */
    PetscReal softness = 4e-25,
      n = 3.0,
      hardness = pow(softness,-1.0/n);
    ctx->nonlinear.viscosity(ctx, hardness, second_invariant, eta, deta);
  }
}


/*! \brief Evaluate the residual.
 *
 *
 * FIXME: I need to document this.
 */
static PetscErrorCode BlatterQ1_residual_local(DMDALocalInfo *info, Node ***velocity, Node ***residual,
                                               BlatterQ1Ctx *ctx)
{
  PetscInt       xs, ys, xm, ym, zm, i, j, k, q, l;
  PetscReal      dx, dy;
  PrmNode        **prm;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  xs = info->zs;
  ys = info->ys;
  xm = info->zm;
  ym = info->ym;
  zm = info->xm;
  dx = ctx->Lx / info->mz;      /* grid spacing in the x direction */
  dy = ctx->Ly / info->my;      /* grid spacing in the y direction */
  ierr = BlatterQ1_begin_2D_parameter_access(info->da, &prm); CHKERRQ(ierr);

  for (i = xs; i < xs + xm; i++) {
    for (j = ys; j < ys + ym; j++) {
      PrmNode parameters[4];            /* 4 nodes (in 2D) */
      PetscReal ds[4][2];       /* 4 quadrature points, 2 dimensions */
      get_nodal_values_2d(prm, i, j, parameters);

      compute_surface_gradient(ctx->Q12D.dchi, parameters, dx, dy, ds);

      for (k = 0; k < zm - 1; k++) {
        PetscInt ls = 0;	/* starting index */
        Node element_velocity[8], *element_residual[8];
        PetscReal zn[8], etabase = 0;

        /* Compute z-coordinates of element nodes: */
        compute_nodal_z_coordinates(parameters, k, zm, zn);

        /* Get nodal values of velocity components: */
        get_nodal_values_3d(velocity, i, j, k, element_velocity);

        /* Get pointers to residual components corresponding to element nodes: */
        get_pointers_to_nodal_values_3d(residual, i, j, k, element_residual);

        if (ctx->no_slip && k == 0) {
          for (l = 0; l < 4; l++) element_velocity[l].u = element_velocity[l].v = 0;
          /* The first 4 basis functions lie on the bottom layer, so their contribution is exactly 0, hence we can skip them */
          ls = 4;
        }

        for (q = 0; q < 8; q++) {   /* for all quadrature points... */
          PetscReal grad_z[3], phi[8], dphi[8][3], jw, eta, deta;
          PetscScalar du[3], dv[3], u, v;

          /* Compute various quantities at this quadrature point in the *physical* element */

          compute_z_gradient(ctx->Q13D.dchi[q], zn, /* inputs (derivatives of shape functions, geometry) */
                             grad_z); /* output (partial derivatives of z with respect to xi,eta,zeta) */

          compute_element_info(ctx->Q13D.chi, ctx->Q13D.dchi, q, dx, dy, grad_z, /* inputs */
                               phi, dphi, &jw); /* outputs (values of shape functions, their derivatives, det(J)*weight */

          compute_nonlinearity(ctx, element_velocity, phi, dphi, /* inputs */
                               &u, &v, du, dv, &eta, &deta); /* outputs u,v, partial derivatives of u,v,
                                                                effective viscosity (eta), derivative of eta with respect to gamma */

          jw /= ctx->rhog;      /* scales residuals to be O(1) */

          if (q == 0) etabase = eta;

          for (l = ls; l < 8; l++) { /* test functions */
            const PetscReal pp = phi[l], *dp = dphi[l];
            element_residual[l]->u += dp[0]*jw*eta*(4.0*du[0] + 2.0*dv[1]) + dp[1]*jw*eta*(du[1] + dv[0]) + dp[2]*jw*eta*du[2] + pp*jw*ctx->rhog*ds[q % 4][0];
            element_residual[l]->v += dp[1]*jw*eta*(2.0*du[0] + 4.0*dv[1]) + dp[0]*jw*eta*(du[1] + dv[0]) + dp[2]*jw*eta*dv[2] + pp*jw*ctx->rhog*ds[q % 4][1];
          }
        }             /* q-loop */

        if (k == 0) { /* we are on a bottom face */
          if (ctx->no_slip) {
            /* Note: Non-Galerkin coarse grid operators are very sensitive to the scaling of Dirichlet boundary
             * conditions.  After shenanigans above, etabase contains the effective viscosity at the closest quadrature
             * point to the bed.  We want the diagonal entry in the Dirichlet condition to have similar magnitude to the
             * diagonal entry corresponding to the adjacent node.  The fundamental scaling of the viscous part is in
             * diagu, diagv below.  This scaling is easy to recognize by considering the finite difference operator after
             * scaling by element size.  The no-slip Dirichlet condition is scaled by this factor, and also in the
             * assembled matrix (see the similar block in BlatterQ1_Jacobian_local).
             *
             * Note that the residual at this Dirichlet node is linear in the state at this node, but also depends
             * (nonlinearly in general) on the neighboring interior nodes through the local viscosity.  This will make
             * a matrix-free Jacobian have extra entries in the corresponding row.  We assemble only the diagonal part,
             * so the solution will exactly satisfy the boundary condition after the first linear iteration.
             */
            const PetscReal dz = PetscRealPart(parameters[0].thickness) / (zm - 1.0);
            const PetscScalar
              diagu = 2*etabase / ctx->rhog*(dx*dy / dz + dx*dz / dy + 4*dy*dz / dx),
              diagv = 2*etabase / ctx->rhog*(dx*dy / dz + 4*dx*dz / dy + dy*dz / dx);
            element_residual[0]->u = ctx->dirichlet_scale*diagu*velocity[i][j][k].u;
            element_residual[0]->v = ctx->dirichlet_scale*diagv*velocity[i][j][k].v;
          } else {              /* Integrate over bottom face to apply boundary condition */

            for (q = 0; q < 4; q++) { /* for each quadrature point on the basal face... */
              const PetscReal jw = 0.25*dx*dy / ctx->rhog, /* FIXME: this Jacobian is incorrect */
                *phi = ctx->Q12D.chi[q];

              PetscScalar u = 0, v = 0, tauc = 0;
              PetscReal taub, dtaub;
              for (l = 0; l < 4; l++) {
                u += phi[l]*element_velocity[l].u;
                v += phi[l]*element_velocity[l].v;
                tauc += phi[l]*parameters[l].tauc;
              }

              ctx->nonlinear.drag(ctx, PetscRealPart(tauc), u, v,
                                  &taub, &dtaub);

              for (l = 0; l < 4; l++) {
                const PetscReal pp = phi[l];
                element_residual[ls + l]->u += pp*jw*taub*u;
                element_residual[ls + l]->v += pp*jw*taub*v;
              }
            } /* end of the quadrature loop */

          } /* end of the generic basal boundary condition */

        } /* end of "if (k == 0)" */
      } /* k-loop */
    } /* j-loop */
  } /* i-loop */

  ierr = BlatterQ1_end_2D_parameter_access(info->da, &prm); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*! \brief Evaluate the Jacobian.
 *
 * FIXME: I need to document this.
 */
static PetscErrorCode BlatterQ1_Jacobian_local(DMDALocalInfo *info, Node ***velocity, Mat A,
                                               Mat B, MatStructure *mstr, BlatterQ1Ctx *ctx)
{
  PetscInt       xs, ys, xm, ym, zm, i, j, k, q, l, ll;
  PetscReal      dx, dy;
  PrmNode        **prm;
  PetscErrorCode ierr;

  (void)A;

  PetscFunctionBegin;
  xs = info->zs;
  ys = info->ys;
  xm = info->zm;
  ym = info->ym;
  zm = info->xm;
  dx = ctx->Lx / info->mz;
  dy = ctx->Ly / info->my;

  ierr = MatZeroEntries(B); CHKERRQ(ierr);

  ierr = BlatterQ1_begin_2D_parameter_access(info->da, &prm); CHKERRQ(ierr);

  for (i = xs; i < xs + xm; i++) {
    for (j = ys; j < ys + ym; j++) {
      PrmNode parameters[4];            /* 4 nodes */
      get_nodal_values_2d(prm, i, j, parameters);

      for (k = 0; k < zm - 1; k++) {
        PetscInt ls = 0;
        Node element_velocity[8];
        PetscReal zn[8], etabase = 0;
        PetscScalar Ke[8*2][8*2];

        PetscMemzero(Ke, sizeof(Ke));

        compute_nodal_z_coordinates(parameters, k, zm, zn);

        get_nodal_values_3d(velocity, i, j, k, element_velocity);

        /* Ensure that values of velocity components at Dirichlet
           (no-slip) nodes are exactly equal to Dirichlet B.C. values. */
        if (ctx->no_slip && k == 0) {
          for (l = 0; l < 4; l++) {
            element_velocity[l].u = 0.0;
            element_velocity[l].v = 0.0;
          }
          ls = 4;
        }

        for (q = 0; q < 8; q++) {
          PetscReal grad_z[3], phi[8], dphi[8][3], jw, eta, deta;
          PetscScalar du[3], dv[3], u, v;

          /* Compute the gradient of z */
          compute_z_gradient(ctx->Q13D.dchi[q], zn, /* inputs */
                             grad_z);               /* output */

          /* compute values of shape functions, their derivatives, and
             quadrature weights at a quadrature point q. */
          compute_element_info(ctx->Q13D.chi, ctx->Q13D.dchi, q, dx, dy, grad_z, /* inputs */
                               phi, dphi, &jw); /* outputs */

          /* Compute u,v, their derivatives, plus effective viscosity and its
             derivative with respect to gamma. */
          compute_nonlinearity(ctx, element_velocity, phi, dphi, /* inputs */
                               &u, &v, du, dv, &eta, &deta); /* outputs */

          jw /= ctx->rhog;      /* residuals are scaled by this factor */

          if (q == 0) {
            etabase = eta;
          }

          for (l = ls; l < 8; l++) { /* test functions */
            const PetscReal *restrict dp = dphi[l];
#if (USE_SSE2_KERNELS!=0)
            /* gcc (up to my 4.5 snapshot) is really bad at hoisting intrinsics so we do it manually */
            __m128d
              p4 = _mm_set1_pd(4), p2 = _mm_set1_pd(2), p05 = _mm_set1_pd(0.5),
              p42 = _mm_setr_pd(4, 2), p24 = _mm_shuffle_pd(p42, p42, _MM_SHUFFLE2(0, 1)),
              du0 = _mm_set1_pd(du[0]), du1 = _mm_set1_pd(du[1]), du2 = _mm_set1_pd(du[2]),
              dv0 = _mm_set1_pd(dv[0]), dv1 = _mm_set1_pd(dv[1]), dv2 = _mm_set1_pd(dv[2]),
              jweta = _mm_set1_pd(jw*eta), jwdeta = _mm_set1_pd(jw*deta),
              dp0 = _mm_set1_pd(dp[0]), dp1 = _mm_set1_pd(dp[1]), dp2 = _mm_set1_pd(dp[2]),
              dp0jweta = _mm_mul_pd(dp0, jweta), dp1jweta = _mm_mul_pd(dp1, jweta), dp2jweta = _mm_mul_pd(dp2, jweta),
              p4du0p2dv1 = _mm_add_pd(_mm_mul_pd(p4, du0), _mm_mul_pd(p2, dv1)), /* 4 du0 + 2 dv1 */
              p4dv1p2du0 = _mm_add_pd(_mm_mul_pd(p4, dv1), _mm_mul_pd(p2, du0)), /* 4 dv1 + 2 du0 */
              pdu2dv2 = _mm_unpacklo_pd(du2, dv2), /* [du2, dv2] */
              du1pdv0 = _mm_add_pd(du1, dv0), /* du1 + dv0 */
              t1 = _mm_mul_pd(dp0, p4du0p2dv1), /* dp0 (4 du0 + 2 dv1) */
              t2 = _mm_mul_pd(dp1, p4dv1p2du0);                                /* dp1 (4 dv1 + 2 du0) */

#endif
            for (ll = l; ll < 8; ll++) {
              const PetscReal *restrict dpl = dphi[ll];
#if (USE_SSE2_KERNELS==0)
              /* The analytic Jacobian in nice, easy-to-read form */
              {
                PetscScalar dgdu, dgdv;
                /* Compute derivatives of gamma with respect to u and v. */
                dgdu = 2.0*du[0]*dpl[0] + dv[1]*dpl[0] + 0.5*(du[1] + dv[0])*dpl[1] + 0.5*du[2]*dpl[2];
                dgdv = 2.0*dv[1]*dpl[1] + du[0]*dpl[1] + 0.5*(du[1] + dv[0])*dpl[0] + 0.5*dv[2]*dpl[2];
                /* Picard part */
                Ke[l*2 + 0][ll*2 + 0] += dp[0]*jw*eta*4.0*dpl[0] + dp[1]*jw*eta*dpl[1] + dp[2]*jw*eta*dpl[2];
                Ke[l*2 + 0][ll*2 + 1] += dp[0]*jw*eta*2.0*dpl[1] + dp[1]*jw*eta*dpl[0];
                Ke[l*2 + 1][ll*2 + 0] += dp[1]*jw*eta*2.0*dpl[0] + dp[0]*jw*eta*dpl[1];
                Ke[l*2 + 1][ll*2 + 1] += dp[1]*jw*eta*4.0*dpl[1] + dp[0]*jw*eta*dpl[0] + dp[2]*jw*eta*dpl[2];
                /* extra Newton terms */
                Ke[l*2 + 0][ll*2 + 0] += dp[0]*jw*deta*dgdu*(4.0*du[0] + 2.0*dv[1]) + dp[1]*jw*deta*dgdu*(du[1] + dv[0]) + dp[2]*jw*deta*dgdu*du[2];
                Ke[l*2 + 0][ll*2 + 1] += dp[0]*jw*deta*dgdv*(4.0*du[0] + 2.0*dv[1]) + dp[1]*jw*deta*dgdv*(du[1] + dv[0]) + dp[2]*jw*deta*dgdv*du[2];
                Ke[l*2 + 1][ll*2 + 0] += dp[1]*jw*deta*dgdu*(4.0*dv[1] + 2.0*du[0]) + dp[0]*jw*deta*dgdu*(du[1] + dv[0]) + dp[2]*jw*deta*dgdu*dv[2];
                Ke[l*2 + 1][ll*2 + 1] += dp[1]*jw*deta*dgdv*(4.0*dv[1] + 2.0*du[0]) + dp[0]*jw*deta*dgdv*(du[1] + dv[0]) + dp[2]*jw*deta*dgdv*dv[2];
              }
#else
              /* This SSE2 code is an exact replica of above, but uses explicit packed instructions for some speed
               * benefit.  On my hardware, these intrinsics are almost twice as fast as above, reducing total assembly cost
               * by 25 to 30 percent. */
              {
                __m128d
                  keu = _mm_loadu_pd(&Ke[l*2+0][ll*2+0]),
                  kev = _mm_loadu_pd(&Ke[l*2+1][ll*2+0]),
                  dpl01 = _mm_loadu_pd(&dpl[0]), dpl10 = _mm_shuffle_pd(dpl01, dpl01, _MM_SHUFFLE2(0, 1)), dpl2 = _mm_set_sd(dpl[2]),
                  t0, t3, pdgduv;
                keu = _mm_add_pd(keu, _mm_add_pd(_mm_mul_pd(_mm_mul_pd(dp0jweta, p42), dpl01),
                                                 _mm_add_pd(_mm_mul_pd(dp1jweta, dpl10),
                                                            _mm_mul_pd(dp2jweta, dpl2))));
                kev = _mm_add_pd(kev, _mm_add_pd(_mm_mul_pd(_mm_mul_pd(dp1jweta, p24), dpl01),
                                                 _mm_add_pd(_mm_mul_pd(dp0jweta, dpl10),
                                                            _mm_mul_pd(dp2jweta, _mm_shuffle_pd(dpl2, dpl2, _MM_SHUFFLE2(0, 1))))));
                pdgduv = _mm_mul_pd(p05, _mm_add_pd(_mm_add_pd(_mm_mul_pd(p42, _mm_mul_pd(du0, dpl01)),
                                                               _mm_mul_pd(p24, _mm_mul_pd(dv1, dpl01))),
                                                    _mm_add_pd(_mm_mul_pd(du1pdv0, dpl10),
                                                               _mm_mul_pd(pdu2dv2, _mm_set1_pd(dpl[2]))))); /* [dgdu, dgdv] */
                t0 = _mm_mul_pd(jwdeta, pdgduv);  /* jw deta [dgdu, dgdv] */
                t3 = _mm_mul_pd(t0, du1pdv0);     /* t0 (du1 + dv0) */
                _mm_storeu_pd(&Ke[l*2+0][ll*2+0], _mm_add_pd(keu, _mm_add_pd(_mm_mul_pd(t1, t0),
                                                                             _mm_add_pd(_mm_mul_pd(dp1, t3),
                                                                                        _mm_mul_pd(t0, _mm_mul_pd(dp2, du2))))));
                _mm_storeu_pd(&Ke[l*2+1][ll*2+0], _mm_add_pd(kev, _mm_add_pd(_mm_mul_pd(t2, t0),
                                                                             _mm_add_pd(_mm_mul_pd(dp0, t3),
                                                                                        _mm_mul_pd(t0, _mm_mul_pd(dp2, dv2))))));
              }
#endif
            } /* ll-loop (trial functions) */
          }   /* l-loop (test functions) */
        }     /* q-loop */

        /* on a bottom face */
        if (k == 0) {
          if (ctx->no_slip) {
            const PetscReal dz = PetscRealPart(parameters[0].thickness) / (zm - 1);
            const PetscScalar diagu = 2*etabase / ctx->rhog*(dx*dy / dz + dx*dz / dy + 4*dy*dz / dx),
              diagv = 2*etabase / ctx->rhog*(dx*dy / dz + 4*dx*dz / dy + dy*dz / dx);
            Ke[0][0] = ctx->dirichlet_scale*diagu;
            Ke[1][1] = ctx->dirichlet_scale*diagv;
          } else {

            for (q = 0; q < 4; q++) {
              const PetscReal jw = 0.25*dx*dy / ctx->rhog, /* FIXME: det(J)*w is wrong here */
                *phi = ctx->Q12D.chi[q];
              PetscScalar u = 0, v = 0, tauc = 0;
              PetscReal taub, dtaub;

              /* Compute u, v, \tau_c at a quadrature point on the bottom face using basis expansions: */
              for (l = 0; l < 4; l++) {
                u += phi[l]*element_velocity[l].u;
                v += phi[l]*element_velocity[l].v;
                tauc += phi[l]*parameters[l].tauc;
              }

              /* Compute the friction coefficient at this quadrature point: */
              ctx->nonlinear.drag(ctx, PetscRealPart(tauc), u, v,
                                  &taub, &dtaub);

              for (l = 0; l < 4; l++) {
                const PetscReal pp = phi[l];
                for (ll = 0; ll < 4; ll++) {
                  const PetscReal ppl = phi[ll];
                  Ke[l*2 + 0][ll*2 + 0] += pp*jw*taub*ppl + pp*jw*dtaub*u*u*ppl;
                  Ke[l*2 + 0][ll*2 + 1] +=                  pp*jw*dtaub*u*v*ppl;
                  Ke[l*2 + 1][ll*2 + 0] +=                  pp*jw*dtaub*v*u*ppl;
                  Ke[l*2 + 1][ll*2 + 1] += pp*jw*taub*ppl + pp*jw*dtaub*v*v*ppl;
                }
              }
            }

          } /* generic basal boundary condition */
        } /* end of "if (k == 0)" */

        /* Set matrix values: */
        {
          const MatStencil rc[8] = {{i, j, k,     0}, {i + 1, j, k,     0}, {i + 1, j + 1, k,     0}, {i, j + 1, k,     0},
                                    {i, j, k + 1, 0}, {i + 1, j, k + 1, 0}, {i + 1, j + 1, k + 1, 0}, {i, j + 1, k + 1, 0}};
          {
            /* fill in lower-triangular part, this is really cheap compared to computing the entries */
            for (l = 0; l < 8; l++) {
              for (ll = l + 1; ll < 8; ll++) {
                Ke[ll*2 + 0][l*2 + 0] = Ke[l*2 + 0][ll*2 + 0];
                Ke[ll*2 + 1][l*2 + 0] = Ke[l*2 + 0][ll*2 + 1];
                Ke[ll*2 + 0][l*2 + 1] = Ke[l*2 + 1][ll*2 + 0];
                Ke[ll*2 + 1][l*2 + 1] = Ke[l*2 + 1][ll*2 + 1];
              }
            }
            ierr = MatSetValuesBlockedStencil(B, 8, rc, 8, rc, &Ke[0][0], ADD_VALUES); CHKERRQ(ierr);
          }
        }

      } /* k-loop */
    } /* j-loop */
  } /* i-loop */
  ierr = BlatterQ1_end_2D_parameter_access(info->da, &prm); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatSetOption(B, MAT_SYMMETRIC, PETSC_TRUE); CHKERRQ(ierr);
  *mstr = SAME_NONZERO_PATTERN;

  PetscFunctionReturn(0);
}

/*! \brief Allocate the fine grid 3D DMDA and set up the SNES object.

  \param[in] com communicator
  \param[in] pism_da PISM-side 2D DMDA (used to get grid info)
  \param[in] Mz number of vertical levels
  \param[in] ctx the application context
  \param[out] result SNES object that will be used with SNESSolve later
*/
PetscErrorCode BlatterQ1_create(MPI_Comm com, DM pism_da,
                                PetscInt Mz,
                                BlatterQ1Ctx *ctx, SNES *result) {
  PetscErrorCode ierr;
  PetscInt dim, Mx, My, Nx, Ny;
  DM da;
  const PetscInt *lx, *ly;

  ierr = DMDAGetInfo(pism_da,
                     &dim,
                     &My,
                     &Mx,
                     NULL, /* Mz */
                     &Ny,  /* number of processors in y-direction */
                     &Nx,  /* number of processors in x-direction */
                     NULL, /* ditto, z-direction */
                     NULL, /* number of degrees of freedom per node */
                     NULL, /* stencil width */
                     NULL, NULL, NULL, /* types of ghost nodes at the boundary */
                     NULL); CHKERRQ(ierr); /* stencil type */
  assert(dim == 2);

  ierr = DMDAGetOwnershipRanges(pism_da, &ly, &lx, NULL);

  ierr = DMDACreate3d(com,
                      DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_PERIODIC, DMDA_BOUNDARY_PERIODIC,
                      DMDA_STENCIL_BOX,
                      Mz, My, Mx,
                      1, Ny, Nx, /* number of processors in z, y, x directions. (Always *one* in the z-direction.) */
                      sizeof(Node)/sizeof(PetscScalar),
                      1, /* stencil width */
                      NULL, ly, lx, /* number of nodes per processor */
                      &da); CHKERRQ(ierr);

  ierr = DMDASetFieldName(da, 0, "x-velocity"); CHKERRQ(ierr);
  ierr = DMDASetFieldName(da, 1, "y-velocity"); CHKERRQ(ierr);

  ierr = BlatterQ1_setup_level(ctx, da); CHKERRQ(ierr);

  /* ADD_VALUES, because BlatterQ1_residual_local contributes to ghosted values. */
  ierr = DMDASNESSetFunctionLocal(da, ADD_VALUES,
                                  (DMDASNESFunctionLocal)BlatterQ1_residual_local,
                                  ctx); CHKERRQ(ierr);
  ierr = DMDASNESSetJacobianLocal(da,
                                  (DMDASNESJacobianLocal)BlatterQ1_Jacobian_local,
                                  ctx); CHKERRQ(ierr);

  ierr = DMCoarsenHookAdd(da, BlatterQ1_coarsening_hook,
                          BlatterQ1_restriction_hook, ctx); CHKERRQ(ierr);
  /* FIXME: do we need a refinement hook? */

  ierr = DMSetApplicationContext(da, ctx); CHKERRQ(ierr);

  ierr = SNESCreate(com, result); CHKERRQ(ierr);
  ierr = SNESSetDM(*result, da); CHKERRQ(ierr);
  ierr = DMDestroy(&da); CHKERRQ(ierr);

  ierr = SNESSetComputeInitialGuess(*result, BlatterQ1_initial_guess, ctx); CHKERRQ(ierr);
  ierr = SNESSetFromOptions(*result); CHKERRQ(ierr);

  return 0;
}

/*! \brief Compute the surface gradient at all 4 quadrature points (in the map plane).

  Note that there are 8 quadrature points in a 3D hexahedral element, but since
  surface elevation does not depend on \f$z\f$, values at the first 4 points (bottom
  "level") are all we need.

  Using the 2D \f$Q_1\f$ basis expansion for the function \f$s(x,y)\f$,

  \f{eqnarray}{
  \frac{\partial s}{\partial x}(q_j) &=& \sum_{i=1}^4 \frac{\partial \phi_i}{\partial x}(q_j) s_i,\\
  \frac{\partial s}{\partial y}(q_j) &=& \sum_{i=1}^4 \frac{\partial \phi_i}{\partial y}(q_j) s_i.
  \f}

  Here \f$s_i\f$ are nodal values of the surface elevation, \f$\phi_i\f$ are 2D
  \f$Q_1\f$ global basis functions and \f$q_j\f$ are quadrature points.

  Now, if the grid has equal spacing in \f$x\f$ and \f$y\f$ directions, then

  \f{eqnarray}{
  \frac{\partial \phi_i}{\partial x} &=& \frac{2}{\Delta x} \frac{\partial \chi_i}{\partial \xi}\\
  \frac{\partial \phi_i}{\partial y} &=& \frac{2}{\Delta y} \frac{\partial \chi_i}{\partial \eta}
  \f}

  This gives
  \f{eqnarray}{
  \frac{\partial s}{\partial x}(q_j) &=& \frac{2}{\Delta x} \sum_{i=1}^4 \frac{\partial \chi_i}{\partial x}(q_j) s_i,\\
  \frac{\partial s}{\partial y}(q_j) &=& \frac{2}{\Delta y} \sum_{i=1}^4 \frac{\partial \chi_i}{\partial y}(q_j) s_i.
  \f}

  Note also that 3D \f$Q_1\f$ element basis functions evaluated at an element face
  reduce to 2D element basis functions, and so this is equivalent to assuming
  that \f$s(x,y,z) = s(x,y)\f$ for all \f$z\f$ and using 3D \f$Q_1\f$ basis expansion
  in this computation.

  \param[in]  dchi values of derivatives of 2D element basis functions \f$\chi\f$ at quadrature points
  \param[in]  parameters 2D parameters at element nodes
  \param[in]  dx,dy grid spacing in x and y directions
  \param[out] ds values of the surface gradient
*/
static void compute_surface_gradient(PetscReal dchi[4][4][2],
                                     const PrmNode parameters[], PetscReal dx, PetscReal dy,
                                     PetscReal ds[4][2])
{
  PetscInt i, q;

  for (q = 0; q < 4; ++q) {
    ds[q][0] = 0.0;
    ds[q][1] = 0.0;

    for (i = 0; i < 4; ++i) {
      ds[q][0] += dchi[q][i][0] * (parameters[i].ice_bottom + parameters[i].thickness);
      ds[q][1] += dchi[q][i][1] * (parameters[i].ice_bottom + parameters[i].thickness);
    }

    ds[q][0] *= 2.0/dx;
    ds[q][1] *= 2.0/dy;
  }
}

/*! \brief Compute z-coordinates of the nodes of a hexahedral element.

  The bottom (top) mesh surface follows the bottom (top) surface of the ice.
  Each "column" of nodes is equally-spaced, independently from others.

  \param[in]  parameters 2D parameters at element nodes
  \param[in]  k  index of the vertical level
  \param[in]  zm total number of vertical levels
  \param[out] zn z-coordinates of element nodes
*/
static void compute_nodal_z_coordinates(const PrmNode parameters[], PetscInt k, PetscInt zm, PetscReal zn[])
{
  const PetscScalar zm1 = zm - 1,
    znl[8] = {parameters[0].ice_bottom + parameters[0].thickness*(PetscScalar)k / zm1,
              parameters[1].ice_bottom + parameters[1].thickness*(PetscScalar)k / zm1,
              parameters[2].ice_bottom + parameters[2].thickness*(PetscScalar)k / zm1,
              parameters[3].ice_bottom + parameters[3].thickness*(PetscScalar)k / zm1,
              parameters[0].ice_bottom + parameters[0].thickness*(PetscScalar)(k + 1) / zm1,
              parameters[1].ice_bottom + parameters[1].thickness*(PetscScalar)(k + 1) / zm1,
              parameters[2].ice_bottom + parameters[2].thickness*(PetscScalar)(k + 1) / zm1,
              parameters[3].ice_bottom + parameters[3].thickness*(PetscScalar)(k + 1) / zm1};
  PetscInt i;
  for (i = 0; i < 8; i++) {
    zn[i] = PetscRealPart(znl[i]);
  }
}
