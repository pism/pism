/* Copyright (C) 2013, 2014 Jed Brown and the PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef _FE3DTools_H_
#define _FE3DTools_H_

#include <petscsys.h>

#if !defined __STDC_VERSION__ || __STDC_VERSION__ < 199901L
#  if defined __cplusplus       /* C++ restrict is nonstandard and compilers have inconsistent rules about where it can be used */
#    define restrict
#  else
#    define restrict PETSC_RESTRICT
#  endif
#endif

/*! \file FE3DTools.h

  This file contains helper macros and functions for setting up 3D
  @f$ Q_1 @f$ finite element systems.

  In the finite element formulation, let @f$ \phi_i @f$ be a trial function,
  i.e. an element of a basis for the space @f$ S_h @f$ approximating the
  solution space and let @f$ \psi @f$ be a test function.

  Then we write

  \f[
  u = \sum_i U_i \phi_i.
  \f]

  Trial functions @f$ \phi_i @f$ are piecewise-trilinear (linear in each of
  @f$ x @f$ , @f$ y @f$ , @f$ z @f$ separately) on elements.

  Instead of defining @f$ \phi_i @f$ for each element, we define them once
  on the **reference element**, the cube
  @f$ [-1,1]\times[-1,1]\times[-1,1] @f$ and use an invertible map from this
  element to a physical hexahedral element.

  **On** the reference element, pre-images of trial functions are
  called @f$ \chi_i @f$ . These @f$ \chi_i @f$ are also known as *element basis
  functions*, as opposed to *global* basis functions @f$ \phi_i @f$ .

  In the Galerkin formulation, test functions @f$ \psi @f$ are the same as
  trial functions @f$ \phi_i @f$ , but it is still helpful to use different
  letters for test and trial functions (this makes notation simpler).
 */


/*! Extract values of the 3D field x at nodes of an element (i,j,k) and store
   them in the array n.

   Use this to *get* values of a field in a FEM assembly loop.
*/
#define get_nodal_values_3d(x,i,j,k,n) do {     \
    (n)[0] = (x)[i][j][k];                      \
    (n)[1] = (x)[i+1][j][k];                    \
    (n)[2] = (x)[i+1][j+1][k];                  \
    (n)[3] = (x)[i][j+1][k];                    \
    (n)[4] = (x)[i][j][k+1];                    \
    (n)[5] = (x)[i+1][j][k+1];                  \
    (n)[6] = (x)[i+1][j+1][k+1];                \
    (n)[7] = (x)[i][j+1][k+1];                  \
  } while (0)

/*! Extract memory addresses corresponding to nodes of an element (i,j,k) in
   the field x and store them in the array n.

   Use this to *set* values of a field in a FEM assembly loop.

   This is used in the resudual evaluation code.
*/
#define get_pointers_to_nodal_values_3d(x,i,j,k,n) do { \
    (n)[0] = &(x)[i][j][k];                     \
    (n)[1] = &(x)[i+1][j][k];                   \
    (n)[2] = &(x)[i+1][j+1][k];                 \
    (n)[3] = &(x)[i][j+1][k];                   \
    (n)[4] = &(x)[i][j][k+1];                   \
    (n)[5] = &(x)[i+1][j][k+1];                 \
    (n)[6] = &(x)[i+1][j+1][k+1];               \
    (n)[7] = &(x)[i][j+1][k+1];                 \
  } while (0)

/*! Extract values of a 2D field x at nodes of an element (i,j) and store them
   in the array n.

   Use this to *get* values of a field in a FEM assembly loop.

   This is used in the code computing boundary integrals.
*/
#define get_nodal_values_2d(x,i,j,n) do {       \
    (n)[0] = (x)[i][j];                         \
    (n)[1] = (x)[i+1][j];                       \
    (n)[2] = (x)[i+1][j+1];                     \
    (n)[3] = (x)[i][j+1];                       \
  } while (0)

/*! \brief Compute partial derivatives of @f$ z(\xi,\eta,\zeta) @f$ with
    respect to @f$ \xi @f$ , @f$ \eta @f$ , and @f$ \zeta @f$ .

  The function @f$ z(\xi,\eta,\zeta) @f$ is defined by the map from the reference
  element to a physical element, i.e.
  \f[
  z (\xi, \eta, \zeta)  =  \sum_{j = 1}^8 z_j \cdot \chi_j (\xi,\eta,\zeta) .
  \f]

  Here @f$ \xi @f$ , @f$ \eta @f$ and @f$ \zeta @f$ are range from @f$ -1 @f$
  to @f$ 1 @f$ , so that the triple @f$ (\xi,\eta,\zeta) @f$ describes an
  arbitrary point in the reference element.

  These derivatives are used to compute the Jacobian of this map (they
  appear in the third column).

  The output of this computation is used in compute_element_info() and nowhere
  else.

  \param[in] dchi derivatives of element basis functions @f$ \phi @f$ with respect to
  @f$ \xi @f$ , @f$ \eta @f$ , @f$ \zeta @f$ .
  \param[in]  zn[] z-coordinates of the nodes of the current element
  \param[out] dz[] partial derivatives of z
 */
void compute_z_gradient(PetscReal dchi[][3], const PetscReal zn[], PetscReal dz[])
{
  PetscInt i;
  dz[0] = dz[1] = dz[2] = 0;
  for (i = 0; i < 8; i++) {
    dz[0] += dchi[i][0] * zn[i];
    dz[1] += dchi[i][1] * zn[i];
    dz[2] += dchi[i][2] * zn[i];
  }
}

/*! Compute temporaries at a quadrature point. */
/*! Compute the following:
   - values of shape functions @f$ \phi_i @f$ at a given quadrature point
   - values of partial derivatives of shape functions @f$ \frac{\partial \phi_i}{\partial x_j} @f$
     at a given quadrature point
   - @f$ det(J)\cdot w @f$ factor (product of the determinant of the Jacobian of the map
     from the reference element and the quadrature weight), at a given quadrature point.
     (These are "modified" quadrature weights.)

   Let $J$ be the jacobian of the map from the reference element.

   This function computes @f$ J @f$ , @f$ det(J) @f$ , and @f$ J^{-1} @f$ . The determinant is then
   multiplied by weights corresponding to the @f$ 2\times2\times2 @f$ Gaussian quadrature
   (weights are hard-wired).

   The inverse @f$ J^{-1} @f$ is used to compute partial derivatives of
   shape functions (with respect to @f$ x @f$ , @f$ y @f$ , @f$ z @f$ ) using
   partial derivatives of reference element basis functions @f$ \chi_i @f$ (with
   respect to @f$ \xi @f$ , @f$ \eta @f$ , @f$ \zeta @f$ ).

   In particular, @f$ \nabla \phi = J^{-1} (\nabla \chi) @f$ .

   **Note:** both the Jacobian and its inverse are stored *transposed*
   here. (For no apparent reason.)

   \todo Quadrature weights are hard-wired!

   \param[in]  chi  values of element basis functions @f$ \chi_i @f$ at quadrature points
   \param[in]  dchi values of partial derivatives of @f$ \chi_i @f$ at quadrature points
   \param[in]  q    index of the quadrature point
   \param[in]  dx,dy horizontal grid spacing (x- and y-direction)
   \param[in]  dz   partial derivatives of z computed by calling compute_z_gradient()
   \param[out] phi  values of shape functions @f$ \phi_i @f$ at the quadrature point q
   \param[out] dphi values of partial derivatives of shape functions @f$ \phi_i @f$ at q
   \param[out] jw @f$ det(J)\cdot w @f$
 */
void compute_element_info(PetscReal chi[8][8],PetscReal dchi[8][8][3],
			  PetscInt q, PetscReal dx, PetscReal dy, const PetscReal dz[restrict],
			  PetscReal phi[restrict],
			  PetscReal dphi[restrict][3],
			  PetscReal *restrict jw)
{
  const PetscReal jac[3][3] = {{dx / 2, 0,      0},
                               {0,      dy / 2, 0},
                               {dz[0],  dz[1],  dz[2]}},
    ijac[3][3] = {{1 / jac[0][0], 0, 0},
                  {0, 1 / jac[1][1], 0},
                  { - jac[2][0] / (jac[0][0]*jac[2][2]), - jac[2][1] / (jac[1][1]*jac[2][2]), 1 / jac[2][2]}},
      jdet = jac[0][0]*jac[1][1]*jac[2][2];

  for (int i = 0; i < 8; i++) {
    const PetscReal *dchi_q = dchi[q][i];
    phi[i] = chi[q][i];
    dphi[i][0] = dchi_q[0]*ijac[0][0] + dchi_q[1]*ijac[1][0] + dchi_q[2]*ijac[2][0]; /* x-derivative */
    dphi[i][1] = dchi_q[0]*ijac[0][1] + dchi_q[1]*ijac[1][1] + dchi_q[2]*ijac[2][1]; /* y-derivative */
    dphi[i][2] = dchi_q[0]*ijac[0][2] + dchi_q[1]*ijac[1][2] + dchi_q[2]*ijac[2][2]; /* z-derivative */
  }
  *jw = 1.0 * jdet;		/* hard-wired quadrature weights */
}


/*! Compute values of shape functions and their derivatives at quadrature points.
 *
 * Given coordinated of the nodes of the 2D @f$ Q_1 @f$ reference element
 * \f{align*}{
 * \xi  &= (-1, 1, 1, -1)\\
 * \eta &= (-1, -1, 1, 1)
 * \f}
 * we define 2D element basis functions
 * \f{align*}{
 * \chi_i(\xi,\eta) &=\frac 1 4 (1 + \xi_i \xi)(1 + \eta_i \eta)
 * \f}
 * 
 * This function pre-computes values of these element basis functions
 * and values of their derivatives at quadrature points for the
 * @f$ 2\times2 @f$ Gaussian quadrature on the reference element.
 * 
 * Note that with this choice of the reference element quadrature points are
 * \f[
 * (\xi^*_i, \eta^*_i) = \left( \frac{\xi_i}{\sqrt3}, \frac{\eta_i}{\sqrt3} \right).
 * \f]
 *
 * The partial derivatives are
 * \f{align*}{
 * \frac{\partial \chi_i}{\partial \xi} &= \frac 1 4 \xi_i (1 + \eta_i \eta)\\
 * \frac{\partial \chi_i}{\partial \eta} &= \frac 1 4 \eta_i (1 + \xi_i \xi)
 * \f}
 *
 * These 2D basis functions are used to approximate boundary integrals.
 */
void initialize_Q12D(PetscReal chi[4][4], PetscReal dchi[4][4][2])
{
  /* Coordinates of the nodes of the reference element: */
  PetscReal xis[4]  = {-1.0,  1.0,  1.0, -1.0};
  PetscReal etas[4] = {-1.0, -1.0,  1.0,  1.0};
  int q, j;

  for (q = 0; q < 4; ++q) {	/* for all quadrature points... */
    /* compute coordinates of this quadrature point */
    PetscReal xi_q = xis[q]/sqrt(3.0), eta_q = etas[q]/sqrt(3.0);

    for (j = 0; j < 4; ++j) {
      /* compute values of shape functions */
      chi[q][j] = 0.25 * (1.0 + xis[j]*xi_q) * (1.0 + etas[j]*eta_q);

      /* compute derivatives of shape functions */
      dchi[q][j][0] = 0.25 *  xis[j] * (1.0 + etas[j] * eta_q);
      dchi[q][j][1] = 0.25 * etas[j] * (1.0 +  xis[j] *  xi_q);
    }
  }
}

/*! Compute values of shape functions and their derivatives at quadrature points.
 *
 * Given coordinated of the nodes of the 3D @f$ Q_1 @f$ reference element
 * \f{align*}{
 * \xi   &= (-1,  1,  1, -1, -1,  1, 1, -1)\\
 * \eta  &= (-1, -1,  1,  1, -1, -1, 1,  1)\\
 * \zeta &= (-1, -1, -1, -1,  1,  1, 1,  1)
 * \f}
 * we define 3D element basis functions
 * \f{align*}{
 * \chi_i(\xi,\eta,\zeta) &=\frac 1 8 (1 + \xi_i \xi)(1 + \eta_i \eta)(1 + \zeta_i \zeta)
 * \f}
 * 
 * This function pre-computes values of these element basis functions
 * and values of their derivatives at quadrature points for the
 * @f$ 2\times2\times2 @f$ Gaussian quadrature on the reference element.
 * 
 * Note that with this choice of the reference element quadrature points are
 * \f[
 * (\xi^*_i, \eta^*_i, \zeta^*_i) = \left( \frac{\xi_i}{\sqrt3}, \frac{\eta_i}{\sqrt3}, \frac{\zeta_i}{\sqrt3} \right).
 * \f]
 *
 * The partial derivatives are
 * \f{align*}{
 * \frac{\partial \chi_i}{\partial \xi} &= \frac 1 8 \xi_i (1 + \eta_i \eta)(1 + \zeta_i \zeta)\\
 * \frac{\partial \chi_i}{\partial \eta} &= \frac 1 8 \eta_i (1 + \xi_i \xi)(1 + \zeta_i \zeta)\\
 * \frac{\partial \chi_i}{\partial \zeta} &= \frac 1 8 \zeta_i (1 + \xi_i \xi)(1 + \eta_i \eta)
 * \f}
 */
void initialize_Q13D(PetscReal chi[8][8], PetscReal dchi[8][8][3])
{
  /* Coordinated of the nodes of the reference element: */
  PetscReal xis[8]   = {-1.0,  1.0,  1.0, -1.0, -1.0,  1.0, 1.0, -1.0};
  PetscReal etas[8]  = {-1.0, -1.0,  1.0,  1.0, -1.0, -1.0, 1.0,  1.0};
  PetscReal zetas[8] = {-1.0, -1.0, -1.0, -1.0,  1.0,  1.0, 1.0,  1.0};
  int q, j;

  for (q = 0; q < 8; ++q) {	/* for all quadrature points... */
    /* compute coordinates of this quadrature point */
    PetscReal xi_q = xis[q]/sqrt(3.0), eta_q = etas[q]/sqrt(3.0), zeta_q = zetas[q]/sqrt(3.0);

    for (j = 0; j < 8; ++j) {
      /* compute values of shape functions */
      chi[q][j] = 0.125 * (1.0 + xis[j]*xi_q) * (1.0 + etas[j]*eta_q) * (1.0 + zetas[j]*zeta_q);

      /* compute derivatives of shape functions */
      dchi[q][j][0] = 0.125 *   xis[j] * (1.0 + etas[j] * eta_q) * (1.0 + zetas[j] * zeta_q);
      dchi[q][j][1] = 0.125 *  etas[j] * (1.0 +  xis[j] *  xi_q) * (1.0 + zetas[j] * zeta_q);
      dchi[q][j][2] = 0.125 * zetas[j] * (1.0 +  xis[j] *  xi_q) * (1.0 +  etas[j] *  eta_q);
    }
  }
}

#endif /* _FE3DTools_H_ */
