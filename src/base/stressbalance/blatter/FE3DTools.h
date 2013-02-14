/* Copyright (C) 2013 Jed Brown and the PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
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

/*! Extract values of the 3D field x at nodes of an element (i,j,k) and store
   them in the array n.

   Use this to *get* values of a field in a FEM assembly loop.
*/
#define get_nodal_values_3d(x,i,j,k,n) do {              \
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
*/
#define get_pointers_to_nodal_values_3d(x,i,j,k,n) do {           \
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
*/
#define get_nodal_values_2d(x,i,j,n) do {               \
    (n)[0] = (x)[i][j];                         \
    (n)[1] = (x)[i+1][j];                       \
    (n)[2] = (x)[i+1][j+1];                     \
    (n)[3] = (x)[i][j+1];                       \
  } while (0)

/*! Compute partial derivatives of z(xi,eta,zeta) with respect to xi, eta, and zeta.

  The function z(xi,eta,zeta) is defined by the map from the reference
  element to a physical element.

  These derivatives are used to compute the Jacobian of this map (they
  appear in the third column).

  The output of this computation is used in compute_element_info() (and nowhere
  else, as far as I can tell).

  \param[in]  dphi derivatives of element basis functions phi with respect to xi, eta, zeta.
  \param[in]  zn[] z-coordinates of the nodes of the current element
  \param[out] dz[] partial derivatives of z
 */
void compute_z_gradient(const PetscReal dphi[][3], const PetscReal zn[], PetscReal dz[])
{
  PetscInt i;
  dz[0] = dz[1] = dz[2] = 0;
  for (i = 0; i < 8; i++) {
    dz[0] += dphi[i][0] * zn[i];
    dz[1] += dphi[i][1] * zn[i];
    dz[2] += dphi[i][2] * zn[i];
  }
}

/*! Compute temporaries at a quadrature point. */
/*! Compute the following:
   - values of shape functions at a given quadrature point
   - values of partial derivatives of shape functions at a given quadrature point
   - det(J)*w factor (product of the determinant of the Jacobian of the map
     from the reference element and the quadrature weight), at a given quadrature point

   Let J be the jacobian of the map from the reference element.

   This function computes \f$J\f$, \f$det(J)\f$, and \f$J^{-1}\f$. The determinant is then
   multiplied by weights corresponding to the 2x2x2 Gaussian quadrature
   (weights are hard-wired).

   The inverse \f$J^{-1}\f$ is used to compute partial derivatives of shape functions
   (with respect to x,y,z) using partial derivatives of reference element basis
   functions (with respect to xi,eta,zeta).

   We have \f$\nabla \phi = J^{-1} (\nabla \chi)\f$.

   Note: both the Jacobian and its inverse are stored *transposed* here. (For
   no apparent reason.)

   FIXME: Quadrature weight are hard-wired!

   \param[in]  chi  values of element basis functions at quadrature points
   \param[in]  dchi values of partial derivatives of chi at quadrature points
   \param[in]  q    index of the quadrature point
   \param[in]  dx   horizontal grid spacing (x-direction)
   \param[in]  dy   horizontal grid spacing (y-direction)
   \param[in]  dz   partial derivatives of z computed by calling compute_z_gradient()
   \param[out] phi  values of shape functions at the quadrature point q
   \param[out] dphi values of partial derivatives of shape functions at q
   \param[out] jw   \f$det(J)*w\f$
 */
void compute_element_info(PetscReal chi[8][8],PetscReal dchi[8][8][3],
			  PetscInt q, PetscReal dx, PetscReal dy, const PetscReal dz[restrict],
			  PetscReal phi[restrict],
			  PetscReal dphi[restrict][3],
			  PetscReal *restrict jw)
{
  const PetscReal jac[3][3] = {{dx / 2, 0, 0}, {0, dy / 2, 0}, {dz[0], dz[1], dz[2]}},
    ijac[3][3] = {{1 / jac[0][0], 0, 0},
                  {0, 1 / jac[1][1], 0},
                  { - jac[2][0] / (jac[0][0]*jac[2][2]), - jac[2][1] / (jac[1][1]*jac[2][2]), 1 / jac[2][2]}},
      jdet = jac[0][0]*jac[1][1]*jac[2][2];
  PetscInt i;

  for (i = 0; i < 8; i++) {
    const PetscReal *dphir = dchi[q][i];
    phi[i] = chi[q][i];
    dphi[i][0] = dphir[0]*ijac[0][0] + dphir[1]*ijac[1][0] + dphir[2]*ijac[2][0];
    dphi[i][1] = dphir[0]*ijac[0][1] + dphir[1]*ijac[1][1] + dphir[2]*ijac[2][1];
    dphi[i][2] = dphir[0]*ijac[0][2] + dphir[1]*ijac[1][2] + dphir[2]*ijac[2][2];
  }
  *jw = 1.0 * jdet;		/* hard-wired quadrature weights */
}


/*! Compute values of shape functions and their derivatives at quadrature points.
 *
 * This corresponds to 2D Q1 elements and the 2x2 Gaussian quadrature.
 */
void initialize_Q12D(PetscReal chi[4][4], PetscReal dchi[4][4][2])
{
  /* Coordinated of the nodes of the reference element: */
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
 * This corresponds to 3D Q1 elements and the 2x2x2 Gaussian quadrature.
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
