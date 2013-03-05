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

#ifndef _BLATTER_IMPLEMENTATION_H_
#define _BLATTER_IMPLEMENTATION_H_

/*! This file contains declarations that need to be visible to the BlatterStressBalance C++ class.
 *
 */

#ifdef __cplusplus
extern "C" {
#endif

/*! Storage for components of the horizontal ice velocity (at each node of the
   3D grid).*/
typedef struct {
  PetscScalar u, v;
} Node;

/*! Storage for model parameters (at each node of the 2D grid). */
typedef struct {
  PetscScalar ice_bottom; /*!< elevation of the bottom surface of the ice */
  PetscScalar thickness;  /*!< thickness */
  PetscScalar tauc;	  /*!< till yield stress */
} PrmNode;

typedef struct {
  PetscReal Lx;			/*!< full width of the model domain in the x-direction */
  PetscReal Ly;			/*!< full width of the model domain in the y-direction */
  PetscReal dirichlet_scale;	/*!< scaling for Dirichlet boundary conditions */
  PetscReal rhog;		/*!< ice density times standard gravity */
  PetscBool no_slip;		/*!< whether to use the no slip condition at the base */
  struct {
    /*! Function evaluating effective viscosity as a function of ice
        hardness and the second invariant \f$ \gamma \f$

	\param[in] ctx pointer to this context (BlatterQ1Ctx struct)
	\param[in] hardness ice hardness
	\param[in] gamma second invariant
	\param[out] eta effective viscosity \f$ \eta \f$.
	\param[out] derivative of the effective viscosity with respect
	  to the second invariant, \f$ \frac{\partial \eta}{\partial \gamma} \f$
    */
    void (*viscosity)(void* ctx, PetscReal hardness, PetscReal gamma,
		      PetscReal *eta, PetscReal *deta);

    /*! Function evaluating the basal drag coefficient \f$ \tau_b \f$ as a
        function of the basal yield stress \f$ \tau_c \f$ and
	\f$ \gamma_b = u_b^2 + v_b^2 \f$.

	\param[in] ctx pointer to this context (BlatterQ1Ctx struct)
	\param[in] tauc basal yield stress \f$ \tau_c \f$
	\param[in] gamma_b \f$ \gamma_b = u_b^2 + v_b^2 \f$
	\param[out] taub \f$ \tau_b \f$
	\param[out] dtaub derivative of \f$ \tau_b \f$ with respect to \f$ \gamma_b \f$
     */
    void (*drag)(void *ctx, PetscReal tauc, PetscReal u, PetscReal v,
		 PetscReal *taub, PetscReal *dtaub);
  } nonlinear;

  /*! 3D Q1 elements with the 2x2x2 (8-point) Gaussian quadrature */
  struct {
    /*! Values of shape functions at quadrature points.
     *
     * chi[q][i] corresponds to quadrature point q, shape function i.
     */
    PetscReal chi[8][8];
    /*! Derivatives of shape functions at quadrature points.
     *
     * dchi[q][i][d] corresponds to quadrature point q, shape function i, direction d.
     */
    PetscReal dchi[8][8][3];
  } Q13D;

  /*! 2D Q1 elements with the 2x2 (4-point) Gaussian quadrature */
  struct {
    /*! Values of shape functions at quadrature points.
     *
     * chi[q][i] corresponds to quadrature point q, shape function i.
     */
    PetscReal chi[4][4];
    /*! Derivatives of shape functions at quadrature points.
     *
     * dchi[q][i][d] corresponds to quadrature point q, shape function i, direction d.
     */
    PetscReal dchi[4][4][2];
  } Q12D;

  /*! Pointer to a PISM-side class that might contain parameters used
      by viscosity() and drag().  */
  void *extra;
} BlatterQ1Ctx;

PetscErrorCode BlatterQ1_begin_2D_parameter_access(DM da, PrmNode ***prm);
PetscErrorCode BlatterQ1_end_2D_parameter_access(DM da, PrmNode ***prm);

PetscErrorCode BlatterQ1_begin_hardness_access(DM da, PetscScalar ****hardness);
PetscErrorCode BlatterQ1_end_hardness_access(DM da, PetscScalar ****hardness);

PetscErrorCode BlatterQ1_create(MPI_Comm com, DM pism_da, PetscInt Mz,
				BlatterQ1Ctx *ctx, SNES *result);
#ifdef __cplusplus
}
#endif


#endif /* _BLATTER_IMPLEMENTATION_H_ */
