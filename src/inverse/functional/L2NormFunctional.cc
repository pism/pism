// Copyright (C) 2012  David Maxwell
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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

#include "L2NormFunctional.hh"

PetscErrorCode L2NormFunctional2S::valueAt(IceModelVec2S &x, PetscReal *OUTPUT) {

  PetscErrorCode   ierr;

  // The value of the objective
  PetscReal value = 0;

  PetscReal **x_a;
  PetscReal x_e[FEQuadrature::Nk];
  PetscReal x_q[FEQuadrature::Nq];
  ierr = x.get_array(x_a); CHKERRQ(ierr);

  // Jacobian times weights for quadrature.
  PetscScalar JxW[FEQuadrature::Nq];
  m_quadrature.getWeightedJacobian(JxW);

  // Loop through all LOCAL elements.
  PetscInt xs = m_element_index.lxs, xm = m_element_index.lxm,
           ys = m_element_index.lys, ym = m_element_index.lym;
  for (PetscInt i=xs; i<xs+xm; i++) {
    for (PetscInt j=ys; j<ys+ym; j++) {

      // Obtain values of x at the quadrature points for the element.
      m_dofmap.extractLocalDOFs(i,j,x_a,x_e);
      m_quadrature.computeTrialFunctionValues(x_e,x_q);

      for (PetscInt q=0; q<FEQuadrature::Nq; q++) {
        const PetscReal x_qq = x_q[q];
        value += JxW[q]*x_qq*x_qq;
      } // q
    } // j
  } // i

  ierr = PISMGlobalSum(&value, OUTPUT, m_grid.com); CHKERRQ(ierr);

  ierr = x.end_access(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode L2NormFunctional2S::dot(IceModelVec2S &a, IceModelVec2S &b, PetscReal *OUTPUT) {

  PetscErrorCode   ierr;

  // The value of the objective
  PetscReal value = 0;

  PetscReal **a_a;
  PetscReal a_q[FEQuadrature::Nq];

  PetscReal **b_a;
  PetscReal b_q[FEQuadrature::Nq];

  ierr = a.get_array(a_a); CHKERRQ(ierr);
  ierr = b.get_array(b_a); CHKERRQ(ierr);

  // Jacobian times weights for quadrature.
  PetscScalar JxW[FEQuadrature::Nq];
  m_quadrature.getWeightedJacobian(JxW);

  // Loop through all LOCAL elements.
  PetscInt xs = m_element_index.lxs, xm = m_element_index.lxm,
           ys = m_element_index.lys, ym = m_element_index.lym;
  for (PetscInt i=xs; i<xs+xm; i++) {
    for (PetscInt j=ys; j<ys+ym; j++) {

      // Obtain values of x at the quadrature points for the element.
      m_quadrature.computeTrialFunctionValues(i,j,m_dofmap,a_a,a_q);
      m_quadrature.computeTrialFunctionValues(i,j,m_dofmap,b_a,b_q);

      for (PetscInt q=0; q<FEQuadrature::Nq; q++) {
        value += JxW[q]*a_q[q]*b_q[q];
      } // q
    } // j
  } // i

  ierr = PISMGlobalSum(&value, OUTPUT, m_grid.com); CHKERRQ(ierr);

  ierr = a.end_access(); CHKERRQ(ierr);
  ierr = b.end_access(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode L2NormFunctional2S::gradientAt(IceModelVec2S &x, IceModelVec2S &gradient) {

  PetscErrorCode   ierr;

  // Clear the gradient before doing anything with it!
  ierr = gradient.set(0); CHKERRQ(ierr);

  PetscReal **x_a;
  PetscReal x_q[FEQuadrature::Nq];
  ierr = x.get_array(x_a); CHKERRQ(ierr);

  PetscReal **gradient_a;
  PetscReal gradient_e[FEQuadrature::Nk];
  ierr = gradient.get_array(gradient_a); CHKERRQ(ierr);

  // An Nq by Nk array of test function values.
  const FEFunctionGerm (*test)[FEQuadrature::Nk] = m_quadrature.testFunctionValues();

  // Jacobian times weights for quadrature.
  PetscScalar JxW[FEQuadrature::Nq];
  m_quadrature.getWeightedJacobian(JxW);

  // Loop through all local and ghosted elements.
  PetscInt xs = m_element_index.xs, xm = m_element_index.xm,
           ys = m_element_index.ys, ym = m_element_index.ym;
  for (PetscInt i=xs; i<xs+xm; i++) {
    for (PetscInt j=ys; j<ys+ym; j++) {

      // Reset the DOF map for this element.
      m_dofmap.reset(i,j,m_grid);

      // Obtain values of x at the quadrature points for the element.
      m_quadrature.computeTrialFunctionValues(i,j,m_dofmap,x_a,x_q);

      // Zero out the element-local residual in prep for updating it.
      for(PetscInt k=0;k<FEQuadrature::Nk;k++){
        gradient_e[k] = 0;
      }

      for (PetscInt q=0; q<FEQuadrature::Nq; q++) {
        const PetscReal x_qq = x_q[q];
        for(PetscInt k=0; k<FEQuadrature::Nk; k++ ) {
          gradient_e[k] += 2*JxW[q]*x_qq*test[q][k].val;
        } // k
      } // q
      m_dofmap.addLocalResidualBlock(gradient_e,gradient_a);
    } // j
  } // i

  ierr = x.end_access(); CHKERRQ(ierr);
  ierr = gradient.end_access(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode L2NormFunctional2V::valueAt(IceModelVec2V &x, PetscReal *OUTPUT) {

  PetscErrorCode   ierr;

  // The value of the objective
  PetscReal value = 0;

  PISMVector2 **x_a;
  PISMVector2 x_e[FEQuadrature::Nk];
  PISMVector2 x_q[FEQuadrature::Nq];

  ierr = x.get_array(x_a); CHKERRQ(ierr);

  // Jacobian times weights for quadrature.
  PetscScalar JxW[FEQuadrature::Nq];
  m_quadrature.getWeightedJacobian(JxW);

  // Loop through all local and ghosted elements.
  PetscInt xs = m_element_index.lxs, xm = m_element_index.lxm,
           ys = m_element_index.lys, ym = m_element_index.lym;
  for (PetscInt i=xs; i<xs+xm; i++) {
    for (PetscInt j=ys; j<ys+ym; j++) {

      // Obtain values of x at the quadrature points for the element.
      m_dofmap.extractLocalDOFs(i,j,x_a,x_e);
      m_quadrature.computeTrialFunctionValues(x_e,x_q);

      for (PetscInt q=0; q<FEQuadrature::Nq; q++) {
        const PISMVector2 &x_qq = x_q[q];
        value += JxW[q]*(x_qq.u*x_qq.u+x_qq.v*x_qq.v);
      } // q
    } // j
  } // i

  ierr = PISMGlobalSum(&value, OUTPUT, m_grid.com); CHKERRQ(ierr);

  ierr = x.end_access(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode L2NormFunctional2V::dot(IceModelVec2V &a, IceModelVec2V &b, PetscReal *OUTPUT) {

  PetscErrorCode   ierr;

  // The value of the objective
  PetscReal value = 0;

  PISMVector2 **a_a;
  PISMVector2 a_q[FEQuadrature::Nq];

  PISMVector2 **b_a;
  PISMVector2 b_q[FEQuadrature::Nq];

  ierr = a.get_array(a_a); CHKERRQ(ierr);
  ierr = b.get_array(b_a); CHKERRQ(ierr);

  // Jacobian times weights for quadrature.
  PetscScalar JxW[FEQuadrature::Nq];
  m_quadrature.getWeightedJacobian(JxW);

  // Loop through all LOCAL elements.
  PetscInt xs = m_element_index.lxs, xm = m_element_index.lxm,
           ys = m_element_index.lys, ym = m_element_index.lym;
  for (PetscInt i=xs; i<xs+xm; i++) {
    for (PetscInt j=ys; j<ys+ym; j++) {

      // Obtain values of x at the quadrature points for the element.
      m_quadrature.computeTrialFunctionValues(i,j,m_dofmap,a_a,a_q);
      m_quadrature.computeTrialFunctionValues(i,j,m_dofmap,b_a,b_q);

      for (PetscInt q=0; q<FEQuadrature::Nq; q++) {
        value += JxW[q]*(a_q[q].u*b_q[q].u+a_q[q].v*b_q[q].v);
      } // q
    } // j
  } // i

  ierr = PISMGlobalSum(&value, OUTPUT, m_grid.com); CHKERRQ(ierr);

  ierr = a.end_access(); CHKERRQ(ierr);
  ierr = b.end_access(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode L2NormFunctional2V::gradientAt(IceModelVec2V &x, IceModelVec2V &gradient) {

  PetscErrorCode   ierr;

  // Clear the gradient before doing anything with it!
  ierr = gradient.set(0); CHKERRQ(ierr);

  PISMVector2 **x_a;
  PISMVector2 x_q[FEQuadrature::Nq];
  ierr = x.get_array(x_a); CHKERRQ(ierr);

  PISMVector2 **gradient_a;
  PISMVector2 gradient_e[FEQuadrature::Nk];
  ierr = gradient.get_array(gradient_a); CHKERRQ(ierr);

  // An Nq by Nk array of test function values.
  const FEFunctionGerm (*test)[FEQuadrature::Nk] = m_quadrature.testFunctionValues();

  // Jacobian times weights for quadrature.
  PetscScalar JxW[FEQuadrature::Nq];
  m_quadrature.getWeightedJacobian(JxW);

  // Loop through all local and ghosted elements.
  PetscInt xs = m_element_index.xs, xm = m_element_index.xm,
           ys = m_element_index.ys, ym = m_element_index.ym;
  for (PetscInt i=xs; i<xs+xm; i++) {
    for (PetscInt j=ys; j<ys+ym; j++) {

      // Reset the DOF map for this element.
      m_dofmap.reset(i,j,m_grid);

      // Obtain values of x at the quadrature points for the element.
      m_quadrature.computeTrialFunctionValues(i,j,m_dofmap,x_a,x_q);

      // Zero out the element-local residual in prep for updating it.
      for(PetscInt k=0;k<FEQuadrature::Nk;k++){
        gradient_e[k].u = 0;
        gradient_e[k].v = 0;
      }

      for (PetscInt q=0; q<FEQuadrature::Nq; q++) {
        const PISMVector2 &x_qq = x_q[q];
        for(PetscInt k=0; k<FEQuadrature::Nk; k++ ) {
          PetscReal gcommon =2*JxW[q]*test[q][k].val;
          gradient_e[k].u += gcommon*x_qq.u;
          gradient_e[k].v += gcommon*x_qq.v;
        } // k
      } // q
      m_dofmap.addLocalResidualBlock(gradient_e,gradient_a);
    } // j
  } // i

  ierr = x.end_access(); CHKERRQ(ierr);
  ierr = gradient.end_access(); CHKERRQ(ierr);
  return 0;
}
