// Copyright (C) 2012, 2014  David Maxwell
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

#include "IP_L2NormFunctional.hh"

namespace pism {

PetscErrorCode IP_L2NormFunctional2S::valueAt(IceModelVec2S &x, double *OUTPUT) {

  PetscErrorCode   ierr;

  // The value of the objective
  double value = 0;

  double x_e[FEQuadrature::Nk];
  double x_q[FEQuadrature::Nq];

  IceModelVec::AccessList list(x);

  // Jacobian times weights for quadrature.
  const double* JxW = m_quadrature.getWeightedJacobian();

  // Loop through all LOCAL elements.
  int xs = m_element_index.lxs, xm = m_element_index.lxm,
    ys = m_element_index.lys, ym = m_element_index.lym;
  for (int i=xs; i<xs+xm; i++) {
    for (int j=ys; j<ys+ym; j++) {

      // Obtain values of x at the quadrature points for the element.
      m_dofmap.extractLocalDOFs(i, j, x, x_e);
      m_quadrature.computeTrialFunctionValues(x_e, x_q);

      for (unsigned int q=0; q<FEQuadrature::Nq; q++) {
        const double x_qq = x_q[q];
        value += JxW[q]*x_qq*x_qq;
      } // q
    } // j
  } // i

  ierr = GlobalSum(m_grid.com, &value,  OUTPUT); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IP_L2NormFunctional2S::dot(IceModelVec2S &a, IceModelVec2S &b, double *OUTPUT) {

  PetscErrorCode   ierr;

  // The value of the objective
  double value = 0;

  double a_q[FEQuadrature::Nq];

  double b_q[FEQuadrature::Nq];

  IceModelVec::AccessList list(a);
  list.add(b);

  // Jacobian times weights for quadrature.
  const double* JxW = m_quadrature.getWeightedJacobian();

  // Loop through all LOCAL elements.
  int xs = m_element_index.lxs, xm = m_element_index.lxm,
           ys = m_element_index.lys, ym = m_element_index.lym;
  for (int i=xs; i<xs+xm; i++) {
    for (int j=ys; j<ys+ym; j++) {

      // Obtain values of x at the quadrature points for the element.
      m_quadrature.computeTrialFunctionValues(i, j, m_dofmap, a, a_q);
      m_quadrature.computeTrialFunctionValues(i, j, m_dofmap, b, b_q);

      for (unsigned int q=0; q<FEQuadrature::Nq; q++) {
        value += JxW[q]*a_q[q]*b_q[q];
      } // q
    } // j
  } // i

  ierr = GlobalSum(m_grid.com, &value,  OUTPUT); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IP_L2NormFunctional2S::gradientAt(IceModelVec2S &x, IceModelVec2S &gradient) {

  // Clear the gradient before doing anything with it!
  gradient.set(0);

  double x_q[FEQuadrature::Nq];
  double gradient_e[FEQuadrature::Nk];

  IceModelVec::AccessList list(x);
  list.add(gradient);

  // An Nq by Nk array of test function values.
  const FEFunctionGerm (*test)[FEQuadrature::Nk] = m_quadrature.testFunctionValues();

  // Jacobian times weights for quadrature.
  const double* JxW = m_quadrature.getWeightedJacobian();

  // Loop through all local and ghosted elements.
  int xs = m_element_index.xs, xm = m_element_index.xm,
           ys = m_element_index.ys, ym = m_element_index.ym;
  for (int i=xs; i<xs+xm; i++) {
    for (int j=ys; j<ys+ym; j++) {

      // Reset the DOF map for this element.
      m_dofmap.reset(i, j, m_grid);

      // Obtain values of x at the quadrature points for the element.
      m_quadrature.computeTrialFunctionValues(i, j, m_dofmap, x, x_q);

      // Zero out the element-local residual in prep for updating it.
      for (unsigned int k=0; k<FEQuadrature::Nk; k++) {
        gradient_e[k] = 0;
      }

      for (unsigned int q=0; q<FEQuadrature::Nq; q++) {
        const double x_qq = x_q[q];
        for (unsigned int k=0; k<FEQuadrature::Nk; k++) {
          gradient_e[k] += 2*JxW[q]*x_qq*test[q][k].val;
        } // k
      } // q
      m_dofmap.addLocalResidualBlock(gradient_e, gradient);
    } // j
  } // i

  return 0;
}

PetscErrorCode IP_L2NormFunctional2V::valueAt(IceModelVec2V &x, double *OUTPUT) {

  // The value of the objective
  double value = 0;

  Vector2 x_e[FEQuadrature::Nk];
  Vector2 x_q[FEQuadrature::Nq];

  IceModelVec::AccessList list(x);

  // Jacobian times weights for quadrature.
  const double* JxW = m_quadrature.getWeightedJacobian();

  // Loop through all local and ghosted elements.
  int xs = m_element_index.lxs, xm = m_element_index.lxm,
           ys = m_element_index.lys, ym = m_element_index.lym;
  for (int i=xs; i<xs+xm; i++) {
    for (int j=ys; j<ys+ym; j++) {

      // Obtain values of x at the quadrature points for the element.
      m_dofmap.extractLocalDOFs(i, j, x, x_e);
      m_quadrature_vector.computeTrialFunctionValues(x_e, x_q);

      for (unsigned int q=0; q<FEQuadrature::Nq; q++) {
        const Vector2 &x_qq = x_q[q];
        value += JxW[q]*(x_qq.u*x_qq.u+x_qq.v*x_qq.v);
      } // q
    } // j
  } // i

  GlobalSum(m_grid.com, &value,  OUTPUT);

  return 0;
}

PetscErrorCode IP_L2NormFunctional2V::dot(IceModelVec2V &a, IceModelVec2V &b, double *OUTPUT) {

  // The value of the objective
  double value = 0;

  Vector2 a_q[FEQuadrature::Nq];

  Vector2 b_q[FEQuadrature::Nq];

  IceModelVec::AccessList list(a);
  list.add(b);

  // Jacobian times weights for quadrature.
  const double* JxW = m_quadrature.getWeightedJacobian();

  // Loop through all LOCAL elements.
  int xs = m_element_index.lxs, xm = m_element_index.lxm,
           ys = m_element_index.lys, ym = m_element_index.lym;
  for (int i=xs; i<xs+xm; i++) {
    for (int j=ys; j<ys+ym; j++) {

      // Obtain values of x at the quadrature points for the element.
      m_quadrature_vector.computeTrialFunctionValues(i, j, m_dofmap, a, a_q);
      m_quadrature_vector.computeTrialFunctionValues(i, j, m_dofmap, b, b_q);

      for (unsigned int q=0; q<FEQuadrature::Nq; q++) {
        value += JxW[q]*(a_q[q].u*b_q[q].u+a_q[q].v*b_q[q].v);
      } // q
    } // j
  } // i

  GlobalSum(m_grid.com, &value,  OUTPUT);

  return 0;
}

PetscErrorCode IP_L2NormFunctional2V::gradientAt(IceModelVec2V &x, IceModelVec2V &gradient) {

  // Clear the gradient before doing anything with it!
  gradient.set(0);

  Vector2 x_q[FEQuadrature::Nq];
  Vector2 gradient_e[FEQuadrature::Nk];

  IceModelVec::AccessList list(x);
  list.add(gradient);

  // An Nq by Nk array of test function values.
  const FEFunctionGerm (*test)[FEQuadrature::Nk] = m_quadrature.testFunctionValues();

  // Jacobian times weights for quadrature.
  const double* JxW = m_quadrature.getWeightedJacobian();

  // Loop through all local and ghosted elements.
  int xs = m_element_index.xs, xm = m_element_index.xm,
           ys = m_element_index.ys, ym = m_element_index.ym;
  for (int i=xs; i<xs+xm; i++) {
    for (int j=ys; j<ys+ym; j++) {

      // Reset the DOF map for this element.
      m_dofmap.reset(i, j, m_grid);

      // Obtain values of x at the quadrature points for the element.
      m_quadrature_vector.computeTrialFunctionValues(i, j, m_dofmap, x, x_q);

      // Zero out the element-local residual in prep for updating it.
      for (unsigned int k=0; k<FEQuadrature::Nk; k++) {
        gradient_e[k].u = 0;
        gradient_e[k].v = 0;
      }

      for (unsigned int q=0; q<FEQuadrature::Nq; q++) {
        const Vector2 &x_qq = x_q[q];
        for (unsigned int k=0; k<FEQuadrature::Nk; k++) {
          double gcommon =2*JxW[q]*test[q][k].val;
          gradient_e[k].u += gcommon*x_qq.u;
          gradient_e[k].v += gcommon*x_qq.v;
        } // k
      } // q
      m_dofmap.addLocalResidualBlock(gradient_e, gradient);
    } // j
  } // i

  return 0;
}

} // end of namespace pism
