// Copyright (C) 2013, 2014, 2015  David Maxwell and Constantine Khroulev
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

#include "IPGroundedIceH1NormFunctional.hh"
#include "base/util/error_handling.hh"
#include "base/util/IceGrid.hh"
#include "base/util/pism_const.hh"

namespace pism {
namespace inverse {

void IPGroundedIceH1NormFunctional2S::valueAt(IceModelVec2S &x, double *OUTPUT) {

  using fem::Quadrature;

  // The value of the objective
  double value = 0;

  double x_e[Quadrature::Nk];
  double x_q[Quadrature::Nq], dxdx_q[Quadrature::Nq], dxdy_q[Quadrature::Nq];

  IceModelVec::AccessList list;
  list.add(x);

  // Jacobian times weights for quadrature.
  const double* JxW = m_quadrature.getWeightedJacobian();

  fem::DirichletData_Scalar dirichletBC;
  dirichletBC.init(m_dirichletIndices, NULL);

  list.add(m_ice_mask);
  MaskQuery iceQuery(m_ice_mask);

  // Loop through all LOCAL elements.
  int xs = m_element_index.lxs, xm = m_element_index.lxm,
           ys = m_element_index.lys, ym = m_element_index.lym;
  for (int j=ys; j<ys+ym; j++) {
    for (int i=xs; i<xs+xm; i++) {
      bool all_grounded_ice = iceQuery.grounded_ice(i, j) & iceQuery.grounded_ice(i+1, j) &
        iceQuery.grounded_ice(i, j+1) & iceQuery.grounded_ice(i+1, j+1);
      if (! all_grounded_ice) {
        continue;
      }

      m_dofmap.reset(i, j, *m_grid);

      // Obtain values of x at the quadrature points for the element.
      m_dofmap.extractLocalDOFs(x, x_e);
      if (dirichletBC) {
        dirichletBC.update_homogeneous(m_dofmap, x_e);
      }
      m_quadrature.computeTrialFunctionValues(x_e, x_q, dxdx_q, dxdy_q);

      for (unsigned int q=0; q<Quadrature::Nq; q++) {
        value += JxW[q]*(m_cL2*x_q[q]*x_q[q]+ m_cH1*(dxdx_q[q]*dxdx_q[q]+dxdy_q[q]*dxdy_q[q]));
      } // q
    } // j
  } // i

  GlobalSum(m_grid->com, &value, OUTPUT, 1);

  dirichletBC.finish();
}

void IPGroundedIceH1NormFunctional2S::dot(IceModelVec2S &a, IceModelVec2S &b, double *OUTPUT) {

  using fem::Quadrature;

  // The value of the objective
  double value = 0;

  double a_e[Quadrature::Nk];
  double a_q[Quadrature::Nq], dadx_q[Quadrature::Nq], dady_q[Quadrature::Nq];

  IceModelVec::AccessList list;
  list.add(a);

  double b_e[Quadrature::Nk];
  double b_q[Quadrature::Nq], dbdx_q[Quadrature::Nq], dbdy_q[Quadrature::Nq];
  list.add(b);

  // Jacobian times weights for quadrature.
  const double* JxW = m_quadrature.getWeightedJacobian();

  fem::DirichletData_Scalar dirichletBC;
  dirichletBC.init(m_dirichletIndices, NULL);

  list.add(m_ice_mask);
  MaskQuery iceQuery(m_ice_mask);

  // Loop through all LOCAL elements.
  int xs = m_element_index.lxs, xm = m_element_index.lxm,
           ys = m_element_index.lys, ym = m_element_index.lym;
  for (int j=ys; j<ys+ym; j++) {
    for (int i=xs; i<xs+xm; i++) {
      bool all_grounded_ice = iceQuery.grounded_ice(i, j) & iceQuery.grounded_ice(i+1, j) &
        iceQuery.grounded_ice(i, j+1) & iceQuery.grounded_ice(i+1, j+1);
      if (! all_grounded_ice) {
        continue;
      }

      m_dofmap.reset(i, j, *m_grid);

      // Obtain values of x at the quadrature points for the element.
      m_dofmap.extractLocalDOFs(a, a_e);
      if (dirichletBC) {
        dirichletBC.update_homogeneous(m_dofmap, a_e);
      }
      m_quadrature.computeTrialFunctionValues(a_e, a_q, dadx_q, dady_q);

      m_dofmap.extractLocalDOFs(b, b_e);
      if (dirichletBC) {
        dirichletBC.update_homogeneous(m_dofmap, b_e);
      }
      m_quadrature.computeTrialFunctionValues(b_e, b_q, dbdx_q, dbdy_q);

      for (unsigned int q=0; q<Quadrature::Nq; q++) {
        value += JxW[q]*(m_cL2*a_q[q]*b_q[q]+ m_cH1*(dadx_q[q]*dbdx_q[q]+dady_q[q]*dbdy_q[q]));
      } // q
    } // j
  } // i

  GlobalSum(m_grid->com, &value, OUTPUT, 1);

  dirichletBC.finish();
}


void IPGroundedIceH1NormFunctional2S::gradientAt(IceModelVec2S &x, IceModelVec2S &gradient) {

  using fem::Quadrature;

  // Clear the gradient before doing anything with it!
  gradient.set(0);

  double x_e[Quadrature::Nk];
  double x_q[Quadrature::Nq], dxdx_q[Quadrature::Nq], dxdy_q[Quadrature::Nq];

  IceModelVec::AccessList list;
  list.add(x);

  double gradient_e[Quadrature::Nk];
  list.add(gradient);

  // An Nq by Nk array of test function values.
  const fem::FunctionGerm (*test)[Quadrature::Nk] = m_quadrature.testFunctionValues();

  // Jacobian times weights for quadrature.
  const double* JxW = m_quadrature.getWeightedJacobian();

  fem::DirichletData_Scalar dirichletBC;
  dirichletBC.init(m_dirichletIndices, NULL);

  list.add(m_ice_mask);
  MaskQuery iceQuery(m_ice_mask);

  // Loop through all local and ghosted elements.
  int xs = m_element_index.xs, xm = m_element_index.xm,
           ys = m_element_index.ys, ym = m_element_index.ym;
  for (int j=ys; j<ys+ym; j++) {
    for (int i=xs; i<xs+xm; i++) {
      bool all_grounded_ice = iceQuery.grounded_ice(i, j) & iceQuery.grounded_ice(i+1, j) &
        iceQuery.grounded_ice(i, j+1) & iceQuery.grounded_ice(i+1, j+1);
      if (! all_grounded_ice) {
        continue;
      }

      // Reset the DOF map for this element.
      m_dofmap.reset(i, j, *m_grid);

      // Obtain values of x at the quadrature points for the element.
      m_dofmap.extractLocalDOFs(i, j, x, x_e);
      if (dirichletBC) {
        dirichletBC.constrain(m_dofmap);
        dirichletBC.update_homogeneous(m_dofmap, x_e);
      }
      m_quadrature.computeTrialFunctionValues(x_e, x_q, dxdx_q, dxdy_q);

      // Zero out the element-local residual in prep for updating it.
      for (unsigned int k=0; k<Quadrature::Nk; k++) {
        gradient_e[k] = 0;
      }

      for (unsigned int q=0; q<Quadrature::Nq; q++) {
        const double &x_qq=x_q[q];
        const double &dxdx_qq=dxdx_q[q], &dxdy_qq=dxdy_q[q];
        for (unsigned int k=0; k<Quadrature::Nk; k++) {
          gradient_e[k] += 2*JxW[q]*(m_cL2*x_qq*test[q][k].val +
            m_cH1*(dxdx_qq*test[q][k].dx + dxdy_qq*test[q][k].dy));
        } // k
      } // q
      m_dofmap.addLocalResidualBlock(gradient_e, gradient);
    } // j
  } // i

  dirichletBC.finish();
}

void IPGroundedIceH1NormFunctional2S::assemble_form(Mat form) {

  using fem::Quadrature;

  PetscErrorCode ierr;
  int         i, j;

  // Zero out the Jacobian in preparation for updating it.
  ierr = MatZeroEntries(form);
  PISM_CHK(ierr, "MatZeroEntries");

  // Jacobian times weights for quadrature.
  const double* JxW = m_quadrature.getWeightedJacobian();

  fem::DirichletData_Scalar zeroLocs;
  zeroLocs.init(m_dirichletIndices, NULL);

  IceModelVec::AccessList list;
  list.add(m_ice_mask);
  MaskQuery iceQuery(m_ice_mask);

  // Values of the finite element test functions at the quadrature points.
  // This is an Nq by Nk array of function germs (Nq=#of quad pts, Nk=#of test functions).
  const fem::FunctionGerm (*test)[Quadrature::Nk] = m_quadrature.testFunctionValues();

  // Loop through all the elements.
  int xs = m_element_index.xs, xm = m_element_index.xm,
           ys = m_element_index.ys, ym = m_element_index.ym;
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      bool all_grounded_ice = iceQuery.grounded_ice(i, j) & iceQuery.grounded_ice(i+1, j) &
         iceQuery.grounded_ice(i, j+1) & iceQuery.grounded_ice(i+1, j+1);
      if (! all_grounded_ice) {
        continue;
      }

      // Element-local Jacobian matrix (there are Quadrature::Nk vector valued degrees
      // of freedom per elment, for a total of (2*Quadrature::Nk)*(2*Quadrature::Nk) = 16
      // entries in the local Jacobian.
      double K[Quadrature::Nk][Quadrature::Nk];

      // Initialize the map from global to local degrees of freedom for this element.
      m_dofmap.reset(i, j, *m_grid);

      // Don't update rows/cols where we project to zero.
      if (zeroLocs) {
        zeroLocs.constrain(m_dofmap);
      }

      // Build the element-local Jacobian.
      ierr = PetscMemzero(K, sizeof(K));
      PISM_CHK(ierr, "PetscMemzero");

      for (unsigned int q=0; q<Quadrature::Nq; q++) {
        for (unsigned int k = 0; k < Quadrature::Nk; k++) {   // Test functions
          for (unsigned int l = 0; l < Quadrature::Nk; l++) { // Trial functions
            const fem::FunctionGerm &test_qk=test[q][k];
            const fem::FunctionGerm &test_ql=test[q][l];
            K[k][l]     += JxW[q]*(m_cL2*test_qk.val*test_ql.val
              +  m_cH1*(test_qk.dx*test_ql.dx + test_qk.dy*test_ql.dy));
          } // l
        } // k
      } // q
      m_dofmap.addLocalJacobianBlock(&K[0][0], form);
    } // j
  } // i

  if (zeroLocs) {
    zeroLocs.fix_jacobian(form);
  }
  zeroLocs.finish();


  ierr = MatAssemblyBegin(form, MAT_FINAL_ASSEMBLY);
  PISM_CHK(ierr, "MatAssemblyBegin");

  ierr = MatAssemblyEnd(form, MAT_FINAL_ASSEMBLY);
  PISM_CHK(ierr, "MatAssemblyEnd");
}

} // end of namespace inverse
} // end of namespace pism
