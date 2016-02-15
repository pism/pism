// Copyright (C) 2012, 2014, 2015, 2016  David Maxwell and Constantine Khroulev
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

#include "IP_H1NormFunctional.hh"
#include "base/util/error_handling.hh"
#include "base/util/IceGrid.hh"
#include "base/util/pism_const.hh"
#include "base/util/pism_utilities.hh"

namespace pism {
namespace inverse {

void IP_H1NormFunctional2S::valueAt(IceModelVec2S &x, double *OUTPUT) {

  using fem::Quadrature2x2;
  const unsigned int Nk = fem::ShapeQ1::Nk;
  const unsigned int Nq = Quadrature2x2::Nq;

  // The value of the objective
  double value = 0;

  double x_e[Nk];
  double x_q[Nq], dxdx_q[Nq], dxdy_q[Nq];
  IceModelVec::AccessList list(x);

  // Jacobian times weights for quadrature.
  const double* JxW = m_quadrature.weighted_jacobian();

  fem::DirichletData_Scalar dirichletBC(m_dirichletIndices, NULL);

  // Loop through all LOCAL elements.
  const int
    xs = m_element_index.lxs,
    xm = m_element_index.lxm,
    ys = m_element_index.lys,
    ym = m_element_index.lym;

  for (int j=ys; j<ys+ym; j++) {
    for (int i=xs; i<xs+xm; i++) {
      m_element_map.reset(i, j, *m_grid);

      // Obtain values of x at the quadrature points for the element.
      m_element_map.nodal_values(x, x_e);
      if (dirichletBC) {
        dirichletBC.update_homogeneous(m_element_map, x_e);
      }
      m_quadrature.quadrature_point_values(x_e, x_q, dxdx_q, dxdy_q);

      for (unsigned int q=0; q<Nq; q++) {
        value += JxW[q]*(m_cL2*x_q[q]*x_q[q]+ m_cH1*(dxdx_q[q]*dxdx_q[q]+dxdy_q[q]*dxdy_q[q]));
      } // q
    } // j
  } // i

  GlobalSum(m_grid->com, &value, OUTPUT, 1);

  dirichletBC.finish();
}

void IP_H1NormFunctional2S::dot(IceModelVec2S &a, IceModelVec2S &b, double *OUTPUT) {

  using fem::Quadrature2x2;
  const unsigned int Nk = fem::ShapeQ1::Nk;
  const unsigned int Nq = Quadrature2x2::Nq;

  // The value of the objective
  double value = 0;

  double a_e[Nk];
  double a_q[Nq], dadx_q[Nq], dady_q[Nq];

  double b_e[Nk];
  double b_q[Nq], dbdx_q[Nq], dbdy_q[Nq];

  IceModelVec::AccessList list(a);
  list.add(b);

  // Jacobian times weights for quadrature.
  const double* JxW = m_quadrature.weighted_jacobian();

  fem::DirichletData_Scalar dirichletBC(m_dirichletIndices, NULL);

  // Loop through all LOCAL elements.
  const int
    xs = m_element_index.lxs,
    xm = m_element_index.lxm,
    ys = m_element_index.lys,
    ym = m_element_index.lym;

  for (int j=ys; j<ys+ym; j++) {
    for (int i=xs; i<xs+xm; i++) {
      m_element_map.reset(i, j, *m_grid);

      // Obtain values of x at the quadrature points for the element.
      m_element_map.nodal_values(a, a_e);
      if (dirichletBC) {
        dirichletBC.update_homogeneous(m_element_map, a_e);
      }
      m_quadrature.quadrature_point_values(a_e, a_q, dadx_q, dady_q);

      m_element_map.nodal_values(b, b_e);
      if (dirichletBC) {
        dirichletBC.update_homogeneous(m_element_map, b_e);
      }
      m_quadrature.quadrature_point_values(b_e, b_q, dbdx_q, dbdy_q);

      for (unsigned int q=0; q<Nq; q++) {
        value += JxW[q]*(m_cL2*a_q[q]*b_q[q]+ m_cH1*(dadx_q[q]*dbdx_q[q]+dady_q[q]*dbdy_q[q]));
      } // q
    } // j
  } // i

  GlobalSum(m_grid->com, &value, OUTPUT, 1);

  dirichletBC.finish();
}


void IP_H1NormFunctional2S::gradientAt(IceModelVec2S &x, IceModelVec2S &gradient) {

  using fem::Quadrature2x2;
  const unsigned int Nk = fem::ShapeQ1::Nk;
  const unsigned int Nq = Quadrature2x2::Nq;

  // Clear the gradient before doing anything with it!
  gradient.set(0);

  double x_e[Nk];
  double x_q[Nq], dxdx_q[Nq], dxdy_q[Nq];

  double gradient_e[Nk];

  IceModelVec::AccessList list(x);
  list.add(gradient);

  // An Nq by Nk array of test function values.
  const fem::Germ<double> (*test)[Nk] = m_quadrature.testFunctionValues();

  // Jacobian times weights for quadrature.
  const double* JxW = m_quadrature.weighted_jacobian();

  fem::DirichletData_Scalar dirichletBC(m_dirichletIndices, NULL);

  // Loop through all local and ghosted elements.
  const int
    xs = m_element_index.xs,
    xm = m_element_index.xm,
    ys = m_element_index.ys,
    ym = m_element_index.ym;

  for (int j=ys; j<ys+ym; j++) {
    for (int i=xs; i<xs+xm; i++) {

      // Reset the DOF map for this element.
      m_element_map.reset(i, j, *m_grid);

      // Obtain values of x at the quadrature points for the element.
      m_element_map.nodal_values(x, x_e);
      if (dirichletBC) {
        dirichletBC.constrain(m_element_map);
        dirichletBC.update_homogeneous(m_element_map, x_e);
      }
      m_quadrature.quadrature_point_values(x_e, x_q, dxdx_q, dxdy_q);

      // Zero out the element-local residual in prep for updating it.
      for (unsigned int k=0; k<Nk; k++) {
        gradient_e[k] = 0;
      }

      for (unsigned int q=0; q<Nq; q++) {
        const double &x_qq=x_q[q];
        const double &dxdx_qq=dxdx_q[q], &dxdy_qq=dxdy_q[q];
        for (unsigned int k=0; k<Nk; k++) {
          gradient_e[k] += 2*JxW[q]*(m_cL2*x_qq*test[q][k].val +
            m_cH1*(dxdx_qq*test[q][k].dx + dxdy_qq*test[q][k].dy));
        } // k
      } // q
      m_element_map.add_residual_contribution(gradient_e, gradient);
    } // j
  } // i

  dirichletBC.finish();
}

void IP_H1NormFunctional2S::assemble_form(Mat form) {

  using fem::Quadrature2x2;
  const unsigned int Nk = fem::ShapeQ1::Nk;
  const unsigned int Nq = Quadrature2x2::Nq;

  PetscErrorCode ierr;

  // Zero out the Jacobian in preparation for updating it.
  ierr = MatZeroEntries(form);
  PISM_CHK(ierr, "MatZeroEntries");

  // Jacobian times weights for quadrature.
  const double* JxW = m_quadrature.weighted_jacobian();

  fem::DirichletData_Scalar zeroLocs(m_dirichletIndices, NULL);

  // Values of the finite element test functions at the quadrature points.
  // This is an Nq by Nk array of function germs (Nq=#of quad pts, Nk=#of test functions).
  const fem::Germ<double> (*test)[Nk] = m_quadrature.testFunctionValues();

  // Loop through all the elements.
  const int
    xs = m_element_index.xs,
    xm = m_element_index.xm,
    ys = m_element_index.ys,
    ym = m_element_index.ym;

  ParallelSection loop(m_grid->com);
  try {
    for (int j=ys; j<ys+ym; j++) {
      for (int i=xs; i<xs+xm; i++) {
        // Element-local Jacobian matrix (there are Nk vector valued degrees
        // of freedom per elment, for a total of (2*Nk)*(2*Nk) = 16
        // entries in the local Jacobian.
        double K[Nk][Nk];

        // Initialize the map from global to local degrees of freedom for this element.
        m_element_map.reset(i, j, *m_grid);

        // Don't update rows/cols where we project to zero.
        if (zeroLocs) {
          zeroLocs.constrain(m_element_map);
        }

        // Build the element-local Jacobian.
        ierr = PetscMemzero(K, sizeof(K));
        PISM_CHK(ierr, "PetscMemzero");

        for (unsigned int q=0; q<Nq; q++) {
          for (unsigned int k = 0; k < Nk; k++) {   // Test functions
            for (unsigned int l = 0; l < Nk; l++) { // Trial functions
              const fem::Germ<double> &test_qk=test[q][k];
              const fem::Germ<double> &test_ql=test[q][l];
              K[k][l] += JxW[q]*(m_cL2*test_qk.val*test_ql.val +
                                 m_cH1*(test_qk.dx*test_ql.dx +
                                        test_qk.dy*test_ql.dy));
            } // l
          } // k
        } // q
        m_element_map.add_jacobian_contribution(&K[0][0], form);
      } // j
    } // i
  } catch (...) {
    loop.failed();
  }
  loop.check();

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
