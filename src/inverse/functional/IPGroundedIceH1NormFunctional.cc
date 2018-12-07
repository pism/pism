// Copyright (C) 2013, 2014, 2015, 2016, 2017  David Maxwell and Constantine Khroulev
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
#include "pism/util/error_handling.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/IceModelVec2CellType.hh"

namespace pism {
namespace inverse {

void IPGroundedIceH1NormFunctional2S::valueAt(IceModelVec2S &x, double *OUTPUT) {

  const unsigned int Nk     = fem::q1::n_chi;
  const unsigned int Nq     = m_quadrature.n();
  const unsigned int Nq_max = fem::MAX_QUADRATURE_SIZE;

  // The value of the objective
  double value = 0;

  double x_e[Nk];
  double x_q[Nq_max], dxdx_q[Nq_max], dxdy_q[Nq_max];

  IceModelVec::AccessList list{&x, &m_ice_mask};

  // Jacobian times weights for quadrature.
  const double* W = m_quadrature.weights();

  fem::DirichletData_Scalar dirichletBC(m_dirichletIndices, NULL);

  // Loop through all LOCAL elements.
  const int
    xs = m_element_index.lxs,
    xm = m_element_index.lxm,
    ys = m_element_index.lys,
    ym = m_element_index.lym;

  for (int j=ys; j<ys+ym; j++) {
    for (int i=xs; i<xs+xm; i++) {
      bool all_grounded_ice = (m_ice_mask.grounded_ice(i, j) and
                               m_ice_mask.grounded_ice(i+1, j) and
                               m_ice_mask.grounded_ice(i, j+1) and
                               m_ice_mask.grounded_ice(i+1, j+1));

      if (! all_grounded_ice) {
        continue;
      }

      m_element.reset(i, j);

      // Obtain values of x at the quadrature points for the element.
      m_element.nodal_values(x, x_e);
      if (dirichletBC) {
        dirichletBC.enforce_homogeneous(m_element, x_e);
      }
      quadrature_point_values(m_quadrature, x_e, x_q, dxdx_q, dxdy_q);

      for (unsigned int q=0; q<Nq; q++) {
        value += W[q]*(m_cL2*x_q[q]*x_q[q]+ m_cH1*(dxdx_q[q]*dxdx_q[q]+dxdy_q[q]*dxdy_q[q]));
      } // q
    } // j
  } // i

  GlobalSum(m_grid->com, &value, OUTPUT, 1);
}

void IPGroundedIceH1NormFunctional2S::dot(IceModelVec2S &a, IceModelVec2S &b, double *OUTPUT) {

  const unsigned int Nk     = fem::q1::n_chi;
  const unsigned int Nq     = m_quadrature.n();
  const unsigned int Nq_max = fem::MAX_QUADRATURE_SIZE;

  // The value of the objective
  double value = 0;

  double a_e[Nk];
  double a_q[Nq_max], dadx_q[Nq_max], dady_q[Nq_max];

  IceModelVec::AccessList list{&a, &b, &m_ice_mask};

  double b_e[Nk];
  double b_q[Nq_max], dbdx_q[Nq_max], dbdy_q[Nq_max];

  // Jacobian times weights for quadrature.
  const double* W = m_quadrature.weights();

  fem::DirichletData_Scalar dirichletBC(m_dirichletIndices, NULL);

  // Loop through all LOCAL elements.
  const int
    xs = m_element_index.lxs,
    xm = m_element_index.lxm,
    ys = m_element_index.lys,
    ym = m_element_index.lym;

  for (int j=ys; j<ys+ym; j++) {
    for (int i=xs; i<xs+xm; i++) {
      bool all_grounded_ice = (m_ice_mask.grounded_ice(i, j) and
                               m_ice_mask.grounded_ice(i+1, j) and
                               m_ice_mask.grounded_ice(i, j+1) and
                               m_ice_mask.grounded_ice(i+1, j+1));

      if (! all_grounded_ice) {
        continue;
      }

      m_element.reset(i, j);

      // Obtain values of x at the quadrature points for the element.
      m_element.nodal_values(a, a_e);
      if (dirichletBC) {
        dirichletBC.enforce_homogeneous(m_element, a_e);
      }
      quadrature_point_values(m_quadrature, a_e, a_q, dadx_q, dady_q);

      m_element.nodal_values(b, b_e);
      if (dirichletBC) {
        dirichletBC.enforce_homogeneous(m_element, b_e);
      }
      quadrature_point_values(m_quadrature, b_e, b_q, dbdx_q, dbdy_q);

      for (unsigned int q=0; q<Nq; q++) {
        value += W[q]*(m_cL2*a_q[q]*b_q[q]+ m_cH1*(dadx_q[q]*dbdx_q[q]+dady_q[q]*dbdy_q[q]));
      } // q
    } // j
  } // i

  GlobalSum(m_grid->com, &value, OUTPUT, 1);
}


void IPGroundedIceH1NormFunctional2S::gradientAt(IceModelVec2S &x, IceModelVec2S &gradient) {

  const unsigned int Nk     = fem::q1::n_chi;
  const unsigned int Nq     = m_quadrature.n();
  const unsigned int Nq_max = fem::MAX_QUADRATURE_SIZE;

  // Clear the gradient before doing anything with it!
  gradient.set(0);

  double x_e[Nk];
  double x_q[Nq_max], dxdx_q[Nq_max], dxdy_q[Nq_max];

  IceModelVec::AccessList list{&x, &gradient, &m_ice_mask};

  double gradient_e[Nk];

  // An Nq by Nk array of test function values.
  const fem::Germs *test = m_quadrature.test_function_values();

  // Jacobian times weights for quadrature.
  const double* W = m_quadrature.weights();

  fem::DirichletData_Scalar dirichletBC(m_dirichletIndices, NULL);

  // Loop through all local and ghosted elements.
  const int
    xs = m_element_index.xs,
    xm = m_element_index.xm,
    ys = m_element_index.ys,
    ym = m_element_index.ym;

  for (int j=ys; j<ys+ym; j++) {
    for (int i=xs; i<xs+xm; i++) {
      bool all_grounded_ice = (m_ice_mask.grounded_ice(i, j) and
                               m_ice_mask.grounded_ice(i+1, j) and
                               m_ice_mask.grounded_ice(i, j+1) and
                               m_ice_mask.grounded_ice(i+1, j+1));

      if (! all_grounded_ice) {
        continue;
      }

      // Reset the DOF map for this element.
      m_element.reset(i, j);

      // Obtain values of x at the quadrature points for the element.
      m_element.nodal_values(x, x_e);
      if (dirichletBC) {
        dirichletBC.constrain(m_element);
        dirichletBC.enforce_homogeneous(m_element, x_e);
      }
      quadrature_point_values(m_quadrature, x_e, x_q, dxdx_q, dxdy_q);

      // Zero out the element-local residual in prep for updating it.
      for (unsigned int k=0; k<Nk; k++) {
        gradient_e[k] = 0;
      }

      for (unsigned int q=0; q<Nq; q++) {
        const double &x_qq=x_q[q];
        const double &dxdx_qq=dxdx_q[q], &dxdy_qq=dxdy_q[q];
        for (unsigned int k=0; k<Nk; k++) {
          gradient_e[k] += 2*W[q]*(m_cL2*x_qq*test[q][k].val +
                                   m_cH1*(dxdx_qq*test[q][k].dx + dxdy_qq*test[q][k].dy));
        } // k
      } // q
      m_element.add_contribution(gradient_e, gradient);
    } // j
  } // i
}

void IPGroundedIceH1NormFunctional2S::assemble_form(Mat form) {

  const unsigned int Nk = fem::q1::n_chi;
  const unsigned int Nq = m_quadrature.n();

  PetscErrorCode ierr;
  int         i, j;

  // Zero out the Jacobian in preparation for updating it.
  ierr = MatZeroEntries(form);
  PISM_CHK(ierr, "MatZeroEntries");

  // Jacobian times weights for quadrature.
  const double* W = m_quadrature.weights();

  fem::DirichletData_Scalar zeroLocs(m_dirichletIndices, NULL);

  IceModelVec::AccessList list(m_ice_mask);

  // Values of the finite element test functions at the quadrature points.
  // This is an Nq by Nk array of function germs (Nq=#of quad pts, Nk=#of test functions).
  const fem::Germs *test = m_quadrature.test_function_values();

  // Loop through all the elements.
  const int
    xs = m_element_index.xs,
    xm = m_element_index.xm,
    ys = m_element_index.ys,
    ym = m_element_index.ym;

  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      bool all_grounded_ice = (m_ice_mask.grounded_ice(i, j) and
                               m_ice_mask.grounded_ice(i+1, j) and
                               m_ice_mask.grounded_ice(i, j+1) and
                               m_ice_mask.grounded_ice(i+1, j+1));
      if (! all_grounded_ice) {
        continue;
      }

      // Element-local Jacobian matrix (there are Nk vector valued degrees
      // of freedom per elment, for a total of (2*Nk)*(2*Nk) = 16
      // entries in the local Jacobian.
      double K[Nk][Nk];

      // Initialize the map from global to local degrees of freedom for this element.
      m_element.reset(i, j);

      // Don't update rows/cols where we project to zero.
      if (zeroLocs) {
        zeroLocs.constrain(m_element);
      }

      // Build the element-local Jacobian.
      ierr = PetscMemzero(K, sizeof(K));
      PISM_CHK(ierr, "PetscMemzero");

      for (unsigned int q=0; q<Nq; q++) {
        for (unsigned int k = 0; k < Nk; k++) {   // Test functions
          const fem::Germ &test_qk=test[q][k];
          for (unsigned int l = 0; l < Nk; l++) { // Trial functions
            const fem::Germ &test_ql=test[q][l];
            K[k][l]     += W[q]*(m_cL2*test_qk.val*test_ql.val
              +  m_cH1*(test_qk.dx*test_ql.dx + test_qk.dy*test_ql.dy));
          } // l
        } // k
      } // q
      m_element.add_contribution(&K[0][0], form);
    } // j
  } // i

  if (zeroLocs) {
    zeroLocs.fix_jacobian(form);
  }

  ierr = MatAssemblyBegin(form, MAT_FINAL_ASSEMBLY);
  PISM_CHK(ierr, "MatAssemblyBegin");

  ierr = MatAssemblyEnd(form, MAT_FINAL_ASSEMBLY);
  PISM_CHK(ierr, "MatAssemblyEnd");
}

} // end of namespace inverse
} // end of namespace pism
