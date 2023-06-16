// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2020, 2022  David Maxwell and Constantine Khroulev
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

#include "IPTotalVariationFunctional.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/array/Scalar.hh"

namespace pism {
namespace inverse {

IPTotalVariationFunctional2S::IPTotalVariationFunctional2S(std::shared_ptr<const IceGrid> grid,
                                                           double c, double exponent, double eps,
                                                           array::Scalar *dirichletLocations) :
    IPFunctional<array::Scalar>(grid), m_dirichletIndices(dirichletLocations),
    m_c(c), m_lebesgue_exp(exponent), m_epsilon_sq(eps*eps) {
}

void IPTotalVariationFunctional2S::valueAt(array::Scalar &x, double *OUTPUT) {

  const unsigned int Nk     = fem::q1::n_chi;
  const unsigned int Nq     = m_element.n_pts();
  const unsigned int Nq_max = fem::MAX_QUADRATURE_SIZE;

  // The value of the objective
  double value = 0;

  double x_e[Nk];
  double x_q[Nq_max], dxdx_q[Nq_max], dxdy_q[Nq_max];

  array::AccessScope list(x);

  fem::DirichletData_Scalar dirichletBC(m_dirichletIndices, NULL);

  // Loop through all LOCAL elements.
  const int
    xs = m_element_index.lxs,
    xm = m_element_index.lxm,
    ys = m_element_index.lys,
    ym = m_element_index.lym;

  for (int j = ys; j < ys + ym; j++) {
    for (int i = xs; i < xs + xm; i++) {
      m_element.reset(i, j);

      // Obtain values of x at the quadrature points for the element.
      m_element.nodal_values(x.array(), x_e);
      if (dirichletBC) {
        dirichletBC.enforce_homogeneous(m_element, x_e);
      }
      m_element.evaluate(x_e, x_q, dxdx_q, dxdy_q);

      for (unsigned int q = 0; q < Nq; q++) {
        auto W = m_element.weight(q);
        value += m_c*W*pow(m_epsilon_sq + dxdx_q[q]*dxdx_q[q] + dxdy_q[q]*dxdy_q[q], m_lebesgue_exp / 2);
      } // q
    } // j
  } // i

  GlobalSum(m_grid->com, &value, OUTPUT, 1);
}

void IPTotalVariationFunctional2S::gradientAt(array::Scalar &x, array::Scalar &gradient) {

  const unsigned int Nk     = fem::q1::n_chi;
  const unsigned int Nq     = m_element.n_pts();
  const unsigned int Nq_max = fem::MAX_QUADRATURE_SIZE;

  // Clear the gradient before doing anything with it.
  gradient.set(0);

  double x_e[Nk];
  double x_q[Nq_max], dxdx_q[Nq_max], dxdy_q[Nq_max];

  double gradient_e[Nk];

  array::AccessScope list{&x, &gradient};

  fem::DirichletData_Scalar dirichletBC(m_dirichletIndices, NULL);

  // Loop through all local and ghosted elements.
  const int
    xs = m_element_index.xs,
    xm = m_element_index.xm,
    ys = m_element_index.ys,
    ym = m_element_index.ym;

  for (int j = ys; j < ys + ym; j++) {
    for (int i = xs; i < xs + xm; i++) {

      // Reset the DOF map for this element.
      m_element.reset(i, j);

      // Obtain values of x at the quadrature points for the element.
      m_element.nodal_values(x.array(), x_e);
      if (dirichletBC) {
        dirichletBC.constrain(m_element);
        dirichletBC.enforce_homogeneous(m_element, x_e);
      }
      m_element.evaluate(x_e, x_q, dxdx_q, dxdy_q);

      // Zero out the element-local residual in preparation for updating it.
      for (unsigned int k = 0; k < Nk; k++) {
        gradient_e[k] = 0;
      }

      for (unsigned int q = 0; q < Nq; q++) {
        auto W = m_element.weight(q);
        const double &dxdx_qq = dxdx_q[q], &dxdy_qq = dxdy_q[q];
        for (unsigned int k = 0; k < Nk; k++) {
          gradient_e[k] += m_c*W*(m_lebesgue_exp)*pow(m_epsilon_sq + dxdx_q[q]*dxdx_q[q] + dxdy_q[q]*dxdy_q[q], m_lebesgue_exp / 2 - 1)
            *(dxdx_qq*m_element.chi(q, k).dx + dxdy_qq*m_element.chi(q, k).dy);
        } // k
      } // q
      m_element.add_contribution(gradient_e, gradient.array());
    } // j
  } // i
}

} // end of namespace inverse
} // end of namespace pism
