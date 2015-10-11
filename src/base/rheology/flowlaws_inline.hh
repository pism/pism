/* Copyright (C) 2015 PISM Authors
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

#ifndef _FLOWLAWS_INLINE_H_
#define _FLOWLAWS_INLINE_H_

// This header is included by flowlaws.hh. Do not include it
// manually.

namespace pism {
namespace rheology {

//! \brief Computes the regularized effective viscosity and its derivative with respect to the
//! second invariant \f$ \gamma \f$.
/*!
 *
 * @f{align*}{
 * \nu &= \frac{1}{2} B \left( \epsilon + \gamma \right)^{(1-n)/(2n)},\\
 * \diff{\nu}{\gamma} &= \frac{1}{2} B \cdot \frac{1-n}{2n} \cdot \left(\epsilon + \gamma \right)^{(1-n)/(2n) - 1}, \\
 * &= \frac{1-n}{2n} \cdot \frac{1}{2} B \left( \epsilon + \gamma \right)^{(1-n)/(2n)} \cdot \frac{1}{\epsilon + \gamma}, \\
 * &= \frac{1-n}{2n} \cdot \frac{\nu}{\epsilon + \gamma}.
 * @f}
 * Here @f$ \gamma @f$ is the second invariant
 * @f{align*}{
 * \gamma &= \frac{1}{2} D_{ij} D_{ij}\\
 * &= \frac{1}{2}\, ((u_x)^2 + (v_y)^2 + (u_x + v_y)^2 + \frac{1}{2}\, (u_y + v_x)^2) \\
 * @f}
 * and
 * @f[ D_{ij}(\mathbf{u}) = \frac{1}{2}\left(\diff{u_{i}}{x_{j}} + \diff{u_{j}}{x_{i}}\right). @f]
 *
 * Either one of \c nu and \c dnu can be NULL if the corresponding output is not needed.
 *
 * \param[in] hardness ice hardness
 * \param[in] gamma the second invariant
 * \param[out] nu effective viscosity
 * \param[out] dnu derivative of \f$ \nu \f$ with respect to \f$ \gamma \f$
 */
inline void FlowLaw::effective_viscosity(double hardness, double gamma,
                                         double *nu, double *dnu) const {
  const double
    my_nu = 0.5 * hardness * pow(m_schoofReg + gamma, m_viscosity_power);

  if (PetscLikely(nu != NULL)) {
    *nu = my_nu;
  }

  if (PetscLikely(dnu != NULL)) {
    *dnu = m_viscosity_power * my_nu / (m_schoofReg + gamma);
  }
}

inline double FlowLaw::exponent() const {
  return m_n;
}

inline double FlowLaw::enhancement_factor() const {
  return m_e;
}


} // end of namespace rheology
} // end of namespace pism

#endif /* _FLOWLAWS_INLINE_H_ */
