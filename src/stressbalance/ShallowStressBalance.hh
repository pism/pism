// Copyright (C) 2010--2019, 2021, 2022, 2025 Constantine Khroulev and Ed Bueler
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

#ifndef _SHALLOWSTRESSBALANCE_H_
#define _SHALLOWSTRESSBALANCE_H_

#include "pism/util/Component.hh"
#include "pism/util/array/Vector.hh"
#include "pism/util/EnthalpyConverter.hh"

namespace pism {
namespace rheology {
class FlowLaw;
}

class Grid;
class IceBasalResistancePlasticLaw;

namespace stressbalance {

class Inputs;

//! Shallow stress balance (such as the SSA).
class ShallowStressBalance : public Component {
public:
  ShallowStressBalance(std::shared_ptr<const Grid> g);
  virtual ~ShallowStressBalance();

  //  initialization and I/O:

  void init();

  virtual void update(const Inputs &inputs, bool full_update) = 0;

  //! \brief Get the thickness-advective 2D velocity.
  const array::Vector1& velocity() const;

  //! \brief Get the basal frictional heating (for the adaptive energy time-stepping).
  const array::Scalar& basal_frictional_heating();

  void compute_basal_frictional_heating(const array::Vector &velocity,
                                        const array::Scalar &tauc,
                                        const array::CellType &mask,
                                        array::Scalar &result) const;
  // helpers:

  //! \brief Produce a report string for the standard output.
  virtual std::string stdout_report() const;

  std::shared_ptr<const rheology::FlowLaw> flow_law() const;

  std::shared_ptr<EnthalpyConverter> enthalpy_converter() const;

  const IceBasalResistancePlasticLaw* sliding_law() const;

  double flow_enhancement_factor() const;
protected:
  virtual void init_impl();

  virtual DiagnosticList diagnostics_impl() const;

  IceBasalResistancePlasticLaw *m_basal_sliding_law;
  std::shared_ptr<rheology::FlowLaw> m_flow_law;
  std::shared_ptr<EnthalpyConverter> m_EC;

  array::Vector2 m_velocity;
  array::Scalar m_basal_frictional_heating;

  //! flow enhancement factor
  double m_e_factor;
};

//! Returns zero velocity field, zero friction heating, and zero for D^2.
/*!
  This derived class is used in the non-sliding SIA approximation. This
  implementation ignores any basal resistance fields (e.g. yield stress from
  the IceModel or other user of this class).
*/
class ZeroSliding : public ShallowStressBalance {
public:
  ZeroSliding(std::shared_ptr<const Grid> g);
  virtual ~ZeroSliding() = default;

  virtual void update(const Inputs &inputs, bool full_update);

protected:
};

class PrescribedSliding : public ZeroSliding {
public:
  PrescribedSliding(std::shared_ptr<const Grid> g);
  virtual ~PrescribedSliding() = default;
  virtual void update(const Inputs &inputs, bool full_update);
protected:
  virtual void init_impl();
};

} // end of namespace stressbalance
} // end of namespace pism

#endif /* _SHALLOWSTRESSBALANCE_H_ */
