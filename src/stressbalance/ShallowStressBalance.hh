// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 Constantine Khroulev and Ed Bueler
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
#include "pism/util/iceModelVec.hh"
#include "pism/util/EnthalpyConverter.hh"

namespace pism {
namespace rheology {
class FlowLaw;
}

class IceGrid;
class IceBasalResistancePlasticLaw;
class IceModelVec2CellType;

namespace stressbalance {

/*!
 * Vertically-averaged ocean pressure difference at the calving front, used in the implementation of
 * the stress boundary condition at the calving front in SSA stress balance solvers.
 */
double ocean_pressure_difference(bool shelf, bool dry_mode, double H, double bed, double sea_level,
                                 double rho_ice, double rho_ocean, double g);

class Inputs;

//! Shallow stress balance (such as the SSA).
class ShallowStressBalance : public Component {
public:
  ShallowStressBalance(IceGrid::ConstPtr g);
  virtual ~ShallowStressBalance();

  //  initialization and I/O:

  void init();

  virtual void update(const Inputs &inputs, bool full_update) = 0;

  //! \brief Get the thickness-advective 2D velocity.
  const IceModelVec2V& velocity() const;

  //! \brief Get the basal frictional heating (for the adaptive energy time-stepping).
  const IceModelVec2S& basal_frictional_heating();

  void compute_2D_stresses(const IceModelVec2V &velocity,
                           const IceModelVec2S &hardness,
                           const IceModelVec2CellType &cell_type,
                           IceModelVec2 &result) const;

  void compute_basal_frictional_heating(const IceModelVec2V &velocity,
                                        const IceModelVec2S &tauc,
                                        const IceModelVec2CellType &mask,
                                        IceModelVec2S &result) const;
  // helpers:

  //! \brief Produce a report string for the standard output.
  virtual std::string stdout_report() const;

  std::shared_ptr<const rheology::FlowLaw> flow_law() const;

  EnthalpyConverter::Ptr enthalpy_converter() const;

  const IceBasalResistancePlasticLaw* sliding_law() const;
protected:
  virtual void init_impl();
  
  virtual DiagnosticList diagnostics_impl() const;

  IceBasalResistancePlasticLaw *m_basal_sliding_law;
  std::shared_ptr<rheology::FlowLaw> m_flow_law;
  EnthalpyConverter::Ptr m_EC;

  IceModelVec2V m_velocity;
  IceModelVec2S m_basal_frictional_heating;
};

//! Returns zero velocity field, zero friction heating, and zero for D^2.
/*!
  This derived class is used in the non-sliding SIA approximation. This
  implementation ignores any basal resistance fields (e.g. yield stress from
  the IceModel or other user of this class).
*/
class ZeroSliding : public ShallowStressBalance {
public:
  ZeroSliding(IceGrid::ConstPtr g);
  virtual ~ZeroSliding();
  
  virtual void update(const Inputs &inputs, bool full_update);

protected:
};

class PrescribedSliding : public ZeroSliding {
public:
  PrescribedSliding(IceGrid::ConstPtr g);
  virtual ~PrescribedSliding();
  virtual void update(const Inputs &inputs, bool full_update);
protected:
  virtual void init_impl();
};

} // end of namespace stressbalance
} // end of namespace pism

#endif /* _SHALLOWSTRESSBALANCE_H_ */
