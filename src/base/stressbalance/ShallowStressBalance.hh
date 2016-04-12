// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016 Constantine Khroulev and Ed Bueler
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

#include "base/util/PISMComponent.hh"
#include "base/util/iceModelVec.hh"
#include "base/enthalpyConverter.hh"

namespace pism {
namespace rheology {
class FlowLaw;
}

class IceGrid;
class IceBasalResistancePlasticLaw;
class IceModelVec2CellType;

namespace stressbalance {

//! Shallow stress balance (such as the SSA).
class ShallowStressBalance : public Component {
public:
  ShallowStressBalance(IceGrid::ConstPtr g, EnthalpyConverter::Ptr e);
  virtual ~ShallowStressBalance();

  //  initialization and I/O:

  void init();
  void set_boundary_conditions(const IceModelVec2Int &locations,
                               const IceModelVec2V &velocities);

  virtual void update(bool fast,
                      double sea_level,
                      const IceModelVec2S &melange_back_pressure) = 0;

  //! \brief Get the thickness-advective 2D velocity.
  const IceModelVec2V& velocity();

  //! \brief Get the basal frictional heating (for the adaptive energy time-stepping).
  const IceModelVec2S& basal_frictional_heating();

  void compute_2D_principal_strain_rates(const IceModelVec2V &velocity,
                                         const IceModelVec2CellType &mask,
                                         IceModelVec2 &result);

  void compute_2D_stresses(const IceModelVec2V &velocity,
                           const IceModelVec2CellType &mask,
                           IceModelVec2 &result);

  void compute_basal_frictional_heating(const IceModelVec2V &velocity,
                                        const IceModelVec2S &tauc,
                                        const IceModelVec2CellType &mask,
                                        IceModelVec2S &result);
  // helpers:

  //! \brief Produce a report string for the standard output.
  virtual std::string stdout_report();

  const rheology::FlowLaw* flow_law();

  EnthalpyConverter::Ptr enthalpy_converter();

  const IceBasalResistancePlasticLaw* sliding_law();
protected:
  virtual void init_impl();
  
  virtual void get_diagnostics_impl(std::map<std::string, Diagnostic::Ptr> &dict,
                                    std::map<std::string, TSDiagnostic::Ptr> &ts_dict);

  double m_sea_level;
  IceBasalResistancePlasticLaw *m_basal_sliding_law;
  rheology::FlowLaw *m_flow_law;
  EnthalpyConverter::Ptr m_EC;

  IceModelVec2V m_velocity;
  const IceModelVec2V *m_bc_values;
  const IceModelVec2Int *m_bc_mask;
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
  ZeroSliding(IceGrid::ConstPtr g, EnthalpyConverter::Ptr e);
  virtual ~ZeroSliding();
  
  virtual void update(bool fast, double sea_level, const IceModelVec2S &melange_back_pressure);

  //! Writes requested couplings fields to file and/or asks an attached
  //! model to do so.
protected:
  virtual void write_variables_impl(const std::set<std::string> &/*vars*/, const PIO &/*nc*/);
  virtual void add_vars_to_output_impl(const std::string &keyword,
                                       std::set<std::string> &result);
  virtual void define_variables_impl(const std::set<std::string> &/*vars*/, const PIO &/*nc*/,
                                     IO_Type /*nctype*/);
};

class PrescribedSliding : public ZeroSliding {
public:
  PrescribedSliding(IceGrid::ConstPtr g, EnthalpyConverter::Ptr e);
  virtual ~PrescribedSliding();
  virtual void update(bool fast, double sea_level, const IceModelVec2S &melange_back_pressure);
  virtual void init();
};

} // end of namespace stressbalance
} // end of namespace pism

#endif /* _SHALLOWSTRESSBALANCE_H_ */
