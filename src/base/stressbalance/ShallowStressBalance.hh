// Copyright (C) 2010, 2011, 2012, 2013, 2014 Constantine Khroulev and Ed Bueler
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

#include "PISMComponent.hh"
#include "iceModelVec.hh"
#include "IceGrid.hh"
#include "flowlaws.hh"
#include "flowlaw_factory.hh"
#include <PISMDiagnostic.hh>

namespace pism {

class IceFlowLaw;
class EnthalpyConverter;
class IceBasalResistancePlasticLaw;

//! Shallow stress balance (such as the SSA).
class ShallowStressBalance : public Component
{
public:
  ShallowStressBalance(const IceGrid &g, EnthalpyConverter &e);
  virtual ~ShallowStressBalance();

  //  initialization and I/O:

  virtual void init() {
  }

  virtual void set_boundary_conditions(IceModelVec2Int &locations,
                                                 IceModelVec2V &velocities) {
    m_vel_bc = &velocities;
    bc_locations = &locations;
  }

  //! \brief Set the sea level used to check for floatation. (Units: meters,
  //! relative to the geoid.)
  void set_sea_level_elevation(double new_sea_level) {
    sea_level = new_sea_level;
  }

  virtual void update(bool fast, IceModelVec2S &melange_back_pressure) = 0;

  // interface to the data provided by the stress balance object:
  virtual void get_diagnostics(std::map<std::string, Diagnostic*> &dict,
                               std::map<std::string, TSDiagnostic*> &/*ts_dict*/);

  //! \brief Get the thickness-advective (SSA) 2D velocity.
  virtual void get_2D_advective_velocity(IceModelVec2V* &result) {
    result = &m_velocity;
  }

  //! \brief Get the basal frictional heating (for the adaptive energy time-stepping).
  virtual void get_basal_frictional_heating(IceModelVec2S* &result) {
    result = &basal_frictional_heating;
  }

  virtual void compute_2D_principal_strain_rates(IceModelVec2V &velocity,
                                                           IceModelVec2Int &mask,
                                                           IceModelVec2 &result);

  virtual void compute_2D_stresses(IceModelVec2V &velocity, IceModelVec2Int &mask,
                                             IceModelVec2 &result);

  virtual void compute_basal_frictional_heating(IceModelVec2V &velocity,
                                                IceModelVec2S &tauc,
                                                IceModelVec2Int &mask,
                                                IceModelVec2S &result);
  // helpers:

  //! \brief Produce a report string for the standard output.
  virtual void stdout_report(std::string &result) {
    result = "";
  }

  const IceFlowLaw* get_flow_law() {
    return flow_law;
  }

  EnthalpyConverter& get_enthalpy_converter() {
    return EC;
  }

  const IceBasalResistancePlasticLaw* get_sliding_law() {
    return basal_sliding_law;
  }
protected:
  double sea_level;
  IceBasalResistancePlasticLaw *basal_sliding_law;
  IceFlowLaw *flow_law;
  EnthalpyConverter &EC;

  IceModelVec2V m_velocity, *m_vel_bc;
  IceModelVec2Int *bc_locations;
  IceModelVec2S basal_frictional_heating;
};

//! Returns zero velocity field, zero friction heating, and zero for D^2.
/*!
  This derived class is used in the non-sliding SIA approximation. This
  implementation ignores any basal resistance fields (e.g. yield stress from
  the IceModel or other user of this class).
*/
class ZeroSliding : public ShallowStressBalance
{
public:
  ZeroSliding(const IceGrid &g, EnthalpyConverter &e);
  virtual ~ZeroSliding();
  
  virtual void update(bool fast, IceModelVec2S &melange_back_pressure);

  virtual void add_vars_to_output(const std::string &/*keyword*/, std::set<std::string> &/*result*/);

  //! Defines requested couplings fields and/or asks an attached model
  //! to do so.
  virtual void define_variables(const std::set<std::string> &/*vars*/, const PIO &/*nc*/,
                                          IO_Type /*nctype*/);

  //! Writes requested couplings fields to file and/or asks an attached
  //! model to do so.
  virtual void write_variables(const std::set<std::string> &/*vars*/, const PIO &/*nc*/);
};

class PrescribedSliding : public ZeroSliding {
public:
  PrescribedSliding(const IceGrid &g, EnthalpyConverter &e);
  virtual ~PrescribedSliding();
  virtual void update(bool fast, IceModelVec2S &melange_back_pressure);
  virtual void init();
};

} // end of namespace pism

#endif /* _SHALLOWSTRESSBALANCE_H_ */
