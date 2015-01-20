// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015 Constantine Khroulev and Ed Bueler
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

#ifndef _SSB_MODIFIER_H_
#define _SSB_MODIFIER_H_

#include "iceModelVec.hh"
#include "PISMComponent.hh"

namespace pism {

class Vars;
class IceFlowLaw;
class EnthalpyConverter;

//! Shallow stress balance modifier (such as the non-sliding SIA).
class SSB_Modifier : public Component
{
public:
  SSB_Modifier(const IceGrid &g, EnthalpyConverter &e);
  virtual ~SSB_Modifier();

  virtual void init() {
  }

  virtual void update(IceModelVec2V *vel_input, bool fast) = 0;

  //! \brief Get the diffusive (SIA) vertically-averaged flux on the staggered grid.
  virtual void get_diffusive_flux(IceModelVec2Stag* &result) {
    result = &m_diffusive_flux;
  }

  //! \brief Get the max diffusivity (for the adaptive time-stepping).
  virtual void get_max_diffusivity(double &result) {
    result = m_D_max;
  }

  virtual void get_horizontal_3d_velocity(IceModelVec3* &u_result,
                                          IceModelVec3* &v_result) {
    u_result = &m_u;
    v_result = &m_v;
  }

  virtual void get_volumetric_strain_heating(IceModelVec3* &result) {
    result = &m_strain_heating;
  }

  virtual void stdout_report(std::string &result) {
    result = "";
  }

  IceFlowLaw* get_flow_law() {
    return m_flow_law;
  }
protected:
  IceFlowLaw *m_flow_law;
  EnthalpyConverter &m_EC;
  double m_D_max;
  IceModelVec2Stag m_diffusive_flux;
  IceModelVec3 m_u, m_v, m_strain_heating;
};


//! The trivial Shallow Stress Balance modifier.
class ConstantInColumn : public SSB_Modifier
{
public:
  ConstantInColumn(const IceGrid &g, EnthalpyConverter &e);
  virtual ~ConstantInColumn();

  virtual void init();

  virtual void update(IceModelVec2V *vel_input, bool fast);
  virtual void add_vars_to_output(const std::string &/*keyword*/, std::set<std::string> &/*result*/) {
  }

  //! Defines requested couplings fields to file and/or asks an attached
  //! model to do so.
  virtual void define_variables(const std::set<std::string> &/*vars*/, const PIO &/*nc*/,
                                          IO_Type /*nctype*/) {
    // empty
  }

  //! Writes requested couplings fields to file and/or asks an attached
  //! model to do so.
  virtual void write_variables(const std::set<std::string> &/*vars*/, const PIO &/*nc*/) {
    // empty
  }

};

} // end of namespace pism

#endif /* _SSB_MODIFIER_H_ */
