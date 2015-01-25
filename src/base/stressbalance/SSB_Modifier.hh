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
class SSB_Modifier : public Component {
public:
  SSB_Modifier(const IceGrid &g, const EnthalpyConverter &e);
  virtual ~SSB_Modifier();

  virtual void init();

  virtual void update(const IceModelVec2V &vel_input, bool fast) = 0;

  //! \brief Get the diffusive (SIA) vertically-averaged flux on the staggered grid.
  virtual const IceModelVec2Stag& diffusive_flux();

  //! \brief Get the max diffusivity (for the adaptive time-stepping).
  virtual double max_diffusivity();

  const IceModelVec3& velocity_u();

  const IceModelVec3& velocity_v();

  const IceModelVec3& volumetric_strain_heating();

  virtual std::string stdout_report();

  IceFlowLaw* flow_law();
protected:
  IceFlowLaw *m_flow_law;
  const EnthalpyConverter &m_EC;
  double m_D_max;
  IceModelVec2Stag m_diffusive_flux;
  IceModelVec3 m_u, m_v, m_strain_heating;
};


//! The trivial Shallow Stress Balance modifier.
class ConstantInColumn : public SSB_Modifier {
public:
  ConstantInColumn(const IceGrid &g, const EnthalpyConverter &e);
  virtual ~ConstantInColumn();

  virtual void init();

  virtual void update(const IceModelVec2V &vel_input, bool fast);

  //! Writes requested couplings fields to file and/or asks an attached
  //! model to do so.
  virtual void write_variables(const std::set<std::string> &vars, const PIO &nc);
protected:
  virtual void add_vars_to_output_impl(const std::string &keyword, std::set<std::string> &result);
  virtual void define_variables_impl(const std::set<std::string> &vars, const PIO &nc,
                                IO_Type nctype);
};

} // end of namespace pism

#endif /* _SSB_MODIFIER_H_ */
