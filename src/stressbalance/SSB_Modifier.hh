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

#ifndef _SSB_MODIFIER_H_
#define _SSB_MODIFIER_H_

#include "pism/util/iceModelVec.hh"
#include "pism/util/Component.hh"
#include "pism/util/EnthalpyConverter.hh"

namespace pism {

class Vars;

namespace rheology {
class FlowLaw;
}

namespace stressbalance {

class Inputs;

//! Shallow stress balance modifier (such as the non-sliding SIA).
class SSB_Modifier : public Component {
public:
  SSB_Modifier(IceGrid::ConstPtr g);
  virtual ~SSB_Modifier();

  virtual void init();

  virtual void update(const IceModelVec2V &sliding_velocity,
                      const Inputs &inputs,
                      bool full_update) = 0;

  //! \brief Get the diffusive (SIA) vertically-averaged flux on the staggered grid.
  virtual const IceModelVec2Stag& diffusive_flux();

  //! \brief Get the max diffusivity (for the adaptive time-stepping).
  virtual double max_diffusivity() const;

  const IceModelVec3& velocity_u() const;

  const IceModelVec3& velocity_v() const;

  const IceModelVec3& volumetric_strain_heating() const;

  virtual std::string stdout_report() const;

  std::shared_ptr<const rheology::FlowLaw> flow_law() const;

protected:
  std::shared_ptr<rheology::FlowLaw> m_flow_law;
  EnthalpyConverter::Ptr m_EC;
  double m_D_max;
  IceModelVec2Stag m_diffusive_flux;
  IceModelVec3 m_u, m_v, m_strain_heating;
};


//! The trivial Shallow Stress Balance modifier.
class ConstantInColumn : public SSB_Modifier {
public:
  ConstantInColumn(IceGrid::ConstPtr g);
  virtual ~ConstantInColumn();

  virtual void init();

  virtual void update(const IceModelVec2V &sliding_velocity,
                      const Inputs &inputs,
                      bool full_update);
};

} // end of namespace stressbalance
} // end of namespace pism

#endif /* _SSB_MODIFIER_H_ */
