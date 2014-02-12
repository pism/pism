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

#ifndef _SSB_MODIFIER_H_
#define _SSB_MODIFIER_H_

#include "iceModelVec.hh"
#include "PISMComponent.hh"

class PISMVars;
class IceFlowLaw;
class EnthalpyConverter;

//! Shallow stress balance modifier (such as the non-sliding SIA).
class SSB_Modifier : public PISMComponent
{
public:
  SSB_Modifier(IceGrid &g, EnthalpyConverter &e, const PISMConfig &c)
    : PISMComponent(g, c), EC(e)
  { D_max = 0.0; variables = NULL; allocate(); }
  virtual ~SSB_Modifier() {}

  virtual PetscErrorCode init(PISMVars &vars) { variables = &vars; return 0; }

  virtual PetscErrorCode update(IceModelVec2V *vel_input, bool fast) = 0;

  //! \brief Get the diffusive (SIA) vertically-averaged flux on the staggered grid.
  virtual PetscErrorCode get_diffusive_flux(IceModelVec2Stag* &result)
  { result = &diffusive_flux; return 0; }

  //! \brief Get the max diffusivity (for the adaptive time-stepping).
  virtual PetscErrorCode get_max_diffusivity(double &result)
  { result = D_max; return 0; }

  virtual PetscErrorCode get_horizontal_3d_velocity(IceModelVec3* &u_result,
                                                    IceModelVec3* &v_result)
  { u_result = &u; v_result = &v; return 0; }

  virtual PetscErrorCode get_volumetric_strain_heating(IceModelVec3* &result)
  { result = &strain_heating; return 0; }

  //! \brief Extends the computational grid (vertically).
  virtual PetscErrorCode extend_the_grid(int old_Mz);

  virtual PetscErrorCode stdout_report(std::string &result)
  { result = ""; return 0; }

  IceFlowLaw* get_flow_law()
  { return flow_law; }
protected:
  virtual PetscErrorCode allocate();

  IceFlowLaw *flow_law;
  EnthalpyConverter &EC;
  double D_max;
  IceModelVec2Stag diffusive_flux;
  IceModelVec3 u, v, strain_heating;

  PISMVars *variables;
};


//! The trivial Shallow Stress Balance modifier.
class ConstantInColumn : public SSB_Modifier
{
public:
  ConstantInColumn(IceGrid &g, EnthalpyConverter &e, const PISMConfig &c);
  virtual ~ConstantInColumn();

  virtual PetscErrorCode init(PISMVars &vars);

  virtual PetscErrorCode update(IceModelVec2V *vel_input, bool fast);
  virtual void add_vars_to_output(std::string /*keyword*/, std::set<std::string> &/*result*/)
  { }

  //! Defines requested couplings fields to file and/or asks an attached
  //! model to do so.
  virtual PetscErrorCode define_variables(std::set<std::string> /*vars*/, const PIO &/*nc*/,
                                          PISM_IO_Type /*nctype*/)
  { return 0; }

  //! Writes requested couplings fields to file and/or asks an attached
  //! model to do so.
  virtual PetscErrorCode write_variables(std::set<std::string> /*vars*/, const PIO &/*nc*/)
  { return 0; }
};
#endif /* _SSB_MODIFIER_H_ */
