// Copyright (C) 2010, 2011, 2012 Constantine Khroulev and Ed Bueler
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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
class SSB_Modifier : public PISMComponent_Diag
{
public:
  SSB_Modifier(IceGrid &g, EnthalpyConverter &e, const NCConfigVariable &c)
    : PISMComponent_Diag(g, c), EC(e)
  { D_max = u_max = v_max = 0.0; variables = NULL; allocate(); }
  virtual ~SSB_Modifier() {}

  virtual PetscErrorCode init(PISMVars &vars) { variables = &vars; return 0; }

  using PISMComponent_Diag::update;
  virtual PetscErrorCode update(IceModelVec2V *vel_input,
                                IceModelVec2S *D2_input,
                                bool fast) = 0;

  //! \brief Get the diffusive (SIA) vertically-averaged flux on the staggered grid.
  virtual PetscErrorCode get_diffusive_flux(IceModelVec2Stag* &result)
  { result = &diffusive_flux; return 0; }

  //! \brief Get the max diffusivity (for the adaptive time-stepping).
  virtual PetscErrorCode get_max_diffusivity(PetscReal &result)
  { result = D_max; return 0; }

  virtual PetscErrorCode get_horizontal_3d_velocity(IceModelVec3* &u_result, IceModelVec3* &v_result)
  { u_result = &u; v_result = &v; return 0; }

  virtual PetscErrorCode get_max_horizontal_velocity(PetscReal &max_u, PetscReal &max_v)
  { max_u = u_max; max_v = v_max; return 0; }

  virtual PetscErrorCode get_volumetric_strain_heating(IceModelVec3* &result)
  { result = &Sigma; return 0; }

  //! \brief Extends the computational grid (vertically).
  virtual PetscErrorCode extend_the_grid(PetscInt old_Mz);

  virtual PetscErrorCode stdout_report(string &result)
  { result = ""; return 0; }

  IceFlowLaw* get_flow_law()
  { return flow_law; }
protected:
  virtual PetscErrorCode allocate();

  IceFlowLaw *flow_law;
  EnthalpyConverter &EC;
  PetscReal D_max, u_max, v_max;
  IceModelVec2Stag diffusive_flux;
  IceModelVec3 u, v, Sigma;

  PISMVars *variables;
};


//! The trivial Shallow Stress Balance modifier.
class SSBM_Trivial : public SSB_Modifier
{
public:
  SSBM_Trivial(IceGrid &g, EnthalpyConverter &e, const NCConfigVariable &c);
  virtual ~SSBM_Trivial();

  virtual PetscErrorCode init(PISMVars &vars);

  using PISMComponent_Diag::update;
  virtual PetscErrorCode update(IceModelVec2V *vel_input, IceModelVec2S *D2_input,
                                bool fast);
  virtual void add_vars_to_output(string /*keyword*/,
                                  map<string,NCSpatialVariable> &/*result*/)
  { }

  //! Defines requested couplings fields to file and/or asks an attached
  //! model to do so.
  virtual PetscErrorCode define_variables(set<string> /*vars*/, const PIO &/*nc*/,
                                          PISM_IO_Type /*nctype*/)
  { return 0; }

  //! Writes requested couplings fields to file and/or asks an attached
  //! model to do so.
  virtual PetscErrorCode write_variables(set<string> /*vars*/, const PIO &/*nc*/)
  { return 0; }

protected:
  PetscErrorCode compute_sigma(IceModelVec2S *D2_input, IceModelVec3 &result);
  IceModelVec3 *enthalpy;
  IceModelVec2S *thickness;
};
#endif /* _SSB_MODIFIER_H_ */
