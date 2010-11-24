// Copyright (C) 2010 Constantine Khroulev and Ed Bueler
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
#include "PISMVars.hh"
#include "flowlaws.hh"
#include "enthalpyConverter.hh"

//! Shallow stress balance modifier (such as the non-sliding SIA).
class SSB_Modifier
{
public:
  SSB_Modifier(IceGrid &g, IceFlowLaw &i, EnthalpyConverter &e, const NCConfigVariable &c)
    : grid(g), ice(i), EC(e), config(c)
  { D_max = u_max = v_max = 0.0; allocate(); }
  virtual ~SSB_Modifier() {}

  virtual PetscErrorCode init(PISMVars &/*vars*/) { return 0; }

  //! \brief Adds more variable names to result (to respect -o_size and
  //! -save_size).
  /*!
    Keyword can be one of "small", "medium" or "big".
   */
  virtual void add_to_output(string /*keyword*/, set<string> &/*result*/) {}

  //! Writes requested fields to a file.
  virtual PetscErrorCode write_fields(set<string> /*vars*/, string /*filename*/)
  { return 0; }

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
protected:
  virtual PetscErrorCode allocate();

  IceGrid &grid;
  IceFlowLaw &ice;
  EnthalpyConverter &EC;
  const NCConfigVariable &config;
  PetscReal D_max, u_max, v_max;
  IceModelVec2Stag diffusive_flux;
  IceModelVec3 u, v, Sigma;
};


//! The trivial Shallow Stress Balance modifier.
class SSBM_Trivial : public SSB_Modifier
{
public:
  SSBM_Trivial(IceGrid &g, IceFlowLaw &i, EnthalpyConverter &e, const NCConfigVariable &c)
    : SSB_Modifier(g, i, e, c) {}
  virtual ~SSBM_Trivial() {}

  virtual PetscErrorCode init(PISMVars &vars);

  virtual PetscErrorCode update(IceModelVec2V *vel_input, IceModelVec2S *D2_input,
                                bool fast);
protected:
  PetscErrorCode compute_sigma(IceModelVec2S *D2_input, IceModelVec3 &result);
  IceModelVec3 *enthalpy;
  IceModelVec2S *thickness;
};
#endif /* _SSB_MODIFIER_H_ */
