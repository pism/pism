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

#ifndef _SHALLOWSTRESSBALANCE_H_
#define _SHALLOWSTRESSBALANCE_H_

#include "PISMComponent.hh"
#include "iceModelVec.hh"
#include "IceGrid.hh"
#include "flowlaws.hh"

class PISMVars;
class IceFlowLaw;
class EnthalpyConverter;
class IceBasalResistancePlasticLaw;

//! Shallow stress balance (such as the SSA).
class ShallowStressBalance : public PISMComponent_Diag
{
public:
  ShallowStressBalance(IceGrid &g, IceBasalResistancePlasticLaw &b,
                       EnthalpyConverter &e, const NCConfigVariable &conf)
    : PISMComponent_Diag(g, conf), basal(b), flow_law(NULL), EC(e)
  {
    vel_bc = NULL; bc_locations = NULL; variables = NULL;
    max_u = max_v = 0.0;
    sea_level = 0;
    allocate();
  }

  virtual ~ShallowStressBalance() {}

  //  initialization and I/O:

  virtual PetscErrorCode init(PISMVars &vars)
  { variables = &vars; return 0; }

  virtual PetscErrorCode set_boundary_conditions(IceModelVec2Int &locations,
                                                 IceModelVec2V &velocities)
  {
    vel_bc = &velocities;
    bc_locations = &locations;
    return 0;
  }

  //! \brief Set the sea level used to check for floatation. (Units: meters,
  //! relative to the geoid.)
  void set_sea_level_elevation(PetscReal new_sea_level)
  { sea_level = new_sea_level; }

  // interface to the data provided by the stress balance object:

  //! \brief Get the thickness-advective (SSA) 2D velocity.
  virtual PetscErrorCode get_advective_2d_velocity(IceModelVec2V* &result)
  { result = &velocity; return 0; }

  //! \brief Get the max advective velocity (for the adaptive mass-continuity time-stepping).
  virtual PetscErrorCode get_max_2d_velocity(PetscReal &u_max, PetscReal &v_max)
  { u_max = max_u; v_max = max_v; return 0; }

  //! \brief Get the basal frictional heating (for the adaptive energy time-stepping).
  virtual PetscErrorCode get_basal_frictional_heating(IceModelVec2S* &result)
  { result = &basal_frictional_heating; return 0; }

  virtual PetscErrorCode get_D2(IceModelVec2S* &result)
  { result = &D2; return 0; }

  virtual PetscErrorCode compute_principal_strain_rates(IceModelVec2 &/*result*/)
  { SETERRQ(grid.com, 1,"not implemented in base class"); }

  virtual PetscErrorCode compute_2D_stresses(IceModelVec2 &/*result*/)
  { SETERRQ(grid.com,1,"not implemented in base class"); }


  // helpers:

  //! \brief Extends the computational grid (vertically).
  virtual PetscErrorCode extend_the_grid(PetscInt /*old_Mz*/)
  { return 0; }
  //! \brief Produce a report string for the standard output.
  virtual PetscErrorCode stdout_report(string &result)
  { result = ""; return 0; }

  IceFlowLaw* get_flow_law()
  { return flow_law; }
protected:
  virtual PetscErrorCode allocate();

  PetscReal sea_level;
  PISMVars *variables;
  IceBasalResistancePlasticLaw &basal;
  IceFlowLaw *flow_law;
  EnthalpyConverter &EC;

  IceModelVec2V velocity, *vel_bc;
  IceModelVec2Int *bc_locations;
  IceModelVec2S basal_frictional_heating, D2;
  PetscReal max_u, max_v;
};


//! Returns zero velocity field, zero friction heating, and zero for D^2.
/*!
  This derived class is used in the non-sliding SIA approximation. This
  implementation ignores any basal resistance fields (e.g. yield stress from
  the IceModel or other user of this class).
 */
class SSB_Trivial : public ShallowStressBalance
{
public:
  SSB_Trivial(IceGrid &g, IceBasalResistancePlasticLaw &b,
              EnthalpyConverter &e, const NCConfigVariable &conf)
    : ShallowStressBalance(g, b, e, conf) {}
  virtual ~SSB_Trivial() {}
  virtual PetscErrorCode update(bool fast);

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
};

#endif /* _SHALLOWSTRESSBALANCE_H_ */

