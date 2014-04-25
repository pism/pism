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
#include "PISMConfig.hh"

namespace pism {

class Vars;
class IceFlowLaw;
class EnthalpyConverter;
class IceBasalResistancePlasticLaw;

//! Shallow stress balance (such as the SSA).
class ShallowStressBalance : public Component
{
public:
  ShallowStressBalance(IceGrid &g, EnthalpyConverter &e, const Config &conf);
  virtual ~ShallowStressBalance();

  //  initialization and I/O:

  virtual PetscErrorCode init(Vars &vars)
  { variables = &vars; return 0; }

  virtual PetscErrorCode set_boundary_conditions(IceModelVec2Int &locations,
                                                 IceModelVec2V &velocities)
  {
    m_vel_bc = &velocities;
    bc_locations = &locations;
    return 0;
  }

  //! \brief Set the sea level used to check for floatation. (Units: meters,
  //! relative to the geoid.)
  void set_sea_level_elevation(double new_sea_level)
  { sea_level = new_sea_level; }

  virtual PetscErrorCode update(bool fast,
                                IceModelVec2S &melange_back_pressure) = 0;

  // interface to the data provided by the stress balance object:
  virtual void get_diagnostics(std::map<std::string, Diagnostic*> &dict,
                               std::map<std::string, TSDiagnostic*> &/*ts_dict*/);

  //! \brief Get the thickness-advective (SSA) 2D velocity.
  virtual PetscErrorCode get_2D_advective_velocity(IceModelVec2V* &result)
  { result = &m_velocity; return 0; }

  //! \brief Get the basal frictional heating (for the adaptive energy time-stepping).
  virtual PetscErrorCode get_basal_frictional_heating(IceModelVec2S* &result)
  { result = &basal_frictional_heating; return 0; }

  virtual PetscErrorCode compute_2D_principal_strain_rates(IceModelVec2V &velocity,
                                                           IceModelVec2Int &mask,
                                                           IceModelVec2 &result);

  virtual PetscErrorCode compute_2D_stresses(IceModelVec2V &velocity, IceModelVec2Int &mask,
                                             IceModelVec2 &result);

  virtual PetscErrorCode compute_basal_frictional_heating(IceModelVec2V &velocity,
                                                          IceModelVec2S &tauc,
                                                          IceModelVec2Int &mask,
                                                          IceModelVec2S &result);
  // helpers:

  //! \brief Extends the computational grid (vertically).
  virtual PetscErrorCode extend_the_grid(int /*old_Mz*/)
  { return 0; }
  //! \brief Produce a report string for the standard output.
  virtual PetscErrorCode stdout_report(std::string &result)
  { result = ""; return 0; }

  const IceFlowLaw* get_flow_law()
  { return flow_law; }

  EnthalpyConverter& get_enthalpy_converter()
  { return EC; }

  const IceBasalResistancePlasticLaw* get_sliding_law()
  { return basal_sliding_law; }
protected:
  virtual PetscErrorCode allocate();

  double sea_level;
  Vars *variables;
  IceBasalResistancePlasticLaw *basal_sliding_law;
  IceFlowLaw *flow_law;
  EnthalpyConverter &EC;

  IceModelVec2V m_velocity, *m_vel_bc;
  IceModelVec2Int *bc_locations;
  IceModelVec2S basal_frictional_heating;
};

class SSB_beta : public Diag<ShallowStressBalance>
{
public:
  SSB_beta(ShallowStressBalance *m, IceGrid &g, Vars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes the gravitational driving stress (diagnostically).
class SSB_taud : public Diag<ShallowStressBalance>
{
public:
  SSB_taud(ShallowStressBalance *m, IceGrid &g, Vars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes the magnitude of the gravitational driving stress
//! (diagnostically).
class SSB_taud_mag : public Diag<ShallowStressBalance>
{
public:
  SSB_taud_mag(ShallowStressBalance *m, IceGrid &g, Vars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! @brief Computes the basal shear stress @f$ \tau_b @f$.
class SSB_taub : public Diag<ShallowStressBalance>
{
public:
  SSB_taub(ShallowStressBalance *m, IceGrid &g, Vars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes the magnitude of the basal shear stress
//! (diagnostically).
class SSB_taub_mag : public Diag<ShallowStressBalance>
{
public:
  SSB_taub_mag(ShallowStressBalance *m, IceGrid &g, Vars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
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
  ZeroSliding(IceGrid &g, EnthalpyConverter &e, const Config &conf);
  virtual ~ZeroSliding();
  
  virtual PetscErrorCode update(bool fast, IceModelVec2S &melange_back_pressure);

  virtual void add_vars_to_output(std::string /*keyword*/, std::set<std::string> &/*result*/);

  //! Defines requested couplings fields and/or asks an attached model
  //! to do so.
  virtual PetscErrorCode define_variables(std::set<std::string> /*vars*/, const PIO &/*nc*/,
                                          IO_Type /*nctype*/);

  //! Writes requested couplings fields to file and/or asks an attached
  //! model to do so.
  virtual PetscErrorCode write_variables(std::set<std::string> /*vars*/, const PIO &/*nc*/);
};

class PrescribedSliding : public ZeroSliding {
public:
  PrescribedSliding(IceGrid &g, EnthalpyConverter &e, const Config &conf);
  virtual ~PrescribedSliding();
  virtual PetscErrorCode update(bool fast, IceModelVec2S &melange_back_pressure);
  virtual PetscErrorCode init(Vars &vars);
};

} // end of namespace pism

#endif /* _SHALLOWSTRESSBALANCE_H_ */
