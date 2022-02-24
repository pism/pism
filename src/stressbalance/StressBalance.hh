// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2021, 2022 Constantine Khroulev and Ed Bueler
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

#ifndef _PISMSTRESSBALANCE_H_
#define _PISMSTRESSBALANCE_H_

#include <memory>               // std::shared_ptr

#include "pism/util/Component.hh"     // derives from Component
#include "pism/util/IceModelVec3.hh"
#include "pism/stressbalance/timestepping.hh"

namespace pism {

class IceModelVec2S;
class IceModelVec2Stag;
class Geometry;

namespace rheology {
class FlowLaw;
} // end of namespace rheology

//! Stress balance models and related diagnostics.
namespace stressbalance {

class ShallowStressBalance;
class SSB_Modifier;

class Inputs {
public:
  Inputs();

  const Geometry *geometry;
  bool new_bed_elevation;

  const IceModelVec2S *basal_melt_rate;
  const IceModelVec2S *basal_yield_stress;
  const IceModelVec2S *water_column_pressure;
  const IceModelVec2S *fracture_density;

  const IceModelVec3  *enthalpy;
  const IceModelVec3  *age;

  const IceModelVec2S *bc_mask;
  const IceModelVec2V *bc_values;

  // inputs used by regional stress balance models
  const IceModelVec2S *no_model_mask;
  const IceModelVec2S *no_model_ice_thickness;
  const IceModelVec2S *no_model_surface_elevation;

  void dump(const char *filename) const;
};

//! The class defining PISM's interface to the shallow stress balance code.
/*!
  Generally all the nontrivial fields are updated by a call to update().  The rest
  of the methods generally provide access to precomputed results.  The following
  diagram shows where these results are generally used in the rest of PISM.  (It
  does not show the call graph, as would doxygen.)

  \image html stressbalance-out.png "\b Methods of StressBalance, and the uses of their results.  Dotted edges show scalars and dashed edges show fields.  Dashed boxes inside the StressBalance object are important methods which may be present in shallow cases.  The age time step has inputs which are a strict subset of the inputs of the energy time step."

  this command fails: \dotfile stressbalance-out.dot
*/
class StressBalance : public Component
{
public:
  StressBalance(IceGrid::ConstPtr g,
                std::shared_ptr<ShallowStressBalance> sb,
                std::shared_ptr<SSB_Modifier> ssb_mod);
  virtual ~StressBalance();

  //! \brief Initialize the StressBalance object.
  void init();

  //! \brief Update all the fields if (full_update), only update diffusive flux
  //! and max. diffusivity otherwise.
  void update(const Inputs &inputs, bool full_update);

  //! \brief Get the thickness-advective (SSA) 2D velocity.
  const IceModelVec2V& advective_velocity() const;

  //! \brief Get the diffusive (SIA) vertically-averaged flux on the staggered grid.
  const IceModelVec2Stag& diffusive_flux() const;

  //! \brief Get the max diffusivity (for the adaptive time-stepping).
  double max_diffusivity() const;

  CFLData max_timestep_cfl_2d() const;
  CFLData max_timestep_cfl_3d() const;

  // for the energy/age time step:

  //! \brief Get components of the the 3D velocity field.
  const IceModelVec3& velocity_u() const;
  const IceModelVec3& velocity_v() const;
  const IceModelVec3& velocity_w() const;

  //! \brief Get the basal frictional heating.
  const IceModelVec2S& basal_frictional_heating() const;

  const IceModelVec3& volumetric_strain_heating() const;

  //! \brief Produce a report string for the standard output.
  std::string stdout_report() const;

  //! \brief Returns a pointer to a shallow stress balance solver implementation.
  const ShallowStressBalance* shallow() const;

  //! \brief Returns a pointer to a stress balance modifier implementation.
  const SSB_Modifier* modifier() const;
protected:
  virtual DiagnosticList diagnostics_impl() const;
  virtual TSDiagnosticList ts_diagnostics_impl() const;

  virtual void define_model_state_impl(const File &output) const;
  virtual void write_model_state_impl(const File &output) const;

  virtual void compute_vertical_velocity(const CellTypeArray1 &mask,
                                         const IceModelVec3 &u,
                                         const IceModelVec3 &v,
                                         const IceModelVec2S *bmr,
                                         IceModelVec3 &result);
  virtual void compute_volumetric_strain_heating(const Inputs &inputs);

  CFLData m_cfl_2d, m_cfl_3d;

  IceModelVec3 m_w, m_strain_heating;

  std::shared_ptr<ShallowStressBalance> m_shallow_stress_balance;
  std::shared_ptr<SSB_Modifier> m_modifier;
};

std::shared_ptr<StressBalance> create(const std::string &model_name,
                                      IceGrid::ConstPtr grid,
                                      bool regional);

void compute_2D_principal_strain_rates(const IceModelVec2V &velocity,
                                       const CellTypeArray1 &mask,
                                       IceModelVec3 &result);

void compute_2D_stresses(const rheology::FlowLaw &flow_law,
                         const IceModelVec2V &velocity,
                         const IceModelVec2S &hardness,
                         const CellTypeArray1 &cell_type,
                         IceModelVec3 &result);

} // end of namespace stressbalance
} // end of namespace pism

#endif /* _PISMSTRESSBALANCE_H_ */
