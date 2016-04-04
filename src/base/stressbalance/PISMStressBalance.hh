// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016 Constantine Khroulev and Ed Bueler
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

#include "base/util/PISMComponent.hh"     // derives from Component
#include "base/util/iceModelVec.hh"

namespace pism {

class IceModelVec2CellType;

//! Stress balance models and related diagnostics.
namespace stressbalance {

class ShallowStressBalance;
class SSB_Modifier;

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
  StressBalance(IceGrid::ConstPtr g, ShallowStressBalance *sb, SSB_Modifier *ssb_mod);
  virtual ~StressBalance();

  //! \brief Initialize the StressBalance object.
  void init();

  //! \brief Set the vertically-averaged ice velocity boundary condition.
  /*!
   * Does not affect the SIA computation.
   */
  void set_boundary_conditions(const IceModelVec2Int &locations,
                               const IceModelVec2V &velocities);

  void set_basal_melt_rate(const IceModelVec2S &bmr);

  //! \brief Update all the fields if (not fast), only update diffusive flux
  //! and max. diffusivity otherwise.
  void update(bool fast, double sea_level,
              const IceModelVec2S &melange_back_pressure);

  //! \brief Get the thickness-advective (SSA) 2D velocity.
  const IceModelVec2V& advective_velocity();

  //! \brief Get the diffusive (SIA) vertically-averaged flux on the staggered grid.
  const IceModelVec2Stag& diffusive_flux();

  //! \brief Get the max diffusivity (for the adaptive time-stepping).
  double max_diffusivity();

  // for the energy/age time step:

  //! \brief Get the 3D velocity (for the energy/age time-stepping).
  const IceModelVec3& velocity_u();
  const IceModelVec3& velocity_v();
  const IceModelVec3& velocity_w();

  //! \brief Get the basal frictional heating (for the energy time-stepping).
  const IceModelVec2S& basal_frictional_heating();

  const IceModelVec3& volumetric_strain_heating();

  // for the calving, etc.:

  //! \brief Get the largest and smallest eigenvalues of the strain rate tensor.
  void compute_2D_principal_strain_rates(const IceModelVec2V &velocity,
                                         const IceModelVec2CellType &mask,
                                         IceModelVec2 &result);

  //! \brief Get the components of the 2D deviatoric stress tensor.
  void compute_2D_stresses(const IceModelVec2V &velocity,
                           const IceModelVec2CellType &mask,
                           IceModelVec2 &result);

  //! \brief Produce a report string for the standard output.
  std::string stdout_report();

  //! \brief Returns a pointer to a stress balance solver implementation.
  ShallowStressBalance* get_stressbalance();

  //! \brief Returns a pointer to a stress balance modifier implementation.
  SSB_Modifier* get_ssb_modifier();
protected:
  virtual void get_diagnostics_impl(std::map<std::string, Diagnostic::Ptr> &dict,
                                    std::map<std::string, TSDiagnostic::Ptr> &ts_dict);

  virtual void write_variables_impl(const std::set<std::string> &vars, const PIO &nc);
  virtual void add_vars_to_output_impl(const std::string &keyword, std::set<std::string> &result);
  virtual void define_variables_impl(const std::set<std::string> &vars, const PIO &nc,
                                     IO_Type nctype);

  virtual void compute_vertical_velocity(const IceModelVec3 &u,
                                         const IceModelVec3 &v,
                                         const IceModelVec2S *bmr,
                                         IceModelVec3 &result);
  virtual void compute_volumetric_strain_heating();

  IceModelVec3 m_w, m_strain_heating;
  const IceModelVec2S *m_basal_melt_rate;

  ShallowStressBalance *m_shallow_stress_balance;
  SSB_Modifier *m_modifier;
};

} // end of namespace stressbalance
} // end of namespace pism

#endif /* _PISMSTRESSBALANCE_H_ */

