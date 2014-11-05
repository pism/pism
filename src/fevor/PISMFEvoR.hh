/* Copyright (C) 2014 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "PISMComponent.hh"
#include "iceModelVec.hh"
#include <fevor_distribution.hh>

namespace pism {

class EnthalpyConverter;
class StressBalance;
class PSB_pressure;
class PSB_tauxz;
class PSB_tauyz;

/*! PISM-side wrapper around the FEvoR code. Provides the
 *  spatially-variable enhancement factor field.
 *
 * Terminology:
 *
 *   particles exist in PISM and contain one or more distributions of
 *   crystals that are tracked through time. They essentially are
 *   infinitesimally small. Distributions exist in FEvoR and contain
 *   sets of independent crystals (or in the case of NNI weakly
 *   dependent). In PISM-FEvoR you will likely not need to access the
 *   crystals directly. Methods should be provided through FEvoR's
 *   distribution class FEvoR::Distribution.
 */
class PISMFEvoR : public Component_TS {
public:
  PISMFEvoR(IceGrid &g, const Config &conf,
            EnthalpyConverter *EC, StressBalance *stress_balance);
  virtual ~PISMFEvoR();

  PetscErrorCode init(Vars &vars);

  virtual PetscErrorCode max_timestep(double t, double &dt, bool &restrict);
  virtual PetscErrorCode update(double t, double dt);

  virtual void add_vars_to_output(const std::string &keyword, std::set<std::string> &result);

  virtual PetscErrorCode define_variables(const std::set<std::string> &vars, const PIO &nc,
                                          IO_Type nctype);

  virtual PetscErrorCode write_variables(const std::set<std::string> &vars, const PIO& nc);
private:
  StressBalance *m_stress_balance;
  EnthalpyConverter *m_EC;

  std::vector<unsigned int> m_packing_dimensions;
  /* An isotropic distribution for calculating enhancement factor. The
   * enhancement factor is defined as ratio of ice's response relative
   * to isotropic ice. Since we need isotropic ice's response to any
   * input stress, this is the easiest way to provide it but it may be
   * the most computationally heavy. Possible efficiency improvement
   * here.
   */
  FEvoR::Distribution m_d_iso;
  std::vector<FEvoR::Distribution> m_distributions;
  std::vector<double> m_p_x, m_p_y, m_p_z, m_p_e;

  // Diagnostics -- total number of recrystallization events in time step
  std::vector<unsigned int> m_n_migration_recrystallizations,
    m_n_polygonization_recrystallizations;

  PetscErrorCode allocate();

  PetscErrorCode set_initial_distribution_parameters();

  PetscErrorCode load_distributions(const std::string &input_file);
  PetscErrorCode save_distributions(const PIO &nc);

  PetscErrorCode save_diagnostics(const PIO &nc);

  PetscErrorCode update_particle_position(double &x, double &y, double &z,
                                          double u, double v, double w,
                                          double m_dt);

  PetscErrorCode evaluate_at_point(IceModelVec3 &input,
                                   double x, double y, double z,
                                   double &result);

  PetscErrorCode pointcloud_to_grid(const std::vector<double> &x,
                                    const std::vector<double> &z,
                                    const std::vector<double> &values,
                                    IceModelVec3 &result);

  IceModelVec3 m_enhancement_factor;
  IceModelVec3 *m_enthalpy;
  PSB_pressure *m_pressure;
  PSB_tauxz *m_tauxz;
  PSB_tauyz *m_tauyz;
};

} // end of namespace pism
