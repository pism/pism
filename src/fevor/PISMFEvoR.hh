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

namespace pism {

class EnthalpyConverter;
class StressBalance;
class PSB_pressure;
class PSB_tauxz;
class PSB_tauyz;

/*! PISM-side wrapper around the FEvoR code. Provides the
 *  spatially-variable enhancement factor field.
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
  PetscErrorCode allocate();
  
  PetscErrorCode update_particle_position(double &x, double &y, double &z,
                                          const double &u,
                                          const double &v, 
                                          const double &w,
                                          const double &m_dt);
  
  PetscErrorCode evaluate_at_point(IceModelVec3 &input,
                                   double x, double y, double z, 
                                   double &result);

  PetscErrorCode pointcloud_to_grid(const std::vector<double> &x,
                                    const std::vector<double> &z,
                                    const std::vector<double> &values,
                                    IceModelVec3 &result);

  IceModelVec3 m_enhancement_factor;
  IceModelVec3 *m_enthalpy;
  StressBalance *m_stress_balance;
  PSB_pressure *m_pressure;
  PSB_tauxz *m_tauxz;
  PSB_tauyz *m_tauyz;
  EnthalpyConverter *m_EC;
};

} // end of namespace pism
