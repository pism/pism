/* Copyright (C) 2013, 2014 PISM Authors
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
#ifndef _PISMEIGENCALVING_H_
#define _PISMEIGENCALVING_H_

#include "iceModelVec.hh"
#include "PISMComponent.hh"

class PISMStressBalance;

class PISMEigenCalving : public PISMComponent
{
public:
  PISMEigenCalving(IceGrid &g, const PISMConfig &conf,
                   PISMStressBalance *stress_balance);
  virtual ~PISMEigenCalving();

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(double dt,
                                IceModelVec2Int &pism_mask,
                                IceModelVec2S &Href,
                                IceModelVec2S &ice_thickness);

  virtual PetscErrorCode max_timestep(double my_t, double &my_dt, bool &restrict);

  // empty methods that we're required to implement:
  virtual void add_vars_to_output(std::string keyword, std::set<std::string> &result);
  virtual PetscErrorCode define_variables(std::set<std::string> vars, const PIO &nc,
                                          PISM_IO_Type nctype);
  virtual PetscErrorCode write_variables(std::set<std::string> vars, const PIO& nc);
protected:
  IceModelVec2 m_strain_rates;
  IceModelVec2S m_thk_loss;
  const int m_stencil_width;
  IceModelVec2Int *m_mask;
  PISMStressBalance *m_stress_balance;
  double m_K;
  bool m_restrict_timestep;

  PetscErrorCode update_strain_rates();
  PetscErrorCode remove_narrow_tongues(IceModelVec2Int &pism_mask,
                                       IceModelVec2S &ice_thickness);
};


#endif /* _PISMEIGENCALVING_H_ */
