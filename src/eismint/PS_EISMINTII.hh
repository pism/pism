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

#ifndef _PS_EISMINTII_H_
#define _PS_EISMINTII_H_

#include "PISMSurface.hh"
#include "iceModelVec.hh"

namespace pism {

/** EISMINT II climate inputs.
 *
 * This class should be removed together with the pisms executable
 * (once I get to that).
 */
class PS_EISMINTII : public PISMSurfaceModel {
public:
  PS_EISMINTII(IceGrid &g, const PISMConfig &conf, int experiment);
  ~PS_EISMINTII();

  // the interface:
  PetscErrorCode init(PISMVars &vars);
  void attach_atmosphere_model(PISMAtmosphereModel *input);
  PetscErrorCode ice_surface_mass_flux(IceModelVec2S &result);
  PetscErrorCode ice_surface_temperature(IceModelVec2S &result);
  PetscErrorCode update(PetscReal t, PetscReal dt);
  void add_vars_to_output(std::string keyword, std::set<std::string> &result);
  PetscErrorCode define_variables(std::set<std::string> vars, const PIO &nc,
                                  PISM_IO_Type nctype);
  PetscErrorCode write_variables(std::set<std::string> vars, const PIO &nc);
private:
  IceModelVec2S m_climatic_mass_balance, m_ice_surface_temp;
  PetscErrorCode allocate();
  PetscErrorCode initialize_using_formulas();
  int m_experiment;
  double m_M_max, m_R_el, m_S_T, m_S_b, m_T_min;
};


} // end of namespace pism

#endif /* _PS_EISMINTII_H_ */
