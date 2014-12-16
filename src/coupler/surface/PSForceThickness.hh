// Copyright (C) 2011, 2012, 2013, 2014 PISM Authors
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

#ifndef _PSFORCETHICKNESS_H_
#define _PSFORCETHICKNESS_H_

#include "PSModifier.hh"
#include "iceModelVec.hh"
#include "NCVariable.hh"

namespace pism {

//! A class implementing a modified surface mass balance which forces
//! ice thickness to a given target by the end of the run.
class PSForceThickness : public PSModifier {
public:
  PSForceThickness(IceGrid &g, SurfaceModel *input);

  virtual ~PSForceThickness();
  virtual void init(Vars &vars);
  virtual void attach_atmosphere_model(AtmosphereModel *input);
  virtual void ice_surface_mass_flux(IceModelVec2S &result);
  virtual void ice_surface_temperature(IceModelVec2S &result);
  virtual void max_timestep(double my_t, double &my_dt, bool &restrict);
  virtual void add_vars_to_output(const std::string &keyword, std::set<std::string> &result);
  virtual void define_variables(const std::set<std::string> &vars, const PIO &nc, IO_Type nctype);
  virtual void write_variables(const std::set<std::string> &vars, const PIO &nc);
private:
  std::string m_input_file;
  double m_alpha, m_alpha_ice_free_factor,  m_ice_free_thickness_threshold;
  IceModelVec2S *m_ice_thickness; //!< current ice thickness produced by IceModel.
  IceModelVec2S m_target_thickness, m_ftt_mask;
  IceModelVec2Int *m_pism_mask;
  NCSpatialVariable m_climatic_mass_balance, m_climatic_mass_balance_original, m_ice_surface_temp;
};

} // end of namespace pism

#endif /* _PSFORCETHICKNESS_H_ */
