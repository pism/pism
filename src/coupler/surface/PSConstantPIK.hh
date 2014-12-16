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

#ifndef _PSCONSTANTPIK_H_
#define _PSCONSTANTPIK_H_

#include "PISMSurface.hh"
#include "iceModelVec.hh"
#include "PISMAtmosphere.hh"

namespace pism {

//! \brief A class implementing a constant-in-time surface model for the surface mass balance.
//!
//! Reads data from a PISM input file.
//!
//! Ice surface temperature is parameterized as in PISM-PIK, using a latitude
//! and surface elevation-dependent formula.

class PSConstantPIK : public SurfaceModel {
public:
  PSConstantPIK(IceGrid &g);

  virtual void init(Vars &vars);

  virtual void attach_atmosphere_model(AtmosphereModel *input);

  virtual void get_diagnostics(std::map<std::string, Diagnostic*> &dict,
                               std::map<std::string, TSDiagnostic*> &ts_dict);
  virtual void update(double my_t, double my_dt);
  virtual void ice_surface_mass_flux(IceModelVec2S &result);
  virtual void ice_surface_temperature(IceModelVec2S &result);
  virtual void define_variables(const std::set<std::string> &vars, const PIO &nc, IO_Type nctype);
  virtual void write_variables(const std::set<std::string> &vars, const PIO &nc);
  virtual void add_vars_to_output(const std::string &keyword, std::set<std::string> &result);
protected:
  std::string input_file;
  IceModelVec2S climatic_mass_balance, ice_surface_temp;
  IceModelVec2S *lat, *usurf;
};

} // end of namespace pism

#endif /* _PSCONSTANTPIK_H_ */
