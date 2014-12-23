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

#ifndef _PSANOMALY_H_
#define _PSANOMALY_H_

#include "PGivenClimate.hh"
#include "PSModifier.hh"

namespace pism {

//! \brief Reads and uses climatic_mass_balance and ice_surface_temp \b anomalies from a file.
class PSAnomaly : public PGivenClimate<PSModifier,SurfaceModel>
{
public:
  PSAnomaly(const IceGrid &g, SurfaceModel* in);
  virtual ~PSAnomaly();

  virtual void init();
  virtual void update(double my_t, double my_dt);

  virtual void ice_surface_mass_flux(IceModelVec2S &result);
  virtual void ice_surface_temperature(IceModelVec2S &result);

  virtual void define_variables(const std::set<std::string> &vars, const PIO &nc, IO_Type nctype);
  virtual void write_variables(const std::set<std::string> &vars, const PIO &nc);
  virtual void add_vars_to_output(const std::string &keyword, std::set<std::string> &result);
protected:
  NCSpatialVariable climatic_mass_balance, ice_surface_temp;
  IceModelVec2T *climatic_mass_balance_anomaly, *ice_surface_temp_anomaly;
};

} // end of namespace pism

#endif /* _PSANOMALY_H_ */
