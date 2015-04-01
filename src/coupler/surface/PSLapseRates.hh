// Copyright (C) 2011, 2012, 2013, 2014, 2015 PISM Authors
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

#ifndef _PSLAPSERATES_H_
#define _PSLAPSERATES_H_

#include "coupler/util/PLapseRates.hh"
#include "coupler/PISMSurface.hh"
#include "PSModifier.hh"

namespace pism {
namespace surface {

class LapseRates : public PLapseRates<SurfaceModel,SurfaceModifier>
{
public:
  LapseRates(const IceGrid &g, SurfaceModel* in);
  virtual ~LapseRates();
protected:
  virtual void init_impl();

  virtual void ice_surface_mass_flux_impl(IceModelVec2S &result);
  virtual void ice_surface_temperature_impl(IceModelVec2S &result);

  virtual void write_variables_impl(const std::set<std::string> &vars, const PIO &nc);
  virtual void add_vars_to_output_impl(const std::string &keyword, std::set<std::string> &result);
  virtual void define_variables_impl(const std::set<std::string> &vars,
                                     const PIO &nc, IO_Type nctype);
protected:
  double m_smb_lapse_rate;
  SpatialVariableMetadata m_climatic_mass_balance, m_ice_surface_temp;
};

} // end of namespace surface
} // end of namespace pism

#endif /* _PSLAPSERATES_H_ */
