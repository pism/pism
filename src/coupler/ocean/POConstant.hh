// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 PISM Authors
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

#ifndef _POCONSTANT_H_
#define _POCONSTANT_H_

#include "coupler/PISMOcean.hh"

namespace pism {
namespace ocean {
//! \brief A class implementing a constant (in terms of the ocean inputs) ocean
//! model. Uses configuration parameters for the sea level elevation and
//! sub-shelf heat flux.
class Constant : public OceanModel {
public:
  Constant(IceGrid::ConstPtr g);
  virtual ~Constant();
protected:
  virtual MaxTimestep max_timestep_impl(double t) const;
  virtual void update_impl(double my_t, double my_dt);
  virtual void init_impl();
  virtual void sea_level_elevation_impl(double &result) const;
  virtual void shelf_base_temperature_impl(IceModelVec2S &result) const;
  virtual void shelf_base_mass_flux_impl(IceModelVec2S &result) const;
protected:
  double m_meltrate;
};

} // end of namespace ocean
} // end of namespace pism
#endif /* _POCONSTANT_H_ */
