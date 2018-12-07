// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#ifndef _IBSURFACEMODEL_H_
#define _IBSURFACEMODEL_H_

#include "pism/util/iceModelVec.hh"
#include "pism/coupler/AtmosphereModel.hh"
#include "pism/coupler/SurfaceModel.hh"

namespace pism {
namespace icebin {

//! \brief A class implementing a constant-in-time surface model for the surface mass balance.
//!
//! Reads data from a PISM input file.
//!
//! Ice surface temperature is parameterized as in PISM-IBSurfaceModel, using a latitude
//! and surface elevation-dependent formula.

class IBSurfaceModel : public pism::surface::SurfaceModel {
public:
  IBSurfaceModel(IceGrid::ConstPtr grid);

protected:
  virtual void init_impl(const Geometry &geometry);
  virtual void update_impl(const Geometry &geometry, double t, double dt);
  virtual MaxTimestep max_timestep_impl(double t) const;

  virtual void define_model_state_impl(const PIO &output) const;
  virtual void write_model_state_impl(const PIO &output) const;

  virtual const IceModelVec2S& mass_flux_impl() const;
  virtual const IceModelVec2S& temperature_impl() const;

protected:
  bool _initialized;
  std::string m_input_file;

public:
  // Inputs from IceBin
  pism::IceModelVec2S icebin_wflux;
  pism::IceModelVec2S icebin_deltah;
  pism::IceModelVec2S icebin_massxfer; // [kg m-2 s-1]
  pism::IceModelVec2S icebin_enthxfer; // [J m-2 s-1]
  // Calculated
  pism::IceModelVec2S surface_temp;
};

} // end of namespace surface
} // end of namespace pism

#endif /* _IBSURFACEMODEL_H_ */
