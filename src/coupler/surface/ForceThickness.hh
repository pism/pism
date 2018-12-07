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

#ifndef _PSFORCETHICKNESS_H_
#define _PSFORCETHICKNESS_H_

#include "pism/coupler/SurfaceModel.hh"

namespace pism {

class IceModelVec2CellType;

namespace surface {

//! A class implementing a modified surface mass balance which forces
//! ice thickness to a given target by the end of the run.
class ForceThickness : public SurfaceModel {
public:
  ForceThickness(IceGrid::ConstPtr g, std::shared_ptr<SurfaceModel> input);
  virtual ~ForceThickness();
protected:
  virtual void init_impl(const Geometry &geometry);
  virtual void update_impl(const Geometry &geometry, double t, double dt);

  virtual void define_model_state_impl(const PIO &output) const;
  virtual void write_model_state_impl(const PIO &output) const;

  virtual const IceModelVec2S& mass_flux_impl() const;

  virtual MaxTimestep max_timestep_impl(double t) const;
private:
  void adjust_mass_flux(double time,
                        const IceModelVec2S &ice_thickness,
                        const IceModelVec2CellType &cell_type,
                        IceModelVec2S &result) const;

  double m_alpha, m_alpha_ice_free_factor,  m_ice_free_thickness_threshold;
  double m_start_time;
  IceModelVec2S m_target_thickness;
  IceModelVec2Int m_ftt_mask;

  IceModelVec2S::Ptr m_mass_flux;
};

} // end of namespace surface
} // end of namespace pism

#endif /* _PSFORCETHICKNESS_H_ */
