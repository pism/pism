/* Copyright (C) 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#ifndef _PBLINGLECLARK_H_
#define _PBLINGLECLARK_H_

#include <memory>               // std::unique_ptr

#include "BedDef.hh"

namespace pism {
namespace bed {

class LingleClarkSerial;

//! A wrapper class around LingleClarkSerial.
class LingleClark : public BedDef {
public:
  LingleClark(IceGrid::ConstPtr g);
  virtual ~LingleClark();

  const IceModelVec2S& total_displacement() const;

  const IceModelVec2S& viscous_displacement() const;

  const IceModelVec2S& relief() const;

  void step(const IceModelVec2S &ice_thickness,
            const IceModelVec2S &sea_level_elevation,
            double dt);

protected:
  virtual void define_model_state_impl(const PIO &output) const;
  virtual void write_model_state_impl(const PIO &output) const;

  DiagnosticList diagnostics_impl() const;

  MaxTimestep max_timestep_impl(double t) const;
  void init_impl(const InputOptions &opts, const IceModelVec2S &ice_thickness,
                 const IceModelVec2S &sea_level_elevation);
  void bootstrap_impl(const IceModelVec2S &bed_elevation,
                      const IceModelVec2S &bed_uplift,
                      const IceModelVec2S &ice_thickness,
                      const IceModelVec2S &sea_level_elevation);
  void update_impl(const IceModelVec2S &ice_thickness,
                   const IceModelVec2S &sea_level_elevation,
                   double t, double dt);

  //! Total (viscous and elastic) bed displacement.
  IceModelVec2S m_bed_displacement;

  //! Storage on rank zero. Used to pass the load to the serial deformation model and get
  //! bed displacement back.
  petsc::Vec::Ptr m_work0;

  //! Bed relief relative to the bed displacement.
  IceModelVec2S m_relief;

  //! Ice-equivalent load thickness.
  IceModelVec2S m_load_thickness;

  //! Serial viscoelastic bed deformation model.
  std::unique_ptr<LingleClarkSerial> m_serial_model;

  //! extended grid for the viscous plate displacement
  IceGrid::Ptr m_extended_grid;

  //! Viscous displacement on the extended grid (part of the model state).
  IceModelVec2S m_viscous_bed_displacement;
  //! rank 0 storage using the extended grid
  petsc::Vec::Ptr m_viscous_bed_displacement0;
};

} // end of namespace bed
} // end of namespace pism

#endif /* _PBLINGLECLARK_H_ */
