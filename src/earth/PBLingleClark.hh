/* Copyright (C) 2013, 2014, 2015, 2016, 2017 PISM Authors
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

#include "PISMBedDef.hh"

namespace pism {
namespace bed {

class BedDeformLC;

//! A wrapper class around BedDeformLC.
class PBLingleClark : public BedDef {
public:
  PBLingleClark(IceGrid::ConstPtr g);
  virtual ~PBLingleClark();

  void uplift_problem(const IceModelVec2S &ice_thickness,
                      const IceModelVec2S &bed_uplift);

protected:
  MaxTimestep max_timestep_impl(double t) const;
  void init_impl();
  void init_with_inputs_impl(const IceModelVec2S &bed_elevation,
                             const IceModelVec2S &bed_uplift,
                             const IceModelVec2S &ice_thickness);
  void update_with_thickness_impl(const IceModelVec2S &ice_thickness,
                                  double my_t, double my_dt);

  void allocate();

  //! Temporary storage
  IceModelVec2S m_work;
  //! Ice thickness at the beginning of the run.
  IceModelVec2S m_H_start;
  //! Bed elevation at the beginning of the run.
  IceModelVec2S m_topg_start;

  //! Storage on rank zero. Used to pass the load to the serial deformation model and get plate
  //! displacement back.
  petsc::Vec::Ptr m_work_0;
  BedDeformLC *m_bdLC;

  //! extended grid for the viscous plate displacement
  IceGrid::Ptr m_extended_grid;
  //! viscous plate displacement on the extended grid (part of the model state)
  IceModelVec2S m_plate_displacement;
};

} // end of namespace bed
} // end of namespace pism

#endif /* _PBLINGLECLARK_H_ */
