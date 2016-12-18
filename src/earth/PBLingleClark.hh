/* Copyright (C) 2013, 2014, 2015, 2016 PISM Authors
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

#include <fftw3.h>

#include "PISMBedDef.hh"
#include "deformation.hh"

namespace pism {
namespace bed {

//! A wrapper class around BedDeformLC.
class PBLingleClark : public BedDef {
public:
  PBLingleClark(IceGrid::ConstPtr g);
  virtual ~PBLingleClark();

protected:
  MaxTimestep max_timestep_impl(double t) const;
  void init_impl();
  void init_with_inputs_impl(const IceModelVec2S &bed_elevation,
                             const IceModelVec2S &bed_uplift,
                             const IceModelVec2S &ice_thickness);
  void update_with_thickness_impl(const IceModelVec2S &ice_thickness,
                                  double my_t, double my_dt);
  void correct_topg();
  void allocate();

  // Vecs on processor 0:
  //! ice thickness
  petsc::Vec::Ptr m_Hp0;
  //! bed elevation
  petsc::Vec::Ptr m_bedp0;
  //! initial (start-of-the-run) thickness
  petsc::Vec::Ptr m_Hstartp0;
  //! initial bed elevation
  petsc::Vec::Ptr m_bedstartp0;
  //! bed uplift
  petsc::Vec::Ptr m_upliftp0;
  BedDeformLC *m_bdLC;
};

} // end of namespace bed
} // end of namespace pism

#endif /* _PBLINGLECLARK_H_ */
