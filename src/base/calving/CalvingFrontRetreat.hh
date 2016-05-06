/* Copyright (C) 2016 PISM Authors
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

#ifndef CALVINGFRONTRETREAT_H
#define CALVINGFRONTRETREAT_H

#include "base/util/PISMComponent.hh"
#include "base/util/iceModelVec.hh"


namespace pism {

class IceModelVec2CellType;

//! An abstract class implementing calving front retreat resulting from application of a
//! spatially-variable horizontal retreat rate.
/*! The retreat rate may correspond to frontal melting or calving. Requires the "part_grid"
    mechanism.
 */
class CalvingFrontRetreat : public Component {
public:
  CalvingFrontRetreat(IceGrid::ConstPtr g);
  virtual ~CalvingFrontRetreat();

  void update(double dt,
              double sea_level,
              const IceModelVec2S &bed_topography,
              IceModelVec2CellType &pism_mask,
              IceModelVec2S &Href,
              IceModelVec2S &ice_thickness);

  MaxTimestep max_timestep();

protected:
  virtual void compute_calving_rate(const IceModelVec2CellType &mask,
                                    IceModelVec2S &result) = 0;

  IceModelVec2S m_thk_loss, m_horizontal_calving_rate;
  bool m_restrict_timestep;
};

} // end of namespace pism


#endif /* CALVINGFRONTRETREAT_H */
