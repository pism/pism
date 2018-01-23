/* Copyright (C) 2016, 2017, 2018 PISM Authors
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

#include "pism/util/Component.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/IceModelVec2CellType.hh"


namespace pism {

//! An abstract class implementing calving front retreat resulting from application of a
//! spatially-variable horizontal retreat rate.
/*! The retreat rate may correspond to frontal melting or calving. Requires the "part_grid"
    mechanism.
 */
class CalvingFrontRetreat : public Component {
public:
  CalvingFrontRetreat(IceGrid::ConstPtr g, unsigned int mask_stencil_width);
  virtual ~CalvingFrontRetreat();

  void update(double dt,
              const IceModelVec2S &sea_level,
              const IceModelVec2Int &ice_thickness_bc_mask,
              const IceModelVec2S &bed_topography,
              IceModelVec2CellType &pism_mask,
              IceModelVec2S &Href,
              IceModelVec2S &ice_thickness);

  const IceModelVec2S& calving_rate() const;

protected:
  MaxTimestep max_timestep_impl(double t) const ;

  virtual void compute_calving_rate(const IceModelVec2CellType &mask,
                                    IceModelVec2S &result) const = 0;

  mutable IceModelVec2CellType m_mask;
  mutable IceModelVec2S m_tmp;
  IceModelVec2S m_horizontal_calving_rate, m_surface_topography;
  bool m_restrict_timestep;
};

/*! @brief Calving (or frontal melt) rate diagnostic. */
class CalvingRate : public Diag<CalvingFrontRetreat>
{
public:
  CalvingRate(const CalvingFrontRetreat *m,
              const std::string &name,
              const std::string &long_name);
  IceModelVec::Ptr compute_impl() const;
};

} // end of namespace pism


#endif /* CALVINGFRONTRETREAT_H */
