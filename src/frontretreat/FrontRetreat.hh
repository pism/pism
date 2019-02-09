/* Copyright (C) 2016, 2017, 2018, 2019 PISM Authors
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

#ifndef FRONTRETREAT_H
#define FRONTRETREAT_H

#include "pism/util/Component.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/IceModelVec2CellType.hh"
#include "pism/util/MaxTimestep.hh"

namespace pism {

//! An abstract class implementing calving front retreat resulting from application of a
//! spatially-variable horizontal retreat rate.
/*! The retreat rate may correspond to frontal melting or calving. Requires the "part_grid"
    mechanism.
 */
class FrontRetreat : public Component {
public:
  FrontRetreat(IceGrid::ConstPtr g, unsigned int mask_stencil_width);
  virtual ~FrontRetreat();

  const IceModelVec2S& retreat_rate() const;

protected:
  void update_geometry(double dt,
                       const IceModelVec2S &sea_level,
                       const IceModelVec2S &bed_topography,
                       const IceModelVec2Int &bc_mask,
                       const IceModelVec2S &horizontal_retreat_rate,
                       IceModelVec2CellType &cell_type,
                       IceModelVec2S &Href,
                       IceModelVec2S &ice_thickness);

  void prepare_mask(const IceModelVec2CellType &input,
                    IceModelVec2CellType &output) const;

  /*!
   * Combines information about maximum time step length computed using given front
   * retreat rate.
   */
  struct Timestep {
    MaxTimestep dt;
    double rate_max;
    double rate_mean;
    int N_cells;
  };

  Timestep max_timestep(const IceModelVec2S &horizontal_retreat_rate) const;

  mutable IceModelVec2CellType m_mask;
  mutable IceModelVec2S m_tmp;
  IceModelVec2S m_horizontal_retreat_rate, m_surface_topography;
  bool m_restrict_timestep;
};

/*! @brief Retreat rate due to calving (or frontal melt). */
class FrontRetreatRate : public Diag<FrontRetreat>
{
public:
  FrontRetreatRate(const FrontRetreat *m,
                   const std::string &name,
                   const std::string &long_name);
  IceModelVec::Ptr compute_impl() const;
};

} // end of namespace pism


#endif /* FRONTRETREAT_H */
