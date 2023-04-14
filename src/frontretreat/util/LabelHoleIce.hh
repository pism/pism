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

#ifndef _PISMLABELHOLEICE_H_
#define _PISMLABELHOLEICE_H_

#include "pism/util/Component.hh"
#include "pism/util/iceModelVec.hh"

namespace pism {

class IceModelVec2CellType;

namespace calving {

/*! \brief PISM ice shelf hole labeling */
/*!
 * Identifies and labels holes in floating ice shelves, which could
 * cause unrealistic mass loss if calving and the iceberg remover
 * work together.
 *
 * Holes may occur in tiny ice shelves by extreme atmospheric or
 * oceanographic conditions. Since these holes are considered as
 * ocean points, may occur, which releases icebergs into the hole,
 * which are afterwards removed by the "IcebergRemover". Continuing
 * calving could disintegrate ice shelves from the interior, while
 * the stabilizing effect of the created ice melange is not taken into
 * account. Ultimately, a catastrophic ice shelf disintegration is
 * triggered and PISM losses a lot of mass. It is similar to mass falling
 * into a _black hole_  from where it can not escape, hence, you may call
 * it "black hole calving".
 *
 * Holes in ice shelves are ocean points that are surrounded by floating
 * ice, grounded ice, or ice free bedrock.
 *
 * Please note; we do not distinguish between ocean and lake. A lake would
 * be a water body entirely surrounded by grounded ice or ice free bedrock.
 *
 * Like the IcebergRemover, this class uses a serial connected component
 * labeling algorithm to identify and label holes in ice shelves as
 * "enclosed ocean" (MASK_ICE_FREE_ENCLOSED_OCEAN).
 *
 */
class LabelHoleIce : public Component
{
public:
  LabelHoleIce(IceGrid::ConstPtr g);
  virtual ~LabelHoleIce();

  virtual void init();
  void update(IceModelVec2CellType &pism_mask);

  //todo:rm?;void open_ocean_margin_retreat(const IceModelVec2T &retreat_mask,
  void open_ocean_mask_margin_retreat(const IceModelVec2S &bed,
				      const IceModelVec2S &sea_level,
				      const IceModelVec2S &ice_area_specific_volume,
				      const IceModelVec2S &ice_thickness);
				      //todo:alternative;const IceModelVec2S &retreat_mask);

  void open_ocean_mask_margin(const IceModelVec2S &bed,
			      const IceModelVec2S &sea_level);

protected:
  IceModelVec2Int m_bc_open_ocean_mask;
  IceModelVec2S m_enclosed_ocean_mask;
  petsc::Vec::Ptr m_mask_enclose_ocean_p0;
};

} // end of namespace calving
} // end of namespace pism

#endif /* _PISMLABELHOLEICE_H_ */
