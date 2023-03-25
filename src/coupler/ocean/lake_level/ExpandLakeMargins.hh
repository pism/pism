/* Copyright (C) 2018, 2023 PISM Authors
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

#ifndef PISM_EXPAND_LAKE_MARGINS_H
#define PISM_EXPAND_LAKE_MARGINS_H

#include "pism/coupler/LakeLevel.hh"

namespace pism {
namespace ocean {
namespace lake_level {

class ExpandLakeMargins : public LakeLevel {
public:
  ExpandLakeMargins(IceGrid::ConstPtr g, std::shared_ptr<LakeLevel> in);
  ~ExpandLakeMargins();

private:
  void update_impl(const Geometry &geometry, double t, double dt);
  void init_impl(const Geometry &geometry);

};

} // end of namespace lake_level
} // end of namespace ocean
} // end of namespace pism

#endif /* PISM_EXPAND_LAKE_MARGINS_H */
