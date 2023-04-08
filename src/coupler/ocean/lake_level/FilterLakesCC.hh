/* Copyright (C) 2023 PISM Authors
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

#ifndef PISM_FILTERLAKESCC_H
#define PISM_FILTERLAKESCC_H

#include "pism/util/connected_components_lakecc.hh"

namespace pism {

/*! Removes narrow lakes
 *
 * The FilterLakesCC class checks the lakes’ geometry. Only lakes containing one (or more)
 * cells that have at least a given number of neighbors that also are part of that lake
 * are retained.
 *
 * We do this to remove narrow lakes which are often related to under-resolved topography.
 */
class FilterLakesCC : public ValidCC<ConnectedComponents> {
public:
  FilterLakesCC(IceGrid::ConstPtr g, double fill_value);
  virtual ~FilterLakesCC() = default;
  void filter_map(int n_filter, IceModelVec2S &lake_level);

private:
  bool ForegroundCond(int i, int j) const;
};

} // end of namespace pism

#endif /* PISM_FILTERLAKESCC_H */
