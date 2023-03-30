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

#ifndef PISM_MASKCC_H
#define PISM_MASKCC_H

#include "connected_components_lakecc.hh"

namespace pism {

class MaskCC : public SinkCC {
public:
  MaskCC(IceGrid::ConstPtr g);
  ~MaskCC();
  void compute_mask(IceModelVec2Int &mask);

protected:
  virtual bool ForegroundCond(int i, int j) const;

private:
  void labelOutMask(int run_number, const VecList &lists, IceModelVec2Int &result);
};

} // end of namespace pism

#endif /* PISM_MASKCC_H */
