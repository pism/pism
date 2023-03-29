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

#ifndef PISM_ISOLATIONCC_H
#define PISM_ISOLATIONCC_H

#include "pism/util/connected_components_lakecc.hh"

namespace pism {

/*!
 * This class marks cells as invalid that are either ice-covered or not connected by an
 * ice-free corridor within the domain margin. It is used to restrict the formation of
 * lakes in the ice sheet interior and of subglacial lakes where thin ice covers a deep
 * basin.
 */
class IsolationCC : public SinkCC {
public:
  IsolationCC(IceGrid::ConstPtr g, const IceModelVec2S &thk,
              double thk_theshold);
  virtual ~IsolationCC() = default;
  void find_isolated_spots(IceModelVec2Int &result);

protected:
  virtual bool ForegroundCond(int i, int j) const;
  void labelIsolatedSpots(int run_number, const VecList &lists, IceModelVec2Int &result);
  void prepare_mask();

private:
  const double m_thk_threshold;
  const IceModelVec2S *m_thk;
};

} // end of namespace pism

#endif /* PISM_ISOLATIONCC_H */
