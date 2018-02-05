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

#ifndef _PISMOCEANKILL_H_
#define _PISMOCEANKILL_H_

#include "pism/util/Component.hh"
#include "pism/util/iceModelVec.hh"

namespace pism {
namespace calving {

/**
 * This class implements the "ocean_kill" mechanism: calving at a
 * fixed calving front determined using ice thickness.
 *
 * FIXME: the ice extent computation should depend on the current sea
 * level elevation (I suppose).
 */
class OceanKill : public Component {
public:
  OceanKill(IceGrid::ConstPtr g);
  virtual ~OceanKill();

  void init();
  void update(IceModelVec2Int &pism_mask, IceModelVec2S &ice_thickness);
  const IceModelVec2Int& mask() const;

protected:
  virtual DiagnosticList diagnostics_impl() const;

  IceModelVec2Int m_mask;
};

} // end of namespace calving
} // end of namespace pism

#endif /* _PISMOCEANKILL_H_ */
