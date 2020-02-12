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

#ifndef _PISMCALVINGATTHICKNESS_H_
#define _PISMCALVINGATTHICKNESS_H_

#include "pism/util/Component.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/IceModelVec2CellType.hh"

namespace pism {
namespace calving {

/*! \brief Calving mechanism removing the ice at the shelf front that
  has thickness below a given threshold. */
class CalvingAtThickness : public Component
{
public:
  CalvingAtThickness(IceGrid::ConstPtr g);
  virtual ~CalvingAtThickness();

  virtual void init();
  void update(IceModelVec2CellType &pism_mask, IceModelVec2S &ice_thickness);
  const IceModelVec2S& threshold() const;

protected:
  virtual DiagnosticList diagnostics_impl() const;
  IceModelVec2S m_calving_threshold;
  IceModelVec2CellType m_old_mask;
};

} // end of namespace calving
} // end of namespace pism

#endif /* _PISMCALVINGATTHICKNESS_H_ */
