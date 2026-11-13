/* Copyright (C) 2013, 2014, 2015, 2016, 2017, 2021, 2022, 2023 PISM Authors
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

#ifndef PISM_FLOAT_KILL_H
#define PISM_FLOAT_KILL_H

#include "pism/util/Component.hh"
#include "pism/util/array/CellType.hh"

namespace pism {

namespace calving {

/*! \brief Calving mechanism removing floating ice. */
class FloatKill : public Component
{
public:
  FloatKill(std::shared_ptr<const Grid> g);
  virtual ~FloatKill() = default;

  virtual void init();
  void update(array::Scalar &cell_type, array::Scalar &ice_thickness);

protected:
  bool m_margin_only, m_calve_near_grounding_line, m_calve_cliff_front;

  array::CellType1 m_old_mask;
};


} // end of namespace calving
} // end of namespace pism

#endif /* PISM_FLOAT_KILL_H */
