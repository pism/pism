/* Copyright (C) 2016--2026 PISM Authors
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

#ifndef PISM_CLIFFCALVINGTENSILE_H
#define PISM_CLIFFCALVINGTENSILE_H

#include "pism/util/Component.hh"
#include "pism/util/array/Scalar.hh"

namespace pism {

namespace array {
class CellType1;
}

namespace calving {

class CliffCalvingTensile : public Component {
public:
  CliffCalvingTensile(std::shared_ptr<const Grid> grid);
  virtual ~CliffCalvingTensile() = default;

  void init();

  void update(const array::CellType1 &cell_type, const array::Scalar1 &ice_thickness,
              const array::Scalar &sea_level, const array::Scalar &bed_elevation);

  const array::Scalar &calving_rate() const;

protected:
  DiagnosticList diagnostics_impl() const;
  
  array::Scalar1 m_calving_rate; 

  double m_I;
  double m_alpha;

  //! if true, use floatation thickness instead of modelled ice thickness
  bool m_use_floatation_thickness;
};

} // end of namespace calving
} // end of namespace pism

#endif /* PISM_CLIFFCALVINGTENSILE_H */
