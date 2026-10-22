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

#ifndef PISM_CLIFFCALVINGSHEAR_H
#define PISM_CLIFFCALVINGSHEAR_H

#include "pism/util/Component.hh"
#include "pism/util/array/Scalar.hh"

namespace pism {

namespace array {
class CellType1;
}

namespace calving {

class CliffCalvingShear : public Component {
public:
  CliffCalvingShear(std::shared_ptr<const Grid> grid);
  virtual ~CliffCalvingShear() = default;

  void init();

  void update(const array::CellType1 &cell_type, const array::Scalar1 &ice_thickness,
              const array::Scalar &sea_level, const array::Scalar &bed_elevation);

  const array::Scalar &calving_rate() const;

protected:
  DiagnosticList diagnostics_impl() const;
  
  //! calving rate field
  array::Scalar1 m_calving_rate;

  //! scaling factor for calving rate
  double m_C0;

  //! maximum cliff calving rate due to mélange buttressing
  double m_max_cliff_calving_rate;

  //! if true, use floatation thickness instead of modelled ice thickness
  bool m_use_floatation_thickness;
};

} // end of namespace calving
} // end of namespace pism

#endif /* PISM_CLIFFCALVINGSHEAR_H */
