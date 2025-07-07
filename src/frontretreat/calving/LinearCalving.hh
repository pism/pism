/* Copyright (C) 2016, 2017, 2018, 2019, 2021, 2022, 2025 PISM Authors
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

#ifndef LINEARCALVING_H
#define LINEARCALVING_H

#include "pism/util/Component.hh"
#include "pism/util/array/Scalar.hh"
#include "pism/util/array/CellType.hh"

namespace pism {

class Geometry;

namespace calving {

/*!
 * @brief Linear calving mechanism
 *
 * This class implements a linear calving parameterization based on
 * the cliff height above water level. The calving rate is computed as:
 *
 * c = a * H_c + b
 *
 * where:
 * - c is the calving rate
 * - a is the linear coefficient (year^-1)
 * - H_c is the cliff height above water level
 * - b is the constant offset (m/year)
 *
 * The calving rate is set to zero if the calculated value is negative.
 * This mechanism is designed for grounded ice cliffs that form when
 * ice shelves collapse, exposing grounded ice to the ocean.
 *
 * Reference: Parsons et al. (2025)
 */
class LinearCalving : public Component {
public:
  LinearCalving(std::shared_ptr<const Grid> grid);
  virtual ~LinearCalving() = default;

  void init();

  void update(const array::CellType1 &cell_type, const array::Scalar &ice_thickness,
              const array::Scalar &sea_level, const array::Scalar &bed_elevation);

  const array::Scalar & calving_rate() const;

protected:
  DiagnosticList diagnostics_impl() const;
  
protected:
  array::Scalar1 m_calving_rate; 

  double m_a,                    
         m_b;
};

} // end of namespace calving
} // end of namespace pism

#endif /* LINEARCALVING_H */  
