/* Copyright (C) 2016, 2018, 2019, 2021, 2022 PISM Authors
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

#ifndef STRESSCALVING_H
#define STRESSCALVING_H

#include "pism/util/Component.hh"
#include "pism/util/array/Scalar.hh"
#include "pism/util/array/Array3D.hh"
#include "pism/util/array/CellType.hh"

namespace pism {

namespace calving {

/*! @brief An abstract class containing fields used by all stress-based calving methods. */
class StressCalving : public Component {
public:
  StressCalving(IceGrid::ConstPtr grid, unsigned int stencil_width);
  virtual ~StressCalving() = default;

  const array::Scalar &calving_rate() const;

protected:
  const int m_stencil_width;

  array::Array3D m_strain_rates;

  array::Scalar m_calving_rate;

  array::CellType1 m_cell_type;
};


} // end of namespace calving
} // end of namespace pism


#endif /* STRESSCALVING_H */
