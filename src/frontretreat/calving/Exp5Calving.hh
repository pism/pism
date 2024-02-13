/* Copyright (C) 2013, 2014, 2015, 2016, 2017, 2018, 2021, 2022, 2024 PISM Authors
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

#ifndef _PISMEXP5CALVING_H_
#define _PISMEXP5CALVING_H_

#include "pism/util/Component.hh"
#include "pism/util/array/CellType.hh"
#include "pism/util/array/Scalar.hh"

#include "pism/util/array/Array2D.hh"



namespace pism {

class Geometry;

namespace calving {

/*! \brief Calving mechanism Exp5 from https://github.com/JRowanJordan/CalvingMIP/wiki/Experiment-5. */
class Exp5Calving : public Component
{
public:
  Exp5Calving(IceGrid::ConstPtr grid);
  virtual ~Exp5Calving() = default;

  void init();

  void update(const array::CellType1 &cell_type,
              const array::Vector1 &ice_velocity, 
              const array::Scalar &ice_thickness);

  const array::Scalar &calving_rate() const;


protected:
  DiagnosticList diagnostics_impl() const;
  double calving_threshold;
  array::Scalar1 m_calving_rate;

};

} // end of namespace calving
} // end of namespace pism

#endif /* _PISMEXP5CALVING_H_ */
