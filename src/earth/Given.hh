/* Copyright (C) 2020, 2021, 2022 PISM Authors
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

#ifndef PISM_BEDDEF_GIVEN
#define PISM_BEDDEF_GIVEN

#include "BedDef.hh"
#include "pism/util/array/Forcing.hh"

namespace pism {
namespace bed {

/*!
 * Read given bed deformation history from input files.
 */
class Given : public BedDef {
public:
  Given(std::shared_ptr<const IceGrid> grid);
  virtual ~Given() = default;
protected:
  void init_impl(const InputOptions &opts, const array::Scalar &ice_thickness,
                 const array::Scalar &sea_level_elevation);

  void update_impl(const array::Scalar &ice_thickness,
                   const array::Scalar &sea_level_elevation,
                   double t, double dt);

  array::Scalar m_topg_reference;

  std::shared_ptr<array::Forcing> m_topg_delta;
};

} // end of namespace bed
} // end of namespace pism

#endif /* PISM_BEDDEF_GIVEN */
