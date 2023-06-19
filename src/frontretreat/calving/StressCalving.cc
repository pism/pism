/* Copyright (C) 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023 PISM Authors
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

#include "pism/frontretreat/calving/StressCalving.hh"

namespace pism {
namespace calving {

StressCalving::StressCalving(std::shared_ptr<const Grid> grid,
                             unsigned int stencil_width)
  : Component(grid),
    m_stencil_width(stencil_width),
    m_strain_rates(m_grid, "strain_rates", array::WITH_GHOSTS, 2),
    m_calving_rate(m_grid, "calving_rate"),
    m_cell_type(m_grid, "cell_type")
{

  m_strain_rates.metadata(0).set_name("eigen1");
  m_strain_rates.set_attrs("internal",
                           "major principal component of horizontal strain-rate",
                           "second-1", "second-1", "", 0);

  m_strain_rates.metadata(1).set_name("eigen2");
  m_strain_rates.set_attrs("internal",
                           "minor principal component of horizontal strain-rate",
                           "second-1", "second-1", "", 1);

  m_calving_rate.set_attrs("internal", "horizontal calving rate", "m s-1", "m year-1", "", 0);

  m_cell_type.set_attrs("internal", "cell type mask", "", "", "", 0);
}

const array::Scalar &StressCalving::calving_rate() const {
  return m_calving_rate;
}

} // end of namespace calving
} // end of namespace pism
