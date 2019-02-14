/* Copyright (C) 2018 PISM Authors
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

#include "CompleteFrontalMelt.hh"

namespace pism {
namespace frontalmelt {

// "modifier" constructor
CompleteFrontalMelt::CompleteFrontalMelt(IceGrid::ConstPtr g, std::shared_ptr<FrontalMelt> input)
  : FrontalMelt(g, input) {

  m_frontal_melt_rate = allocate_frontal_melt_rate(g);
}

// "model" constructor
CompleteFrontalMelt::CompleteFrontalMelt(IceGrid::ConstPtr g)
  : CompleteFrontalMelt(g, nullptr) {
  // empty
}

CompleteFrontalMelt::~CompleteFrontalMelt() {
  // empty
}

const IceModelVec2S& CompleteFrontalMelt::frontal_melt_rate_impl() const {
  return *m_frontal_melt_rate;
}

} // end of namespace ocean
} // end of namespace pism
