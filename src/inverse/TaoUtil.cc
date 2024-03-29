// Copyright (C) 2012, 2014, 2015, 2023  David Maxwell and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include "pism/inverse/TaoUtil.hh"

namespace pism {
namespace taoutil {

TAOTerminationReason::TAOTerminationReason(TaoConvergedReason r)  {
  m_reason = r;
}

void TAOTerminationReason::get_description(std::ostream &desc, int indent_level) {
  for (int i=0; i < indent_level; i++) {
    desc << sm_indent;
  }
  desc << TaoConvergedReasons[m_reason];
}

} // end of namespace taoutil
} // end of namespace pism
