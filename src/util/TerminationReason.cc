// Copyright (C) 2012, 2014, 2023  David Maxwell
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

#include "pism/util/TerminationReason.hh"

namespace pism {

const char *TerminationReason::sm_indent = "  ";

KSPTerminationReason::KSPTerminationReason(KSPConvergedReason r)  {
  m_reason = r;
}
void KSPTerminationReason::get_description(std::ostream &desc, int indent_level) {
  for (int i=0; i < indent_level; i++) {
    desc << sm_indent;
  }
  desc << KSPConvergedReasons[m_reason];
}

SNESTerminationReason::SNESTerminationReason(SNESConvergedReason r) {
  m_reason = r;
}
void SNESTerminationReason::get_description(std::ostream &desc,int indent_level) {
  for (int i=0; i < indent_level; i++) {
    desc << sm_indent;
  }
  desc << SNESConvergedReasons[m_reason];
}

void GenericTerminationReason::get_description(std::ostream &desc,int indent_level) {
  for (int i=0; i < indent_level; i++) {
    desc << sm_indent;
  }
  desc << m_description;
}

} // end of namespace pism
