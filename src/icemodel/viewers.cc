// Copyright (C) 2004-2011, 2013, 2014, 2015, 2016, 2017, 2020, 2021, 2022 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <cstring>
#include <cmath>

#include "IceModel.hh"

#include "pism/util/ConfigInterface.hh"
#include "pism/util/Diagnostic.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/Vars.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/petscwrappers/Viewer.hh"

namespace pism {

void IceModel::view_field(const array::Array *field) {
  unsigned int viewer_size = (unsigned int)m_config->get_number("output.runtime.viewer.size");

  if (field->ndims() != 2) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "map-plane views of 3D quantities are not supported.");
  }

  auto name = field->get_name();

  if (m_viewers[name].empty()) {
    for (size_t k = 0; k < field->ndof(); ++k) {
      auto dof_name = field->metadata(k).get_string("short_name");
      auto v = std::make_shared<petsc::Viewer>(m_grid->com,
                                               dof_name, viewer_size,
                                               m_grid->Lx(), m_grid->Ly());
      m_viewers[name].emplace_back(v);
    }
  }

  field->view(m_viewers[name]);
}

//! Update the runtime graphical viewers.
/*!
Most viewers are updated by this routine, but some other are updated elsewhere.
 */
void IceModel::update_viewers() {

  auto viewers = set_split(m_config->get_string("output.runtime.viewer.variables"), ',');

  // map-plane viewers
  for (auto v : viewers) {
    if (m_grid->variables().is_available(v)) {
      this->view_field(m_grid->variables().get(v));
    } else {
      // if not found, try to compute:
      auto diag = m_diagnostics.find(v);

      if (diag != m_diagnostics.end()) {
        this->view_field(diag->second->compute().get());
      }
    }
  }
}

} // end of namespace pism
