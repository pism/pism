// Copyright (C) 2004-2011, 2013, 2014, 2015, 2016, 2017 Jed Brown, Ed Bueler and Constantine Khroulev
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

namespace pism {

void IceModel::view_field(const IceModelVec *field) {
  unsigned int viewer_size = (unsigned int)m_config->get_double("output.runtime.viewer.size");

  unsigned int dims = field->ndims();

  if (dims != 2) {
    throw RuntimeError(PISM_ERROR_LOCATION, "map-plane views of 3D quantities are not supported.");
  }

  if (field->ndof() == 1) {    // scalar fields
    std::string name = field->metadata().get_string("short_name");
    petsc::Viewer::Ptr viewer = m_viewers[name];

    if (not viewer) {
      m_viewers[name].reset(new petsc::Viewer(m_grid->com, name, viewer_size, m_grid->Lx(), m_grid->Ly()));
      viewer = m_viewers[name];
    }

    const IceModelVec2S *v2d = dynamic_cast<const IceModelVec2S*>(field);
    if (v2d == NULL) {
      throw RuntimeError(PISM_ERROR_LOCATION, "get_ndims() returns GRID_2D but dynamic_cast gives a NULL");
    }

    v2d->view(viewer, petsc::Viewer::Ptr());

  } else if (field->ndof() == 2) { // vector fields
    std::string
      name_1 = field->metadata(0).get_string("short_name"),
      name_2 = field->metadata(1).get_string("short_name");
    petsc::Viewer::Ptr
      v1 = m_viewers[name_1],
      v2 = m_viewers[name_2];

    if (not v1) {
      m_viewers[name_1].reset(new petsc::Viewer(m_grid->com, name_1, viewer_size, m_grid->Lx(), m_grid->Ly()));
      v1 = m_viewers[name_1];
    }

    if (not v2) {
      m_viewers[name_2].reset(new petsc::Viewer(m_grid->com, name_2, viewer_size, m_grid->Lx(), m_grid->Ly()));
      v2 = m_viewers[name_2];
    }

    const IceModelVec2 *v2d = dynamic_cast<const IceModelVec2*>(field);
    if (v2d == NULL) {
      throw RuntimeError(PISM_ERROR_LOCATION, "get_ndims() returns GRID_2D but dynamic_cast gives a NULL");
    }

    v2d->view(v1, v2);
  }
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
