/* Copyright (C) 2015 PISM Authors
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

#include "Profiling.hh"
#include "error_handling.hh"

namespace pism {

// PETSc profiling events

Profiling::Profiling() {
  PetscErrorCode ierr = PetscClassIdRegister("PISM", &m_classid);
  PISM_CHK(ierr, "PetscClassIdRegister");
}

void Profiling::begin(const char * name) const {
  PetscLogEvent event = 0;
  PetscErrorCode ierr;

  if (m_events.find(name) == m_events.end()) {
    // not registered yet
    ierr = PetscLogEventRegister(name, m_classid, &event);
    PISM_CHK(ierr, "PetscLogEventRegister");
    m_events[name] = event;
  } else {
    event = m_events[name];
  }
  ierr = PetscLogEventBegin(event, 0, 0, 0, 0);
  PISM_CHK(ierr, "PetscLogEventBegin");
}

void Profiling::end(const char * name) const {
  PetscLogEvent event = 0;
  if (m_events.find(name) == m_events.end()) {
    throw RuntimeError::formatted("cannot end event \"%s\" because it was not started",
                                  name);
  } else {
    event = m_events[name];
  }
  PetscErrorCode ierr = PetscLogEventEnd(event, 0, 0, 0, 0);
  PISM_CHK(ierr, "PetscLogEventEnd");
}

void Profiling::stage_begin(const char * name) const {
  PetscLogStage stage = 0;
  PetscErrorCode ierr;

  if (m_stages.find(name) == m_stages.end()) {
    // not registered yet
    ierr = PetscLogStageRegister(name, &stage);
    PISM_CHK(ierr, "PetscLogStageRegister");
    m_stages[name] = stage;
  } else {
    stage = m_stages[name];
  }
  ierr = PetscLogStagePush(stage);
  PISM_CHK(ierr, "PetscLogStagePush");
}

void Profiling::stage_end(const char * name) const {
  (void) name;
  PetscErrorCode ierr = PetscLogStagePop();
  PISM_CHK(ierr, "PetscLogStagePop");
}

} // end of namespace pism
