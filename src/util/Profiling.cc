/* Copyright (C) 2015, 2016, 2021 PISM Authors
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

#include <petscviewer.h>

#include "Profiling.hh"
#include "error_handling.hh"

namespace pism {

// PETSc profiling events

Profiling::Profiling() {
  PetscErrorCode ierr = PetscClassIdRegister("PISM", &m_classid);
  PISM_CHK(ierr, "PetscClassIdRegister");
}

//! Enable PETSc logging.
void Profiling::start() const {
#if PETSC_VERSION_LE(3,6,3)
  PetscErrorCode ierr = PetscLogBegin(); PISM_CHK(ierr, "PetscLogBegin");
#else
  PetscErrorCode ierr = PetscLogAllBegin(); PISM_CHK(ierr, "PetscLogAllBegin");
#endif
}

//! Save detailed profiling data to a Python script.
void Profiling::report(const std::string &filename) const {
  PetscErrorCode ierr;
  PetscViewer log_viewer;

  ierr = PetscViewerCreate(PETSC_COMM_WORLD, &log_viewer);
  PISM_CHK(ierr, "PetscViewerCreate");

  ierr = PetscViewerSetType(log_viewer, PETSCVIEWERASCII);
  PISM_CHK(ierr, "PetscViewerSetType");

  ierr = PetscViewerFileSetName(log_viewer, filename.c_str());
  PISM_CHK(ierr, "PetscViewerFileSetName");

  ierr = PetscViewerPushFormat(log_viewer, PETSC_VIEWER_ASCII_INFO_DETAIL);
  PISM_CHK(ierr, "PetscViewerPushFormat");

  ierr = PetscViewerSetUp(log_viewer);
  PISM_CHK(ierr, "PetscViewerSetUp");

  ierr = PetscLogView(log_viewer);
  PISM_CHK(ierr, "PetscLogView"); PISM_CHK(ierr, "PetscLogView");

  ierr = PetscViewerPopFormat(log_viewer);
  PISM_CHK(ierr, "PetscViewerPopFormat");

  ierr = PetscViewerDestroy(&log_viewer);
  PISM_CHK(ierr, "PetscViewerDestroy");
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

  if (m_events.find(name) == m_events.end()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "cannot end event \"%s\" because it was not started",
                                  name);
  }

  PetscErrorCode ierr = PetscLogEventEnd(m_events[name], 0, 0, 0, 0);
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
