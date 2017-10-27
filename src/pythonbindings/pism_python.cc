// Copyright (C) 2011, 2014, 2015, 2016, 2017 David Maxwell and Constantine Khroulev
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

#include <signal.h>
#include <errno.h>
#include <petscsys.h>

#include "pism_python.hh"
#include "util/pism_utilities.hh"
#include "util/error_handling.hh"

namespace pism {
namespace python {

void set_abort_on_sigint(bool abort) {
  gSIGINT_is_fatal = abort;
}

SigInstaller::SigInstaller(int sig, void (*new_handler)(int) ) {
  m_old_handler = signal(sig, new_handler);
  m_sig = sig;
}

void SigInstaller::release() {
  if (m_old_handler) {
    signal(m_sig, m_old_handler);
  }
  m_old_handler = NULL;
}

SigInstaller::~SigInstaller() {
  release();
}

int gSignalSet = 0;
bool gSIGINT_is_fatal;

int check_signal() {
  int rv = gSignalSet;
  if (rv) {
    gSignalSet = 0;
  }
  return 0;
}

void sigint_handler(int sig) {
  if (sig == SIGINT) {
    if (gSIGINT_is_fatal) {
      throw pism::RuntimeError(ErrorLocation(), "caught signal SIGTERM.");
    } else {
      PetscPrintf(PETSC_COMM_WORLD, "\nCaught signal SIGTERM, waiting to exit.\n");
      gSignalSet = SIGINT;
    }
  }
}

} // end of namespace python
} // end of namespace pism
