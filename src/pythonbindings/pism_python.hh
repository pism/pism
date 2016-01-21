// Copyright (C) 2011, 2014, 2015 David Maxwell and Constantine Khroulev
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


#ifndef _PISM_PYTHON_
#define _PISM_PYTHON_

namespace pism {

//! @brief Code added for use in Python wrappers.
namespace python {

void set_abort_on_sigint(bool abort);

int check_signal();
void sigint_handler(int sig);

extern bool gSIGINT_is_fatal;

//! Installs a signal handler on construction; deinstalls on destruction.
class SigInstaller
{
public:
  //! Installs handle \a new_handler for signal \a sig.
  SigInstaller(int sig, void (*new_handler)(int));
  //! Restores the signal handler to its previous value.
  void release();

  ~SigInstaller();
private:
  void (*m_old_handler)(int);
  int m_sig;
};

} // end of namespace python
} // end of namespace pism

#endif
