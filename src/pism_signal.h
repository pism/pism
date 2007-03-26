// Copyright (C) 2007 Jed Brown
//
// This file is part of Pism.
//
// Pism is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// Pism is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with Pism; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef __pism_signal_h
#define __pism_signal_h

#ifdef __cplusplus
extern "C" 
{
#endif

#include <signal.h>

  extern volatile sig_atomic_t pism_signal;

  void pism_signal_handler(int sig);

#ifdef __cplusplus
}
#endif

#endif
