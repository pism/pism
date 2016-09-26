// Copyright (C) 2007--2016 Jed Brown, Ed Bueler and Constantine Khroulev
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

#ifndef __pism_const_hh
#define __pism_const_hh

#include <string>
#include <petsclog.h>

namespace pism {

#ifndef __GNUC__
#  define  __attribute__(x)  /* nothing */
#endif

extern const char *PISM_Revision;
extern const char *PISM_DefaultConfigFile;

const int TEMPORARY_STRING_LENGTH = 32768; // 32KiB ought to be enough.

std::string pism_timestamp(MPI_Comm com);
std::string pism_username_prefix(MPI_Comm com);
std::string pism_args_string();
std::string pism_filename_add_suffix(const std::string &filename,
                                     const std::string &separator,
                                     const std::string &suffix);

double wall_clock_hours(MPI_Comm com, double start_time);
PetscLogDouble GetTime();

} // end of namespace pism

#endif
