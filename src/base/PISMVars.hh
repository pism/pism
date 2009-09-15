// Copyright (C) 2009 Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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

#ifndef __PISMVars_hh
#define __PISMVars_hh

#include <map>
#include <string>
#include "iceModelVec.hh"


//! \brief A class for passing PISM variables from the core to other parts of
//! the code (such as climate couplers).
class PISMVars {
public:
  PetscErrorCode add(IceModelVec &);
  void remove(string);
  IceModelVec* get(string) const;

protected:
  mutable map<string, IceModelVec*> variables;
};

#endif // __PISMVars_hh
