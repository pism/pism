// Copyright (C) 2009, 2010, 2013, 2014 Constantine Khroulev
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

#ifndef __PISMVars_hh
#define __PISMVars_hh

#include <map>
#include <set>
#include <string>
#include "iceModelVec.hh"

namespace pism {


//! \brief A class for passing PISM variables from the core to other parts of
//! the code (such as climate couplers).
class Vars {
public:
  PetscErrorCode add(IceModelVec &);
  PetscErrorCode add(IceModelVec &, const std::string &name);
  void remove(const std::string &);
  IceModelVec* get(const std::string &) const;
  std::set<std::string> keys() const;
  PetscErrorCode check_for_nan() const;

protected:
  mutable std::map<std::string, IceModelVec*> variables,
    standard_names;             //!< stores standard names of variables that
  //! have standard names, allowing looking them
  //! up using either short or standard names and
  //! preserving the one-to-one map from keys
  //! (strings) to pointers (represented by
  //! "variables").
};

} // end of namespace pism

#endif // __PISMVars_hh
