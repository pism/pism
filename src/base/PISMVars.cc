// Copyright (C) 2009, 2010 Constantine Khroulev
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

#include "PISMVars.hh"

//!Add an IceModelVec to the dictionary.
/*!
  Adds both short_name and standard_name (if present); should be called \b after setting the metadata.

  This code will only work for IceModelVecs with dof == 1.
 */
PetscErrorCode PISMVars::add(IceModelVec &v) {

  string short_name = v.string_attr("short_name");

  // an IceModelVec always has a short_name:
  if (variables[short_name] == NULL)
    variables[short_name] = &v;
  else
    SETERRQ1(1, "PISMVars::add(): an IceModelVec with the short_name '%s' was added already.",
	     short_name.c_str());

  // if it doesn't have a standard name, we're done
  if (!v.has_attr("standard_name"))
    return 0;

  string standard_name = v.string_attr("standard_name");
  if (variables[standard_name] == NULL)
    variables[standard_name] = &v;
  else
    SETERRQ1(1, "PISMVars::add(): an IceModelVec with the standard_name '%s' was added already.",
	     standard_name.c_str());

  return 0;
}

//! Removes a variable with the key \c name from the dictionary.
/*! If an IceModelVec has short_name "foo" and the standard_name "bar",
  remove_variable("foo") will not remove the entry corresponding to "bar".
 */
void PISMVars::remove(string name) {
  variables.erase(name);
}

//! \brief Returns a pointer to an IceModelVec containing variable \c name or
//! NULL if that variable was not found.
IceModelVec* PISMVars::get(string name) const {
  return variables[name];
}

//! Returns a set of pointers to all the variables in the dictionary.
set<IceModelVec*> PISMVars::get_variables() const {
  set <IceModelVec*> result;

  map<string,IceModelVec*>::iterator i = variables.begin();
  while (i != variables.end()) {
    result.insert((*i).second);
    ++i;
  }

  return result;
}

//! Debugging: checks if IceModelVecs in the dictionary have NANs.
PetscErrorCode PISMVars::check_for_nan() const {
  PetscErrorCode ierr;
  set<IceModelVec*> vars = get_variables();

  set<IceModelVec*>::iterator j = vars.begin();
  while (j != vars.end()) {
    ierr = (*j)->has_nan(); CHKERRQ(ierr);

    j++;
  }

  return 0;
}
