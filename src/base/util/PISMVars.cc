// Copyright (C) 2009--2011 Constantine Khroulev
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

//! \brief Add an IceModelVec v using the name \c name.
PetscErrorCode PISMVars::add(IceModelVec &v, string name) {

  variables[name] = &v;

  return 0;
}

//!Add an IceModelVec to the dictionary.
/*!
  Adds standard_name if present, otherwise uses short_name.

  This code will only work for IceModelVecs with dof == 1.
 */
PetscErrorCode PISMVars::add(IceModelVec &v) {

  string short_name = v.string_attr("name");

  if (v.has_attr("standard_name")) {

    string standard_name = v.string_attr("standard_name");
    if (standard_names[standard_name] == NULL)
      standard_names[standard_name] = &v;
    else
      SETERRQ1(PETSC_COMM_SELF, 1, "PISMVars::add(): an IceModelVec with the standard_name '%s' was added already.",
               standard_name.c_str());

  }

  if (variables[short_name] == NULL)
    variables[short_name] = &v;
  else
    SETERRQ1(PETSC_COMM_SELF, 1, "PISMVars::add(): an IceModelVec with the short_name '%s' was added already.",
             short_name.c_str());

  return 0;
}

//! Removes a variable with the key \c name from the dictionary.
void PISMVars::remove(string name) {
  IceModelVec *v = variables[name];

  if (v != NULL) {              // the argument is a "short" name
    if (v->has_attr("standard_name")) {
      string std_name = v->string_attr("standard_name");

      variables.erase(name);
      standard_names.erase(std_name);
    }
  } else {
    v = standard_names[name];

    if (v != NULL) {            // the argument is a standard_name
      string short_name = v->string_attr("name");

      variables.erase(short_name);
      standard_names.erase(name);
    }
  }

}

//! \brief Returns a pointer to an IceModelVec containing variable \c name or
//! NULL if that variable was not found.
/*!
 * Checks standard_name first, then short name
 */
IceModelVec* PISMVars::get(string name) const {
  map<string, IceModelVec* >::const_iterator j = standard_names.find(name);
  if (j != standard_names.end())
    return (j->second);

  j = variables.find(name);
  if (j != variables.end())
    return (j->second);

  return NULL;
}

//! \brief Returns the set of keys (variable names) in the dictionary.
/*!
 * Provides one (short) name per variable.
 *
 * This means that one can safely iterate over these keys, reading, writing,
 * displaying or de-allocating variables (without worrying that a variable will
 * get written or de-allocated twice).
 */
set<string> PISMVars::keys() const {
  set<string> result;

  map<string,IceModelVec*>::iterator i = variables.begin();
  while (i != variables.end())
    result.insert((*i++).first);

  return result;
}

//! Debugging: checks if IceModelVecs in the dictionary have NANs.
PetscErrorCode PISMVars::check_for_nan() const {
  PetscErrorCode ierr;
  set<string> names = keys();

  set<string>::iterator i = names.begin();
  while (i != names.end()) {
    IceModelVec *var = get(*i);
    ierr = var->has_nan(); CHKERRQ(ierr);
    i++;
  }

  return 0;
}
