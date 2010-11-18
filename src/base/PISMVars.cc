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
    if (variables[standard_name] == NULL)
      variables[standard_name] = &v;
    else
      SETERRQ1(1, "PISMVars::add(): an IceModelVec with the standard_name '%s' was added already.",
               standard_name.c_str());

  } else {

    if (variables[short_name] == NULL)
      variables[short_name] = &v;
    else
      SETERRQ1(1, "PISMVars::add(): an IceModelVec with the short_name '%s' was added already.",
               short_name.c_str());
  }

  if (variables_short[short_name] == NULL)
    variables_short[short_name] = &v;
  else
    SETERRQ1(1, "PISMVars::add(): an IceModelVec with the short_name '%s' was added already.",
             short_name.c_str());

  return 0;
}

//! Removes a variable with the key \c name from the dictionary.
/*! If an IceModelVec has short_name "foo" and the standard_name "bar",
  remove_variable("foo") will not remove the entry corresponding to "bar".
 */
void PISMVars::remove(string name) {
  variables.erase(name);
  variables_short.erase(name);
}

//! \brief Returns a pointer to an IceModelVec containing variable \c name or
//! NULL if that variable was not found.
/*!
 * Checks standard_name first, then short name
 */
IceModelVec* PISMVars::get(string name) const {
  map<string, IceModelVec* >::const_iterator j = variables.find(name);
  if (j != variables.end())
    return (j->second);

  j = variables_short.find(name);
  if (j != variables_short.end())
    return (j->second);

  return NULL;
}

//! \brief Returns the set of keys (variable names) in the dictionary.
/*!
 * Provides one name per variable (the standard_name if it was set, otherwise
 * the short name).
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
