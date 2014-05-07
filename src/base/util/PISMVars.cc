// Copyright (C) 2009--2011, 2013, 2014 Constantine Khroulev
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

#include "PISMVars.hh"

namespace pism {

//! \brief Add an IceModelVec v using the name `name`.
PetscErrorCode Vars::add(IceModelVec &v, const std::string &name) {

  variables[name] = &v;

  return 0;
}

//!Add an IceModelVec to the dictionary.
/*!
  Adds standard_name if present, otherwise uses short_name.

  This code will only work for IceModelVecs with dof == 1.
 */
PetscErrorCode Vars::add(IceModelVec &v) {

  NCSpatialVariable &m = v.metadata();
  std::string name = v.name();

  if (m.has_attribute("standard_name")) {

    std::string standard_name = m.get_string("standard_name");
    if (standard_names[standard_name] == NULL)
      standard_names[standard_name] = &v;
    else
      SETERRQ1(PETSC_COMM_SELF, 1, "Vars::add(): an IceModelVec with the standard_name '%s' was added already.",
               standard_name.c_str());

  }

  if (variables[name] == NULL)
    variables[name] = &v;
  else
    SETERRQ1(PETSC_COMM_SELF, 1, "Vars::add(): an IceModelVec with the short_name '%s' was added already.",
             name.c_str());

  return 0;
}

//! Removes a variable with the key `name` from the dictionary.
void Vars::remove(const std::string &name) {
  IceModelVec *v = variables[name];
  NCSpatialVariable &m = v->metadata();

  if (v != NULL) {              // the argument is a "short" name
    if (m.has_attribute("standard_name")) {
      std::string std_name = m.get_string("standard_name");

      variables.erase(name);
      standard_names.erase(std_name);
    }
  } else {
    v = standard_names[name];

    if (v != NULL) {            // the argument is a standard_name
      std::string short_name = v->name();

      variables.erase(short_name);
      standard_names.erase(name);
    }
  }

}

//! \brief Returns a pointer to an IceModelVec containing variable `name` or
//! NULL if that variable was not found.
/*!
 * Checks standard_name first, then short name
 */
IceModelVec* Vars::get(const std::string &name) const {
  std::map<std::string, IceModelVec* >::const_iterator j = standard_names.find(name);
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
std::set<std::string> Vars::keys() const {
  std::set<std::string> result;

  std::map<std::string,IceModelVec*>::iterator i = variables.begin();
  while (i != variables.end())
    result.insert((*i++).first);

  return result;
}

//! Debugging: checks if IceModelVecs in the dictionary have NANs.
PetscErrorCode Vars::check_for_nan() const {
  PetscErrorCode ierr;
  std::set<std::string> names = keys();

  std::set<std::string>::iterator i = names.begin();
  while (i != names.end()) {
    IceModelVec *var = get(*i);
    ierr = var->has_nan(); CHKERRQ(ierr);
    i++;
  }

  return 0;
}

} // end of namespace pism
