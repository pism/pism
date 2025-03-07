// Copyright (C) 2009--2011, 2013, 2014, 2015, 2016, 2017, 2020, 2021, 2022, 2023, 2024 Constantine Khroulev
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

#include "array/Scalar.hh"
#include <memory>
using std::dynamic_pointer_cast;

#include "pism/util/Vars.hh"
#include "pism/util/VariableMetadata.hh"
#include "pism/util/array/CellType.hh"
#include "pism/util/array/Vector.hh"
#include "pism/util/array/Array3D.hh"
#include "pism/util/error_handling.hh"

namespace pism {

Vars::Vars() {
}

bool Vars::is_available(const std::string &name) const {
  // check if "name" is a standard name
  if (m_standard_names.find(name) != m_standard_names.end()) {
    return true;
  }
  // check if "name" is a short name of a read-only variable
  if (m_variables.find(name) != m_variables.end()) {
    return true;
  }
  // check if "name" is a short name of a "shared" variable
  if (m_variables_shared.find(name) != m_variables_shared.end()) {
    return true;
  }
  return false;
}

//! \brief Add an array::Array v using the name `name`.
void Vars::add(const array::Array &v, const std::string &name) {

  if (m_variables.find(name) != m_variables.end()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Vars::add(): an array::Array with the name '%s'"
                                  " was added already.",
                                  name.c_str());
  }
  m_variables[name] = &v;
}

//!Add an array::Array to the dictionary.
/*!
  Adds standard_name if present, otherwise uses short_name.

  This code will only work for array::Arrays with dof == 1.
*/
void Vars::add(const array::Array &v) {

  const SpatialVariableMetadata &m = v.metadata();
  std::string name = v.get_name();

  if (m.has_attribute("standard_name")) {

    std::string standard_name = m["standard_name"];
    if (m_standard_names[standard_name].empty()) {
      m_standard_names[standard_name] = name;
    } else {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Vars::add(): an array::Array with the standard_name '%s'"
                                    " was added already.",
                                    standard_name.c_str());
    }
  }

  if (m_variables[name] == NULL) {
    m_variables[name] = &v;
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Vars::add(): an array::Array with the name '%s'"
                                  " was added already.",
                                  name.c_str());
  }
}

//! Removes a variable with the key `name` from the dictionary.
void Vars::remove(const std::string &name) {

  const array::Array *v = m_variables[name];
  const SpatialVariableMetadata &m = v->metadata();

  if (v != NULL) {              // the argument is a "short" name
    m_variables.erase(name);
    if (m.has_attribute("standard_name")) {
      std::string std_name = m["standard_name"];

      m_standard_names.erase(std_name);
    }
  } else {
    std::string &short_name = m_standard_names[name];
    v = m_variables[short_name];

    if (v != NULL) {            // the argument is a standard_name
      m_variables.erase(short_name);
      m_standard_names.erase(name);
    }
  }
}

//! \brief Returns a pointer to an array::Array containing variable `name` or
//! NULL if that variable was not found.
/*!
 * Checks standard_name first, then short name
 */
const array::Array* Vars::get(const std::string &name) const {
  const array::Array *tmp = get_internal(name);
  if (tmp == NULL) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "variable '%s' is not available", name.c_str());
  }
  return tmp;
}

const array::Array* Vars::get_internal(const std::string &name) const {

  auto j = m_standard_names.find(name);
  if (j != m_standard_names.end()) {
    std::string short_name = j->second;

    auto k = m_variables.find(short_name);
    if (k != m_variables.end()) {
      return k->second;
    }
  }

  auto k = m_variables.find(name);
  if (k != m_variables.end()) {
    return (k->second);
  }

  std::shared_ptr<array::Array> shared = get_internal_shared(name);
  if ((bool)shared) {
    return shared.get();
  }

  return NULL;
}

template<class A>
const A* get_(const Vars *vars, const std::string &name, const std::string &kind) {
  auto tmp = dynamic_cast<const A*>(vars->get(name));

  if (tmp == NULL) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "%s variable '%s' is not available",
                                  kind.c_str(), name.c_str());
  }

  return tmp;
}

const array::Scalar* Vars::get_2d_scalar(const std::string &name) const {
  return get_<array::Scalar>(this, name, "2D scalar");
}

const array::Scalar1* Vars::get_2d_scalar1(const std::string &name) const {
  return get_<array::Scalar1>(this, name, "2D scalar (ghosted)");
}

const array::Scalar2* Vars::get_2d_scalar2(const std::string &name) const {
  return get_<array::Scalar2>(this, name, "2D scalar (ghosted, stencil width=2)");
}

const array::Vector* Vars::get_2d_vector(const std::string &name) const {
  return get_<array::Vector>(this, name, "2D vector");
}

const array::CellType* Vars::get_2d_cell_type(const std::string &name) const {
  return get_<array::CellType>(this, name, "2D cell type");
}

const array::Array3D* Vars::get_3d_scalar(const std::string &name) const {
  return get_<array::Array3D>(this, name, "3D scalar");
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

  for (const auto &v : m_variables) {
    result.insert(v.first);
  }

  for (const auto &v : m_variables_shared) {
    result.insert(v.first);
  }

  return result;
}

void Vars::add_shared(std::shared_ptr<array::Array> variable) {

  const SpatialVariableMetadata &m = variable->metadata();
  std::string name = variable->get_name();

  if (m.has_attribute("standard_name")) {

    std::string standard_name = m["standard_name"];
    if (m_standard_names[standard_name].empty()) {
      m_standard_names[standard_name] = name;
    } else {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Vars::add_shared(): an array::Array with the standard_name '%s'"
                                    " was added already.",
                                    standard_name.c_str());
    }
  }

  if (m_variables_shared.find(name) == m_variables_shared.end()) {
    m_variables_shared[name] = variable;
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Vars::add_shared(): an array::Array with the name '%s'"
                                  " was added already.",
                                  name.c_str());
  }
}


void Vars::add_shared(std::shared_ptr<array::Array> variable, const std::string &name) {

  if (m_variables_shared.find(name) != m_variables_shared.end()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Vars::add_shared(): an array::Array with the name '%s'"
                                  " was added already.",
                                  name.c_str());
  }
  m_variables_shared[name] = std::move(variable);
}


bool Vars::is_available_shared(const std::string &name) const {

  // check the standard name
  if (m_standard_names.find(name) != m_standard_names.end()) {
    std::string short_name = m_standard_names[name];
    // return true if the corresponding short name is one of a
    // "shared" variable
    return (m_variables_shared.find(short_name) != m_variables_shared.end());
  }

  // check if "name" is a short name of a "shared" variable
  return (m_variables_shared.find(name) != m_variables_shared.end());
}

template <class A>
std::shared_ptr<A> get_shared_(const Vars *vars, const std::string &name, const std::string &kind) {
  auto tmp = dynamic_pointer_cast<A, array::Array>(vars->get_shared(name));

  if (not(bool) tmp) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "shared %s variable '%s' is not available",
                                  kind.c_str(), name.c_str());
  }

  return tmp;
}

std::shared_ptr<array::Array> Vars::get_shared(const std::string &name) const {
  auto tmp = get_internal_shared(name);
  if (not (bool)tmp) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "shared variable '%s' is not available", name.c_str());
  }
  return tmp;
}

std::shared_ptr<array::Scalar> Vars::get_2d_scalar_shared(const std::string &name) const {
  return get_shared_<array::Scalar>(this, name, "2D scalar");
}

std::shared_ptr<array::Scalar1> Vars::get_2d_scalar1_shared(const std::string &name) const {
  return get_shared_<array::Scalar1>(this, name, "2D scalar (ghosted)");
}

std::shared_ptr<array::Scalar2> Vars::get_2d_scalar2_shared(const std::string &name) const {
  return get_shared_<array::Scalar2>(this, name, "2D scalar (ghosted, stencil width=2)");
}


std::shared_ptr<array::Vector> Vars::get_2d_vector_shared(const std::string &name) const {
  return get_shared_<array::Vector>(this, name, "2D vector");
}

std::shared_ptr<array::CellType> Vars::get_2d_cell_type_shared(const std::string &name) const {
  return get_shared_<array::CellType>(this, name, "2D cell type");
}


std::shared_ptr<array::Array3D> Vars::get_3d_scalar_shared(const std::string &name) const {
  return get_shared_<array::Array3D>(this, name, "3D scalar");
}

std::set<std::string> Vars::keys_shared() const {

  std::set<std::string> result;
  for (const auto &v : m_variables_shared) {
    result.insert(v.first);
  }

  return result;
}

std::shared_ptr<array::Array> Vars::get_internal_shared(const std::string &name) const {

  auto j = m_standard_names.find(name);
  if (j != m_standard_names.end()) {
    std::string short_name = j->second;

    auto k = m_variables_shared.find(short_name);
    if (k != m_variables_shared.end()) {
      return k->second;
    }

    return std::shared_ptr<array::Array>();
  }

  auto k = m_variables_shared.find(name);
  if (k != m_variables_shared.end()) {
    return (k->second);
  }

  return std::shared_ptr<array::Array>();
}

} // end of namespace pism
