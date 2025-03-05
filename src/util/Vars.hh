// Copyright (C) 2009--2024 PISM Authors
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

#ifndef PISM_VARS_H
#define PISM_VARS_H

#include <map>
#include <set>
#include <string>
#include <memory>

namespace pism {


namespace array {
class Array3D;
class Array;
class CellType;
class Scalar;
class Scalar1;
class Scalar2;
class Vector;
} // end of namespace array

//! \brief A class for passing PISM variables from the core to other parts of
//! the code (such as climate couplers).
class Vars {
public:
  Vars();
  void add(const array::Array &);
  void add(const array::Array &, const std::string &name);
  void remove(const std::string &name);
  bool is_available(const std::string &name) const;

  const array::Array* get(const std::string &name) const;
  const array::Scalar* get_2d_scalar(const std::string &name) const;
  const array::Scalar1* get_2d_scalar1(const std::string &name) const;
  const array::Scalar2* get_2d_scalar2(const std::string &name) const;
  const array::Vector* get_2d_vector(const std::string &name) const;
  const array::CellType* get_2d_cell_type(const std::string &name) const;
  const array::Array3D* get_3d_scalar(const std::string &name) const;

  std::set<std::string> keys() const;

  void add_shared(std::shared_ptr<array::Array>);
  void add_shared(std::shared_ptr<array::Array>, const std::string &name);

  bool is_available_shared(const std::string &name) const;

  std::shared_ptr<array::Array> get_shared(const std::string &name) const;
  std::shared_ptr<array::Scalar> get_2d_scalar_shared(const std::string &name) const;
  std::shared_ptr<array::Scalar1> get_2d_scalar1_shared(const std::string &name) const;
  std::shared_ptr<array::Scalar2> get_2d_scalar2_shared(const std::string &name) const;
  std::shared_ptr<array::Vector> get_2d_vector_shared(const std::string &name) const;
  std::shared_ptr<array::CellType> get_2d_cell_type_shared(const std::string &name) const;
  std::shared_ptr<array::Array3D> get_3d_scalar_shared(const std::string &name) const;

  std::set<std::string> keys_shared() const;
private:
  const array::Array* get_internal(const std::string &name) const;
  mutable std::map<std::string, const array::Array*> m_variables;
  //! stores standard names of variables that
  //! have standard names, allowing looking them
  //! up using either short or standard names and
  //! preserving the one-to-one map from keys
  //! (strings) to pointers (represented by
  //! "variables").
  mutable std::map<std::string, std::string> m_standard_names;

  //! variables in *shared ownership*
  mutable std::map<std::string, std::shared_ptr<array::Array>> m_variables_shared;

  std::shared_ptr<array::Array> get_internal_shared(const std::string &name) const;

  // Hide copy constructor / assignment operator.
  Vars(Vars const &);
  Vars & operator=(Vars const &);
};

} // end of namespace pism

#endif // PISM_VARS_H
