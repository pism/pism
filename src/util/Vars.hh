// Copyright (C) 2009, 2010, 2013, 2014, 2015, 2016, 2017, 2022 Constantine Khroulev
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

#ifndef __Vars_hh
#define __Vars_hh

#include <map>
#include <set>
#include <string>
#include <memory>

namespace pism {

class IceModelVec;
class IceModelVec2S;
class IceModelVec2V;
class IceModelVec2S;
class IceModelVec3;
namespace array {
using CellType0 = class CellType;
} // end of namespace array

//! \brief A class for passing PISM variables from the core to other parts of
//! the code (such as climate couplers).
class Vars {
public:
  Vars();
  void add(const IceModelVec &);
  void add(const IceModelVec &, const std::string &name);
  void remove(const std::string &name);
  bool is_available(const std::string &name) const;

  const IceModelVec* get(const std::string &name) const;
  const IceModelVec2S* get_2d_scalar(const std::string &name) const;
  const IceModelVec2S* get_2d_mask(const std::string &name) const;
  const IceModelVec2V* get_2d_vector(const std::string &name) const;
  const array::CellType0* get_2d_cell_type(const std::string &name) const;
  const IceModelVec3* get_3d_scalar(const std::string &name) const;

  std::set<std::string> keys() const;

  typedef std::shared_ptr<IceModelVec> VecPtr;
  typedef std::shared_ptr<IceModelVec2S> Vec2SPtr;
  typedef std::shared_ptr<IceModelVec2S> Vec2IntPtr;
  typedef std::shared_ptr<IceModelVec2V> Vec2VPtr;
  typedef std::shared_ptr<array::CellType0> Vec2CellTypePtr;
  typedef std::shared_ptr<IceModelVec3> Vec3Ptr;

  void add_shared(VecPtr);
  void add_shared(VecPtr, const std::string &name);

  bool is_available_shared(const std::string &name) const;

  VecPtr get_shared(const std::string &name) const;
  Vec2SPtr get_2d_scalar_shared(const std::string &name) const;
  Vec2VPtr get_2d_vector_shared(const std::string &name) const;
  Vec2IntPtr get_2d_mask_shared(const std::string &name) const;
  Vec2CellTypePtr get_2d_cell_type_shared(const std::string &name) const;
  Vec3Ptr get_3d_scalar_shared(const std::string &name) const;

  std::set<std::string> keys_shared() const;
private:
  const IceModelVec* get_internal(const std::string &name) const;
  mutable std::map<std::string, const IceModelVec*> m_variables;
  //! stores standard names of variables that
  //! have standard names, allowing looking them
  //! up using either short or standard names and
  //! preserving the one-to-one map from keys
  //! (strings) to pointers (represented by
  //! "variables").
  mutable std::map<std::string, std::string> m_standard_names;

  //! variables in *shared ownership*
  mutable std::map<std::string, VecPtr> m_variables_shared;

  VecPtr get_internal_shared(const std::string &name) const;

  // Hide copy constructor / assignment operator.
  Vars(Vars const &);
  Vars & operator=(Vars const &);
};

} // end of namespace pism

#endif // __Vars_hh
