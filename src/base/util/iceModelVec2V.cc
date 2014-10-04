// Copyright (C) 2009--2014 Constantine Khroulev
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

#include "iceModelVec.hh"
#include "pism_const.hh"
#include "IceGrid.hh"
#include "iceModelVec_helpers.hh"

namespace pism {

IceModelVec2V::IceModelVec2V() : IceModelVec2() {
  m_dof = 2;
  begin_end_access_use_dof = false;
}

PetscErrorCode  IceModelVec2V::create(IceGrid &my_grid, const std::string &my_short_name, IceModelVecKind ghostedp,
                                      unsigned int stencil_width) {

  PetscErrorCode ierr = IceModelVec2::create(my_grid, my_short_name, ghostedp,
                                             stencil_width, m_dof); CHKERRQ(ierr);

  m_metadata[0].init_2d("u" + my_short_name, my_grid);
  m_metadata[1].init_2d("v" + my_short_name, my_grid);

  m_name = "vel" + my_short_name;

  return 0;
}

PetscErrorCode IceModelVec2V::get_array(Vector2** &a) {
  PetscErrorCode ierr;
  ierr = begin_access(); CHKERRQ(ierr);
  a = static_cast<Vector2**>(array);
  return 0;
}

PetscErrorCode IceModelVec2V::magnitude(IceModelVec2S &result) const {
  IceModelVec::AccessList list;
  list.add(*this);
  list.add(result);

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    result(i, j) = (*this)(i, j).magnitude();
  }

  return 0;
}

PetscErrorCode IceModelVec2V::set_name(const std::string &new_name, int component) {
  (void) component;

  std::string tmp = new_name;
  reset_attrs(0);
  reset_attrs(1);
  
  m_name = "vel" + tmp;

  m_metadata[0].set_name("u" + tmp);
  m_metadata[1].set_name("v" + tmp);

  return 0;
}

//! Sets the variable's various names without changing any other metadata
PetscErrorCode IceModelVec2V::rename(const std::string &short_name, const std::string &long_name, 
                                     const std::string &standard_name, int component)
{
  (void) component;

  if (!short_name.empty())
  {
    std::string tmp = short_name;
    m_name = "vel" + tmp;

    m_metadata[0].set_name("u" + tmp);
    m_metadata[1].set_name("v" + tmp);
  }

  if (!long_name.empty()) {
    std::string xprefix = "X component of ";
    std::string yprefix = "Y component of ";
    m_metadata[0].set_string("long_name", xprefix + long_name);
    m_metadata[1].set_string("long_name", yprefix + long_name);
  }

  if (!standard_name.empty()) {
    m_metadata[0].set_string("standard_name", standard_name);
    m_metadata[1].set_string("standard_name", standard_name);
  }

  return 0;
}  

//! Sets the variable's various names without changing any other metadata
PetscErrorCode IceModelVec2V::rename(const std::string & short_name,
                                     const std::vector<std::string> &long_names, 
                                     const std::string & standard_name)
{
  if (!short_name.empty())
  {
    std::string tmp = short_name;

    m_name = "vel" + tmp;

    m_metadata[0].set_name("u" + tmp);
    m_metadata[1].set_name("v" + tmp);
  }

  m_metadata[0].set_string("long_name", long_names[0]);
  m_metadata[1].set_string("long_name", long_names[1]);

  if (!standard_name.empty()) {
    m_metadata[0].set_string("standard_name", standard_name);
    m_metadata[1].set_string("standard_name", standard_name);
  }

  return 0;
}

PetscErrorCode IceModelVec2V::add(double alpha, IceModelVec &x) {
  return add_2d<IceModelVec2V>(this, alpha, &x, this);
}

PetscErrorCode IceModelVec2V::add(double alpha, IceModelVec &x, IceModelVec &result) const {
  return add_2d<IceModelVec2V>(this, alpha, &x, &result);
}

PetscErrorCode IceModelVec2V::copy_to(IceModelVec &destination) const {
  return copy_2d<IceModelVec2V>(this, &destination);
}

} // end of namespace pism
