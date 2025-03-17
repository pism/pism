/* Copyright (C) 2015, 2016, 2025 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef _VEC_H_
#define _VEC_H_

#include <petscvec.h>
#include <memory>

#include "pism/util/Wrapper.hh"

namespace pism {
namespace petsc {

class DM;

/** Wrapper around PETSc's Vec. Simplifies memory management.
 *
 * The constructor takes ownership of the Vec argument passed to it.
 *
 * The destructor call VecDestroy().
 */
class Vec : public Wrapper< ::Vec > {
public:
  Vec();
  Vec(::Vec v);
  ~Vec();
};

//! Wrapper around VecGetArray and VecRestoreArray.
class VecArray {
public:
  VecArray(::Vec v);
  ~VecArray();
  double* get();
private:
  ::Vec m_v;
  double *m_array;
};

//! Wrapper around VecGetArray2d and VecRestoreArray2d.
class VecArray2D {
public:
  VecArray2D(::Vec vec, int my_Mx, int my_My);
  VecArray2D(::Vec vec, int my_Mx, int my_My, int i0, int j0);
  ~VecArray2D();

  inline double& operator()(int i, int j) {
    return m_array[j + m_j_offset][i + m_i_offset];
  }
private:
  int m_Mx, m_My, m_i_offset, m_j_offset;
  ::Vec m_v;
  double **m_array;
};

class DMDAVecArray {
public:
  DMDAVecArray(std::shared_ptr<DM> dm, ::Vec v);
  ~DMDAVecArray();
  void* get();
private:
  std::shared_ptr<DM> m_dm;
  ::Vec m_v;
  void *m_array;
};

class DMDAVecArrayDOF {
public:
  DMDAVecArrayDOF(std::shared_ptr<DM> dm, ::Vec v);
  ~DMDAVecArrayDOF();
  void* get();
private:
  std::shared_ptr<DM> m_dm;
  ::Vec m_v;
  void *m_array;
};

class TemporaryGlobalVec : public Vec {
public:
  TemporaryGlobalVec(std::shared_ptr<DM> dm);
  ~TemporaryGlobalVec();
private:
  std::shared_ptr<DM> m_dm;
};

} // end of namespace petsc
} // end of namespace pism


#endif /* _VEC_H_ */
