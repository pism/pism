// Copyright (C) 2012, 2014, 2015, 2017, 2021 David Maxwell and Constantine Khroulev
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

#ifndef IPTWOBLOCKVEC_HH
#define IPTWOBLOCKVEC_HH

#include <petscis.h>
#include <petscvec.h>

#include "pism/util/petscwrappers/Vec.hh"
#include "pism/util/petscwrappers/IS.hh"
#include "pism/util/petscwrappers/VecScatter.hh"

namespace pism {
namespace inverse {

class IPTwoBlockVec {
public:
  IPTwoBlockVec(Vec a, Vec b);
  ~IPTwoBlockVec() = default;

  IS blockAIndexSet();
  IS blockBIndexSet();

  void scatter(Vec a, Vec b);
  void scatterToA(Vec a);
  void scatterToB(Vec b);

  void scatter(Vec ab, Vec a, Vec b);
  void scatterToA(Vec ab, Vec a);
  void scatterToB(Vec ab, Vec b);

  void gather(Vec a, Vec b);
  void gatherFromA(Vec a);
  void gatherFromB(Vec b);

  void gather(Vec a, Vec b, Vec ab);
  void gatherFromA(Vec a, Vec ab);
  void gatherFromB(Vec b, Vec ab);

  operator Vec () {
    return m_ab;
  }

protected:
  void scatter_begin_end(VecScatter s, Vec a, Vec b, ScatterMode m);
  petsc::Vec m_ab;

  PetscInt m_na_local, m_na_global, m_nb_local, m_nb_global;

  petsc::IS m_a_in_ab;
  petsc::IS m_b_in_ab;

  petsc::VecScatter m_scatter_a;
  petsc::VecScatter m_scatter_b;
};

} // end of namespace inverse
} // end of namespace pism

#endif /* end of include guard: IPTWOBLOCKVEC_HH */
