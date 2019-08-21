// Copyright (C) 2009-2011, 2013, 2014, 2015, 2016, 2019 PISM Authors
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

#ifndef __columnSystem_hh
#define __columnSystem_hh

#include <string>
#include <ostream>
#include <vector>

namespace pism {

//! Virtual base class.  Abstracts a tridiagonal system to solve in a column of ice and/or bedrock.
/*!
  Because both the age evolution and conservation of energy equations require us to set up
  and solve a tridiagonal system of equations, this is structure is worth abstracting.

  This base class just holds the tridiagonal system and the ability to
  solve it, but does not insert entries into the relevant matrix locations.
  Derived classes will actually set up instances of the system.

  The sequence requires setting the column-independent (public) data members,
  calling the initAllColumns() routine, and then setting up and solving
  the system in each column.

  The tridiagonal algorithm here comes from Numerical Recipes in C
  Section 2.4, page 50.  It solves the system:

@verbatim
  [b_1  c_1   0   ...                           ] [  u_1  ]   [   r_1   ]
  [a_2  b_2  c_2  ...                           ] [  u_2  ]   [   r_2   ]
  [               ...                           ].[  ...  ] = [   ...   ]
  [               ...  a_{N-1}  b_{N-1}  c_{N-1}] [u_{N-1}]   [ r_{N-1} ]
  [               ...     0     a_N      b_N    ] [  u_N  ]   [   r_N   ]
@endverbatim

  HOWEVER... the version in this code is different from Numerical
  Recipes in two ways:

  - Indexing is zero-based
  -  Variables have been renamed.

@verbatim
  NR      PISM
  ==================
  a       L       "Lower Diagonal" (L doesn't use index 0)
  b       D       "Diagonal"
  c       U       "Upper Diagonal"
  u       x
  r       rhs
  bet     b
  j       k
  n       n
  gam     work
@endverbatim

  Therefore... this version of the code solves the following problem:

@verbatim
  [D_0  U_0   0   ...                           ] [  x_0  ]   [   r_0   ]
  [L_1  D_1  U_1  ...                           ] [  x_1  ]   [   r_1   ]
  [               ...                           ].[  ...  ] = [   ...   ]
  [               ...  L_{N-2}  D_{N-2}  U_{N-2}] [x_{N-2}]   [ r_{N-2} ]
  [               ...     0     L_{N-1}  D_{N-1}] [x_{N-1}]   [ r_{N-1} ]
@endverbatim
*/
class TridiagonalSystem {
public:
  TridiagonalSystem(unsigned int max_size, const std::string &prefix);

  double norm1(unsigned int system_size) const;
  double ddratio(unsigned int system_size) const;
  void reset();

  // uses an output argument to allow re-using storage and avoid
  // copying
  void solve(unsigned int system_size, std::vector<double> &result);
  void solve(unsigned int system_size, double *result);

  void save_system_with_solution(const std::string &filename,
                                 unsigned int system_size,
                                 const std::vector<double> &solution);

  //! Save the system to a stream using the ASCII MATLAB (Octave)
  //! format. Virtual to allow saving more info in derived classes.
  void save_system(std::ostream &output,
                   unsigned int system_size) const;


  void save_matrix(std::ostream &output,
                   unsigned int system_size,
                   const std::string &info) const;

  void save_vector(std::ostream &output,
                   const std::vector<double> &v,
                   unsigned int system_size, const std::string &info) const;

  std::string prefix() const;

  double& L(size_t i) {
    return m_L[i];
  }
  double& D(size_t i) {
    return m_D[i];
  }
  double& U(size_t i) {
    return m_U[i];
  }
  double& RHS(size_t i) {
    return m_rhs[i];
  }
private:
  unsigned int m_max_system_size;         // maximum system size
  std::vector<double> m_L, m_D, m_U, m_rhs, m_work; // vectors for tridiagonal system

  std::string m_prefix;
};

class IceModelVec3;
class ColumnInterpolation;

//! Base class for tridiagonal systems in the ice.
/*! Adds data members used in time-dependent systems with advection
  (dx, dy, dz, dt, velocity components).
 */
class columnSystemCtx {
public:
  columnSystemCtx(const std::vector<double>& storage_grid, const std::string &prefix,
                  double dx, double dy, double dt,
                  const IceModelVec3 &u3, const IceModelVec3 &v3, const IceModelVec3 &w3);
  ~columnSystemCtx();

  void save_to_file(const std::vector<double> &x);
  void save_to_file(const std::string &filename, const std::vector<double> &x);

  unsigned int ks() const;
  double dz() const;
  const std::vector<double>& z() const;
  void fine_to_coarse(const std::vector<double> &fine, int i, int j,
                      IceModelVec3& coarse) const;
protected:
  TridiagonalSystem *m_solver;

  ColumnInterpolation *m_interp;

  //! current system size; corresponds to the highest vertical level within the ice
  unsigned int m_ks;
  //! current column indexes
  int m_i, m_j;

  double m_dx, m_dy, m_dz, m_dt;

  //! u-component of the ice velocity
  std::vector<double> m_u;
  //! v-component of the ice velocity
  std::vector<double> m_v;
  //! w-component of the ice velocity
  std::vector<double> m_w;
  //! levels of the fine vertical grid
  std::vector<double> m_z;

  //! pointers to 3D velocity components
  const IceModelVec3 &m_u3, &m_v3, &m_w3;

  void init_column(int i, int j, double ice_thickness);

  void reportColumnZeroPivotErrorMFile(unsigned int M);

  void init_fine_grid(const std::vector<double>& storage_grid);

  void coarse_to_fine(const IceModelVec3 &coarse, int i, int j, double* fine) const;
};

} // end of namespace pism

#endif  /* __columnSystem_hh */

