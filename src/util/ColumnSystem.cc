// Copyright (C) 2004-2019 PISM Authors
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

#include <cassert>
#include <fstream>
#include <iostream>

#include "pism/util/pism_utilities.hh"
#include "pism/util/iceModelVec.hh"
#include "ColumnSystem.hh"

#include "pism/util/error_handling.hh"
#include "pism/util/ColumnInterpolation.hh"

namespace pism {


//! Allocate a tridiagonal system of maximum size max_system_size.
/*!
Let N = `max_system_size`.  Then allocated locations are like this:
\verbatim
D[0]   U[0]    0      0      0    ...
L[1]   D[1]   U[1]    0      0    ...
 0     L[2]   D[2]   U[2]    0    ...
 0      0     L[3]   D[3]   U[3]  ...
\endverbatim
with the last row
\verbatim
0       0     ...     0  L[N-1]  D[N-1]
\endverbatim
Thus the index into the arrays L, D, U is always the row number.
 */
TridiagonalSystem::TridiagonalSystem(unsigned int max_size,
                                     const std::string &prefix)
  : m_max_system_size(max_size), m_prefix(prefix) {
  assert(m_max_system_size >= 1 && m_max_system_size < 1e6);

  m_L.resize(m_max_system_size);
  m_D.resize(m_max_system_size);
  m_U.resize(m_max_system_size);
  m_rhs.resize(m_max_system_size);
  m_work.resize(m_max_system_size);
}

//! Zero all entries.
void TridiagonalSystem::reset() {
#if PISM_DEBUG==1
  memset(&m_L[0],    0, (m_max_system_size)*sizeof(double));
  memset(&m_U[0],    0, (m_max_system_size)*sizeof(double));
  memset(&m_D[0],    0, (m_max_system_size)*sizeof(double));
  memset(&m_rhs[0],  0, (m_max_system_size)*sizeof(double));
  memset(&m_work[0], 0, (m_max_system_size)*sizeof(double));
#endif
}

//! Compute 1-norm, which is max sum of absolute values of columns.
double TridiagonalSystem::norm1(unsigned int system_size) const {
  assert(system_size <= m_max_system_size);
  if (system_size == 1) {
    return fabs(m_D[0]);   // only 1x1 case is special
  }
  double z = fabs(m_D[0]) + fabs(m_L[1]);
  for (unsigned int k = 1; k < system_size; k++) {  // k is column index (zero-based)
    z = std::max(z, fabs(m_U[k-1])) + fabs(m_D[k]) + fabs(m_L[k+1]);
  }
  z = std::max(z, fabs(m_U[system_size-2]) + fabs(m_D[system_size-1]));
  return z;
}


//! Compute diagonal-dominance ratio.  If this is less than one then the matrix is strictly diagonally-dominant.
/*!
Let \f$A = (a_{ij})\f$ be the tridiagonal matrix
described by L, D, U for row indices 0 through `n`.  The computed ratio is
  \f[ \max_{j=1, \dots, n} \frac{|a_{j, j-1}|+|a_{j, j+1}|}{|a_{jj}|}, \f]
where \f$a_{1, 0}\f$ and \f$a_{n, n+1}\f$ are interpreted as zero.

If this is smaller than one then it is a theorem that the tridiagonal solve will
succeed.

We return -1.0 if the absolute value of any diagonal element is less than
1e-12 of the 1-norm of the matrix.
 */
double TridiagonalSystem::ddratio(unsigned int system_size) const {
  assert(system_size <= m_max_system_size);

  const double scale = norm1(system_size);

  if ((fabs(m_D[0]) / scale) < 1.0e-12) {
    return -1.0;
  }
  double z = fabs(m_U[0]) / fabs(m_D[0]);

  for (unsigned int k = 1; k < system_size - 1; k++) {  // k is row index (zero-based)
    if ((fabs(m_D[k]) / scale) < 1.0e-12) {
      return -1.0;
    }
    const double s = fabs(m_L[k]) + fabs(m_U[k]);
    z = std::max(z, s / fabs(m_D[k]));
  }

  if ((fabs(m_D[system_size - 1]) / scale) < 1.0e-12) {
    return -1.0;
  }
  z = std::max(z, fabs(m_L[system_size - 1]) / fabs(m_D[system_size - 1]));

  return z;
}

//! Utility for simple ascii view of a vector (one-dimensional column) quantity.
/*!
Give first argument NULL to get standard out.  No binary viewer.

Give description string as `info` argument.

Result should be executable as a Python (SciPy) script.

Does not stop on non-fatal errors.
 */
void TridiagonalSystem::save_vector(std::ostream &output,
                                    const std::vector<double> &v,
                                    unsigned int system_size,
                                    const std::string &variable) const {
  assert(system_size >= 1);

  output << variable << " = numpy.array([";
  for (unsigned int k = 0; k < system_size; k++) {
    output << v[k] << ", ";
  }
  output << "])" << std::endl;
}


//! View the tridiagonal matrix.
void TridiagonalSystem::save_matrix(std::ostream &output,
                                    unsigned int system_size,
                                    const std::string &variable) const {

  if (system_size < 2) {
    std::cout << "\n\n<nmax >= 2 required to view tri-diagonal matrix " << variable
              << " ... skipping view" << std::endl;
    return;
  }

  save_vector(output, m_U, system_size, variable + "_U");
  save_vector(output, m_D, system_size, variable + "_D");
  save_vector(output, m_L, system_size, variable + "_L");

  // prepare to convert to a sparse matrix
  output << variable << "_diagonals = ["
         << variable << "_L[1:], "
         << variable << "_D, "
         << variable << "_U[:-1]]" << std::endl;
  output << "import scipy.sparse" << std::endl;
  output << variable << " = scipy.sparse.diags(" << variable << "_diagonals, [-1, 0, 1])" << std::endl;
}


//! View the tridiagonal system A x = b to an output stream, both A as a full matrix and b as a vector.
void TridiagonalSystem::save_system(std::ostream &output,
                                    unsigned int system_size) const {
  save_matrix(output, system_size, m_prefix + "_A");
  save_vector(output, m_rhs, system_size, m_prefix + "_rhs");
}

//! Write system matrix, right-hand-side, and (provided) solution into an already-named Python
//! script.
void TridiagonalSystem::save_system_with_solution(const std::string &filename,
                                                  unsigned int M,
                                                  const std::vector<double> &x) {
  std::ofstream output(filename.c_str());
  output << "import numpy" << std::endl;
  output << "# system has 1-norm = " << norm1(M)
         << " and diagonal-dominance ratio = " << ddratio(M) << std::endl;

  save_system(output, M);
  save_vector(output, x, M, m_prefix + "_x");
}


void TridiagonalSystem::solve(unsigned int system_size, std::vector<double> &result) {
  result.resize(m_max_system_size);

  solve(system_size, result.data());
}


//! The actual code for solving a tridiagonal system.
/*!
This is modified slightly from a Numerical Recipes version.

Input size n is size of instance.  Requires n <= TridiagonalSystem::m_max_system_size.

Solution of system in x.
 */
void TridiagonalSystem::solve(unsigned int system_size, double *result) {
  assert(system_size >= 1);
  assert(system_size <= m_max_system_size);

  if (m_D[0] == 0.0) {
    throw RuntimeError(PISM_ERROR_LOCATION, "zero pivot at row 1");
  }

  double b = m_D[0];

  result[0] = m_rhs[0] / b;
  for (unsigned int k = 1; k < system_size; ++k) {
    m_work[k] = m_U[k - 1] / b;

    b = m_D[k] - m_L[k] * m_work[k];

    if (b == 0.0) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "zero pivot at row %d", k + 1);
    }

    result[k] = (m_rhs[k] - m_L[k] * result[k-1]) / b;
  }

  for (int k = system_size - 2; k >= 0; --k) {
    result[k] -= m_work[k + 1] * result[k + 1];
  }
}

std::string TridiagonalSystem::prefix() const {
  return m_prefix;
}

//! A column system is a kind of a tridiagonal system.
columnSystemCtx::columnSystemCtx(const std::vector<double>& storage_grid,
                                 const std::string &prefix,
                                 double dx, double dy, double dt,
                                 const IceModelVec3 &u3,
                                 const IceModelVec3 &v3,
                                 const IceModelVec3 &w3)
  : m_dx(dx), m_dy(dy), m_dt(dt), m_u3(u3), m_v3(v3), m_w3(w3) {
  assert(dx > 0.0);
  assert(dy > 0.0);
  assert(dt > 0.0);

  init_fine_grid(storage_grid);

  m_solver = new TridiagonalSystem(m_z.size(), prefix);

  m_interp = new ColumnInterpolation(storage_grid, m_z);

  m_u.resize(m_z.size());
  m_v.resize(m_z.size());
  m_w.resize(m_z.size());
}

columnSystemCtx::~columnSystemCtx() {
  delete m_solver;
  delete m_interp;
}

unsigned int columnSystemCtx::ks() const {
  return m_ks;
}

double columnSystemCtx::dz() const {
  return m_dz;
}

const std::vector<double>& columnSystemCtx::z() const {
  return m_z;
}

void columnSystemCtx::fine_to_coarse(const std::vector<double> &fine, int i, int j,
                                     IceModelVec3& coarse) const {
  double *array = coarse.get_column(i, j);
  m_interp->fine_to_coarse(&fine[0], array);
}

void columnSystemCtx::coarse_to_fine(const IceModelVec3 &coarse, int i, int j,
                                     double* fine) const {
  const double *array = coarse.get_column(i, j);
  m_interp->coarse_to_fine(array, m_ks, fine);
}

void columnSystemCtx::init_fine_grid(const std::vector<double>& storage_grid) {
  // Compute m_dz as the minimum vertical spacing in the coarse
  // grid:
  unsigned int Mz = storage_grid.size();
  double Lz = storage_grid.back();
  m_dz = Lz;
  for (unsigned int k = 1; k < Mz; ++k) {
    m_dz = std::min(m_dz, storage_grid[k] - storage_grid[k - 1]);
  }

  size_t Mz_fine = static_cast<size_t>(ceil(Lz / m_dz) + 1);
  m_dz = Lz / (Mz_fine - 1);

  m_z.resize(Mz_fine);
  // compute levels of the fine grid:
  for (unsigned int k = 0; k < Mz_fine; ++k) {
    m_z[k] = storage_grid[0] + k * m_dz;
  }
  // Note that it *is* allowed to go over Lz.
}

void columnSystemCtx::init_column(int i, int j,
                                  double ice_thickness) {
  m_i  = i;
  m_j  = j;
  m_ks = static_cast<unsigned int>(floor(ice_thickness / m_dz));

  // Force m_ks to be in the allowed range.
  if (m_ks >= m_z.size()) {
    m_ks = m_z.size() - 1;
  }

  m_solver->reset();

  // post-condition
#if PISM_DEBUG==1
  // check if m_ks is valid
  if (m_ks >= m_z.size()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "ks = %d computed at i = %d, j = %d is invalid,\n"
                                  "possibly because of invalid ice thickness (%f meters) or dz (%f meters).",
                                  m_ks, m_i, m_j, ice_thickness, m_dz);
  }
#endif
}

//! Write system matrix and right-hand-side into an Python script.  The file name contains ZERO_PIVOT_ERROR.
void columnSystemCtx::reportColumnZeroPivotErrorMFile(unsigned int M) {

  auto filename = pism::printf("%s_i%d_j%d_ZERO_PIVOT_ERROR.py",
                               m_solver->prefix().c_str(), m_i, m_j);

  std::ofstream output(filename);
  output << "# system has 1-norm = " << m_solver->norm1(M)
         << " and diagonal-dominance ratio = " << m_solver->ddratio(M) << std::endl;

  m_solver->save_system(output, M);
}


//! @brief Write system matrix, right-hand-side, and (provided)
//! solution into Python script. Constructs file name from m_prefix.
void columnSystemCtx::save_to_file(const std::vector<double> &x) {

  auto filename = pism::printf("%s_i%d_j%d.py", m_solver->prefix().c_str(), m_i, m_j);

  std::cout << "saving "
            << m_solver->prefix() << " column system at (i,j)"
            << " = (" << m_i << "," << m_j << ") to " << filename << " ...\n" << std::endl;

  this->save_to_file(filename, x);
}

void columnSystemCtx::save_to_file(const std::string &filename,
                                   const std::vector<double> &x) {
  m_solver->save_system_with_solution(filename, m_z.size(), x);
}

} // end of namespace pism
