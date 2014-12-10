// Copyright (C) 2004-2014 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "pism_const.hh"
#include "iceModelVec.hh"
#include "columnSystem.hh"
#include <cassert>

#include <fstream>
#include <iostream>

#include "error_handling.hh"
#include "ColumnInterpolation.hh"

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

Result should be executable as part of a Matlab/Octave script.

Does not stop on non-fatal errors.
 */
void TridiagonalSystem::save_vector(std::ostream &output,
                                    const std::vector<double> &v,
                                    unsigned int system_size,
                                    const std::string &variable) const {
  assert(system_size >= 1);

  output << "%% viewing ColumnSystem column object with description '" << variable << "'"
         << " (columns  [k value])" << std::endl;

  output << variable << "_with_index = [..." << std::endl;
  for (unsigned int k = 0; k < system_size; k++) {
    output << "  " << k << " " << v[k];
    if (k == system_size - 1) {
      output << "];" << std::endl;
    } else {
      output << ";" << std::endl;
    }
  }

  output << variable << " = " << variable << "_with_index(:,2);\n" << std::endl;
}


//! View the tridiagonal matrix.  Views as a full matrix if nmax <= 120, otherwise by listing diagonals.
/*!
Give first argument NULL to get standard out.  No binary viewer.

Give description string as `info` argument.
 */
void TridiagonalSystem::save_matrix(std::ostream &output,
                                 unsigned int system_size,
                                 const std::string &variable) const {

  if (system_size < 2) {
    std::cout << "\n\n<nmax >= 2 required to view tri-diagonal matrix " << variable
              << " ... skipping view" << std::endl;
    return;
  }

  if (system_size > 500) {
    std::cout << "\n\n<nmax > 500:" << variable
              << " matrix too big to display as full; viewing tridiagonal matrix diagonals..."
              << std::endl;

    save_vector(output, m_U, system_size + 1, variable + "_super_diagonal_U");
    save_vector(output, m_D, system_size,   variable + "_diagonal_D");

    // discard m_L[0], which is not used
    {
      std::vector<double> L_tmp(system_size - 1);
      for (unsigned int i = 0; i < system_size - 1; ++i) {
        L_tmp[i] = m_L[i + 1];
      }
      save_vector(output, L_tmp, system_size - 1, variable + "_sub_diagonal_L");
    }
  } else {
    output << "\n"
           << variable << " = [..." << std::endl;

    for (unsigned int i = 0; i < system_size; ++i) { // row
      for (unsigned int j = 0; j < system_size; j++) { // column
        double A_ij = 0.0;

        if (j == i - 1) {
          A_ij = m_L[i];
        } else if (j == i) {
          A_ij = m_D[i];
        } else if (j == i + 1) {
          A_ij = m_U[i];
        } else {
          A_ij = 0.0;
        }

        output << A_ij << " ";
      } // column loop

      if (i != system_size - 1) {
        output << ";" << std::endl;
      } else {
        output << "];" << std::endl;
      }
    } // row-loop
  } // end of printing the full matrix
}


//! View the tridiagonal system A x = b to an output stream, both A as a full matrix and b as a vector.
void TridiagonalSystem::save_system(std::ostream &output,
                                    unsigned int system_size) const {
  save_matrix(output, system_size, m_prefix + "_A");
  save_vector(output, m_rhs, system_size, m_prefix + "_rhs");
}

//! Write system matrix, right-hand-side, and (provided) solution into an already-named m-file.
void TridiagonalSystem::save_system_with_solution(const std::string &filename,
                                                  unsigned int M,
                                                  const std::vector<double> &x) {
  std::ofstream output(filename);
  output << "% system has 1-norm = " << norm1(M)
         << " and diagonal-dominance ratio = " << ddratio(M) << std::endl;

  save_system(output, M);
  save_vector(output, x, M, m_prefix + "_x");
}


//! The actual code for solving a tridiagonal system.
/*!
This is modified slightly from a Numerical Recipes version.

Input size n is size of instance.  Requires n <= TridiagonalSystem::m_max_system_size.

Solution of system in x.
 */
void TridiagonalSystem::solve(unsigned int system_size, std::vector<double> &result) {
  assert(system_size >= 1);
  assert(system_size <= m_max_system_size);

  if (m_D[0] == 0.0) {
    throw RuntimeError("zero pivot at row 1");
  }

  result.resize(m_max_system_size);

  double b = m_D[0];

  result[0] = m_rhs[0] / b;
  for (unsigned int k = 1; k < system_size; ++k) {
    m_work[k] = m_U[k - 1] / b;

    b = m_D[k] - m_L[k] * m_work[k];

    if (b == 0.0) {
      throw RuntimeError::formatted("zero pivot at row %d", k + 1);
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
                                 IceModelVec3 *u3, IceModelVec3 *v3, IceModelVec3 *w3)
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
  double *array = NULL;
  coarse.getInternalColumn(i, j, &array);
  m_interp->fine_to_coarse(&fine[0], array);
}

void columnSystemCtx::coarse_to_fine(IceModelVec3 *coarse, int i, int j, int ks,
                                     double* fine) const {
  double *array = NULL;
  coarse->getInternalColumn(i, j, &array);
  m_interp->coarse_to_fine(array, ks, fine);
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
    throw RuntimeError::formatted("ks = %d computed at i = %d, j = %d is invalid,\n"
                                  "possibly because of invalid ice thickness (%f meters) or dz (%f meters).",
                                  m_ks, m_i, m_j, ice_thickness, m_dz);
  }
#endif
}

//! Write system matrix and right-hand-side into an m-file.  The file name contains ZERO_PIVOT_ERROR.
void columnSystemCtx::reportColumnZeroPivotErrorMFile(unsigned int M) {
  char filename[TEMPORARY_STRING_LENGTH];
  snprintf(filename, sizeof(filename), "%s_i%d_j%d_ZERO_PIVOT_ERROR.m",
           m_solver->prefix().c_str(), m_i, m_j);

  std::ofstream output(filename);
  output << "% system has 1-norm = " << m_solver->norm1(M)
         << " and diagonal-dominance ratio = " << m_solver->ddratio(M) << std::endl;

  m_solver->save_system(output, M);
}


//! @brief Write system matrix, right-hand-side, and (provided)
//! solution into an m-file. Constructs file name from m_prefix.
/*!
An example of the use of this procedure is from <c>examples/searise-greenland/</c>
running the enthalpy formulation.  First run spinup.sh in that directory  (FIXME:
which was modified to have equal spacing in z, when I did this example) to
generate `g20km_steady.nc`.  Then:

\code
  $ pismr -calving ocean_kill -e 3 -atmosphere searise_greenland -surface pdd \
          -config_override  config_269.0_0.001_0.80_-0.500_9.7440.nc \
          -no_mass -y 1 -i g20km_steady.nc -view_sys -id 19 -jd 79

    ...

  $ octave -q
  >> enth_i19_j79
  >> whos

    ...

  >> A = enth_A; b = enth_rhs; x = enth_x;
  >> norm(x - (A\b))/norm(x)
  ans =  1.4823e-13
  >> cond(A)
  ans =  2.6190
\endcode

Of course we can also do `spy(A)`, `eig(A)`, and look at individual entries,
and row and column sums, and so on.
 */
void columnSystemCtx::viewColumnInfoMFile(const std::vector<double> &x) {
  std::cout << "saving "
            << m_solver->prefix() << " column system at (i,j)"
            << " = (" << m_i << "," << m_j << ") to m-file...\n" << std::endl;

  char buffer[TEMPORARY_STRING_LENGTH];
  snprintf(buffer, sizeof(buffer), "%s_i%d_j%d.m", m_solver->prefix().c_str(), m_i, m_j);

  m_solver->save_system_with_solution(buffer, m_z.size(), x);
}

} // end of namespace pism
