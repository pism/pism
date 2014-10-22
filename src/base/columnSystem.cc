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

namespace pism {

//! Allocate a tridiagonal system of maximum size nmax.
/*!
Let N = `nmax`.  Then allocated locations are like this:
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
columnSystemCtx::columnSystemCtx(unsigned int nmax, const std::string &my_prefix)
  : m_nmax(nmax), m_prefix(my_prefix) {
  assert(m_nmax >= 1 && m_nmax < 1e6);

  L.resize(m_nmax);
  D.resize(m_nmax);
  U.resize(m_nmax);
  rhs.resize(m_nmax);
  work.resize(m_nmax);

  resetColumn();

  m_indicesValid = false;
}


columnSystemCtx::~columnSystemCtx() {
}

unsigned int columnSystemCtx::ks() const {
  return m_ks;
}

//! Zero all entries.
void columnSystemCtx::resetColumn() {
#if PISM_DEBUG==1
  memset(&L[0],    0, (m_nmax)*sizeof(double));
  memset(&U[0],    0, (m_nmax)*sizeof(double));
  memset(&D[0],    0, (m_nmax)*sizeof(double));
  memset(&rhs[0],  0, (m_nmax)*sizeof(double));
  memset(&work[0], 0, (m_nmax)*sizeof(double)); 
#endif
}


//! Compute 1-norm, which is max sum of absolute values of columns.
double columnSystemCtx::norm1(unsigned int n) const {
  assert(n <= m_nmax);
  if (n == 1)  {
    return fabs(D[0]);   // only 1x1 case is special
  }
  double z = fabs(D[0]) + fabs(L[1]);
  for (unsigned int k = 1; k < n; k++) {  // k is column index (zero-based)
    z = std::max(z, fabs(U[k-1])) + fabs(D[k]) + fabs(L[k+1]);
  }
  z = std::max(z, fabs(U[n-2]) + fabs(D[n-1]));
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
double columnSystemCtx::ddratio(unsigned int n) const {
  assert(n <= m_nmax);

  const double scale = norm1(n);

  if ((fabs(D[0]) / scale) < 1.0e-12) {
    return -1.0;
  }
  double z = fabs(U[0]) / fabs(D[0]);

  for (unsigned int k = 1; k < n - 1; k++) {  // k is row index (zero-based)
    if ((fabs(D[k]) / scale) < 1.0e-12) {
      return -1.0;
    }
    const double s = fabs(L[k]) + fabs(U[k]);
    z = std::max(z, s / fabs(D[k]));
  }

  if ((fabs(D[n - 1]) / scale) < 1.0e-12) {
    return -1.0;
  }
  z = std::max(z, fabs(L[n - 1]) / fabs(D[n - 1]));

  return z;
}


void columnSystemCtx::setIndicesAndClearThisColumn(int my_i, int my_j,
                                                   double ice_thickness,
                                                   double dz,
                                                   unsigned int Mz) {
  // pre-condition
#if PISM_DEBUG==1
  if (m_indicesValid && m_i == my_i && m_j == my_j) {
    throw RuntimeError("setIndicesAndClearThisColumn() called twice in same column");
  }
#endif

  m_i  = my_i;
  m_j  = my_j;
  m_ks = static_cast<unsigned int>(floor(ice_thickness / dz));

  // Force m_ks to be in the allowed range.
  if (m_ks >= Mz) {
    m_ks = Mz - 1;
  }

  resetColumn();

  m_indicesValid = true;

  // post-condition
#if PISM_DEBUG==1
  // check if m_ks is valid
  if (m_ks >= Mz) {
    throw RuntimeError::formatted("ks = %d computed at i = %d, j = %d is invalid,\n"
                                  "possibly because of invalid ice thickness (%f meters) or dz (%f meters).",
                                  m_ks, m_i, m_j, ice_thickness, dz);
  }
#endif
}


//! Utility for simple ascii view of a vector (one-dimensional column) quantity.
/*!
Give first argument NULL to get standard out.  No binary viewer.

Give description string as `info` argument.

Result should be executable as part of a Matlab/Octave script.

Does not stop on non-fatal errors.
 */
void columnSystemCtx::viewVectorValues(std::ostream &output,
                                       const std::vector<double> &v,
                                       unsigned int M,
                                       const std::string &variable) const {
  assert(M >= 1);

  output << "%% viewing ColumnSystem column object with description '" << variable << "'"
         << " (columns  [k value])" << std::endl;

  output << variable << "_with_index = [..." << std::endl;
  for (unsigned int k = 0; k < M; k++) {
    output << "  " << k << " " << v[k];
    if (k == M-1) {
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
void columnSystemCtx::viewMatrix(std::ostream &output,
                                 unsigned int M,
                                 const std::string &variable) const {

  if (M < 2) {
    std::cout << "\n\n<nmax >= 2 required to view tri-diagonal matrix " << variable
              << " ... skipping view" << std::endl;
    return;
  }

  if (M > 500) {
    std::cout << "\n\n<nmax > 500:" << variable
              << " matrix too big to display as full; viewing tridiagonal matrix diagonals..."
              << std::endl;

    viewVectorValues(output, U, M-1, variable + "_super_diagonal_U");
    viewVectorValues(output, D, M,   variable + "_diagonal_D");

    // discard L[0], which is not used
    {
      std::vector<double> L_tmp(M - 1);
      for (unsigned int i = 0; i < M - 1; ++i) {
        L_tmp[i] = L[i + 1];
      }
      viewVectorValues(output, L_tmp, M-1, variable + "_sub_diagonal_L");
    }
  } else {
    output << "\n"
           << variable << " = [..." << std::endl;

    for (unsigned int i = 0; i < M; ++i) { // row
      for (unsigned int j = 0; j < M; j++) { // column
        double A_ij = 0.0;

        if (j == i - 1) {
          A_ij = L[i];
        } else if (j == i) {
          A_ij = D[i];
        } else if (j == i + 1) {
          A_ij = U[i];
        } else {
          A_ij = 0.0;
        }

        output << A_ij << " ";
      } // column loop

      if (i != M-1) {
        output << ";" << std::endl;
      } else {
        output << "];" << std::endl;
      }
    } // row-loop
  } // end of printing the full matrix
}


//! View the tridiagonal system A x = b to an output stream, both A as a full matrix and b as a vector.
void columnSystemCtx::viewSystem(std::ostream &output,
                                 unsigned int M) const {
  viewMatrix(output, M, m_prefix + "_A");
  viewVectorValues(output, rhs, M, m_prefix + "_rhs");
}


//! The actual code for solving a tridiagonal system.  Return code has diagnostic importance.
/*!
This is modified slightly from a Numerical Recipes version.

Input size n is size of instance.  Requires n <= columnSystemCtx::m_nmax.

Solution of system in x.

Success is return code zero.  Positive return code gives location of zero pivot.
Negative return code indicates a software problem.
 */
void columnSystemCtx::solveTridiagonalSystem(unsigned int n, std::vector<double> &x) {
  assert(m_indicesValid == true);
  assert(n >= 1);
  assert(n <= m_nmax);

  if (D[0] == 0.0)
    throw RuntimeError("zero pivot at row 1");

  x.resize(m_nmax);

  double b = D[0];

  x[0] = rhs[0] / b;
  for (unsigned int k = 1; k < n; ++k) {
    work[k] = U[k - 1] / b;

    b = D[k] - L[k] * work[k];

    if (b == 0.0) {
      throw RuntimeError::formatted("zero pivot at row %d", k + 1);
    }

    x[k] = (rhs[k] - L[k] * x[k-1]) / b;
  }

  for (int k = n - 2; k >= 0; --k)
    x[k] -= work[k + 1] * x[k + 1];

  m_indicesValid = false;
}


//! Write system matrix and right-hand-side into an m-file.  The file name contains ZERO_PIVOT_ERROR.
void columnSystemCtx::reportColumnZeroPivotErrorMFile(unsigned int M) {
  char filename[TEMPORARY_STRING_LENGTH];
  snprintf(filename, sizeof(filename), "%s_i%d_j%d_ZERO_PIVOT_ERROR.m",
           m_prefix.c_str(), m_i, m_j);

  std::ofstream output(filename);
  output << "% system has 1-norm = " << norm1(M)
         << " and diagonal-dominance ratio = " << ddratio(M) << std::endl;

  viewSystem(output, M);
}


//! Write system matrix, right-hand-side, and (provided) solution into an m-file.  Constructs file name from m_prefix.
/*!
An example of the use of this procedure is from <c>examples/searise-greenland/</c>
running the enthalpy formulation.  First run spinup.sh in that directory  (FIXME:
which was modified to have equal spacing in z, when I did this example) to
generate `g20km_steady.nc`.  Then:

\code
  $ pismr -calving ocean_kill -e 3 -atmosphere searise_greenland -surface pdd -config_override  config_269.0_0.001_0.80_-0.500_9.7440.nc \
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
void columnSystemCtx::viewColumnInfoMFile(const std::vector<double> &x,
                                          unsigned int M) {
  std::cout << "saving "
            << m_prefix << " column system at (i,j)"
            << " = (" << m_i << "," << m_j << ") to m-file...\n" << std::endl;

  char buffer[TEMPORARY_STRING_LENGTH];
  snprintf(buffer, sizeof(buffer), "%s_i%d_j%d.m", m_prefix.c_str(), m_i, m_j);

  viewColumnInfoMFile(buffer, M, x);
}


//! Write system matrix, right-hand-side, and (provided) solution into an already-named m-file.
void columnSystemCtx::viewColumnInfoMFile(const std::string &filename,
                                          unsigned int M,
                                          const std::vector<double> &x) {
  std::ofstream output(filename);
  output << "% system has 1-norm = " << norm1(M)
         << " and diagonal-dominance ratio = " << ddratio(M) << std::endl;

  viewSystem(output, M);
  viewVectorValues(output, x, M, m_prefix + "_x");
}

} // end of namespace pism
