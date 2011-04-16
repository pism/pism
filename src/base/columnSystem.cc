// Copyright (C) 2004-2011 Jed Brown, Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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

#include <petsc.h>
#include "pism_const.hh"   // e.g. MASK_FLOATING
#include "iceModelVec.hh"
#include "columnSystem.hh"


//! Allocate a tridiagonal system of maximum size my_nmax.
/*!
Let N = \c nmax.  Then allocated locations are like this:
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

Note L[0] is not allocated and U[N-1] is not allocated.
 */
columnSystemCtx::columnSystemCtx(PetscInt my_nmax, string my_prefix) : nmax(my_nmax), prefix(my_prefix) {
  if (nmax < 1) {
    PetscPrintf(PETSC_COMM_WORLD,
      "columnSystemCtx ERROR: nmax of system too small\n");
    PISMEnd();
  }
  if (nmax > 1000000) {
    PetscPrintf(PETSC_COMM_WORLD,
      "columnSystemCtx ERROR: nmax of system unreasonable (> 10^6)\n");
    PISMEnd();
  }
  Lp   = new PetscScalar[nmax-1];
  L    = Lp-1; // ptr arith.; note L[0]=Lp[-1] not allocated
  D    = new PetscScalar[nmax];
  U    = new PetscScalar[nmax-1];
  rhs  = new PetscScalar[nmax];
  work = new PetscScalar[nmax];

  resetColumn();

  indicesValid = false;

}


columnSystemCtx::~columnSystemCtx() {
  delete [] Lp;
  delete [] D;
  delete [] U;
  delete [] rhs;
  delete [] work;
}


//! Zero all entries.
PetscErrorCode columnSystemCtx::resetColumn() {
  PetscErrorCode ierr;
  ierr = PetscMemzero(Lp,   (nmax-1)*sizeof(PetscScalar)); CHKERRQ(ierr);
  ierr = PetscMemzero(U,    (nmax-1)*sizeof(PetscScalar)); CHKERRQ(ierr);
  ierr = PetscMemzero(D,    (nmax)*sizeof(PetscScalar)); CHKERRQ(ierr);
  ierr = PetscMemzero(rhs,  (nmax)*sizeof(PetscScalar)); CHKERRQ(ierr);
  ierr = PetscMemzero(work, (nmax)*sizeof(PetscScalar)); CHKERRQ(ierr);
  return 0;
}


//! Compute 1-norm, which is max sum of absolute values of columns.
PetscScalar columnSystemCtx::norm1(const PetscInt n) const {
  if (n > nmax) {
    PetscPrintf(PETSC_COMM_WORLD,"PISM ERROR:  n > nmax in columnSystemCtx::norm1()\n");
    PISMEnd();
  }
  if (n == 1)  return fabs(D[0]);   // only 1x1 case is special
  PetscScalar z = fabs(D[0]) + fabs(L[1]);
  for (PetscInt k = 1; k < n; k++) {  // k is column index (zero-based)
    z = PetscMax(z, fabs(U[k-1])) + fabs(D[k]) + fabs(L[k+1]);
  }
  z = PetscMax(z, fabs(U[n-2]) + fabs(D[n-1]));
  return z;
}


//! Compute diagonal-dominance ratio.  If this is less than one then the matrix is strictly diagonally-dominant.
/*!
Let \f$A = (a_{ij})\f$ be the tridiagonal matrix
described by L,D,U for row indices 0 through \c n.  The computed ratio is 
  \f[ \max_{j=1,\dots,n} \frac{|a_{j,j-1}|+|a_{j,j+1}|}{|a_{jj}|}, \f]
where \f$a_{1,0}\f$ and \f$a_{n,n+1}\f$ are interpreted as zero.

If this is smaller than one then it is a theorem that the tridiagonal solve will
succeed.

We return -1.0 if the absolute value of any diagonal element is less than
1e-12 of the 1-norm of the matrix.
 */
PetscScalar columnSystemCtx::ddratio(const PetscInt n) const {
  if (n > nmax) {
    PetscPrintf(PETSC_COMM_WORLD,"PISM ERROR:  n > nmax in columnSystemCtx::ddratio()\n");
    PISMEnd();
  }
  const PetscScalar scale = norm1(n);

  if ( (fabs(D[0]) / scale) < 1.0e-12)  return -1.0;
  PetscScalar z = fabs(U[0]) / fabs(D[0]);

  for (PetscInt k = 1; k < n-1; k++) {  // k is row index (zero-based)
    if ( (fabs(D[k]) / scale) < 1.0e-12)  return -1.0;
    const PetscScalar s = fabs(L[k]) + fabs(U[k]);
    z = PetscMax(z, s / fabs(D[k]) );
  } 

  if ( (fabs(D[n-1]) / scale) < 1.0e-12)  return -1.0;
  z = PetscMax(z, fabs(L[n-1]) / fabs(D[n-1]) );

  return z;
}


PetscErrorCode columnSystemCtx::setIndicesAndClearThisColumn(
                  PetscInt my_i, PetscInt my_j, PetscInt my_ks) {
  if (indicesValid) {  SETERRQ(3,
     "setIndicesAndClearThisColumn() called twice in same column"); }
  i = my_i;
  j = my_j;
  ks = my_ks;

  resetColumn();
  
  indicesValid = true;
  return 0;
}


//! Utility for simple ascii view of a vector (one-dimensional column) quantity.
/*!
Give first argument NULL to get standard out.  No binary viewer.

Give description string as \c info argument.

Result should be executable as part of a Matlab/Octave script.
 */
PetscErrorCode columnSystemCtx::viewVectorValues(PetscViewer viewer,
                                                 PetscScalar *v, PetscInt m, const char* info) const {
  PetscErrorCode ierr;

  if (v==NULL) {
    SETERRQ1(2,"columnSystem ERROR: can't view '%s' by v=NULL pointer ... ending ...\n", info);
  }
  if (m<1) {
    SETERRQ1(3,"columnSystem ERROR: can't view '%s' because m<1 ... ending ...\n",info);
  }

  PetscTruth iascii;
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PETSC_COMM_SELF,&viewer); CHKERRQ(ierr);
  }
  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii); CHKERRQ(ierr);
  if (!iascii) { SETERRQ(1,"Only ASCII viewer for ColumnSystem\n"); }

  ierr = PetscViewerASCIIPrintf(viewer,
     "\n%% viewing ColumnSystem column object with description '%s' (columns  [k value])\n",
     info); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
      "%s_with_index = [...\n",info); CHKERRQ(ierr);
  for (PetscInt k=0; k<m; k++) {
    ierr = PetscViewerASCIIPrintf(viewer,
      "  %5d %.12f",k,v[k]); CHKERRQ(ierr);
    if (k == m-1) { 
      ierr = PetscViewerASCIIPrintf(viewer, "];\n"); CHKERRQ(ierr);
    } else {
      ierr = PetscViewerASCIIPrintf(viewer, ";\n"); CHKERRQ(ierr);
    }
  }
  ierr = PetscViewerASCIIPrintf(viewer,
      "%s = %s_with_index(:,2);\n\n",info,info); CHKERRQ(ierr);
  return 0;
}


//! View the tridiagonal matrix.  Views as a full matrix if nmax <= 120, otherwise by listing diagonals.
/*! 
Give first argument NULL to get standard out.  No binary viewer.

Give description string as \c info argument.
 */
PetscErrorCode columnSystemCtx::viewMatrix(PetscViewer viewer, const char* info) const {
  PetscErrorCode ierr;
  PetscTruth iascii;
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PETSC_COMM_SELF,&viewer); CHKERRQ(ierr);
  }
  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii); CHKERRQ(ierr);
  if (!iascii) { SETERRQ(1,"Only ASCII viewer for ColumnSystem\n"); }

  if (L==NULL) {
    SETERRQ1(2,"columnSystemCtx ERROR: can't matrix '%s' because L==NULL ... ending ...\n", info);
  }
  if (D==NULL) {
    SETERRQ1(3,"columnSystemCtx ERROR: can't matrix '%s' because D==NULL ... ending ...\n", info);
  }
  if (U==NULL) {
    SETERRQ1(4,"columnSystemCtx ERROR: can't matrix '%s' because U==NULL ... ending ...\n", info);
  }

  if (nmax < 2) {
    ierr = PetscViewerASCIIPrintf(viewer,
      "\n\n<nmax >= 2 required to view columnSystemCtx tridiagonal matrix '%s' ... skipping view\n",info);
    CHKERRQ(ierr);
    return 0;
  }

  if (nmax > 500) {
    ierr = PetscViewerASCIIPrintf(viewer,
      "\n\n<nmax > 500: columnSystemCtx matrix too big to display as full; viewing tridiagonal matrix '%s' by diagonals ...\n",info); CHKERRQ(ierr);
    char vinfo[PETSC_MAX_PATH_LEN];
    snprintf(vinfo,PETSC_MAX_PATH_LEN, "%s_super_diagonal_U", info);
    ierr = viewVectorValues(viewer,U,nmax-1,vinfo); CHKERRQ(ierr);
    snprintf(vinfo,PETSC_MAX_PATH_LEN, "%s_diagonal_D", info);
    ierr = viewVectorValues(viewer,D,nmax,vinfo); CHKERRQ(ierr);
    snprintf(vinfo,PETSC_MAX_PATH_LEN, "%s_sub_diagonal_L", info);
    ierr = viewVectorValues(viewer,Lp,nmax-1,vinfo); CHKERRQ(ierr);
  } else {
    ierr = PetscViewerASCIIPrintf(viewer,
        "\n%s = [...\n",info); CHKERRQ(ierr);
    for (PetscInt k=0; k<nmax; k++) {    // k+1 is row  (while j+1 is column)
      if (k == 0) {              // viewing first row
        ierr = PetscViewerASCIIPrintf(viewer,"%.12f %.12f ",D[k],U[k]); CHKERRQ(ierr);
        for (PetscInt n=2; n<nmax; n++) {
          ierr = PetscViewerASCIIPrintf(viewer,"%3.1f ",0.0); CHKERRQ(ierr);
        }
      } else if (k < nmax-1) {   // viewing generic row
        for (PetscInt n=0; n<k-1; n++) {
          ierr = PetscViewerASCIIPrintf(viewer,"%3.1f ",0.0); CHKERRQ(ierr);
        }
        ierr = PetscViewerASCIIPrintf(viewer,"%.12f %.12f %.12f ",L[k],D[k],U[k]); CHKERRQ(ierr);
        for (PetscInt n=k+2; n<nmax; n++) {
          ierr = PetscViewerASCIIPrintf(viewer,"%3.1f ",0.0); CHKERRQ(ierr);
        }
      } else {                   // viewing last row
        for (PetscInt n=0; n<k-1; n++) {
          ierr = PetscViewerASCIIPrintf(viewer,"%3.1f ",0.0); CHKERRQ(ierr);
        }
        ierr = PetscViewerASCIIPrintf(viewer,"%.12f %.12f ",L[k],D[k]); CHKERRQ(ierr);
      }
      
      if (k == nmax-1) {
        ierr = PetscViewerASCIIPrintf(viewer,"];\n\n"); CHKERRQ(ierr);  // end final row
      } else {
        ierr = PetscViewerASCIIPrintf(viewer,";\n"); CHKERRQ(ierr);  // end of generic row
      }
    }
  }
  
  return 0;
}


//! View the tridiagonal system A x = b to a PETSc viewer, both A as a full matrix and b as a vector.
PetscErrorCode columnSystemCtx::viewSystem(PetscViewer viewer) const {
  PetscErrorCode ierr;
  char  info[PETSC_MAX_PATH_LEN];
  snprintf(info,PETSC_MAX_PATH_LEN, "%s_A", prefix.c_str());
  ierr = viewMatrix(viewer,info); CHKERRQ(ierr);
  snprintf(info,PETSC_MAX_PATH_LEN, "%s_rhs", prefix.c_str());
  ierr = viewVectorValues(viewer,rhs,nmax,info); CHKERRQ(ierr);
  return 0;
}


//! The actual code for solving a tridiagonal system.  Return code has diagnostic importance.
/*!
This is modified slightly from a Numerical Recipes version.

Input size n is size of instance.  Requires n <= columnSystemCtx::nmax.

Solution of system in x.

Success is return code zero.  Positive return code gives location of zero pivot.
Negative return code indicates a software problem.
 */
PetscErrorCode columnSystemCtx::solveTridiagonalSystem(
                  const PetscInt n, PetscScalar **x) {
#ifdef PISM_DEBUG
  if (x == NULL) { SETERRQ(-999,"x is NULL in columnSystemCtx"); }
  if (*x == NULL) { SETERRQ(-998,"*x is NULL in columnSystemCtx"); }
  if (n < 1) { SETERRQ(-997,"instance size n < 1 in columnSystemCtx"); }
  if (n > nmax) { SETERRQ(-996,"instance size n too large in columnSystemCtx"); }

  if (!indicesValid) { SETERRQ(-995,"column indices not valid in columnSystemCtx"); }
#endif
  PetscScalar b;
  b = D[0];
  if (b == 0.0) { return 1; }
  (*x)[0] = rhs[0]/b;
  for (int k=1; k<n; ++k) {
    work[k] = U[k-1]/b;
    b = D[k] - L[k] * work[k];
    if (b == 0.0) { return k+1; }
    (*x)[k] = (rhs[k] - L[k] * (*x)[k-1]) / b;
  }
  for (int k=n-2; k>=0; --k) {
    (*x)[k] -= work[k+1] * (*x)[k+1];
  }

  indicesValid = false;
  return 0;
}


//! Write system matrix and right-hand-side into an m-file.  The file name contains ZERO_PIVOT_ERROR.
PetscErrorCode columnSystemCtx::reportColumnZeroPivotErrorMFile(const PetscErrorCode errindex) {
  PetscErrorCode ierr;
  char fname[PETSC_MAX_PATH_LEN];
  snprintf(fname, PETSC_MAX_PATH_LEN, "%s_i%d_j%d_ZERO_PIVOT_ERROR_%d.m",
           prefix.c_str(),i,j,errindex);
  ierr = viewColumnInfoMFile(fname, NULL, 0); CHKERRQ(ierr);
  return 0;
}


//! Write system matrix, right-hand-side, and (provided) solution into an m-file.  Constructs file name from prefix.
/*! 
An example of the use of this procedure is from <c>examples/searise-greenland/</c>
running the enthalpy formulation.  First run spinup.sh in that directory  (FIXME:
which was modified to have equal spacing in z, when I did this example) to
generate \c g20km_steady.nc.  Then:

\code
  $ pismr -ocean_kill -e 3 -atmosphere searise_greenland -surface pdd -config_override  config_269.0_0.001_0.80_-0.500_9.7440.nc \
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

Of course we can also do \c spy(A), \c eig(A), and look at individual entries,
and row and column sums, and so on.
 */
PetscErrorCode columnSystemCtx::viewColumnInfoMFile(PetscScalar *x, PetscInt n) {
  PetscErrorCode ierr;
  char fname[PETSC_MAX_PATH_LEN];
  snprintf(fname, PETSC_MAX_PATH_LEN, "%s_i%d_j%d.m", prefix.c_str(), i,j);
  ierr = viewColumnInfoMFile(fname, x, n); CHKERRQ(ierr);
  return 0;
}


//! Write system matrix, right-hand-side, and (provided) solution into an already-named m-file.
/*!
Because this may be called on only one processor, it builds a viewer on MPI
communicator PETSC_COMM_SELF.
 */
PetscErrorCode columnSystemCtx::viewColumnInfoMFile(char *filename, PetscScalar *x, PetscInt n) {
  PetscErrorCode ierr;
  PetscViewer viewer;
  ierr = PetscViewerCreate(PETSC_COMM_SELF, &viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer, PETSC_VIEWER_ASCII);CHKERRQ(ierr);
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, filename);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
        "%%  system has 1-norm = %.3e  and  diagonal-dominance ratio = %.5f\n",
        norm1(n), ddratio(n)); CHKERRQ(ierr);
  ierr = viewSystem(viewer); CHKERRQ(ierr);
  if ((x != NULL) && (n > 0)) {
    char  info[PETSC_MAX_PATH_LEN];
    snprintf(info,PETSC_MAX_PATH_LEN, "%s_x", prefix.c_str());
    ierr = viewVectorValues(viewer, x, n, info); CHKERRQ(ierr);
  }
  ierr = PetscViewerDestroy(viewer); CHKERRQ(ierr);
  return 0;
}

