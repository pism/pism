// Copyright (C) 2009-2011, 2013, 2014 Ed Bueler
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
#include <petsc.h>

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

    [b_1  c_1   0   ...                           ] [  u_1  ]   [   r_1   ]
    [a_2  b_2  c_2  ...                           ] [  u_2  ]   [   r_2   ]
    [               ...                           ].[  ...  ] = [   ...   ]
    [               ...  a_{N-1}  b_{N-1}  c_{N-1}] [u_{N-1}]   [ r_{N-1} ]
    [               ...     0     a_N      b_N    ] [  u_N  ]   [   r_N   ]

HOWEVER... the version in this code is different from Numerical
Recipes in two ways:
   a) Indexing is zero-based
   b) Variables have been renamed.

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


Therefore... this version of the code solves the following problem:

    [D_0  U_0   0   ...                           ] [  x_0  ]   [   r_0   ]
    [L_1  D_1  U_1  ...                           ] [  x_1  ]   [   r_1   ]
    [ 0   L_2  D_2  ...                           ].[  ...  ] = [   ...   ]
    [               ...  L_{N-2}  D_{N-2}  U_{N-2}] [x_{N-2}]   [ r_{N-2} ]
    [               ...     0     L_{N-1}  D_{N-1}] [x_{N-1}]   [ r_{N-1} ]

 */
class columnSystemCtx {

public:
  columnSystemCtx(unsigned int my_nmax, std::string my_prefix);
  virtual ~columnSystemCtx();

  PetscErrorCode setIndicesAndClearThisColumn(int i, int j, int ks);  

  double    norm1(unsigned int n) const;
  double    ddratio(unsigned int n) const;

  PetscErrorCode viewVectorValues(PetscViewer viewer,
                                  const double *v, int m, const char* info) const;
  PetscErrorCode viewMatrix(PetscViewer viewer, const char* info) const;
  virtual PetscErrorCode viewSystem(PetscViewer viewer) const;

  PetscErrorCode reportColumnZeroPivotErrorMFile(const PetscErrorCode errindex);
  PetscErrorCode viewColumnInfoMFile(double *x, unsigned int n);
  PetscErrorCode viewColumnInfoMFile(char *filename, double *x, unsigned int n);

protected:
  unsigned int nmax;
  double *L, *Lp, *D, *U, *rhs, *work; // vectors for tridiagonal system

  int    i, j, ks;

  // deliberately protected so only derived classes can use
  PetscErrorCode solveTridiagonalSystem(unsigned int n, double *x);
  
  std::string      prefix;
private:
  bool        indicesValid;
  PetscErrorCode resetColumn();
};

} // end of namespace pism

#endif  /* __columnSystem_hh */

