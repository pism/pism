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

#endif	/* __columnSystem_hh */

