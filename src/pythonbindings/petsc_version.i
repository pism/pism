// Copyright (C) 2011, 2012 David Maxwell
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



// This python module only contains the name of PETSC_DIR and PETSC_ARCH
// that PISM was compiled against.  These values are used to initialize
// petsc4py with the correct PETSc version when issuing an "import PISM"
// from within python.  We keep these values in a module separate
// from the main python bindings module 'cpp' because loading 'cpp' prior
// to initializing petsc4py was causing crashes upon exiting from python.
// I don't know what the cause was, but presumably some static code in 'cpp'
// was using PETSc, which had not yet been initialized. Now we can safely
// load this module, 'petsc_version' to initialize petsc4py, and then
// bring in the rest of 'cpp'.

%module petsc_version

%{
#include "petsc_version.hh"
%}

%immutable PISM_PETSC_ARCH;
%immutable PISM_PETSC_DIR;

%include "petsc_version.hh"
