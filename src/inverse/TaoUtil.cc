// Copyright (C) 2012  David Maxwell
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

#include "TaoUtil.hh"

TaoInitializer::TaoInitializer(int *argc, char ***argv, char *file, char *help) {
  PetscErrorCode ierr = TaoInitialize(argc,argv,file,help);
  CHKERRCONTINUE(ierr);
  if(ierr) {
    PetscPrintf(PETSC_COMM_WORLD,"FATAL ERROR: Unable to initialize TAO.");
    PISMEnd();
  }
}

TaoInitializer::TaoInitializer(int *argc, char ***argv, char *help) {
  PetscErrorCode ierr = TaoInitialize(argc,argv,NULL,help);
  CHKERRCONTINUE(ierr);
  if(ierr) {
    PetscPrintf(PETSC_COMM_WORLD,"FATAL ERROR: Unable to initialize TAO.");
    PISMEnd();
  }
}

TaoInitializer::TaoInitializer(int *argc, char ***argv) {
  PetscErrorCode ierr = TaoInitialize(argc,argv,NULL,NULL);
  CHKERRCONTINUE(ierr);
  if(ierr) {
    PetscPrintf(PETSC_COMM_WORLD,"FATAL ERROR: Unable to initialize TAO.");
    PISMEnd();
  }
}

TaoInitializer::~TaoInitializer() {
  PetscErrorCode ierr = TaoFinalize();
  CHKERRCONTINUE(ierr);
  if(ierr) {
    PetscPrintf(PETSC_COMM_WORLD,"FATAL ERROR: Unable to finalize TAO.");
    PISMEnd();
  }
}

// TAO failure codes are negative and success codes are positive.  We declare 
// associated description strings as an array and then point to the middle 
// of the array so that we can look up description strings via the error code.
const char *TaoConvergedReasonsShifted[] = {
    " ", " ", 
    "TAO_DIVERGED_USER",
    "TAO_DIVERGED_TR_REDUCTION",
    "TAO_DIVERGED_LS_FAILURE",
    "TAO_DIVERGED_MAXFCN",
    "TAO_DIVERGED_NAN",
    " ",
    "TAO_DIVERGED_MAXITS",
    " ",
    "TAO_CONTINUE_ITERATING",
    "TAO_CONVERGED_FATOL",
    "TAO_CONVERGED_FRTOL",
    "TAO_CONVERGED_GATOL",
    "TAO_CONVERGED_GRTOL",
    "TAO_CONVERGED_GTTOL",
    "TAO_CONVERGED_STEPTOL",
    "TAO_CONVERGED_MINF",
    "TAO_CONVERGED_USER" };
const char *const* TaoConvergedReasons = TaoConvergedReasonsShifted + 10;

TAOTerminationReason::TAOTerminationReason( TaoSolverTerminationReason r)  {
  m_reason = r;
}
void TAOTerminationReason::get_description( std::ostream &desc, int indent_level) {
  for( int i=0; i < indent_level; i++) {
    desc << sm_indent;
  }
  desc << TaoConvergedReasons[m_reason];
}
