// Copyright (C) 2004-2010 Jed Brown, Ed Bueler and Constantine Khroulev
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

#ifndef __iceExactSSAModel_hh
#define __iceExactSSAModel_hh

#include <petscda.h>
#include "grid.hh"
#include "materials.hh"
#include "iceModel.hh"
#include "SSAFD.hh"

class IceExactSSAModel : public IceModel {
public:
  IceExactSSAModel(IceGrid &g, NCConfigVariable &config, NCConfigVariable &overrides, char mytest);
  virtual PetscErrorCode setFromOptions();

  virtual PetscErrorCode misc_setup();
  virtual PetscErrorCode set_grid_defaults();
  virtual PetscErrorCode set_vars_from_options();
  virtual PetscErrorCode init_physics();

  virtual PetscErrorCode run();
  PetscErrorCode         reportErrors();

protected:
    char            test;       // only 'I', 'J', 'M' supported
    PetscTruth      exactOnly;
  
    PetscErrorCode  fillFromExactSolution();
    PetscErrorCode  taucSetI();
    PetscErrorCode  setInitStateAndBoundaryVelsI();
    PetscErrorCode  setInitStateJ();
    PetscErrorCode  setInitStateM();

private:
    // constants for I; "static" o.k. because also const
    static const PetscScalar   m_schoof, L_schoof, aspect_schoof, H0_schoof,
                               B_schoof, p_schoof, DEFAULT_PLASTIC_REGULARIZE;
    // constants for J, M
    static const PetscScalar   LforJ, LforM;

  // we need a pointer to the SSA solver to adjust strength extension parameters
  SSAFD *ssa;
  SSBM_Trivial *modifier;
};

#endif /* __iceExactSSAModel_hh */

