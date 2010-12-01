// Copyright (C) 2010 Constantine Khroulev
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

#include "iceModel.hh"
#include "PISMDiagnostic.hh"


PetscErrorCode IceModel::init_diagnostics() {

  diagnostics["hardav"] = new IceModel_hardav(this, grid, variables);
  diagnostics["rank"]   = new IceModel_rank(this, grid, variables);

  stress_balance->get_diagnostics(diagnostics);

  map<string, PISMDiagnostic*>::iterator j = diagnostics.begin();
  while (j != diagnostics.end()) {
    string name = j->first;
    PISMDiagnostic *diag = j->second;

    int N = diag->get_nvars();
    verbPrintf(2, grid.com, " ** %s [%d variable(s)]\n", name.c_str(), N);

    for (int k = 0; k < N; ++k) {
      NCSpatialVariable *var = diag->get_metadata(k);

      string long_name = var->get_string("long_name");

      verbPrintf(2, grid.com, " * %s\n", long_name.c_str());

      ++j;
    }
  }

  return 0;
}

IceModel_hardav::IceModel_hardav(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {
  
  // set metadata:
  vars[0].init("hardav", grid, GRID_2D);

  const PetscScalar power = 1.0 / model->ice->exponent();
  char unitstr[TEMPORARY_STRING_LENGTH];
  snprintf(unitstr, sizeof(unitstr), "Pa s%f", power);

  set_attrs("vertical average of ice hardness", "",
            unitstr, unitstr, 0);

  vars[0].set("valid_min", 0);
  vars[0].set("_FillValue", -0.01);
}

//! \brief Computes vertically-averaged ice hardness.
PetscErrorCode IceModel_hardav::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  const PetscScalar fillval = -0.01;
  PetscScalar *Eij; // columns of enthalpy values

  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "hardav", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  ierr = model->Enth3.begin_access(); CHKERRQ(ierr);
  ierr = model->vH.begin_access(); CHKERRQ(ierr);
  ierr = result->begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = model->Enth3.getInternalColumn(i,j,&Eij); CHKERRQ(ierr);
      const PetscScalar H = model->vH(i,j);
      if (H > 0.0) {
        (*result)(i,j) = model->ice->averagedHardness_from_enth(H, grid.kBelowHeight(H),
                                                             grid.zlevels, Eij);
      } else { // put negative value below valid range
        (*result)(i,j) = fillval;
      }
    }
  }
  ierr = model->Enth3.end_access(); CHKERRQ(ierr);
  ierr = model->vH.end_access(); CHKERRQ(ierr);
  ierr = result->end_access(); CHKERRQ(ierr);

  output = result;
  return 0;
}


IceModel_rank::IceModel_rank(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {
  
  // set metadata:
  vars[0].init("rank", grid, GRID_2D);
  
  set_attrs("processor rank", "", "", "", 0);
  vars[0].time_independent = true;
}

PetscErrorCode IceModel_rank::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "rank", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  ierr = result->begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i)
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j)
      (*result)(i,j) = grid.rank;
  ierr = result->end_access();

  output = result;
  return 0;
}
