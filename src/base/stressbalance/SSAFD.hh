// Copyright (C) 2004--2011 Jed Brown, Ed Bueler and Constantine Khroulev
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

#ifndef _SSAFD_H_
#define _SSAFD_H_

#include "SSA.hh"
#include <petscksp.h>


//! PISM's SSA solver: the finite difference implementation
class SSAFD : public SSA
{
public:
  SSAFD(IceGrid &g, IceBasalResistancePlasticLaw &b, IceFlowLaw &i, EnthalpyConverter &e,
        const NCConfigVariable &c) :
    SSA(g,b,i,e,c)
  {
    PetscErrorCode ierr = allocate_fd();
    if (ierr != 0) {
      PetscPrintf(grid.com, "FATAL ERROR: SSAFD allocation failed.\n");
      PISMEnd();
    }
  }

  virtual ~SSAFD()
  {
    PetscErrorCode ierr = deallocate_fd();
    if (ierr != 0) {
      PetscPrintf(grid.com, "FATAL ERROR: SSAFD de-allocation failed.\n");
      PISMEnd();
    }
  }

  virtual PetscErrorCode init(PISMVars &vars);

protected:
  virtual PetscErrorCode allocate_fd();

  virtual PetscErrorCode deallocate_fd();

  virtual PetscErrorCode solve();

  virtual PetscErrorCode compute_hardav_staggered(IceModelVec2Stag &result);

  virtual PetscErrorCode compute_nuH_staggered(IceModelVec2Stag &result,
                                               PetscReal epsilon);

  virtual PetscErrorCode compute_nuH_norm(PetscReal &norm,
                                          PetscReal &norm_change);

  virtual PetscErrorCode assemble_matrix(bool include_basal_shear, Mat A);

  virtual PetscErrorCode assemble_rhs(Vec rhs);

  virtual PetscErrorCode writeSSAsystemMatlab();

  virtual PetscErrorCode update_nuH_viewers();

  virtual PetscErrorCode set_diagonal_matrix_entry(Mat A, int i, int j,
                                                   PetscScalar value);

  virtual bool is_marginal(int i, int j);

  // objects used internally
  IceModelVec2Stag hardness, nuH, nuH_old;
  KSP SSAKSP;
  Mat SSAStiffnessMatrix;
  Vec SSARHS;
  bool use_cfbc;
  PetscScalar scaling;

  bool view_nuh;
  PetscViewer nuh_viewer;
  PetscInt nuh_viewer_size;
};

//! Constructs a new SSAFD
SSA * SSAFDFactory(IceGrid &, IceBasalResistancePlasticLaw &, 
                  IceFlowLaw &, EnthalpyConverter &, const NCConfigVariable &);

#endif /* _SSAFD_H_ */

