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

#ifndef _SSAFD_H_
#define _SSAFD_H_

#include "ShallowStressBalance.hh"

//! PISM's SSA solver implementation
class SSAFD : public ShallowStressBalance
{
public:
  SSAFD(IceGrid &g, const NCConfigVariable &config);
  virtual ~SSAFD();
protected:

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode allocate();
  virtual PetscErrorCode deallocate();

  virtual PetscErrorCode update(bool fast);

  virtual PetscErrorCode solve();

  virtual PetscErrorCode compute_nuH_staggered(IceModelVec2Stag &result,
                                               PetscReal epsilon); // done

  virtual PetscErrorCode compute_nuH_norm(IceModelVec2Stag nuH,
                                          IceModelVec2Stag nuH_old,
                                          PetscReal *norm,
                                          PetscReal *norm_change); // done

  virtual PetscErrorCode assemble_matrix(bool include_basal_shear,
                                         IceModelVec2Stag nuH, Mat A);

  virtual PetscErrorCode assemble_rhs(Vec rhs); // done

protected:
  virtual PetscErrorCode compute_driving_stress(IceModelVec2S &taudx,
                                                IceModelVec2S &taudy); // done
  virtual PetscErrorCode compute_hardav_staggered(IceModelVec2Stag &result); // done

  virtual PetscErrorCode compute_basal_frictional_heating(IceModelVec2S &result); // done

  virtual PetscErrorCode compute_D2(IceModelVec2S &result); // done

  IceModelVec2Mask *mask;
  IceModelVec2S *thickness;
  IceModelVec3 *enthalpy;
};


#endif /* _SSAFD_H_ */
