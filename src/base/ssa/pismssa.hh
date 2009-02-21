// Copyright (C) 2009 Jed Brown and Ed Bueler
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

#ifndef _PISM_SSA_HH
#define _PISM_SSA_HH

typedef struct _p_PetscViewer *PetscViewer;
class IceGrid;
class IceModelVec2;
class IceModelVec3;
class IceType;

typedef struct _p_SSA *SSA;

#define SSAType char *

#define SSAFE "fe"

extern PetscCookie SSA_COOKIE;

#if defined(PETSC_USE_DYNAMIC_LIBRARIES)
# define SSARegisterDynamic(s,p,n,f) SSARegister(s,p,n,0)
#else
# define SSARegisterDynamic(s,p,n,f) SSARegister(s,p,n,f)
#endif

extern PetscErrorCode SSACreate(IceGrid*,SSA*);
extern PetscErrorCode SSADestroy(SSA);
extern PetscErrorCode SSASolve(SSA,IceModelVec2 &ubar,IceModelVec2 &vbar);
extern PetscErrorCode SSAView(SSA,PetscViewer);
extern PetscErrorCode SSASetFromOptions(SSA);
extern PetscErrorCode SSASetUp(SSA);
extern PetscErrorCode SSASetType(SSA,const SSAType);
extern PetscErrorCode SSASetIceType(SSA,IceType*);
extern PetscErrorCode SSASetBasalType(SSA,PlasticBasalType*);
extern PetscErrorCode SSASetOceanType(SSA,SeaWaterType*);
extern PetscErrorCode SSASetFields(SSA,IceModelVec2 *mask,IceModelVec2 *sia_uvbar,IceModelVec2 *H,IceModelVec2 *h,IceModelVec2 *tauc,IceModelVec3 *T);
extern PetscErrorCode SSASetFictitiousNuH(SSA,PetscReal);
extern PetscErrorCode SSASetCutoffThickness(SSA,PetscReal);
extern PetscErrorCode SSARegister(const char sname[],const char path[],const char name[],PetscErrorCode (*function)(SSA));
extern PetscErrorCode SSAInitializePackage(const char path[]);
extern PetscErrorCode SSAMapToSplitVecs(SSA ssa,IceModelVec2 &ubar,IceModelVec2 &vbar);

#endif /* _PISM_SSA_HH */

