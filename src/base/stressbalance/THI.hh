/* Copyright (C) 2010, 2011 Jed Brown and the PISM Authors */

/* This file is part of PISM. */

/* PISM is free software; you can redistribute it and/or modify it under the */
/* terms of the GNU General Public License as published by the Free Software */
/* Foundation; either version 2 of the License, or (at your option) any later */
/* version. */

/* PISM is distributed in the hope that it will be useful, but WITHOUT ANY */
/* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS */
/* FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more */
/* details. */

/* You should have received a copy of the GNU General Public License */
/* along with PISM; if not, write to the Free Software */
/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA */

#ifndef _THI_H_
#define _THI_H_

#include "THItools.hh"

typedef struct _p_THI   *THI;

typedef enum {THIASSEMBLY_TRIDIAGONAL,THIASSEMBLY_FULL} THIAssemblyMode;

struct _p_THI {
  PETSCHEADER(int);
  void (*initialize)(THI,PetscReal x,PetscReal y,PrmNode *p);
  PetscInt  nlevels;
  PetscInt  zlevels;
  PetscReal Lx,Ly,Lz;           /* Model domain */
  PetscReal alpha;              /* Bed angle */
  Units     units;
  PetscReal dirichlet_scale;
  PetscReal ssa_friction_scale;
  PRange    eta;
  PRange    beta2;
  struct {
    PetscReal Bd2,eps,exponent;
  } viscosity;
  struct {
    PetscReal irefgam,eps2,exponent;
  } friction;
  PetscReal rhog;
  PetscTruth no_slip;
  PetscTruth tridiagonal;
  PetscTruth coarse2d;
  PetscTruth verbose;
  MatType mattype;
};

PetscErrorCode THICreate(MPI_Comm comm, THI *inthi);
PetscErrorCode THIDestroy(THI thi);

PetscErrorCode DARefineHierarchy_THI(DA dac0,PetscInt nlevels,DA hierarchy[]);
PetscErrorCode DAGetInterpolation_THI(DA dac,DA daf,Mat *A,Vec *scale);
PetscErrorCode DAGetMatrix_THI_Tridiagonal(DA da,const MatType mtype,Mat *J);
PetscErrorCode THIFunctionLocal(DALocalInfo *info,Node ***x,Node ***f,THI thi);
PetscErrorCode THIJacobianLocal_2D(DALocalInfo *info,Node **x,Mat B,THI thi);
PetscErrorCode THIJacobianLocal_3D(DALocalInfo *info,Node ***x,Mat B,THI thi,THIAssemblyMode amode);
PetscErrorCode THIJacobianLocal_3D_Full(DALocalInfo *info,Node ***x,Mat B,THI thi);
PetscErrorCode THIJacobianLocal_3D_Tridiagonal(DALocalInfo *info,Node ***x,Mat B,THI thi);
PetscErrorCode THISetDMMG(THI thi,DMMG *dmmg);
PetscErrorCode THIInitial(DMMG dmmg,Vec X);
PetscErrorCode THISolveStatistics(THI thi,DMMG *dmmg,PetscInt coarsened,const char name[]);

#endif /* _THI_H_ */
