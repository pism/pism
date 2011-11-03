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

#include <petscdmmg.h>

typedef struct _p_THI   *THI;

typedef struct {
  PetscScalar b;                /* bed */
  PetscScalar h;                /* thickness */
  PetscScalar beta2;            /* friction */
} PrmNode;

typedef struct {
  PetscScalar u,v;
} Node;


typedef enum {THIASSEMBLY_TRIDIAGONAL,THIASSEMBLY_FULL} THIAssemblyMode;

PetscErrorCode THICreate(MPI_Comm comm, THI *inthi);
PetscErrorCode THISetup(MPI_Comm comm, DA pism_da2, PetscReal Lx, PetscReal Ly, THI thi, DMMG **dmmg);
PetscErrorCode THIDestroy(THI thi);

PetscErrorCode DAGetPrmNodeArray(DA da, PrmNode ***prm);
PetscErrorCode DARestorePrmNodeArray(DA da, PrmNode ***prm);

PetscErrorCode DAPrmNodeArrayCommBegin(DA da);
PetscErrorCode DAPrmNodeArrayCommEnd(DA da);

#endif /* _THI_H_ */
