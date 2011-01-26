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

// The following are macros (instead of inline functions) so that error handling
// is less cluttered.  They should be replaced with empty macros when in
// optimized mode.

#include <petscmat.h>
#include "iceModelVec.hh"
#include "flowlaws.hh"

static const PetscInt numQuadPoints = 4;
static const PetscReal quadPoints[4][2] = {{ -0.57735026918962573, -0.57735026918962573 },
                                           {  0.57735026918962573, -0.57735026918962573 },
                                           {  0.57735026918962573,  0.57735026918962573 },
                                           { -0.57735026918962573,  0.57735026918962573 }};
static const PetscReal quadWeights[4]  = {1,1,1,1};
static const PetscReal quadWeights1[2] = {1,1};

#undef H
#undef L
#undef M
#undef P
#define H 0.78867513459481287
#define L 0.21132486540518708
#define M (-0.5)
#define P (0.5)
static const PetscReal interp[4*4] = {H*H,H*L,L*L,L*H,  L*H,H*H,H*L,L*L,  L*L,H*L,H*H,L*H,  H*L,L*L,L*H,H*H};
static const PetscReal derivx[4*4] = {M*H,P*H,P*L,M*L,  M*H,P*H,P*L,M*L,  M*L,P*L,P*H,M*H,  M*L,P*L,P*H,M*H};
static const PetscReal derivy[4*4] = {H*M,L*M,L*P,H*P,  L*M,H*M,H*P,L*P,  L*M,H*M,H*P,L*P,  H*M,L*M,L*P,H*P};
static const PetscReal interp1[4]  = {H,L,L,H};
static const PetscReal deriv1[4]   = {M,P,M,P};
#undef H
#undef L
#undef M
#undef P

#define PismValidVelocity(U) do {                               \
    if (!(-   1e+5 < (U).u && (U).u < 1e+5                      \
          && -1e+5 < (U).v && (U).v < 1e+5))                    \
      SETERRQ3(1,"Invalid velocity (%g,%g) not within %g",      \
               (U).u,(U).v,1e+5);                               \
  } while (0)

#define PismValidStrainRate(Du) do {                                    \
    if (!(-   1e+5 < (Du)[0] && (Du)[0] < 1e+5                          \
          && -1e+5 < (Du)[1] && (Du)[1] < 1e+5                          \
          && -1e+5 < (Du)[2] && (Du)[2] < 1e+5))                        \
      SETERRQ4(1,"Invalid Strain Rate (%g,%g,%g) not within %g",        \
               (Du)[0],(Du)[1],(Du)[2],1e+5);                           \
  } while (0)

#define PismValidStress2(f) do {                                        \
    if (!(-   1e4 < (f).u && (f).u < 1e4                                \
          && -1e4 < (f).v && (f).v < 1e4))                              \
      SETERRQ3(1,"Invalid Stress residual (%g,%g) not within %g",       \
               (f).u,(f).v,1e4);                                        \
  } while (0)

#define PismValidFriction(b) do {                                       \
    if (!(0 <= (b) && (b) < 1e25))                                      \
      SETERRQ2(1,"Invalid friction %g not within [0,%g]",(b),1e25);     \
  } while (0)


PetscTruth Floating(const IceFlowLaw &ice, PetscScalar ocean_rho,
                           PetscReal H, PetscReal bed);

void QuadZeroScalar(PetscReal x[]);

void QuadZeroVel(PISMVector2 x[]);

void QuadExtractScalar(PetscInt i,PetscInt j,PetscReal **xg,PetscReal x[]);

void QuadInsertScalar(PetscInt i,PetscInt j,PetscReal x[],PetscReal **xg);

void QuadExtractVel(PetscInt i,PetscInt j,const PISMVector2 **xg,PISMVector2 x[]);

void QuadInsertVel(const MatStencil row[],const PISMVector2 x[],PISMVector2 **xg);

void QuadMatMultScalar(const PetscReal *A,const PetscReal *x,PetscReal *y);

void QuadMatMultVel(const PetscReal *A,const PISMVector2 *x,PISMVector2 *y);

void QuadMatMultTransposeScalar(const PetscReal *A,const PetscReal *x,PetscReal *y);

void QuadMatMultTransposeVel(const PetscReal *A,const PISMVector2 *x,PISMVector2 *y);

int PismIntMask(PetscScalar maskvalue);

PetscErrorCode QuadEvaluateVel(const PISMVector2 *x,PetscInt q,
                               const PetscReal jinvDiag[],
                               PISMVector2 *u,PetscReal Du[]);

PetscErrorCode QuadGetStencils(DALocalInfo *info,PetscInt i,PetscInt j,
                               MatStencil row[],MatStencil col[]);
