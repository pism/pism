/* Copyright (C) 2010-2011 Jed Brown and Ed Bueler */

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

/* ! \file THItools.h Tools used by Jed's Blatter solver ("Toy Hydrostatic Ice"). */

#include <petscdmmg.h>
#include <ctype.h>              // toupper()
#include <private/daimpl.h>     // There is not yet a public interface to manipulate dm.ops

#if !defined __STDC_VERSION__ || __STDC_VERSION__ < 199901L
#  if defined __cplusplus       /* C++ restrict is nonstandard and compilers have inconsistent rules about where it can be used */
#    define restrict
#  else
#    define restrict PETSC_RESTRICT
#  endif
#endif
#if defined __SSE2__
#  include <emmintrin.h>
#endif

/* The SSE2 kernels are only for PetscScalar=double on architectures that support it */
#define USE_SSE2_KERNELS (!defined NO_SSE2                              \
                          && !defined PETSC_USE_COMPLEX                 \
                          && !defined PETSC_USE_SCALAR_SINGLE           \
                          && !defined PETSC_USE_SCALAR_LONG_DOUBLE      \
                          && defined __SSE2__)

typedef enum {QUAD_GAUSS,QUAD_LOBATTO} QuadratureType;
static const char *QuadratureTypes[] = {"gauss","lobatto","QuadratureType","QUAD_",0};
static const PetscReal HexQWeights[8] = {1,1,1,1,1,1,1,1};
static const PetscReal HexQNodes[]    = {-0.57735026918962573, 0.57735026918962573};

#define G 0.57735026918962573
#define H (0.5*(1.+G))
#define L (0.5*(1.-G))
#define M (-0.5)
#define P (0.5)

/* Special quadrature: Lobatto in horizontal, Gauss in vertical */
static const PetscReal HexQInterp_Lobatto[8][8] = {{H,0,0,0,L,0,0,0},
                                                   {0,H,0,0,0,L,0,0},
                                                   {0,0,H,0,0,0,L,0},
                                                   {0,0,0,H,0,0,0,L},
                                                   {L,0,0,0,H,0,0,0},
                                                   {0,L,0,0,0,H,0,0},
                                                   {0,0,L,0,0,0,H,0},
                                                   {0,0,0,L,0,0,0,H}};

static const PetscReal HexQDeriv_Lobatto[8][8][3] = {
  {{M*H,M*H,M},{P*H,0,0}  ,{0,0,0}    ,{0,P*H,0}  ,{M*L,M*L,P},{P*L,0,0}  ,{0,0,0}    ,{0,P*L,0}  },
  {{M*H,0,0}  ,{P*H,M*H,M},{0,P*H,0}  ,{0,0,0}    ,{M*L,0,0}  ,{P*L,M*L,P},{0,P*L,0}  ,{0,0,0}    },
  {{0,0,0}    ,{0,M*H,0}  ,{P*H,P*H,M},{M*H,0,0}  ,{0,0,0}    ,{0,M*L,0}  ,{P*L,P*L,P},{M*L,0,0}  },
  {{0,M*H,0}  ,{0,0,0}    ,{P*H,0,0}  ,{M*H,P*H,M},{0,M*L,0}  ,{0,0,0}    ,{P*L,0,0}  ,{M*L,P*L,P}},
  {{M*L,M*L,M},{P*L,0,0}  ,{0,0,0}    ,{0,P*L,0}  ,{M*H,M*H,P},{P*H,0,0}  ,{0,0,0}    ,{0,P*H,0}  },
  {{M*L,0,0}  ,{P*L,M*L,M},{0,P*L,0}  ,{0,0,0}    ,{M*H,0,0}  ,{P*H,M*H,P},{0,P*H,0}  ,{0,0,0}    },
  {{0,0,0}    ,{0,M*L,0}  ,{P*L,P*L,M},{M*L,0,0}  ,{0,0,0}    ,{0,M*H,0}  ,{P*H,P*H,P},{M*H,0,0}  },
  {{0,M*L,0}  ,{0,0,0}    ,{P*L,0,0}  ,{M*L,P*L,M},{0,M*H,0}  ,{0,0,0}    ,{P*H,0,0}  ,{M*H,P*H,P}}};

/* Standard Gauss */
static const PetscReal HexQInterp_Gauss[8][8] = {{H*H*H,L*H*H,L*L*H,H*L*H, H*H*L,L*H*L,L*L*L,H*L*L},
                                                 {L*H*H,H*H*H,H*L*H,L*L*H, L*H*L,H*H*L,H*L*L,L*L*L},
                                                 {L*L*H,H*L*H,H*H*H,L*H*H, L*L*L,H*L*L,H*H*L,L*H*L},
                                                 {H*L*H,L*L*H,L*H*H,H*H*H, H*L*L,L*L*L,L*H*L,H*H*L},
                                                 {H*H*L,L*H*L,L*L*L,H*L*L, H*H*H,L*H*H,L*L*H,H*L*H},
                                                 {L*H*L,H*H*L,H*L*L,L*L*L, L*H*H,H*H*H,H*L*H,L*L*H},
                                                 {L*L*L,H*L*L,H*H*L,L*H*L, L*L*H,H*L*H,H*H*H,L*H*H},
                                                 {H*L*L,L*L*L,L*H*L,H*H*L, H*L*H,L*L*H,L*H*H,H*H*H}};

static const PetscReal HexQDeriv_Gauss[8][8][3] = {
  {{M*H*H,H*M*H,H*H*M},{P*H*H,L*M*H,L*H*M},{P*L*H,L*P*H,L*L*M},{M*L*H,H*P*H,H*L*M}, {M*H*L,H*M*L,H*H*P},{P*H*L,L*M*L,L*H*P},{P*L*L,L*P*L,L*L*P},{M*L*L,H*P*L,H*L*P}},
  {{M*H*H,L*M*H,L*H*M},{P*H*H,H*M*H,H*H*M},{P*L*H,H*P*H,H*L*M},{M*L*H,L*P*H,L*L*M}, {M*H*L,L*M*L,L*H*P},{P*H*L,H*M*L,H*H*P},{P*L*L,H*P*L,H*L*P},{M*L*L,L*P*L,L*L*P}},
  {{M*L*H,L*M*H,L*L*M},{P*L*H,H*M*H,H*L*M},{P*H*H,H*P*H,H*H*M},{M*H*H,L*P*H,L*H*M}, {M*L*L,L*M*L,L*L*P},{P*L*L,H*M*L,H*L*P},{P*H*L,H*P*L,H*H*P},{M*H*L,L*P*L,L*H*P}},
  {{M*L*H,H*M*H,H*L*M},{P*L*H,L*M*H,L*L*M},{P*H*H,L*P*H,L*H*M},{M*H*H,H*P*H,H*H*M}, {M*L*L,H*M*L,H*L*P},{P*L*L,L*M*L,L*L*P},{P*H*L,L*P*L,L*H*P},{M*H*L,H*P*L,H*H*P}},
  {{M*H*L,H*M*L,H*H*M},{P*H*L,L*M*L,L*H*M},{P*L*L,L*P*L,L*L*M},{M*L*L,H*P*L,H*L*M}, {M*H*H,H*M*H,H*H*P},{P*H*H,L*M*H,L*H*P},{P*L*H,L*P*H,L*L*P},{M*L*H,H*P*H,H*L*P}},
  {{M*H*L,L*M*L,L*H*M},{P*H*L,H*M*L,H*H*M},{P*L*L,H*P*L,H*L*M},{M*L*L,L*P*L,L*L*M}, {M*H*H,L*M*H,L*H*P},{P*H*H,H*M*H,H*H*P},{P*L*H,H*P*H,H*L*P},{M*L*H,L*P*H,L*L*P}},
  {{M*L*L,L*M*L,L*L*M},{P*L*L,H*M*L,H*L*M},{P*H*L,H*P*L,H*H*M},{M*H*L,L*P*L,L*H*M}, {M*L*H,L*M*H,L*L*P},{P*L*H,H*M*H,H*L*P},{P*H*H,H*P*H,H*H*P},{M*H*H,L*P*H,L*H*P}},
  {{M*L*L,H*M*L,H*L*M},{P*L*L,L*M*L,L*L*M},{P*H*L,L*P*L,L*H*M},{M*H*L,H*P*L,H*H*M}, {M*L*H,H*M*H,H*L*P},{P*L*H,L*M*H,L*L*P},{P*H*H,L*P*H,L*H*P},{M*H*H,H*P*H,H*H*P}}};

static const PetscReal (*HexQInterp)[8],(*HexQDeriv)[8][3];

/* Standard 2x2 Gauss quadrature for the bottom layer. */
static const PetscReal QuadQInterp[4][4] = {{H*H,L*H,L*L,H*L},
                                            {L*H,H*H,H*L,L*L},
                                            {L*L,H*L,H*H,L*H},
                                            {H*L,L*L,L*H,H*H}};

static const PetscReal QuadQDeriv[4][4][2] = {
  {{M*H,M*H},{P*H,M*L},{P*L,P*L},{M*L,P*H}},
  {{M*H,M*L},{P*H,M*H},{P*L,P*H},{M*L,P*L}},
  {{M*L,M*L},{P*L,M*H},{P*H,P*H},{M*H,P*L}},
  {{M*L,M*H},{P*L,M*L},{P*H,P*L},{M*H,P*H}}};

#undef G
#undef H
#undef L
#undef M
#undef P

#define HexExtract(x,i,j,k,n) do {              \
    (n)[0] = (x)[i][j][k];                      \
    (n)[1] = (x)[i+1][j][k];                    \
    (n)[2] = (x)[i+1][j+1][k];                  \
    (n)[3] = (x)[i][j+1][k];                    \
    (n)[4] = (x)[i][j][k+1];                    \
    (n)[5] = (x)[i+1][j][k+1];                  \
    (n)[6] = (x)[i+1][j+1][k+1];                \
    (n)[7] = (x)[i][j+1][k+1];                  \
  } while (0)

#define HexExtractRef(x,i,j,k,n) do {           \
    (n)[0] = &(x)[i][j][k];                     \
    (n)[1] = &(x)[i+1][j][k];                   \
    (n)[2] = &(x)[i+1][j+1][k];                 \
    (n)[3] = &(x)[i][j+1][k];                   \
    (n)[4] = &(x)[i][j][k+1];                   \
    (n)[5] = &(x)[i+1][j][k+1];                 \
    (n)[6] = &(x)[i+1][j+1][k+1];               \
    (n)[7] = &(x)[i][j+1][k+1];                 \
  } while (0)

#define QuadExtract(x,i,j,n) do {               \
    (n)[0] = (x)[i][j];                         \
    (n)[1] = (x)[i+1][j];                       \
    (n)[2] = (x)[i+1][j+1];                     \
    (n)[3] = (x)[i][j+1];                       \
  } while (0)

static PetscScalar Sqr(PetscScalar a) {return a*a;}

static void HexGrad(const PetscReal dphi[][3],const PetscReal zn[],PetscReal dz[])
{
  PetscInt i;
  dz[0] = dz[1] = dz[2] = 0;
  for (i=0; i<8; i++) {
    dz[0] += dphi[i][0] * zn[i];
    dz[1] += dphi[i][1] * zn[i];
    dz[2] += dphi[i][2] * zn[i];
  }
}

static void HexComputeGeometry(PetscInt q,PetscReal hx,PetscReal hy,const PetscReal dz[restrict],PetscReal phi[restrict],PetscReal dphi[restrict][3],PetscReal *restrict jw)
{
  const PetscReal
    jac[3][3] = {{hx/2,0,0}, {0,hy/2,0}, {dz[0],dz[1],dz[2]}}
  ,ijac[3][3] = {{1/jac[0][0],0,0}, {0,1/jac[1][1],0}, {-jac[2][0]/(jac[0][0]*jac[2][2]),-jac[2][1]/(jac[1][1]*jac[2][2]),1/jac[2][2]}}
  ,jdet = jac[0][0]*jac[1][1]*jac[2][2];
  PetscInt i;

  for (i=0; i<8; i++) {
    const PetscReal *dphir = HexQDeriv[q][i];
    phi[i] = HexQInterp[q][i];
    dphi[i][0] = dphir[0]*ijac[0][0] + dphir[1]*ijac[1][0] + dphir[2]*ijac[2][0];
    dphi[i][1] = dphir[0]*ijac[0][1] + dphir[1]*ijac[1][1] + dphir[2]*ijac[2][1];
    dphi[i][2] = dphir[0]*ijac[0][2] + dphir[1]*ijac[1][2] + dphir[2]*ijac[2][2];
  }
  *jw = 1.0 * jdet;
}

typedef struct _n_Units *Units;

struct _n_Units {
  /* fundamental */
  PetscReal meter;
  PetscReal kilogram;
  PetscReal second;
  /* derived */
  PetscReal Pascal;
  PetscReal year;
};

static void RangeUpdate(PetscReal *min,PetscReal *max,PetscReal x)
{
  if (x < *min) *min = x;
  if (x > *max) *max = x;
}


typedef struct {
  PetscScalar u,v;
} Node;


/* PRange CLASS */

typedef struct {
  PetscReal min,max,cmin,cmax;
} PRange;

static void PRangeClear(PRange *p)
{
  p->cmin = p->min = 1e100;
  p->cmax = p->max = -1e100;
}

#undef __FUNCT__  
#define __FUNCT__ "PRangeMinMax"
static PetscErrorCode PRangeMinMax(PRange *p,PetscReal min,PetscReal max)
{

  PetscFunctionBegin;
  p->cmin = min;
  p->cmax = max;
  if (min < p->min) p->min = min;
  if (max > p->max) p->max = max;
  PetscFunctionReturn(0);
}


/* PrmNode CLASS */

typedef struct {
  PetscScalar b;                /* bed */
  PetscScalar h;                /* thickness */
  PetscScalar beta2;            /* friction */
} PrmNode;

static void PrmNodeHexGetZ(const PrmNode pn[],PetscInt k,PetscInt zm,PetscReal zn[])
{
  const PetscScalar zm1 = zm-1,
    znl[8] = {pn[0].b + pn[0].h*(PetscScalar)k/zm1,
              pn[1].b + pn[1].h*(PetscScalar)k/zm1,
              pn[2].b + pn[2].h*(PetscScalar)k/zm1,
              pn[3].b + pn[3].h*(PetscScalar)k/zm1,
              pn[0].b + pn[0].h*(PetscScalar)(k+1)/zm1,
              pn[1].b + pn[1].h*(PetscScalar)(k+1)/zm1,
              pn[2].b + pn[2].h*(PetscScalar)(k+1)/zm1,
              pn[3].b + pn[3].h*(PetscScalar)(k+1)/zm1};
  PetscInt i;
  for (i=0; i<8; i++) zn[i] = PetscRealPart(znl[i]);
}

#undef __FUNCT__  
#define __FUNCT__ "DAGetPrmNodeArray"
static PetscErrorCode DAGetPrmNodeArray(DA da,PrmNode ***prm)
{
  PetscErrorCode ierr;
  DA             da2prm;
  Vec            X;

  PetscFunctionBegin;
  ierr = PetscObjectQuery((PetscObject)da,"DA2Prm",(PetscObject*)&da2prm);CHKERRQ(ierr);
  if (!da2prm) SETERRQ(PETSC_ERR_ARG_WRONG,"No DA2Prm composed with given DA");
  ierr = PetscObjectQuery((PetscObject)da,"DA2Prm_Vec",(PetscObject*)&X);CHKERRQ(ierr);
  if (!X) SETERRQ(PETSC_ERR_ARG_WRONG,"No DA2Prm_Vec composed with given DA");
  ierr = DAVecGetArray(da2prm,X,prm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DARestorePrmNodeArray"
static PetscErrorCode DARestorePrmNodeArray(DA da,PrmNode ***prm)
{
  PetscErrorCode ierr;
  DA             da2prm;
  Vec            X;

  PetscFunctionBegin;
  ierr = PetscObjectQuery((PetscObject)da,"DA2Prm",(PetscObject*)&da2prm);CHKERRQ(ierr);
  if (!da2prm) SETERRQ(PETSC_ERR_ARG_WRONG,"No DA2Prm composed with given DA");
  ierr = PetscObjectQuery((PetscObject)da,"DA2Prm_Vec",(PetscObject*)&X);CHKERRQ(ierr);
  if (!X) SETERRQ(PETSC_ERR_ARG_WRONG,"No DA2Prm_Vec composed with given DA");
  ierr = DAVecRestoreArray(da2prm,X,prm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

