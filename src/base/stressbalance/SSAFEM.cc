// Copyright (C) 2009--2011 Jed Brown and Ed Bueler
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

#include "SSAFEM.hh"

// The following are macros (instead of inline functions) so that error handling
// is less cluttered.  They should be replaced with empty macros when in
// optimized mode.

#define PismValidVelocity(U) do {                               \
    if (!(-   1e+5 < (U).x && (U).x < 1e+5                      \
          && -1e+5 < (U).y && (U).y < 1e+5))                    \
      SETERRQ3(1,"Invalid velocity (%g,%g) not within %g",      \
               (U).x,(U).y,1e+5);                               \
  } while (0)

#define PismValidStrainRate(Du) do {                                    \
    if (!(-   1e+5 < (Du)[0] && (Du)[0] < 1e+5                          \
          && -1e+5 < (Du)[1] && (Du)[1] < 1e+5                          \
          && -1e+5 < (Du)[2] && (Du)[2] < 1e+5))                        \
      SETERRQ4(1,"Invalid Strain Rate (%g,%g,%g) not within %g",        \
               (Du)[0],(Du)[1],(Du)[2],1e+5);                           \
  } while (0)

#define PismValidStress2(f) do {                                        \
    if (!(-   1e4 < (f).x && (f).x < 1e4                                \
          && -1e4 < (f).y && (f).y < 1e4))                              \
      SETERRQ3(1,"Invalid Stress residual (%g,%g) not within %g",       \
               (f).x,(f).y,1e4);                                        \
  } while (0)

#define PismValidFriction(b) do {                                       \
    if (!(0 <= (b) && (b) < 1e25))                                      \
      SETERRQ2(1,"Invalid friction %g not within [0,%g]",(b),1e25);     \
  } while (0)

static PetscTruth Floating(const IceFlowLaw *ice, PetscScalar ocean_rho,
                           PetscReal H, PetscReal bed)
{
  return ice->rho*H + ocean_rho*bed < 0 ? PETSC_TRUE : PETSC_FALSE;
}

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

static inline void QuadZeroScalar(PetscReal x[])
{
  PetscMemzero(x,4*sizeof(x[0]));
}
static inline void QuadZeroVel(SSANode x[])
{
  PetscMemzero(x,4*sizeof(x[0]));
}

static inline void QuadExtractScalar(PetscInt i,PetscInt j,PetscReal **xg,PetscReal x[])
{
  x[0] = xg[i][j]; x[1] = xg[i+1][j]; x[2] = xg[i+1][j+1]; x[3] = xg[i][j+1];
}
static inline void QuadInsertScalar(PetscInt i,PetscInt j,PetscReal x[],PetscReal **xg)
{
  xg[i][j] = x[0]; xg[i+1][j] = x[1]; xg[i+1][j+1] = x[2]; xg[i][j+1] = x[3];
}

static inline void QuadExtractVel(PetscInt i,PetscInt j,const SSANode **xg,SSANode x[])
{
  x[0].x = xg[i][j].x;        x[0].y = xg[i][j].y;
  x[1].x = xg[i+1][j].x;      x[1].y = xg[i+1][j].y;
  x[2].x = xg[i+1][j+1].x;    x[2].y = xg[i+1][j+1].y;
  x[3].x = xg[i][j+1].x;      x[3].y = xg[i][j+1].y;
}

static inline void QuadInsertVel(const MatStencil row[],const SSANode x[],SSANode **xg)
{
#if 0
  printf("El(%d,%d): inserting [(%g,%g) (%g,%g) (%g,%g) (%g,%g)]\n",
         i,j,x[0].x,x[0].y,x[1].x,x[1].y,x[2].x,x[2].y,x[3].x,x[3].y);
#endif
  for (int k=0; k<4; k++) {
    if (row[k].i < 0 || row[k].j < 0) continue;
    // The stencil is in PETSc ordering [j][i]
    xg[row[k].j][row[k].i].x += x[k].x;
    xg[row[k].j][row[k].i].y += x[k].y;
  }
}

static inline void QuadMatMultScalar(const PetscReal *A,const PetscReal *x,PetscReal *y)
{
  PetscInt i,j;
  const PetscReal *a;

  for (i=0,a=A; i<4; i++,a+=4) {
    y[i] = 0;
    for (j=0; j<4; j++) {
      y[i] += a[j] * x[j];
    }
  }
}

static inline void QuadMatMultVel(const PetscReal *A,const SSANode *x,SSANode *y)
{
  PetscInt i,j;
  const PetscReal *a;

  for (i=0,a=A; i<4; i++,a+=4) {
    y[i].x = y[i].y = 0;
    for (j=0; j<4; j++) {
      y[i].x += a[j] * x[j].x;
      y[i].y += a[j] * x[j].y;
    }
  }
}

static inline void QuadMatMultTransposeScalar(const PetscReal *A,const PetscReal *x,PetscReal *y)
{
  PetscInt i,j;
  const PetscReal *a;

  QuadZeroScalar(y);
  for (i=0,a=A; i<4; i++,a+=4) {
    for (j=0; j<4; j++) {
      y[j] += a[j] * x[i];
    }
  }
}

static inline void QuadMatMultTransposeVel(const PetscReal *A,const SSANode *x,SSANode *y)
{
  PetscInt i,j;
  const PetscReal *a;

  QuadZeroVel(y);
  for (i=0,a=A; i<4; i++,a+=4) {
    for (j=0; j<4; j++) {
      y[j].x += a[j] * x[i].x;
      y[j].y += a[j] * x[i].y;
    }
  }
}

#undef __FUNCT__
#define __FUNCT__ "QuadEvaluateVel"
static PetscErrorCode QuadEvaluateVel(const SSANode *x,PetscInt q,
                                      const PetscReal jinvDiag[],
                                      SSANode *u,PetscReal Du[])
{

  PetscFunctionBegin;
  u->x = u->y = 0;
  Du[0] = Du[1] = Du[2] = 0;
  for (PetscInt k=0; k<4; k++) {   // loop over trial functions to compute function values at quadrature point
    const PetscInt qk = q*4+k;
    u->x  += interp[qk] * x[k].x; u->y  += interp[qk] * x[k].y; // velocity at q
    Du[0] += jinvDiag[0]*derivx[qk] * x[k].x;                 // u_x
    Du[1] += jinvDiag[1]*derivy[qk] * x[k].y;                 // v_y
    Du[2] += 0.5*(jinvDiag[1]*derivy[qk] * x[k].x + jinvDiag[0]*derivx[qk] * x[k].y); // (u_y + v_x)/2
  }
  PismValidVelocity(*u);
  PismValidStrainRate(Du);
  PetscFunctionReturn(0);
}

static PetscErrorCode QuadGetStencils(DALocalInfo *info,PetscInt i,PetscInt j,
                                      MatStencil row[],MatStencil col[])
{
  PetscErrorCode ierr;
  PetscInt k;

  PetscFunctionBegin;
  // PETSc reverses the meaning of i and j in stencil and info objects
  ierr = PetscMemzero(col,4*sizeof(col[0]));CHKERRQ(ierr); // Probably just paranoia
  col[0].j = i;   col[0].i = j;
  col[1].j = i+1; col[1].i = j;
  col[2].j = i+1; col[2].i = j+1;
  col[3].j = i;   col[3].i = j+1;
  ierr = PetscMemcpy(row,col,4*sizeof(col[0]));CHKERRQ(ierr);
  for (k=0; k<4; k++) {         // We do not sum into rows that are not owned by us
    if (   row[k].i < info->xs || info->xs+info->xm-1 < row[k].i
        || row[k].j < info->ys || info->ys+info->ym-1 < row[k].j) {
      row[k].j = row[k].i = -1;
    }
  }
  PetscFunctionReturn(0);
}


//* Apply a stress along a particular edge.  The canonical element node numbering is
//
// 3 -- 2  with edge numbering   |- 2 -|
// |    |                        3     1
// 0 -- 1                        |_ 0 _|
//
// The applied force will always be directed inwards on the specified edge, with the force at quadrature points
// determined by whether that node is floating and the state of some options (-ssa_boundary_*).  The default is to apply
// a calving front to floating ice and ice that is grounded below sea level and to apply a stress-free condition at
// margins that are grounded above sea level since this probably indicates that it's not intended to be a margin at all.
static PetscErrorCode ApplyBoundaryStress(
              FECTX *fe, PetscReal jacDiag[],
              const IceFlowLaw *ice, PetscReal ocean_rho,
              PetscReal **H,PetscReal **b,const MatStencil col[],PetscInt edge,
              SSANode y[]) {
  static const PetscInt node_list[4][2] = {{0,1},{1,2},{2,3},{3,0}};
  static const SSANode coeff_list[4] = {{0,-1},{-1,0},{0,1},{1,0}};
  static const PetscInt direction[4] = {0,1,0,1};
  const PetscInt *node  = node_list[edge];
  const SSANode  *coeff = &coeff_list[edge];
  PetscInt q,k;
  SSAFEM    *ssa = fe->ssa;
  PetscReal  earth_grav = fe->earth_grav;

  for (q=0; q<2; q++) {
    const PetscReal jw = jacDiag[direction[edge]] * quadWeights1[q];
    PetscReal Hq = 0,bq = 0,f = 0;
    for (k=0; k<2; k++) {
      Hq += interp1[q*2+k] * H[col[node[k]].j][col[node[k]].i];
      bq += interp1[q*2+k] * b[col[node[k]].j][col[node[k]].i];
    }
    if (Floating(ice,ocean_rho,Hq,bq)) {
      // Floating ice can either be stress-free or have a calving face condition
      if (!ssa->boundary.floating_stress_free) {
        f = 0.5 * earth_grav * ice->rho * (1 - ice->rho/ocean_rho) * Hq * Hq;
      }
    } else if (bq < 0) {        // Grounded below sea level at the margin
      if (ssa->boundary.grounded_as_floating) {
        f = 0.5 * earth_grav * ice->rho * (1 - ice->rho/ocean_rho) * Hq * Hq;
      } else if (!ssa->boundary.submarine_stress_free) {
        // A non-floating submarine calving face, this is the default condition
        f = 0.5 * earth_grav * (ice->rho*Hq*Hq - ocean_rho*bq*bq);
      }
    } else {                    // Grounded above sea level at the margin
      // Usually stress-free if the boundary is above sea level, this probably means that the edge of the domain is not
      // really supposed to be a boundary.  If you want it to be a calving face, set
      // -ssa_boundary_calving_face_above_sea_level
      if (ssa->boundary.calving_above_sea_level) {
        f = 0.5 * earth_grav * ice->rho * Hq * Hq;
      }
    }
    f /= fe->ref.Pressure() * fe->ref.Height(); // Nondimensionalize
    for (k=0; k<2; k++) {
      y[node[k]].x += coeff->x * interp1[q*2+k] * jw * f;
      y[node[k]].y += coeff->y * interp1[q*2+k] * jw * f;
    }
  }
  return 0;
}

// This is called when surface and/or temperature have changed.
//
// Since this interfaces with the rest of PISM and doesn't touch velocity or
// stresses, everything is fully dimensional
#undef __FUNCT__
#define __FUNCT__ "SSAFESetUp"
static PetscErrorCode SSAFESetUp(FECTX *fe)
{
  SSAFEM          *ssa = fe->ssa;
  DALocalInfo      info_struct,
                  *info = &info_struct; // DAGetLocalInfo writes into our memory
  PetscReal      **h,
                 **H,
                 **bed,
                 **tauc,
                  *Te[4],
                  *Tq[4];
  PetscReal        jinvDiag[2];
  IceGrid         *grid = fe->grid;
  PetscInt         i,j,k,q,p,
                   Mz = grid->Mz;
  PetscErrorCode   ierr;

  PetscFunctionBegin;
  // We only differentiate dimensional quantities (not velocity) here
  jinvDiag[0] = 2/grid->dx;
  jinvDiag[1] = 2/grid->dy;

  ierr = PetscMalloc4(Mz,PetscReal,&Tq[0],
                      Mz,PetscReal,&Tq[1],
                      Mz,PetscReal,&Tq[2],
                      Mz,PetscReal,&Tq[3]);CHKERRQ(ierr);
  ierr = ssa->T->begin_access();CHKERRQ(ierr);
  ierr = ssa->h->get_array(h);CHKERRQ(ierr);
  ierr = ssa->H->get_array(H);CHKERRQ(ierr);
  ierr = ssa->bed->get_array(bed);CHKERRQ(ierr);
  ierr = ssa->tauc->get_array(tauc);CHKERRQ(ierr);
  ierr = DAGetLocalInfo(ssa->da,info);CHKERRQ(ierr);
  // See SSAFEFunction for discussion of communication
  for (i=info->gys; i<info->gys+info->gym-1; i++) { // Include ghost cells on either end
    for (j=info->gxs; j<info->gxs+info->gxm-1; j++) {
      const PetscInt ij = (i-info->gys)*(info->gxm-1)+(j-info->gxs);  // Index for stored arrays
      PetscReal he[4],He[4],be[4],hq[4],Hq[4],bq[4],taue[4],tauq[4],hxq[4],hyq[4];
      FEStoreNode *feS;
      if ( (i < 0 || info->my <= i || j < 0 || info->mx <= j)
           && (ssa->wrap == DA_NONPERIODIC) )
        SETERRQ(1,"We are not using a periodic grid, this should not happen");
      QuadExtractScalar(i,j,h,he); // surface height
      QuadExtractScalar(i,j,H,He); // thickness
      QuadExtractScalar(i,j,bed,be); // bed elevation
      QuadExtractScalar(i,j,tauc,taue); // basal friction
      QuadMatMultScalar(interp,he,hq);
      QuadMatMultScalar(interp,He,Hq);
      QuadMatMultScalar(interp,be,bq);
      QuadMatMultScalar(interp,taue,tauq);
      QuadMatMultScalar(derivx,he,hxq);
      QuadMatMultScalar(derivy,he,hyq);
      feS = &fe->feStore[4*ij];
      for (q=0; q<4; q++) {
        feS[q].h  = hq[q];
        feS[q].H  = Hq[q];
        feS[q].b  = bq[q];
        feS[q].tauc = tauq[q];
        feS[q].hx = jinvDiag[0] * hxq[q];
        feS[q].hy = jinvDiag[0] * hyq[q];
      }

      // Surface and thickness information is stored, now do the thermal stuff
      ierr = ssa->T->getInternalColumn(i,j,&Te[0]);CHKERRQ(ierr);
      ierr = ssa->T->getInternalColumn(i+1,j,&Te[1]);CHKERRQ(ierr);
      ierr = ssa->T->getInternalColumn(i+1,j+1,&Te[2]);CHKERRQ(ierr);
      ierr = ssa->T->getInternalColumn(i,j+1,&Te[3]);CHKERRQ(ierr);
      // Interpolate to quadrature points at every vertical level
      for (k=0; k<Mz; k++) { // This loop could be cut short at the surface.
        Tq[0][k] = Tq[1][k] = Tq[2][k] = Tq[3][k] = 0;
        for (q=0; q<4; q++) {
          for (p=0; p<4; p++) {
            Tq[q][k] += interp[q*4+p] * Te[p][k];
          }
        }
      }
      for (q=0; q<4; q++) {
        // Evaluate column integrals in flow law at every quadrature point's column
        PetscReal *iS = &fe->integratedStore[fe->sbs*(ij*4+q)]; // Location to put the stored data
        fe->ice->integratedStore(
               feS[q].H,grid->kBelowHeight(feS[q].H),grid->zlevels,Tq[q],iS);
      }
    }
  }
  ierr = ssa->h->end_access();CHKERRQ(ierr);
  ierr = ssa->H->end_access();CHKERRQ(ierr);
  ierr = ssa->bed->end_access();CHKERRQ(ierr);
  ierr = ssa->tauc->end_access();CHKERRQ(ierr);
  ierr = ssa->T->end_access();CHKERRQ(ierr);
  ierr = PetscFree4(Tq[0],Tq[1],Tq[2],Tq[3]);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// The stores are dimensional, everything after \a u (inclusive) is nondimensional
#undef __FUNCT__
#define __FUNCT__ "PointwiseNuHAndBeta"
static inline PetscErrorCode PointwiseNuHAndBeta(
       FECTX *fe,const FEStoreNode *feS,const PetscReal *iS,
       const SSANode *u,const PetscReal Du[],
       PetscReal *nuH,PetscReal *dNuH,PetscReal *beta,PetscReal *dbeta)
{

  PetscFunctionBegin;
  SSAFEM   *ssa = fe->ssa;
  if (feS->H < ssa->shelfExtension->thickness()) {
    *nuH = ssa->shelfExtension->viscosity();
    if (dNuH) *dNuH = 0;
    //SETERRQ(1,"Shold not happen for test I");
  } else {
    PetscReal dimDu[3];
    for (int i=0; i<3; i++) dimDu[i] = fe->ref.StrainRate() * Du[i];
    fe->ice->integratedViscosity(iS,dimDu,nuH,dNuH);
  }
  if (Floating(fe->ice,fe->ocean_rho,feS->H,feS->b)) {
    // The ice is floating here so there is no friction.  Note that the purpose of checking flotation this way is to get
    // subgrid resolution of stress in the vicinity of the grounding line.  According to Goldberg et. al. 2009 (probably
    // will be published in 2009...) this is important to loosen the resolution requirements near the grounding line.
    *beta = 0;
    if (dbeta) *dbeta = 0;
    SETERRQ(1,"Not tested yet");
  } else {
    ssa->basal->dragWithDerivative(feS->tauc,
                                   u->x*fe->ref.Velocity(),u->y*fe->ref.Velocity(),
                                   beta,dbeta);
    if (1) {
      PetscReal good_beta;
      good_beta = ssa->basal->drag(feS->tauc,
                                   u->x*fe->ref.Velocity(),u->y*fe->ref.Velocity());
      if (PetscAbs(*beta - good_beta)/(*beta + good_beta) > 0) { // Use tolerance to not test associativity
        SETERRQ2(1,"`dragWithDerivative' producing different answers from `drag' %e != %e",
                 *beta,good_beta);
      }
    }
  }
  // Return nondimensional values, the factor of 2 comes from a different
  //   definition of integrated viscosity
  *nuH  *=  2 / fe->ref.IntegratedViscosity();
  if (dNuH) *dNuH *= 2 * fe->ref.StrainRate2() / fe->ref.IntegratedViscosity(); // The derivative with respect to second invariant
  *beta /= fe->ref.Drag();
  if (dbeta) *dbeta *= fe->ref.Velocity2() / fe->ref.Drag();
  PismValidFriction(*beta);
  PetscFunctionReturn(0);
}


static void FixDirichletValues(FECTX *fe,PetscReal lmask[],SSANode **siaVel,
                               MatStencil row[],MatStencil col[],SSANode x[])
{
  for (PetscInt k=0; k<4; k++) {
    if (PismIntMask(lmask[k]) == MASK_SHEET) {
      // When in the SIA domain, use the correct values to begin with.  This operation makes the column of the
      // Jacobian corresponding to the Dirichlet dof look like the identity.  Overwriting the output vector (loop
      // near the end of this function) is equivalent to making the row look like the identity.
      //
      // Note that \a row might have some entries already eliminated so we must use \a col (which is more natural
      // anyway).
      x[k].x = siaVel[col[k].j][col[k].i].x / fe->ref.Velocity();
      x[k].y = siaVel[col[k].j][col[k].i].y / fe->ref.Velocity();
      row[k].j = row[k].i = -1;
      col[k].j = col[k].i = -1;
    }
  }
}

#undef __FUNCT__
#define __FUNCT__ "SSAFEFunction"
static PetscErrorCode SSAFEFunction(DALocalInfo *info,
                                    const SSANode **xg, SSANode **yg,
                                    FECTX *fe)
{
  SSAFEM          *ssa = fe->ssa;
  IceGrid         *grid = fe->grid;
  PetscInt         i,j,k,q;
  PetscReal        jacDiag[2],jinvDiag[2],jdet;
  PetscReal      **mask,**H,**bed;
  SSANode        **siaVel;
  PetscTruth       flg;
  PetscErrorCode   ierr;

  PetscFunctionBegin;
  if (!finite(xg[info->ys][info->xs].x)) SETERRQ(1,__FUNCT__ " called with non-finite value");
  ierr = verbPrintf(5,grid->com,"In %s\n",__FUNCT__);CHKERRQ(ierr);

  for (i=grid->xs; i<grid->xs+grid->xm; i++) {
    for (j=grid->ys; j<grid->ys+grid->ym; j++) {
      yg[i][j].x = yg[i][j].y = 0;
    }
  }
  // Since we use uniform cartesian coordinates, the Jacobian is constant and diagonal on every element.
  // Note that the reference element is \f$ [-1,1]^2 \f$ hence the extra factor of 2.
  jacDiag[0] = 0.5*grid->dx/fe->ref.Length();
  jacDiag[1] = 0.5*grid->dy/fe->ref.Length();
  //if (fe->debug.rescale) jacDiag[0] = jacDiag[1] = 1;
  jinvDiag[0] = 1/jacDiag[0];
  jinvDiag[1] = 1/jacDiag[1];
  jdet = jacDiag[0]*jacDiag[1];

  ierr = PetscOptionsHasName(NULL,"-ssa_monitor_solution",&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscPrintf(grid->com,
        "SSA Pointwise solution values (evaluating function at this value)\n");
        CHKERRQ(ierr);
    for (i=grid->xs; i<grid->xs+grid->xm; i++) {
      for (j=grid->ys; j<grid->ys+grid->ym; j++) {
        ierr = PetscPrintf(grid->com,"[%2d,%2d] (%g,%g)\n",i,j,
                           xg[i][j].x,xg[i][j].y);CHKERRQ(ierr);
      }
    }
  }
  // Consider the node layout in 1 dimension for ranks A,B,C
  // Owning rank:     A   |   B   |   C
  // Node index :   0 1 2 | 3 4 5 | 6 7 8
  //
  // Rank B has up-to-date values of the current iterate for [2,3,4,5,6] but we only set function values for the owned
  // portion [3,4,5].  In the finite element context, rank B will evaluate weak forms on elements
  // [(2,3),(3,4),(4,5),(5,6)] but will only insert values for the owned nodes.  This is consistent with the
  // communication pattern natural to finite difference methods for which the DA object is designed.
  //
  // In the finite element context, the usual communication pattern is to have the numbers (0..8) above refer to
  // elements and decide on ownership of interface degrees of freedom.  Then each rank evaluates weak forms over the
  // locally owned elements and performs communication to sum into the global vectors and matrices.  With FD
  // communication, a rank only inserts into locally owned portions of the residual vector and rows of the global matrix
  // while with FE communication each rank adds values into locally owned plus remotely owned interface parts of the
  // residual vector and rows of the global matrix.
  //
  // Note that PETSc uses the opposite meaning of x,y directions in the DA
  ierr = ssa->mask->get_array(mask);CHKERRQ(ierr);
  ierr = ssa->H->get_array(H);CHKERRQ(ierr);
  ierr = ssa->bed->get_array(bed);CHKERRQ(ierr);
  ierr = DAVecGetArray(fe->da,ssa->siaVelLocal,&siaVel);CHKERRQ(ierr);
  for (i=info->gys; i<info->gys+info->gym-1; i++) { // Include ghost cells on either end
    for (j=info->gxs; j<info->gxs+info->gxm-1; j++) {
      const PetscInt ij = (i-info->gys)*(info->gxm-1)+(j-info->gxs); // Index for stored arrays
      SSANode     x[4],y[4];
      PetscReal lmask[4];
      MatStencil row[4],col[4];
      if ((i < 0 || info->my <= i || j < 0 || info->mx <= j) && (ssa->wrap == DA_NONPERIODIC))
        SETERRQ(1,"We are not using a periodic grid, this should not happen");
      QuadExtractVel(i,j,xg,x);
      QuadExtractScalar(i,j,mask,lmask);
      ierr = QuadGetStencils(info,i,j,row,col);CHKERRQ(ierr);
      QuadZeroVel(y);
      // Apply inhomogeneous Neumann conditions.  We do this first because it does not depend on velocity and \a col
      // still contains the indices for all the nodes.  Note that FixDirichletValues changes \a col so that it only has
      // valid (non-negative) indices for the non-Dirichlet nodes.
      //
      // If the left margin of our ghosted patch is equal to the left margin of our owned domain, then we're at the left
      // margin of the global domain.  See ApplyBoundaryStress for the edge numbering.
      if (i==info->gys && i==info->ys) {
        ierr = ApplyBoundaryStress(fe,jacDiag,fe->ice,fe->ocean_rho,H,bed,col,3,y);CHKERRQ(ierr);
      } else if (i==info->gys+info->gym-2 && i==info->ys+info->ym-2) {
        ierr = ApplyBoundaryStress(fe,jacDiag,fe->ice,fe->ocean_rho,H,bed,col,1,y);CHKERRQ(ierr);
      }
      if (j==info->gxs && j==info->xs) {
        ierr = ApplyBoundaryStress(fe,jacDiag,fe->ice,fe->ocean_rho,H,bed,col,0,y);CHKERRQ(ierr);
      } else if (j==info->gxs+info->gxm-2 && j==info->xs+info->xm-2) {
        ierr = ApplyBoundaryStress(fe,jacDiag,fe->ice,fe->ocean_rho,H,bed,col,2,y);CHKERRQ(ierr);
      }
      FixDirichletValues(fe,lmask,siaVel,row,col,x);
      // \a x now contains correct velocities at Dirichlet nodes and \a col has indices for those nodes set to -1
      for (q=0; q<numQuadPoints; q++) {     // loop over quadrature points
        const FEStoreNode *feS = &fe->feStore[ij*4+q];
        const PetscReal *iS  = &fe->integratedStore[(ij*4+q)*fe->sbs];
        const PetscReal    jw  = jdet * quadWeights[q];
        SSANode u,v;
        PetscReal nuH,beta,Du[3],Dv[3];
        ierr = QuadEvaluateVel(x,q,jinvDiag,&u,Du);CHKERRQ(ierr);
        ierr = PointwiseNuHAndBeta(fe,feS,iS,&u,Du,&nuH,NULL,&beta,NULL);CHKERRQ(ierr);
        //if (fe->debug.const_nuH) {nuH = 1;}    // nondimensional
        //if (fe->debug.zero_nuH) {nuH = 0;}
        //if (fe->debug.const_beta) {beta = 1;}  // nondimensional
        //if (fe->debug.zero_beta) {beta = 0;}
        // Coefficients of test function (v) in weak form.
        v.x = beta*u.x + fe->ice->rho * fe->earth_grav * feS->H * feS->hx
                           / fe->ref.DrivingStress();
        v.y = beta*u.y + fe->ice->rho * fe->earth_grav * feS->H * feS->hy
                           / fe->ref.DrivingStress();
        Dv[0] = nuH * Du[0];
        Dv[1] = nuH * Du[1];
        Dv[2] = nuH * Du[2];

        //if (fe->debug.cross) {Dv[2] = 0;}

        // Sum residuals over test functions
        for (k=0; k<4; k++) {   // Loop over test functions
          const PetscInt qk = q*4+k; // this indexing transposes matrices
          y[k].x += interp[qk] * jw * v.x                     // basal drag and driving force
            + derivx[qk] * jinvDiag[0] * jw * (2*Dv[0]+Dv[1]) // diagonal term in colon plus the trace term
            + derivy[qk] * jinvDiag[1] * jw * Dv[2];          // off-diagonal term "half is summed twice"
          y[k].y += interp[qk] * jw * v.y                     // basal drag and driving force
            + derivy[qk] * jinvDiag[1] * jw * (2*Dv[1]+Dv[0]) // diagonal term in colon plus the trace term
            + derivx[qk] * jinvDiag[0] * jw * Dv[2];          // off-diagonal term "half is summed twice"
        }
      }
      QuadInsertVel(row,y,yg);
    }
  }
  ierr = ssa->H->end_access();CHKERRQ(ierr);
  ierr = ssa->bed->end_access();CHKERRQ(ierr);

  // Enforce Dirichlet conditions strongly
  for (i=grid->xs; i<grid->xs+grid->xm; i++) {
    for (j=grid->ys; j<grid->ys+grid->ym; j++) {
      //PismValidStress2(yg[i][j]);
      if (PismIntMask(mask[i][j]) == MASK_SHEET) {
        // Enforce zero sliding strongly
        yg[i][j].x = fe->dirichletScale * (xg[i][j].x - siaVel[i][j].x)
                        / fe->ref.Velocity();
        yg[i][j].y = fe->dirichletScale * (xg[i][j].y - siaVel[i][j].y)
                        / fe->ref.Velocity();
      }
    }
  }
  ierr = ssa->mask->end_access();CHKERRQ(ierr);
  ierr = DAVecRestoreArray(fe->da,ssa->siaVelLocal,&siaVel);CHKERRQ(ierr);

  ierr = PetscOptionsHasName(NULL,"-ssa_monitor_function",&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscPrintf(grid->com,"SSA Solution and Function values (pointwise residuals)\n");CHKERRQ(ierr);
    for (i=grid->xs; i<grid->xs+grid->xm; i++) {
      for (j=grid->ys; j<grid->ys+grid->ym; j++) {
        ierr = PetscSynchronizedPrintf(grid->com,
                 "[%2d,%2d] u=(%12.4e,%12.4e)  f=(%12.4e,%12.4e)\n",
                 i,j,xg[i][j].x,xg[i][j].y,yg[i][j].x,yg[i][j].y);CHKERRQ(ierr);
      }
    }
    ierr = PetscSynchronizedFlush(grid->com);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

// As in function evaluation, velocity is nondimensional
#undef __FUNCT__
#define __FUNCT__ "SSAFEJacobian"
static PetscErrorCode SSAFEJacobian(DALocalInfo *info,
                                    const SSANode **xg, Mat J,
                                    FECTX *fe)
{
  SSAFEM          *ssa  = fe->ssa;
  IceGrid         *grid = fe->grid;
  PetscReal        jacDiag[2],jinvDiag[2],jdet;
  PetscReal      **mask;
  SSANode        **siaVel;
  PetscInt         i,j;
  PetscErrorCode   ierr;

  PetscFunctionBegin;
  ierr = verbPrintf(5,grid->com,"In %s\n",__FUNCT__);CHKERRQ(ierr);

  // Since we use uniform cartesian coordinates, the Jacobian is constant and diagonal on every element.
  // Note that the reference element is \f$ [-1,1]^2 \f$ hence the extra factor of 2.
  jacDiag[0] = 0.5*grid->dx / fe->ref.Length();
  jacDiag[1] = 0.5*grid->dy / fe->ref.Length();
  //if (fe->debug.rescale) jacDiag[0] = jacDiag[1] = 1;
  jinvDiag[0] = 1/jacDiag[0];
  jinvDiag[1] = 1/jacDiag[1];
  jdet = jacDiag[0]*jacDiag[1];

  ierr = MatZeroEntries(J);CHKERRQ(ierr);

  ierr = ssa->mask->get_array(mask);CHKERRQ(ierr);
  ierr = DAVecGetArray(fe->da,ssa->siaVelLocal,&siaVel);CHKERRQ(ierr);
  for (i=info->gys; i<info->gys+info->gym-1; i++) { // Include ghost cells on either end
    for (j=info->gxs; j<info->gxs+info->gxm-1; j++) {
      const PetscInt ij = (i-info->gys)*(info->gxm-1)+(j-info->gxs); // Index for stored arrays
      MatStencil     row[4],col[4];
      SSANode        x[4];
      PetscReal      K[4*4*4],lmask[4];
      if ( (i < 0 || info->my <= i || j < 0 || info->mx <= j) 
           && (ssa->wrap == DA_NONPERIODIC) )
        SETERRQ(1,"We are not using a periodic grid, this should not happen");
      QuadExtractVel(i,j,xg,x);
      QuadExtractScalar(i,j,mask,lmask);
      ierr = QuadGetStencils(info,i,j,row,col);CHKERRQ(ierr);
      FixDirichletValues(fe,lmask,siaVel,row,col,x);
      ierr = PetscMemzero(K,sizeof(K));CHKERRQ(ierr);
      for (PetscInt q=0; q<numQuadPoints; q++) {
        const FEStoreNode *feS = &fe->feStore[ij*4+q];
        const PetscReal   *iS  = &fe->integratedStore[(ij*4+q)*fe->sbs];
        const PetscReal    jw  = jdet*quadWeights[q];
        SSANode w;
        PetscReal nuH,dNuH,beta,dbeta,Dw[3];
        ierr = QuadEvaluateVel(x,q,jinvDiag,&w,Dw);CHKERRQ(ierr);
        ierr = PointwiseNuHAndBeta(fe,feS,iS,&w,Dw,&nuH,&dNuH,&beta,&dbeta);CHKERRQ(ierr);
        //if (fe->debug.const_nuH) {nuH = 1;dNuH = 0;}    // nondimensional
        //if (fe->debug.zero_nuH) {nuH = dNuH = 0;}
        //if (fe->debug.const_beta) {beta = 1;dbeta = 0;} // nondimensional
        //if (fe->debug.zero_beta) {beta = dbeta = 0;}

        for (PetscInt k=0; k<4; k++) {   // Test functions
          for (PetscInt l=0; l<4; l++) { // Trial functions
            const PetscInt qk = q*4+k,ql = q*4+l; // transpose and regular
            const PetscReal ht = interp[qk],h = interp[ql],
              dxt = derivx[qk]*jinvDiag[0],dyt = derivy[qk]*jinvDiag[1],
              dx = jinvDiag[0]*derivx[ql], dy = jinvDiag[1]*derivy[ql],
              // Cross terms appearing with beta'
              bvx = ht*w.x,bvy = ht*w.y,bux = w.x*h,buy = w.y*h,
              // Cross terms appearing with nuH'
              cvx = dxt*(2*Dw[0]+Dw[1]) + dyt*Dw[2],
              cvy = dyt*(2*Dw[1]+Dw[0]) + dxt*Dw[2],
              cux = (2*Dw[0]+Dw[1])*dx + Dw[2]*dy,
              cuy = (2*Dw[1]+Dw[0])*dy + Dw[2]*dx;
            // u-u coupling
            K[k*16+l*2]     += jw*(beta*ht*h + dbeta*bvx*bux + nuH*(2*dxt*dx + dyt*0.5*dy) + dNuH*cvx*cux);
            // u-v coupling
            K[k*16+l*2+1]   += jw*(dbeta*bvx*buy + nuH*(0.5*dyt*dx + dxt*dy) + dNuH*cvx*cuy);
            // v-u coupling
            K[k*16+8+l*2]   += jw*(dbeta*bvy*bux + nuH*(0.5*dxt*dy + dyt*dx) + dNuH*cvy*cux);
            // v-v coupling
            K[k*16+8+l*2+1] += jw*(beta*ht*h + dbeta*bvy*buy + nuH*(2*dyt*dy + dxt*0.5*dx) + dNuH*cvy*cuy);
          }
        }
      }
      ierr = MatSetValuesBlockedStencil(J,4,row,4,col,K,ADD_VALUES);CHKERRQ(ierr);
    }
  }
  for (i=info->ys; i<info->ys+info->ym; i++) { // Include ghost cells on either end
    for (j=info->xs; j<info->xs+info->xm; j++) {
      if (PismIntMask(mask[i][j]) == MASK_SHEET) {
        const PetscReal ident[4] = {fe->dirichletScale,0,0,fe->dirichletScale};
        MatStencil row;
        row.j = i; row.i = j;
        ierr = MatSetValuesBlockedStencil(J,1,&row,1,&row,ident,ADD_VALUES);CHKERRQ(ierr);
      }
    }
  }
  ierr = ssa->mask->end_access();CHKERRQ(ierr);
  ierr = DAVecRestoreArray(fe->da,ssa->siaVelLocal,&siaVel);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatSetOption(J,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SSADestroy_FE"
static PetscErrorCode SSADestroy_FE(FECTX *fe)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = SNESDestroy(fe->snes);CHKERRQ(ierr);
  ierr = PetscFree(fe->integratedStore);CHKERRQ(ierr);
  ierr = PetscFree(fe->feStore);CHKERRQ(ierr);
  ierr = PetscFree(fe);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SSASolve_FE"
static PetscErrorCode SSASolve_FE(FECTX *fe, Vec X)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = SSAFESetUp(fe); CHKERRQ(ierr);
  ierr = SNESSolve(fe->snes,NULL,X);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SSAView_FE"
static PetscErrorCode SSAView_FE(FECTX *fe, PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscTruth iascii;

  PetscFunctionBegin;
  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii);
           CHKERRQ(ierr);
  if (iascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"FE: Store block size %d\n",
                                  fe->sbs);CHKERRQ(ierr);
  } else {
    SETERRQ1(PETSC_ERR_SUP,"Viewer type %s not supported for SSA FE",
             ((PetscObject)viewer)->type_name);
  }
  ierr = SNESView(fe->snes,viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "SSASetFromOptions_FE"
static PetscErrorCode SSASetFromOptions_FE(FECTX *fe)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscOptionsHead("SSA FE options");CHKERRQ(ierr);
  ierr = PetscOptionsReal("-ssa_fe_dirichlet_scale",
                          "Enforce Dirichlet conditions with this additional scaling",
                          "SSAFESetDirichletScale",
                          fe->dirichletScale,
                          &fe->dirichletScale,NULL);CHKERRQ(ierr);
#if 0
  ierr = PetscOptionsTruth("-ssa_debug_rescale",
                           "Arbitrarily choose element Jacobian to be the identity","",
                           fe->debug.rescale,
                           &fe->debug.rescale,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsTruth("-ssa_debug_const_nuH",
                           "Set nuH=1","",
                           fe->debug.const_nuH,
                           &fe->debug.const_nuH,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsTruth("-ssa_debug_const_beta",
                           "Set beta=1","",
                           fe->debug.const_beta,
                           &fe->debug.const_beta,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsTruth("-ssa_debug_zero_nuH",
                           "Set nuH=0","",
                           fe->debug.zero_nuH,
                           &fe->debug.zero_nuH,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsTruth("-ssa_debug_zero_beta",
                           "Set beta=0","",
                           fe->debug.zero_beta,
                           &fe->debug.zero_beta,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsTruth("-ssa_debug_cross",
                           "Drop cross (u_y,v_x) terms","",
                           fe->debug.cross,
                           &fe->debug.cross,NULL);CHKERRQ(ierr);
#endif
  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SSASetUp_FE"
static PetscErrorCode SSASetUp_FE(FECTX *fe)
{
  PetscErrorCode  ierr;
  SSAFEM         *ssa = fe->ssa;

  PetscFunctionBegin;
  ierr = SNESCreate(((PetscObject)ssa)->comm,&fe->snes);CHKERRQ(ierr);
  ierr = SNESSetOptionsPrefix(fe->snes,((PetscObject)ssa)->prefix);CHKERRQ(ierr);

  ierr = DASetLocalFunction(fe->da,(DALocalFunction1)SSAFEFunction);CHKERRQ(ierr);
  ierr = DASetLocalJacobian(fe->da,(DALocalFunction1)SSAFEJacobian);CHKERRQ(ierr);

  ierr = SNESSetFunction(fe->snes,ssa->r,SNESDAFormFunction,fe);CHKERRQ(ierr);
  ierr = SNESSetJacobian(fe->snes,ssa->J,ssa->J,SNESDAComputeJacobian,fe);CHKERRQ(ierr);
  ierr = SNESSetFromOptions(fe->snes);CHKERRQ(ierr);

  {
    DALocalInfo info;
    PetscInt    nElements;
    ierr = DAGetLocalInfo(fe->da,&info);CHKERRQ(ierr);
    nElements = (info.gxm-1)*(info.gym-1); // Includes overlap elements
    // We have a struct for the feStore at each quadrature point
    ierr = PetscMalloc(4*nElements*sizeof(fe->feStore[0]),&fe->feStore);CHKERRQ(ierr);
    // We don't have a struct for each block of the integrated store but we can get its block size
    fe->sbs = fe->ice->integratedStoreSize();
    ierr = PetscMalloc(4*nElements*fe->sbs*sizeof(PetscReal),&fe->integratedStore);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}


PetscErrorCode SSAFEM::allocate_fem() {
  PetscErrorCode ierr;

  ctx = new FECTX;
  // ctx->da not initialized yet
  // ctx->snes not initialized yet
  // ctx->feStore not initialized yet
  // ctx->sbs not initialized yet
  ctx->dirichletScale = 1.0;
  ctx->ice = &(this->ice);
  ctx->ocean_rho = config.get("sea_water_density");
  ctx->earth_grav = config.get("standard_gravity");
  ctx->ref.SetUp();
  ctx->grid = &(this->grid);
  ctx->ssa = this;
  return 0;
}

PetscErrorCode SSAFEM::deallocate_fem() {
  PetscErrorCode ierr;

  delete ctx;
  return 0;
}

PetscErrorCode SSAFEM::solve()
{
  PetscErrorCode ierr;
  PetscViewer    viewer;
  char           filename[PETSC_MAX_PATH_LEN];
  PetscTruth     flg;

  PetscFunctionBegin;
  ierr = PetscOptionsGetString(NULL, "-ssa_view", filename,
                               PETSC_MAX_PATH_LEN, &flg); CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerASCIIOpen(grid.com,filename,&viewer);
             CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"SNES before SSASolve_FE\n");
             CHKERRQ(ierr);
    ierr = SNESView(ctx->snes,viewer);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"solutin vector before SSASolve_FE\n");
             CHKERRQ(ierr);
    ierr = VecView(SSAX,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
  }
  //ierr = PetscLogEventBegin(LOG_SSA_Solve,ssa,0,0,0);CHKERRQ(ierr);
  
  ierr = SSASolve_FE(ctx,SSAX);CHKERRQ(ierr);

  //ierr = PetscLogEventBegin(LOG_SSA_Solve,ssa,0,0,0);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,"-ssa_view_solution",&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerASCIIOpen(grid.com,filename,&viewer);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"solution vector after SSASolve\n");
             CHKERRQ(ierr);
    ierr = VecView(SSAX,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
  }
  //ierr = SSAMapToSplitVecs(ssa,ubar,vbar);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


