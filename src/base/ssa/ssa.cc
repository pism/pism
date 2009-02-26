// Copyright (C) 2009 Jed Brown
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

/**
* \file ssa.cc
*
* \brief The SSA package provides an interface for SSA solvers to be used from IceModel without being part of IceModel.
*
The following is mainly from an email 2009-02-23 from Jed to Ed and Constantine.

The SSA design follows the PETSc creation and usage pattern.  Using this package requires the addition of a single
pointer to IceModel, see <tt>typedef struct _p_SSA *SSA</tt> in pismssa.hh.  All interaction with this object is done
through the interface in pismssa.hh.  Note that the definition of this struct is never provided to IceModel, it's just
an opaque pointer and can only be manipulated through the interface (pismssa.hh).  It starts out NULL which indicates
that it's not being used.  When the -ssa_external option is given (checked in iMssa.hh), an object is created and this
pointer now points to it.  That indicates that it should be used, see IceModel::velocitySSA(PetscInt*), instead of the
solver hard-wired into IceModel.  Note that if -ssa_external is not given, no code in src/base/ssa will be executed.
When this solver is active, running with -help will give all available options.

The SSA object has a fairly small public interface defined in pismssa.hh.  There can be many implementations, but the
user (IceModel) does not need to see them and which implementation is used can be set on the command line.  Like PETSc
objects, the basic design is known as PIMPL (Pointer to IMPLementation) in that the generic functionality (SSA) is a
completely separate object from implementation concerns (SSA_FE).

It follows the usual PETSc creation pattern.  You call SSACreate with an IceGrid which contains the communicator and the
reference DA.  This gives you an object whose type has not been set yet (think of it as a factory at this point, if you
like).  You can set some options and register extra implementations (look at the public interface in pismssa.hh).  Also,
you can call SSASetType() which looks up the type (a string) in the list of registered implementations and sets the
type.  (Now you're working with the actual object, unless you call SSASetType again with different arguments.)

Eventually you call SSASetFromOptions (I should add SSASetOptionsPrefix, but it's probably unlikely anyone will be using
multiple SSA solvers at the same time) which consults the command line and sets the type if it hasn't been set yet or
the user asked for a different type.  Now the object is ready to use, you set the fields (SSA keeps pointers to these so
technically it only needs to be done once, but this operation tells SSA that they may have changed---either the data in
the vector or the location of the vector itself) and call SSASolve.  The last state is saved in the SSA context and used
as the initial iterate the next time SSASolve is called.  Note that you can use -ssa_view to get lots of details about
the solver at the beginning of every solve.

To provide a different implementation, you write something analogous to ssafe.cc and you're set.  No IceModel code needs
to be recompiled.  If someone else is using PISM as a library, they can write another implementation of SSA and use it
with IceModel without changing or recompiling any of PISM.  The implementation can implement some or all of the
implementation's interface in \struct _SSAOps.

I followed this pattern because it serves my purposes well and I'm very familiar with it.  Most of the code in base/ssa
will compile with a standard C compiler.  You can do something similar with an abstract base class and factories, but we
didn't have a convention for how to do this.  On the other hand, PETSc has a proven convention for C-style objects.
Note that I wrote that part before writing IceFactory (in effect creating a convention).

Also, the current design could all be converted to use C++ classes with only cosmetic changes.  SSA would be a class
with no virtual methods and a pointer to a pure virtual class SSAImplBase (which would contain no state and implement no
functions).  Then implementations (SSA_FE) would derive from SSAImplBase and would be required to implement the full
interface (see _SSAOps).

*/

#include "iceModel.hh"
#include "ssaimpl.hh"

PetscCookie SSA_COOKIE;

static PetscFList SSAList = 0;

#if PETSC_VERSION_MAJOR >= 3
static PetscLogEvent LOG_SSA_Solve;
#else
static PetscEvent LOG_SSA_Solve;
#endif

static const char *PismSetupStateName[3] = {"SETUP_GREEN","SETUP_STALE","SETUP_CURRENT"};

class CustomIce : public IceType {
public:
  CustomIce(MPI_Comm c,const char *pre,PismRef &r) : IceType(c,pre), ref(r) {
    // Set defaults
    p = 1.33;
    gamma_reg = 1e-6;
    B0 = 1;
  }
  virtual PetscScalar flow(PetscScalar,PetscScalar,PetscScalar,PetscScalar) const
    {PetscPrintf(comm,"Not implemented"); PetscEnd(); return 0; }
  virtual PetscScalar effectiveViscosityColumn(PetscScalar,PetscInt,const PetscScalar*,PetscScalar,PetscScalar,PetscScalar,PetscScalar,const PetscScalar*,const PetscScalar*) const
    {PetscPrintf(comm,"Not implemented"); PetscEnd(); return 0; }
  #undef __FUNCT__
  #define __FUNCT__ "CustomIce::setFromOptions"
  virtual PetscErrorCode setFromOptions() {
    PetscErrorCode ierr;

    PetscFunctionBegin;
    ierr = PetscOptionsBegin(comm,prefix,"CustomIce options",NULL);CHKERRQ(ierr);
    {
      ierr = PetscOptionsReal("-rheology_p","Exponent in power-law relation",NULL,p,&p,NULL);CHKERRQ(ierr);
      ierr = PetscOptionsReal("-rheology_gamma_reg","Nondimensional regularizing strain rate second invariant",NULL,gamma_reg,&gamma_reg,NULL);CHKERRQ(ierr);
      ierr = PetscOptionsReal("-rheology_B0","Nondimensional viscosity at reference strain rate",NULL,B0,&B0,NULL);CHKERRQ(ierr);
    }
    ierr = PetscOptionsEnd();CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  #undef __FUNCT__
  #define __FUNCT__ "CustomIce::view"
  virtual PetscErrorCode view(PetscViewer viewer) const {
    PetscErrorCode ierr;
    PetscTruth iascii;

    PetscFunctionBegin;
    ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii);CHKERRQ(ierr);
    if (iascii) {
      ierr = PetscViewerASCIIPrintf(viewer,"CustomIce object (%s):\n",prefix);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"Exponent p            %8.1e\n",p);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"Regularization gamma  %8.1e\n",gamma_reg);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"Viscosity at 1 (B0)   %8.1e\n",B0);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
    } else {
      SETERRQ(1,"No binary viewer for this type");
    }
    PetscFunctionReturn(0);
  }
  virtual PetscInt integratedStoreSize() const {return 1;}
  virtual void integratedStore(PetscScalar H, PetscInt /*kbelowH*/, const PetscScalar */*zlevels*/,
                               const PetscScalar /*T*/[], PetscScalar store[]) const
    {store[0] = (H/ref.Height()) * B0 / 2;} // Store nondimensional coefficient
  virtual void integratedViscosity(const PetscScalar store[], const PetscScalar Du[], PetscScalar *nuH, PetscScalar *dNuH) const {
    PetscReal gamma = secondInvariantDu(Du) / ref.StrainRate2(), // nondimensionalize gamma
      gamreg = gamma_reg + gamma,                                // regularized
      nNuH = store[0] * pow(gamreg,(p-2)/2),ndNuH = nNuH / gamreg * (p-2)/2; // nondimensional
    // Return dimensional quantities, SSAFE just undoes this scaling but the rest of Pism is dimensional
    *nuH = nNuH * ref.IntegratedViscosity();
    if (dNuH) *dNuH = ndNuH * ref.IntegratedViscosity() / ref.StrainRate2();
  }
  virtual PetscScalar exponent() const { return 1/(p-1); }
  virtual PetscScalar hardnessParameter(PetscScalar T) const { return B0 * ref.IntegratedViscosity(); } // probably wrong scaling
protected:
  // The second invariant of a symmetric strain rate tensor in compressed form [u_x, v_y, 0.5(u_y+v_x)]
  PetscScalar secondInvariantDu(const PetscScalar Du[]) const
  { return 0.5 * (PetscSqr(Du[0]) + PetscSqr(Du[1]) + PetscSqr(Du[0]+Du[1]) + 2*PetscSqr(Du[2])); }
private:
  PismRef &ref;
  PetscReal p,gamma_reg,B0;
};

#undef __FUNCT__
#define __FUNCT__ "SSARegister"
PetscErrorCode SSARegister(const char sname[],const char path[],const char name[],PetscErrorCode (*function)(SSA))
{
  PetscErrorCode ierr;
  char           fullname[PETSC_MAX_PATH_LEN];

  PetscFunctionBegin;
  ierr = PetscFListConcat(path,name,fullname);CHKERRQ(ierr);
  ierr = PetscFListAdd(&SSAList,sname,fullname,(void (*)(void))function);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// Declarations of the \c create function for all implementations that should be registered by default.
extern PetscErrorCode SSACreate_FE(SSA);

#undef __FUNCT__
#define __FUNCT__ "SSAInitializePackage"
PetscErrorCode SSAInitializePackage(const char path[])
{
  static PetscTruth initialized = PETSC_FALSE;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (initialized) PetscFunctionReturn(0);
  initialized = PETSC_TRUE;
  ierr = SSARegisterDynamic(SSAFE,path,"SSACreate_FE",SSACreate_FE);CHKERRQ(ierr);
#if PETSC_VERSION_MAJOR >= 3
  ierr = PetscCookieRegister("SSA Solver",&SSA_COOKIE);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("SSASolve",SSA_COOKIE,&LOG_SSA_Solve);CHKERRQ(ierr);
#else
  ierr = PetscCookieRegister(&SSA_COOKIE);CHKERRQ(ierr);
  ierr = PetscLogEventRegister(&LOG_SSA_Solve,"SSASolve",SSA_COOKIE);CHKERRQ(ierr);
#endif
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SSAMapToSplitVecs"
// Move velocity from the block_size=2 solution vectors in nondimensional units to the
// split scalar vectors in dimensional (SI) units to be used by the rest of Pism.
PetscErrorCode SSAMapToSplitVecs(SSA ssa,IceModelVec2 &ubar,IceModelVec2 &vbar)
{
  PetscErrorCode   ierr;
  PetscScalar    **u, **v;
  IceGrid         *grid = ssa->grid;
  SSANode        **x;
  PetscScalar      usum=0,vsum=0,c2sum=0,c2max=0;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ssa,SSA_COOKIE,1);
  ierr = ubar.get_array(u);CHKERRQ(ierr);
  ierr = vbar.get_array(v);CHKERRQ(ierr);
  ierr = DAVecGetArray(ssa->da,ssa->x,&x);CHKERRQ(ierr);

  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; i++) {
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; j++) {
      u[i][j] = x[i][j].x * ssa->ref.Velocity();
      v[i][j] = x[i][j].y * ssa->ref.Velocity();
      {
        PetscScalar c2 = PetscSqr(u[i][j]) + PetscSqr(v[i][j]);
        usum += PetscAbs(u[i][j]);
        vsum += PetscAbs(v[i][j]);
        c2sum += c2;
        if (c2max < c2) c2max = c2;
      }
    }
  }

  ierr = DAVecRestoreArray(ssa->da,ssa->x,&x);CHKERRQ(ierr);
  ierr = ubar.end_access();CHKERRQ(ierr);
  ierr = vbar.end_access();CHKERRQ(ierr);

  {
    PetscInt n = grid->xm*grid->ym;
    ierr = verbPrintf(4,grid->com,"SSA solution velocities (m/a): mean|u|=%8.1e  mean|v|=%8.1e  mean(c)=%8.1e  max(c)=%8.1e\n",usum/n*secpera,vsum/n*secpera,sqrt(c2sum)/n*secpera,sqrt(c2max)*secpera);CHKERRQ(ierr);
  }

  // Communicate so that we have stencil width.  I (Jed) am not sure this is needed now, but it will be needed when
  // evaluating fluxes for surface evolution.  This communication happens only once per nonlinear solve so it's not so
  // wasteful anyway ;-)
  ierr = ubar.beginGhostComm();CHKERRQ(ierr);
  ierr = vbar.beginGhostComm();CHKERRQ(ierr);
  ierr = ubar.endGhostComm();CHKERRQ(ierr);
  ierr = vbar.endGhostComm();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "SSACreate"
PetscErrorCode SSACreate(IceGrid *grid,SSA *inssa)
{
  SSA ssa;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidPointer(grid,1);
  PetscValidPointer(inssa,2);
  *inssa = 0;
  ierr = SSAInitializePackage(NULL);CHKERRQ(ierr);
  ierr = PetscHeaderCreate(ssa,_p_SSA,struct _SSAOps,SSA_COOKIE,0,"SSA",grid->com,SSADestroy,SSAView);CHKERRQ(ierr);
  ssa->grid = grid;

  *inssa = ssa;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SSADestroy"
PetscErrorCode SSADestroy(SSA ssa)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ssa,SSA_COOKIE,1);
  if (ssa->ops->Destroy) {ierr = ssa->ops->Destroy(ssa);CHKERRQ(ierr);}
  ierr = VecDestroy(ssa->x);CHKERRQ(ierr);
  ierr = VecDestroy(ssa->r);CHKERRQ(ierr);
  ierr = VecDestroy(ssa->siaVelLocal);CHKERRQ(ierr);
  ierr = MatDestroy(ssa->J);CHKERRQ(ierr);
  ierr = DADestroy(ssa->da);CHKERRQ(ierr);
  {
    // Only delete our ice if we created it.  Since our CustomIce is not public (the one in materials.hh is
    // CustomGlenIce), it's not possible for IceType to give us a CustomIce.
    CustomIce *c = dynamic_cast<CustomIce*>(ssa->ice);
    if (c) delete c;
  }
  ierr = PetscHeaderDestroy(ssa);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SSASetIceType"
PetscErrorCode SSASetIceType(SSA ssa,IceType *ice)
{

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ssa,SSA_COOKIE,1);
  ssa->ice = ice;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SSASetBasalType"
PetscErrorCode SSASetBasalType(SSA ssa,PlasticBasalType *basal)
{

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ssa,SSA_COOKIE,1);
  ssa->basal = basal;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SSASetOceanType"
PetscErrorCode SSASetOceanType(SSA ssa,SeaWaterType *ocean)
{

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ssa,SSA_COOKIE,1);
  ssa->ocean = ocean;
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "SSASetFictitiousNuH"
PetscErrorCode SSASetShelfExtension(SSA ssa,IceShelfExtension *sext)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ssa,SSA_COOKIE,1);
  ssa->shelfExtension = sext;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SSASetFromOptions"
PetscErrorCode SSASetFromOptions(SSA ssa)
{
  PetscErrorCode ierr;
  char type[256];
  PetscTruth optionsSetType,customIce,periodic;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ssa,SSA_COOKIE,1);
  customIce           = PETSC_FALSE;
  periodic            = PETSC_FALSE;
  ierr = PetscOptionsBegin(((PetscObject)ssa)->comm,((PetscObject)ssa)->prefix,"Shallow Stream/Shelf Approximation solver options","SSA");CHKERRQ(ierr);
  {
    ierr = PetscOptionsList("-ssa_type","SSA solver","SSASetType",SSAList,((PetscObject)ssa)->type_name?((PetscObject)ssa)->type_name:SSAFE,type,sizeof(type),&optionsSetType);CHKERRQ(ierr);
    {
      // Change the type if set on the command line
      if (optionsSetType) {ierr = SSASetType(ssa,type);CHKERRQ(ierr);}
      // Use default if still unset (user code did not set it and neither did the command line)
      if (!((PetscObject)ssa)->type_name) {ierr = SSASetType(ssa,SSAFE);CHKERRQ(ierr);}
    }
    ierr = PetscOptionsTruth("-ssa_da_periodic","Use DA_XYPERIODIC for SSA (instead of DA_NONPERIODIC with natural boundary conditions","",periodic,&periodic,NULL);CHKERRQ(ierr);
    ssa->wrap = periodic ? DA_XYPERIODIC : DA_NONPERIODIC;
    ierr = PetscOptionsTruth("-ssa_initial_guess_nonzero","use nonzero inital velocity","",ssa->initialGuessNonzero,&ssa->initialGuessNonzero,NULL);CHKERRQ(ierr);
    // The following boundary condition options are all off by default
    ierr = PetscOptionsTruth("-ssa_boundary_floating_stress_free",
                             "Treat floating margins as stress-free instead of as calving face","",
                             ssa->boundary.floating_stress_free,&ssa->boundary.floating_stress_free,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsTruth("-ssa_boundary_grounded_as_floating",
                             "Treat grounded margins as though they were floating instead of as usual calving face","",
                             ssa->boundary.grounded_as_floating,&ssa->boundary.grounded_as_floating,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsTruth("-ssa_boundary_submarine_stress_free",
                             "Treat submarine grounded margins as stress-free","",
                             ssa->boundary.submarine_stress_free,&ssa->boundary.submarine_stress_free,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsTruth("-ssa_boundary_calving_above_sea_level",
                             "Treat margins above sea level as calving faces","",
                             ssa->boundary.calving_above_sea_level,&ssa->boundary.calving_above_sea_level,NULL);CHKERRQ(ierr);
    // This option can be removed unless we want to enable using completely different ice types in the SSA solver from
    // the rest of the model.
    ierr = PetscOptionsTruth("-ssa_rheology_custom","Use custom ice type","",customIce,&customIce,NULL);CHKERRQ(ierr);
    if (customIce) {
      ssa->ice = new CustomIce(((PetscObject)ssa)->comm,((PetscObject)ssa)->prefix,ssa->ref);
    }
    if (ssa->ops->SetFromOptions) {
      ierr = ssa->ops->SetFromOptions(ssa);CHKERRQ(ierr);
    }
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  {
    CustomIce *c = dynamic_cast<CustomIce*>(ssa->ice);
    // Only call setFromOptions if we created the ice---so it hasn't been called already
    if (c) {
      ierr = c->setFromOptions();CHKERRQ(ierr);
    }
  }

  ierr = SSASetUp(ssa);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SSASetUp"
PetscErrorCode SSASetUp(SSA ssa)
{
  PetscErrorCode  ierr;
  IceGrid        *grid = ssa->grid;
  char            mtype[256] = "baij";

  PetscFunctionBegin;
  // 2 Dofs, nonperiodic.  Otherwise the same as grid.da
  ierr = DACreate2d(grid->com,ssa->wrap,DA_STENCIL_BOX,grid->My,grid->Mx,PETSC_DECIDE,PETSC_DECIDE,2,1,NULL,NULL,&ssa->da);CHKERRQ(ierr);
  ierr = DACreateGlobalVector(ssa->da,&ssa->x);CHKERRQ(ierr);
  ierr = VecDuplicate(ssa->x,&ssa->r);CHKERRQ(ierr);
  ierr = DACreateLocalVector(ssa->da,&ssa->siaVelLocal);CHKERRQ(ierr);
  ierr = PetscOptionsGetString(((PetscObject)ssa)->prefix,"-ssa_mat_type",mtype,sizeof(mtype),NULL);CHKERRQ(ierr);
  ierr = DAGetMatrix(ssa->da,mtype,&ssa->J);CHKERRQ(ierr);
  ierr = MatSetOptionsPrefix(ssa->J,"ssa_");CHKERRQ(ierr);

  ssa->ref.SetUp();

  // The SSA context will normally hold the state from the last solve.  On the first iteration, however, we need an
  // initial guess.
  if (ssa->initialGuessNonzero) {
    SSANode *x;
    PetscInt bs,m;
    ierr = VecGetArray(ssa->x,(PetscScalar**)&x);CHKERRQ(ierr);
    ierr = VecGetBlockSize(ssa->x,&bs);CHKERRQ(ierr);
    ierr = VecGetLocalSize(ssa->x,&m);CHKERRQ(ierr);
    for (PetscInt i=0; i<m/bs; i++) { // Velocity is nondimensional
      x[i].x = 1e-8;
      x[i].y = 0;
    }
    ierr = VecRestoreArray(ssa->x,(PetscScalar**)&x);CHKERRQ(ierr);
  } else {
    ierr = VecZeroEntries(ssa->x);CHKERRQ(ierr);
  }

  if (ssa->ops->SetUp) {
    ierr = ssa->ops->SetUp(ssa);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SSASetType"
PetscErrorCode SSASetType(SSA ssa,const SSAType type)
{
  PetscErrorCode ierr,(*r)(SSA);
  PetscTruth match;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ssa,SSA_COOKIE,1);
  PetscValidCharPointer(type,2);

  ierr = PetscTypeCompare((PetscObject)ssa,type,&match);CHKERRQ(ierr);
  if (match) PetscFunctionReturn(0);
  ierr =  PetscFListFind(SSAList,((PetscObject)ssa)->comm,type,(void(**)(void))&r);CHKERRQ(ierr);
  if (!r) SETERRQ1(PETSC_ERR_ARG_UNKNOWN_TYPE,"Unable to find requested SSA type %s",type);
  /* Destroy the previous private SSA context */
  if (ssa->ops->Destroy) { ierr = (*ssa->ops->Destroy)(ssa);CHKERRQ(ierr); }
  /* Reinitialize function pointers in SSAOps structure */
  ierr = PetscMemzero(ssa->ops,sizeof(struct _SSAOps));CHKERRQ(ierr);
  if (!ssa->ice) SETERRQ(1,"Must set ice type before choosing SSA solver type");
  /* Call the SSACreate_XXX routine for this particular SSA solver */
  ssa->setupcalled = SETUP_GREEN;
  ierr = (*r)(ssa);CHKERRQ(ierr);
  ierr = PetscObjectChangeTypeName((PetscObject)ssa,type);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SSASetFields"
PetscErrorCode SSASetFields(SSA ssa,IceModelVec2 *mask,IceModelVec2 uvbar[],IceModelVec2 *H,IceModelVec2 *h,IceModelVec2 *bed,IceModelVec2 *tauc,IceModelVec3 *T)
{

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ssa,SSA_COOKIE,1);
  ssa->mask    = mask;
  ssa->siaVel  = uvbar;
  ssa->H       = H;
  ssa->h       = h;
  ssa->bed     = bed;
  ssa->tauc    = tauc;
  ssa->T       = T;
  ssa->setupcalled = SETUP_STALE;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SSAUpdateNodalSIAVelocity"
PetscErrorCode SSAUpdateNodalSIAVelocity(SSA ssa)
{
  PetscErrorCode   ierr;
  PetscScalar    **siaVelX,**siaVelY;
  SSANode        **siaVel;
  PetscInt         i,j;
  DALocalInfo      info;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ssa,SSA_COOKIE,1);
  ierr = DAGetLocalInfo(ssa->da,&info);CHKERRQ(ierr);
  ierr = DAVecGetArray(ssa->da,ssa->siaVelLocal,&siaVel);CHKERRQ(ierr);
  ierr = ssa->siaVel[0].get_array(siaVelX);CHKERRQ(ierr);
  ierr = ssa->siaVel[1].get_array(siaVelY);CHKERRQ(ierr);
  // The concept of \a x and \a y is really confusing here.  We want to take inward-looking values off the staggered
  // grid when we are at a boundary.  Petsc-native indexing is [j][i] and we index that way because it is compatible
  // with the grid object, however the \a x and \a y components retain their usual meaning within the Pism code (though
  // not with the command-line and output-file definitions).  Cry, cry.
  for (j=info.ys; j<info.ys+info.ym; j++) {
    for (i=info.xs; i<info.xs+info.xm; i++) {
      // Note that when the grid is periodic we will always use the average of the two staggered points
      if (j == info.gys)        // Left margin, use the staggered point to the right
        siaVel[j][i].x = siaVelX[j][i];
      else if (j == info.gys+info.gym-1) // Right margin, use the staggered point to the left
        siaVel[j][i].x = siaVelX[j-1][i];
      else                      // Interior, take the average
        siaVel[j][i].x = 0.5*(siaVelX[j-1][i] + siaVelX[j][i]);
      if (i == info.gxs)        // Bottom margin, use the staggered point above
        siaVel[j][i].y = siaVelY[j][i];
      else if (i == info.gxs+info.gxm-1) // Top margin, use the staggered point below
        siaVel[j][i].y = siaVelY[j][i-1];
      else                      // Interior, take the average
        siaVel[j][i].y = 0.5*(siaVelY[j-1][i] + siaVelY[j][i]);
    }
  }
  ierr = DAVecRestoreArray(ssa->da,ssa->siaVelLocal,&siaVel);CHKERRQ(ierr);
  ierr = ssa->siaVel[0].end_access();CHKERRQ(ierr);
  ierr = ssa->siaVel[1].end_access();CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(ssa->da,ssa->siaVelLocal,INSERT_VALUES,ssa->siaVelLocal);CHKERRQ(ierr);
  ierr = DALocalToLocalEnd  (ssa->da,ssa->siaVelLocal,INSERT_VALUES,ssa->siaVelLocal);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SSASolve"
PetscErrorCode SSASolve(SSA ssa,IceModelVec2 &ubar,IceModelVec2 &vbar)
{
  PetscErrorCode ierr;
  PetscViewer    viewer;
  char           filename[PETSC_MAX_PATH_LEN];
  PetscTruth     flg;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ssa,SSA_COOKIE,1);
  if (!ssa->ice)            SETERRQ(1,"The ice type for the SSA solver must be set with SSASetIceType\n");
  if (!ssa->basal)          SETERRQ(1,"The basal type for the SSA solver must be set with SSASetBasalType\n");
  if (!ssa->ocean)          SETERRQ(1,"The ocean type for the SSA solver must be set with SSASetOceanType\n");
  if (!ssa->shelfExtension) SETERRQ(1,"The ice shelf extension for the SSA solver must be set with SSASetShelfExtension\n");
  if (ssa->setupcalled < SETUP_CURRENT) {
    ierr = SSAUpdateNodalSIAVelocity(ssa);CHKERRQ(ierr);
    // Debugging discovered that these were not up to date.  Since we (IceModelVec) control access to modification of
    // these vectors, we should be able to cache whether they are current.  Then the operation below would be a no-op if
    // it has already been done.
    ierr = ssa->bed->beginGhostComm();CHKERRQ(ierr);
    ierr = ssa->tauc->beginGhostComm();CHKERRQ(ierr);
    ierr = ssa->bed->endGhostComm();CHKERRQ(ierr);
    ierr = ssa->tauc->endGhostComm();CHKERRQ(ierr);
  }
  ierr = PetscOptionsGetString(((PetscObject)ssa)->prefix,"-ssa_view",filename,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerASCIIOpen(((PetscObject)ssa)->comm,filename,&viewer);CHKERRQ(ierr);
    ierr = SSAView(ssa,viewer);CHKERRQ(ierr);
    ierr = VecView(ssa->x,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
  }
  ierr = PetscLogEventBegin(LOG_SSA_Solve,ssa,0,0,0);CHKERRQ(ierr);
  ierr = ssa->ops->Solve(ssa);CHKERRQ(ierr);
  ierr = PetscLogEventBegin(LOG_SSA_Solve,ssa,0,0,0);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(((PetscObject)ssa)->prefix,"-ssa_view_solution",&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerASCIIOpen(((PetscObject)ssa)->comm,filename,&viewer);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"Solution vector after SSASolve\n");CHKERRQ(ierr);
    ierr = VecView(ssa->x,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
  }
  ierr = SSAMapToSplitVecs(ssa,ubar,vbar);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SSAView"
PetscErrorCode SSAView(SSA ssa,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscTruth iascii;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ssa,SSA_COOKIE,1);
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(((PetscObject)ssa)->comm,&viewer);CHKERRQ(ierr);
  }
  PetscCheckSameComm(ssa,1,viewer,2);
  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii);CHKERRQ(ierr);
  if (iascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"SSA object:(%s)\n",
                                  ((PetscObject)ssa)->prefix ? ((PetscObject)ssa)->prefix : "no prefix");CHKERRQ(ierr);
    ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"type: %s\n",
                                  ((PetscObject)ssa)->type_name ? ((PetscObject)ssa)->type_name : "type not set");CHKERRQ(ierr);
    if (ssa->ops->View) {
      ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
      ierr = ssa->ops->View(ssa,viewer);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer,"setup state: %s\n",PismSetupStateName[ssa->setupcalled]);CHKERRQ(ierr);
    ierr = ssa->ice->view(viewer);CHKERRQ(ierr);
    ierr = ssa->shelfExtension->view(viewer);CHKERRQ(ierr);
    ierr = ssa->ref.View(viewer);CHKERRQ(ierr);
    ierr = DAView(ssa->da,viewer);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
  } else if (ssa->ops->View) {
    ierr = ssa->ops->View(ssa,viewer);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PismRef::View"
PetscErrorCode PismRef::View(PetscViewer viewer) const
{
  PetscErrorCode ierr;
  PetscTruth iascii;

  PetscFunctionBegin;
  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii);CHKERRQ(ierr);
  if (iascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"PismRef object\n");CHKERRQ(ierr);
    ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"Length               %8.1e meters\n",Length());CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"Height               %8.1e meters\n",Height());CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"Time                 %8.1e seconds\n",Time());CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"Velocity             %8.1e meters/second\n",Velocity());CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"VerticalVelocity     %8.1e meters/second\n",VerticalVelocity());CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"StrainRate           %8.1e 1/seconds\n",StrainRate());CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"DrivingStress        %8.1e Pascal\n",DrivingStress());CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"IntegratedViscosity  %8.1e Pascal meter seconds\n",IntegratedViscosity());CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"Drag                 %8.1e Pascal / (meter second)\n",Drag());CHKERRQ(ierr);
    ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
  } else {
    SETERRQ(1,"No binary viewer for this type");
  }
  PetscFunctionReturn(0);
}
