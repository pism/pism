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

#ifndef _SSAFEM_H_
#define _SSAFEM_H_

#include "SSA.hh"
#include <petscsnes.h>

struct FECTX;  // the context we will pass to PETSc is defined below, but we
               // need to refer to it in the SSAFEM class
struct FEStoreNode;
struct SSANode;
struct SSABoundaryOptions {
  PetscTruth floating_stress_free,
             grounded_as_floating,
             submarine_stress_free,
             calving_above_sea_level;
};

//! PISM's SSA solver: the finite element method implementation written by Jed
/*!
Jed's code is in rev 831:
  src/base/ssaJed/*
The following is a wrapper around Jed's code.  The wrapper duplicates the
functionality of SSAFD.
 */
class SSAFEM : public SSA
{
  friend PetscErrorCode SSAFESetUp(FECTX *);
  friend PetscErrorCode PointwiseNuHAndBeta(FECTX *,const FEStoreNode *,const PetscReal *,
                                            const SSANode *,const PetscReal[],
                                            PetscReal *,PetscReal *,PetscReal *,PetscReal *);
  friend PetscErrorCode SSAFEFunction(DALocalInfo *, const SSANode **, SSANode **, FECTX *);
  friend PetscErrorCode SSAFEJacobian(DALocalInfo *, const SSANode **, Mat, FECTX *);
  friend PetscErrorCode ApplyBoundaryStress(FECTX *, PetscReal [], const IceFlowLaw *,
                                            PetscReal, PetscReal **, PetscReal **,
                                            const MatStencil [], PetscInt, SSANode []);
  friend PetscErrorCode SSASetUp_FE(FECTX*);
public:
  SSAFEM(IceGrid &g, IceBasalResistancePlasticLaw &b, IceFlowLaw &i, EnthalpyConverter &e,
         const NCConfigVariable &c) :
    SSA(g,b,i,e,c)
  {
    allocate_fem();  // can't be done by allocate() since constructor is not virtual
    wrap = DA_XYPERIODIC;       // FIXME: DA is always XYPERIODIC; "wrap"
                                // should be removed
  }

  virtual ~SSAFEM()
  {
    deallocate_fem();
  }

  virtual PetscErrorCode init(PISMVars &vars) {
    PetscErrorCode ierr;
    ierr = SSA::init(vars); CHKERRQ(ierr);
    ierr = verbPrintf(2,grid.com,
                      "  [using the finite element method implementation by Jed Brown]\n"); CHKERRQ(ierr);

    mask = dynamic_cast<IceModelVec2Mask*>(vars.get("mask"));
    if (mask == NULL) SETERRQ(1, "mask is not available");

    thickness = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
    if (thickness == NULL) SETERRQ(1, "land_ice_thickness is not available");

    tauc = dynamic_cast<IceModelVec2S*>(vars.get("tauc"));
    if (tauc == NULL) SETERRQ(1, "tauc is not available");

    surface = dynamic_cast<IceModelVec2S*>(vars.get("surface_altitude"));
    if (surface == NULL) SETERRQ(1, "surface_altitude is not available");

    bed = dynamic_cast<IceModelVec2S*>(vars.get("bedrock_altitude"));
    if (bed == NULL) SETERRQ(1, "bedrock_altitude is not available");

    enthalpy = dynamic_cast<IceModelVec3*>(vars.get("enthalpy"));
    if (enthalpy == NULL) SETERRQ(1, "enthalpy is not available");

    return 0;
  }

protected:
  virtual PetscErrorCode allocate_fem();

  virtual PetscErrorCode deallocate_fem();

  virtual PetscErrorCode solve();

  // objects used internally
  FECTX *ctx;
  SSABoundaryOptions boundary;  // FIXME: needs to be initialized somewhere
  DAPeriodicType wrap;          // FIXME: should be removed
  Mat J;                        // FIXME!
  Vec r;                        // FIXME!
};

struct SSANode {
  PetscScalar x,y;
};

struct FEStoreNode {
  PetscReal h,H,tauc,hx,hy,b; // These values are fully dimensional
};

// Manage nondimensionalization [FIXME:  I've forced a degenerate dimensional version]
class PismRef {
public:
  PismRef() { SetUp(); }
  void SetUp() {
    length = 1;
    height = 1;
    time = 1;
    pressure = 1;
  }
  PetscReal Length() const { return length; }
  PetscReal Area() const { return length*length; }
  PetscReal Height() const { return height; }
  PetscReal Time() const { return time; }
  PetscReal Velocity() const { return length/time; }
  PetscReal VerticalVelocity() const { return height/time; }
  PetscReal StrainRate() const { return 1/time; }
  PetscReal Velocity2() const { PetscReal v = Velocity(); return 1.0; }
  PetscReal StrainRate2() const { PetscReal s = StrainRate(); return 1.0; }
  PetscReal Slope() const { return height / length; }
  PetscReal Pressure() const { return pressure; }
  PetscReal DrivingStress() const { return Pressure() * Slope(); }
  PetscReal IntegratedViscosity() const { return DrivingStress() * Length() / StrainRate(); }
  PetscReal Drag() const { return DrivingStress() / Velocity(); }
private:
  PetscReal length,height,time,pressure;
};


/*
enum PismSetupState { SETUP_GREEN, SETUP_STALE, SETUP_CURRENT };
*/



//! Context for Jed's FEM implementation of SSA.
/*!
The first element of this struct *must* be a DA, because of how SNESSetFunction
and SNESSetJacobian use their last arguments.  (See SSASetUp_FE().)
 */
struct FECTX {
  DA           da;
  SNES         snes;
  FEStoreNode *feStore;
  PetscScalar *integratedStore; // Storage for constitutive relation
  PetscInt     sbs;             // Store block size (number of values per quadrature point)
                                // FIXME:  how to initialize correctly?
  PetscReal    dirichletScale;
  IceFlowLaw  *ice;             // constitutive relation; FIXME: initialize!
  PetscReal    ocean_rho;
  PetscReal    earth_grav;
  PismRef      ref;
  IceGrid     *grid;
  SSAFEM      *ssa;
};
  

#endif /* _SSAFEM_H_ */

