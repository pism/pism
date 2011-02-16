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

struct FECTX;
struct FEStoreNode {
  PetscReal h,H,tauc,hx,hy,b; // These values are fully dimensional
};



// Manage nondimensionalization [FIXME:  I've forced a degenerate dimensional version]
class PismRef {
public:
  PismRef() {
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


//! The callbacks from SNES are mediated via SNESDAFormFunction, which has the
//! convention that its context argument is a pointer to a struct argument 
//! having a DA as its first entry.  The FECTX fulfills this requirement, and
//! allows for passing the callback on to an honest SSAFEM object.
class SSAFEM;
struct FECTX {
  DA           da;
  SSAFEM      *ssa;
};

//! SNES callbacks.  These simply forward the call on to an SSAFEM.
PetscErrorCode SSAFEFunction(DALocalInfo *, const PISMVector2 **, 
                                                      PISMVector2 **, FECTX *);
PetscErrorCode SSAFEJacobian(DALocalInfo *, const PISMVector2 **, Mat, FECTX *);


class IceGridElementIndexer 
{
public:
  IceGridElementIndexer(const IceGrid &g);
  
  PetscInt element_count()
  {
    return (xm-xs)*(ym-ys);
  }
  
  PetscInt flatten(PetscInt i, PetscInt j)
  {
    return (i-xs)*ym+(j-ys);
  }
  
  PetscInt xs, xm, ys, ym;
  
};

//! PISM's SSA solver: the finite element method implementation written by Jed
/*!
Jed's code is in rev 831:
  src/base/ssaJed/...
The following is a wrapper around Jed's code.  The wrapper duplicates the
functionality of SSAFD.
 */
class SSAFEM : public SSA
{
  friend PetscErrorCode SSAFEFunction(DALocalInfo *, const PISMVector2 **, PISMVector2 **, FECTX *);
  friend PetscErrorCode SSAFEJacobian(DALocalInfo *, const PISMVector2 **, Mat, FECTX *);
public:
  SSAFEM(IceGrid &g, IceBasalResistancePlasticLaw &b, IceFlowLaw &i, EnthalpyConverter &e,
         const NCConfigVariable &c) :
    SSA(g,b,i,e,c), element_index(g)
  {
    allocate_fem();  // can't be done by allocate() since constructor is not virtual
  }

  virtual ~SSAFEM()
  {
    deallocate_fem();
  }

  virtual PetscErrorCode init(PISMVars &vars);

protected:
  PetscErrorCode setup();

  virtual PetscErrorCode PointwiseNuHAndBeta(const FEStoreNode *,const PetscReal *,
                                             const PISMVector2 *,const PetscReal[],
                                             PetscReal *,PetscReal *,PetscReal *,PetscReal *);

  virtual void FixDirichletValues(PetscReal lmask[],PISMVector2 **BC_vel,
                                  MatStencil row[],MatStencil col[],PISMVector2 x[]);

  virtual PetscErrorCode allocate_fem();

  virtual PetscErrorCode deallocate_fem();

  virtual PetscErrorCode compute_local_function(DALocalInfo *info, const PISMVector2 **xg, PISMVector2 **yg);

  virtual PetscErrorCode compute_local_jacobian(DALocalInfo *info, const PISMVector2 **xg, Mat J);

  virtual PetscErrorCode solve();
  
  virtual PetscErrorCode compute_hardav(IceModelVec2S &result);

  virtual PetscErrorCode view(PetscViewer viewer);

  virtual PetscErrorCode setFromOptions();

  // objects used internally
  IceModelVec2S hardav;         // vertically-averaged ice hardness
  FECTX ctx;
  Mat J;
  Vec r;

  SNES         snes;
  FEStoreNode *feStore;
  PetscScalar *integratedStore; // Storage for constitutive relation
  PetscInt     sbs;             // Store block size (number of values per quadrature point)
                                // FIXME:  how to initialize correctly?
  PetscReal    dirichletScale;
  PetscReal    ocean_rho;
  PetscReal    earth_grav;
  PismRef      ref;

  IceGridElementIndexer element_index;

};


#endif /* _SSAFEM_H_ */

