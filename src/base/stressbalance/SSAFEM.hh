// Copyright (C) 2009--2011 Jed Brown and Ed Bueler and Constantine Khroulev and David Maxwell
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
#include "FETools.hh"
#include <petscsnes.h>


//! Storage for SSA coefficients at a quadrature point.
struct FEStoreNode {
  PetscReal h,H,tauc,hx,hy,b,B;
};


/*
//! Currently unused class that may be used to implement nondimensionalization in SSAFEM.
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
*/

class SSAFEM;

//! Adaptor for gluing SNESDAFormFunction callbacks to an SSAFEM.
/* The callbacks from SNES are mediated via SNESDAFormFunction, which has the
 convention that its context argument is a pointer to a struct 
 having a DA as its first entry.  The SSAFEM_SNESCallbackData fulfills 
 this requirement, and allows for passing the callback on to an honest 
 SSAFEM object. */
struct SSAFEM_SNESCallbackData {
  DA           da;
  SSAFEM      *ssa;
};

//! SNES callbacks.  
/*! These simply forward the call on to the SSAFEM memeber of the SSAFEM_SNESCallbackData */
PetscErrorCode SSAFEFunction(DALocalInfo *, const PISMVector2 **, 
                                                      PISMVector2 **, SSAFEM_SNESCallbackData *);
PetscErrorCode SSAFEJacobian(DALocalInfo *, const PISMVector2 **, Mat, SSAFEM_SNESCallbackData *);

//! Factory function for constructing a new SSAFEM.
SSA * SSAFEMFactory(IceGrid &, IceBasalResistancePlasticLaw &, 
                  IceFlowLaw &, EnthalpyConverter &, const NCConfigVariable &);

//! PISM's SSA solver: the finite element method implementation written by Jed and David
/*!
Jed's original code is in rev 831: src/base/ssaJed/...
The SSAFEM duplicates the functionality of SSAFD, using the finite element method.
*/
class SSAFEM : public SSA
{
  friend PetscErrorCode SSAFEFunction(DALocalInfo *, const PISMVector2 **, PISMVector2 **, SSAFEM_SNESCallbackData *);
  friend PetscErrorCode SSAFEJacobian(DALocalInfo *, const PISMVector2 **, Mat, SSAFEM_SNESCallbackData *);
public:
  SSAFEM(IceGrid &g, IceBasalResistancePlasticLaw &b, IceFlowLaw &i, EnthalpyConverter &e,
         const NCConfigVariable &c) :
    SSA(g,b,i,e,c), element_index(g)
  {
    quadrature.init(grid);
    allocate_fem();  // can't be done by allocate() since constructor is not virtual
  }

  virtual ~SSAFEM()
  {
    deallocate_fem();
  }

  virtual PetscErrorCode init(PISMVars &vars);

protected:
  PetscErrorCode setup();

  virtual PetscErrorCode PointwiseNuHAndBeta(const FEStoreNode *,
                                             const PISMVector2 *,const PetscReal[],
                                             PetscReal *,PetscReal *,PetscReal *,PetscReal *);

  void FixDirichletValues(PetscReal local_treatment_mask[], PetscReal local_bc_mask[],PISMVector2 **BC_vel,
                          PISMVector2 x[], FEDOFMap &my_dofmap);

  virtual PetscErrorCode allocate_fem();

  virtual PetscErrorCode deallocate_fem();

  virtual PetscErrorCode compute_local_function(DALocalInfo *info, const PISMVector2 **xg, PISMVector2 **yg);

  virtual PetscErrorCode compute_local_jacobian(DALocalInfo *info, const PISMVector2 **xg, Mat J);

  virtual PetscErrorCode solve();
  
  virtual PetscErrorCode setFromOptions();

  // objects used internally
  IceModelVec2S hardav;         // vertically-averaged ice hardness
  SSAFEM_SNESCallbackData callback_data;
  Mat J;
  Vec r;

  SNES         snes;
  FEStoreNode *feStore;
  PetscReal    dirichletScale;
  PetscReal    ocean_rho;
  PetscReal    earth_grav;
  // PismRef      ref;

  FEElementMap element_index;
  FEQuadrature quadrature;
  FEDOFMap dofmap;
};


#endif /* _SSAFEM_H_ */

