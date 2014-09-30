// Copyright (C) 2009--2014 Jed Brown and Ed Bueler and Constantine Khroulev and David Maxwell
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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
#include "TerminationReason.hh"

//! Storage for SSA coefficients at a quadrature point.
struct FEStoreNode {
  double H,tauc,b,B;
  PISMVector2 driving_stress;
  int mask;
};


class SSAFEM;

//! Adaptor for gluing SNESDAFormFunction callbacks to an SSAFEM.
/* The callbacks from SNES are mediated via SNESDAFormFunction, which has the
 convention that its context argument is a pointer to a struct 
 having a DA as its first entry.  The SSAFEM_SNESCallbackData fulfills 
 this requirement, and allows for passing the callback on to an honest 
 SSAFEM object. */
struct SSAFEM_SNESCallbackData {
  DM           da;
  SSAFEM      *ssa;
};

//! SNES callbacks.  
/*! These simply forward the call on to the SSAFEM memeber of the SSAFEM_SNESCallbackData */
#if PETSC_VERSION_LT(3,5,0)
PetscErrorCode SSAFEFunction(DMDALocalInfo *, const PISMVector2 **, 
                             PISMVector2 **, SSAFEM_SNESCallbackData *);
PetscErrorCode SSAFEJacobian(DMDALocalInfo *info, const PISMVector2 **xg,
                             Mat A, Mat J,
                             MatStructure *str, SSAFEM_SNESCallbackData *fe);
#else
PetscErrorCode SSAFEFunction(DMDALocalInfo *, const PISMVector2 **,
                             PISMVector2 **, SSAFEM_SNESCallbackData *);
PetscErrorCode SSAFEJacobian(DMDALocalInfo *info, const PISMVector2 **xg,
                             Mat A, Mat J, SSAFEM_SNESCallbackData *fe);
#endif

//! Factory function for constructing a new SSAFEM.
SSA * SSAFEMFactory(IceGrid &, EnthalpyConverter &, const PISMConfig &);

//! PISM's SSA solver: the finite element method implementation written by Jed and David
/*!
Jed's original code is in rev 831: src/base/ssaJed/...
The SSAFEM duplicates most of the functionality of SSAFD, using the finite element method.
*/
class SSAFEM : public SSA
{
#if PETSC_VERSION_LT(3,5,0)
  friend PetscErrorCode SSAFEFunction(DMDALocalInfo *, const PISMVector2 **,
                                      PISMVector2 **, SSAFEM_SNESCallbackData *);
  friend PetscErrorCode SSAFEJacobian(DMDALocalInfo *info, const PISMVector2 **xg,
                                      Mat A, Mat J,
                                      MatStructure *str, SSAFEM_SNESCallbackData *fe);
#else
  friend PetscErrorCode SSAFEFunction(DMDALocalInfo *, const PISMVector2 **,
                                      PISMVector2 **, SSAFEM_SNESCallbackData *);
  friend PetscErrorCode SSAFEJacobian(DMDALocalInfo *info, const PISMVector2 **xg,
                                      Mat A, Mat J, SSAFEM_SNESCallbackData *fe);
#endif
public:
  SSAFEM(IceGrid &g, EnthalpyConverter &e, const PISMConfig &c);

  virtual ~SSAFEM();

  virtual PetscErrorCode init(PISMVars &vars);

  virtual PetscErrorCode cacheQuadPtValues();

protected:

  virtual PetscErrorCode PointwiseNuHAndBeta(const FEStoreNode *,
                                             const PISMVector2 *,const double[],
                                             double *,double *,double *,double *);

  void FixDirichletValues(double local_bc_mask[],PISMVector2 **BC_vel,
                          PISMVector2 x[], FEDOFMap &my_dofmap);

  virtual PetscErrorCode allocate_fem();

  virtual PetscErrorCode deallocate_fem();

  virtual PetscErrorCode compute_local_function(DMDALocalInfo *info, const PISMVector2 **xg, PISMVector2 **yg);

  virtual PetscErrorCode compute_local_jacobian(DMDALocalInfo *info, const PISMVector2 **xg, Mat J);

  virtual PetscErrorCode solve();

  virtual PetscErrorCode solve(TerminationReason::Ptr &reason);

  virtual PetscErrorCode solve_nocache(TerminationReason::Ptr &reason);
  
  virtual PetscErrorCode setFromOptions();


  // objects used internally
  SSAFEM_SNESCallbackData callback_data;

  SNES         snes;
  FEStoreNode *feStore;
  double    dirichletScale;
  double    ocean_rho;
  double    earth_grav;
  double    m_beta_ice_free_bedrock;
  double    m_epsilon_ssa;

  FEElementMap element_index;
  FEQuadrature quadrature;
  FEDOFMap dofmap;
};


#endif /* _SSAFEM_H_ */

