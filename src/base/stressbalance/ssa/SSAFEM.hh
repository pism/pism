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

namespace pism {

//! Storage for SSA coefficients at a quadrature point.
struct SSACoefficients {
  double H,                     //!< ice thickness
    tauc,                       //!< basal yield stress
    b,                          //!< bed elevation
    B;                          //!< ice hardness
  Vector2 driving_stress;
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
PetscErrorCode SSAFEFunction(DMDALocalInfo *, const Vector2 **, 
                             Vector2 **, SSAFEM_SNESCallbackData *);
PetscErrorCode SSAFEJacobian(DMDALocalInfo *info, const Vector2 **xg,
                             Mat A, Mat J,
                             MatStructure *str, SSAFEM_SNESCallbackData *fe);

//! Factory function for constructing a new SSAFEM.
SSA * SSAFEMFactory(IceGrid &, EnthalpyConverter &, const Config &);

//! PISM's SSA solver: the finite element method implementation written by Jed and David
/*!
  Jed's original code is in rev 831: src/base/ssaJed/...
  The SSAFEM duplicates most of the functionality of SSAFD, using the finite element method.
*/
class SSAFEM : public SSA
{
  friend PetscErrorCode pism::SSAFEFunction(DMDALocalInfo *, const Vector2 **, Vector2 **, SSAFEM_SNESCallbackData *);
  friend PetscErrorCode pism::SSAFEJacobian(DMDALocalInfo *info, const Vector2 **xg,
                                            Mat A, Mat J,
                                            MatStructure *str, SSAFEM_SNESCallbackData *fe);
public:
  SSAFEM(IceGrid &g, EnthalpyConverter &e, const Config &c);

  virtual ~SSAFEM();

  virtual PetscErrorCode init(Vars &vars);

  virtual PetscErrorCode cacheQuadPtValues();

protected:

  virtual PetscErrorCode PointwiseNuHAndBeta(const SSACoefficients &,
                                             const Vector2 &, const double[],
                                             double *,double *,double *,double *);

  virtual PetscErrorCode allocate_fem();

  virtual PetscErrorCode deallocate_fem();

  virtual PetscErrorCode compute_local_function(DMDALocalInfo *info, const Vector2 **xg, Vector2 **yg);

  virtual PetscErrorCode compute_local_jacobian(DMDALocalInfo *info, const Vector2 **xg, Mat J);

  virtual PetscErrorCode solve();

  virtual PetscErrorCode solve(TerminationReason::Ptr &reason);

  virtual PetscErrorCode solve_nocache(TerminationReason::Ptr &reason);
  
  virtual PetscErrorCode setFromOptions();


  // objects used internally
  SSAFEM_SNESCallbackData m_callback_data;

  SNES         m_snes;
  SSACoefficients *m_coefficients;
  double    m_dirichletScale;
  double    m_ocean_rho;
  double    m_earth_grav;
  double    m_beta_ice_free_bedrock;
  double    m_epsilon_ssa;

  FEElementMap m_element_index;
  FEQuadrature_Scalar m_quadrature;
  FEQuadrature_Vector m_quadrature_vector;
  FEDOFMap m_dofmap;

private:
  PetscErrorCode monitor_jacobian(Mat Jac);
  PetscErrorCode monitor_function(const Vector2 **velocity_global,
                                  Vector2 **residual_global);
};


} // end of namespace pism

#endif /* _SSAFEM_H_ */

