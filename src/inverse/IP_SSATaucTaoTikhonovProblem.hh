// Copyright (C) 2012, 2014  David Maxwell and Constantine Khroulev
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


#ifndef IP_SSATAUCTAOTIKHONOVPROBLEM_HH_HB8UWICX
#define IP_SSATAUCTAOTIKHONOVPROBLEM_HH_HB8UWICX

#ifdef PISM_USE_TR1
#include <tr1/memory>
#else
#include <memory>
#endif

#include "IPTaoTikhonovProblem.hh"
#include "IP_SSATaucForwardProblem.hh"


#include "TaoUtil.hh"
#include "functional/IPFunctional.hh"

namespace pism {

//! Defines an IPTaoTikhonovProblem for inversion of basal yeild stresses \f$\tau_c\f$ from %SSA velocities.
/*! The forward problem for the inversion is defined by an IP_SSATaucForwardProblem.  The problem itself
  is solved with a TaoBasicSolver as described by the class-level documentation for IPTaoTikhonovProblem.  It
  is a reduced space method, inasmuch as we are performing unconstrained minimization on a Tikhonov functional
  that depends on a design variable (the value of \f$\tau_c\f$). It is compatible with any of the elementary
  TAO minimization algorithms, e.g. tao_cg, tao_lmvm.  If the minimization algorithm tao_blmvm is selected,
  the values of \f$\tau_c\f$ will be constrained by the config variables \a inv_ssa_tauc_min
  and \a inv_ssa_tauc_max.  

  The TAO algorithm tao_lcl is not compatible with IP_SSATaucTaoTikhonovProblem.  Use IP_SSATaucTaoTikhonovProblemLCL
  instead.
*/
class IP_SSATaucTaoTikhonovProblem: public IPTaoTikhonovProblem<IP_SSATaucForwardProblem> {
public:

  IP_SSATaucTaoTikhonovProblem(IP_SSATaucForwardProblem &forward, 
                                IP_SSATaucTaoTikhonovProblem::DesignVec &d0, 
                                IP_SSATaucTaoTikhonovProblem::StateVec &u_obs, double eta, 
                                IPFunctional<IP_SSATaucTaoTikhonovProblem::DesignVec>&designFunctional, 
                                IPFunctional<IP_SSATaucTaoTikhonovProblem::StateVec>&stateFunctional) :
    IPTaoTikhonovProblem<IP_SSATaucForwardProblem>(forward,d0,u_obs,eta,designFunctional,stateFunctional) {};

  virtual ~IP_SSATaucTaoTikhonovProblem() {};

  virtual PetscErrorCode connect(Tao tao);

  //! Callback to TAO to set bounds on \f$\tau_c\f$ for constrained minimization algorithms.
  virtual PetscErrorCode getVariableBounds(Tao tao, Vec lo, Vec hi);

};

} // end of namespace pism

#endif /* end of include guard: IP_SSATAUCTIKHONOVPROBLEM_HH_HB8UWICX */

