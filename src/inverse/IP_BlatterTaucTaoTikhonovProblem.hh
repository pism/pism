// Copyright (C) 2026 Andy Aschwanden and Constantine Khroulev
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

#ifndef IP_BLATTERTAUCTAOTIKHONOVPROBLEM_HH
#define IP_BLATTERTAUCTAOTIKHONOVPROBLEM_HH

#include "pism/inverse/IPTaoTikhonovProblem.hh"
#include "pism/inverse/IP_BlatterTaucForwardProblem.hh"
#include "pism/inverse/TaoUtil.hh"
#include "pism/inverse/functional/IPFunctional.hh"

namespace pism {
namespace inverse {

//! Tikhonov problem for inversion of basal yield stress from Blatter velocities.
/*!
  Analogous to IP_SSATaucTaoTikhonovProblem but using
  IP_BlatterTaucForwardProblem as the forward problem.

  Compatible with TAO minimization algorithms (tao_cg, tao_lmvm, tao_blmvm).
  If tao_blmvm is selected, tauc will be constrained by the config variables
  inverse.stress_balance.tauc_min and inverse.stress_balance.tauc_max.
*/
class IP_BlatterTaucTaoTikhonovProblem
    : public IPTaoTikhonovProblem<IP_BlatterTaucForwardProblem> {
public:
  IP_BlatterTaucTaoTikhonovProblem(
      IP_BlatterTaucForwardProblem &forward,
      DesignVec &d0,
      StateVec &u_obs,
      double eta,
      IPFunctional<DesignVec> &designFunctional,
      IPFunctional<StateVec> &stateFunctional)
      : IPTaoTikhonovProblem<IP_BlatterTaucForwardProblem>(
            forward, d0, u_obs, eta, designFunctional, stateFunctional) {}

  virtual ~IP_BlatterTaucTaoTikhonovProblem() {}

  virtual void connect(Tao tao);

  //! Set bounds on tauc for constrained minimization algorithms.
  virtual void getVariableBounds(Tao tao, Vec lo, Vec hi);
};

} // end of namespace inverse
} // end of namespace pism

#endif /* IP_BLATTERTAUCTAOTIKHONOVPROBLEM_HH */
