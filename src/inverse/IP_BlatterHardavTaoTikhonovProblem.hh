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

#ifndef IP_BLATTERHARDAVTAOTIKHONOVPROBLEM_HH
#define IP_BLATTERHARDAVTAOTIKHONOVPROBLEM_HH

#include "pism/inverse/IPTaoTikhonovProblem.hh"
#include "pism/inverse/IP_BlatterHardavForwardProblem.hh"
#include "pism/inverse/TaoUtil.hh"
#include "pism/inverse/functional/IPFunctional.hh"

namespace pism {
namespace inverse {

//! Tikhonov problem for inversion of ice hardness from Blatter velocities.
/*!
  Analogous to IP_SSAHardavTaoTikhonovProblem but using
  IP_BlatterHardavForwardProblem as the forward problem.

  Compatible with TAO minimization algorithms (tao_cg, tao_lmvm, tao_blmvm).
  If tao_blmvm is selected, hardness will be constrained by the config
  variables inverse.stress_balance.hardav_min and
  inverse.stress_balance.hardav_max.
*/
class IP_BlatterHardavTaoTikhonovProblem
    : public IPTaoTikhonovProblem<IP_BlatterHardavForwardProblem> {
public:
  IP_BlatterHardavTaoTikhonovProblem(
      IP_BlatterHardavForwardProblem &forward,
      DesignVec &d0,
      StateVec &u_obs,
      double eta,
      IPFunctional<DesignVec> &designFunctional,
      IPFunctional<StateVec> &stateFunctional)
      : IPTaoTikhonovProblem<IP_BlatterHardavForwardProblem>(
            forward, d0, u_obs, eta, designFunctional, stateFunctional) {}

  virtual ~IP_BlatterHardavTaoTikhonovProblem() {}

  virtual void connect(Tao tao);

  //! Set bounds on hardness for constrained minimization algorithms.
  virtual void getVariableBounds(Tao tao, Vec lo, Vec hi);
};

} // end of namespace inverse
} // end of namespace pism

#endif /* IP_BLATTERHARDAVTAOTIKHONOVPROBLEM_HH */
