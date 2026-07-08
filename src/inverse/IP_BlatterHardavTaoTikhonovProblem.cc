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

#include "pism/inverse/IP_BlatterHardavTaoTikhonovProblem.hh"

namespace pism {
namespace inverse {

void IP_BlatterHardavTaoTikhonovProblem::connect(Tao tao) {
  PetscErrorCode ierr;

  IPTaoTikhonovProblem<IP_BlatterHardavForwardProblem>::connect(tao);

  const char *type;
  ierr = TaoGetType(tao, &type);
  PISM_CHK(ierr, "TaoGetType");

  if (strcmp(type, "blmvm") == 0) {
    taoutil::TaoGetVariableBoundsCallback<IP_BlatterHardavTaoTikhonovProblem>::connect(
        tao, *this);
  }
}

void IP_BlatterHardavTaoTikhonovProblem::getVariableBounds(Tao /*tao*/, Vec lo, Vec hi) {
  double zeta_min = 0.0, zeta_max = 0.0;
  double hardav_min = 0.0, hardav_max = 0.0;

  hardav_min = m_grid->ctx()->config()->get_number("inverse.stress_balance.hardav_min");
  hardav_max = m_grid->ctx()->config()->get_number("inverse.stress_balance.hardav_max");

  IPDesignVariableParameterization &design_param = m_forward.design_param();
  design_param.fromDesignVariable(hardav_min, &zeta_min);
  design_param.fromDesignVariable(hardav_max, &zeta_max);

  PetscErrorCode ierr = VecSet(lo, zeta_min);
  PISM_CHK(ierr, "VecSet");

  ierr = VecSet(hi, zeta_max);
  PISM_CHK(ierr, "VecSet");
}

} // end of namespace inverse
} // end of namespace pism
