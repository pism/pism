// Copyright (C) 2013, 2014, 2015, 2016, 2023  David Maxwell and Constantine Khroulev
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


#include "pism/inverse/IP_SSATaucTaoTikhonovProblem.hh"

namespace pism {
namespace inverse {

void IP_SSATaucTaoTikhonovProblem::connect(Tao tao) {
  PetscErrorCode ierr;

  IPTaoTikhonovProblem<IP_SSATaucForwardProblem>::connect(tao);

  const char *type;
  ierr = TaoGetType(tao,&type);
  PISM_CHK(ierr, "TaoGetType");

  if (strcmp(type,"blmvm") == 0) {
    taoutil::TaoGetVariableBoundsCallback<IP_SSATaucTaoTikhonovProblem>::connect(tao,*this);
  }  
}


void IP_SSATaucTaoTikhonovProblem::getVariableBounds(Tao /*tao*/, Vec lo, Vec hi) {
  double zeta_min = 0.0, zeta_max = 0.0, tauc_min = 0.0, tauc_max = 0.0;

  tauc_min = m_grid->ctx()->config()->get_number("inverse.ssa.tauc_min");
  tauc_max = m_grid->ctx()->config()->get_number("inverse.ssa.tauc_max");

  IPDesignVariableParameterization &tauc_param = m_forward.tauc_param();
  tauc_param.fromDesignVariable(tauc_min,&zeta_min);
  tauc_param.fromDesignVariable(tauc_max,&zeta_max);

  PetscErrorCode ierr = VecSet(lo,zeta_min);
  PISM_CHK(ierr, "VecSet");

  ierr = VecSet(hi,zeta_max);
  PISM_CHK(ierr, "VecSet");
}


} // end of namespace inverse
} // end of namespace pism
