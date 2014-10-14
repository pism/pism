// Copyright (C) 2013, 2014  David Maxwell and Constantine Khroulev
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


#include "IP_SSATaucTaoTikhonovProblem.hh"

namespace pism {

PetscErrorCode IP_SSATaucTaoTikhonovProblem::connect(Tao tao) {
  PetscErrorCode ierr;

  ierr = IPTaoTikhonovProblem<IP_SSATaucForwardProblem>::connect(tao); CHKERRQ(ierr);

  const char *type;
  ierr = TaoGetType(tao,&type); CHKERRQ(ierr);
  if (strcmp(type,"blmvm") == 0) {
    ierr = TaoGetVariableBoundsCallback<IP_SSATaucTaoTikhonovProblem>::connect(tao,*this); CHKERRQ(ierr);    
  }  
  return 0;
}


PetscErrorCode IP_SSATaucTaoTikhonovProblem::getVariableBounds(Tao /*tao*/, Vec lo, Vec hi) {
  PetscErrorCode ierr;
  double zeta_min, zeta_max, tauc_min, tauc_max;

  tauc_min = m_grid->config.get("inv_ssa_tauc_min");
  tauc_max = m_grid->config.get("inv_ssa_tauc_max");

  IPDesignVariableParameterization &tauc_param = m_forward.tauc_param();
  ierr = tauc_param.fromDesignVariable(tauc_min,&zeta_min); CHKERRQ(ierr);
  ierr = tauc_param.fromDesignVariable(tauc_max,&zeta_max); CHKERRQ(ierr);

  ierr = VecSet(lo,zeta_min); CHKERRQ(ierr);
  ierr = VecSet(hi,zeta_max); CHKERRQ(ierr);
  return 0;
}


} // end of namespace pism
