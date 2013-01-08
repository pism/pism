// Copyright (C) 2011, 2012, 2013 David Maxwell
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

#include "InvTaucParameterization.hh"
#include "pism_options.hh"
#include <cmath>

InvTaucParamIdent g_InvTaucParamIdent;
InvTaucParamSquare g_InvTaucParamSquare;
InvTaucParamExp g_InvTaucParamExp;

PetscErrorCode InvTaucParameterization::init( const NCConfigVariable & config ) { 
  m_tauc_scale = config.get("tauc_param_tauc_scale");
  return 0;
}

PetscErrorCode InvTaucParameterization::convertToTauc( IceModelVec2S &zeta, IceModelVec2S &tauc, bool communicate ) {
  PetscReal **zeta_a, **tauc_a;
  PetscErrorCode ierr;
  
  ierr = zeta.get_array(zeta_a); CHKERRQ(ierr);
  ierr = tauc.get_array(tauc_a); CHKERRQ(ierr);
  IceGrid &grid = *zeta.get_grid();
  for(PetscInt i= grid.xs; i< grid.xs+grid.xm; i++) {
    for(PetscInt j= grid.ys; j< grid.ys+grid.ym; j++) {
      ierr = this->toTauc(zeta_a[i][j], &tauc_a[i][j], NULL); CHKERRQ(ierr);
      if(std::isnan(tauc_a[i][j])) {
        PetscPrintf(PETSC_COMM_WORLD,"made a tauc nan zeta=%g tauc=%g\n",zeta_a[i][j],tauc_a[i][j]);
      } 
    }
  }
  ierr = zeta.end_access(); CHKERRQ(ierr);
  ierr = tauc.end_access(); CHKERRQ(ierr);
  if(communicate) {
    ierr = tauc.update_ghosts(); CHKERRQ(ierr);
  }
  return 0;
}

PetscErrorCode InvTaucParameterization::convertFromTauc( IceModelVec2S &tauc, IceModelVec2S &zeta, bool communicate ) {
  PetscReal **zeta_a, **tauc_a;
  PetscErrorCode ierr;
  ierr = zeta.get_array(zeta_a); CHKERRQ(ierr);
  ierr = tauc.get_array(tauc_a); CHKERRQ(ierr);
  IceGrid &grid = *zeta.get_grid();
  for(PetscInt i= grid.xs; i< grid.xs+grid.xm; i++) {
    for(PetscInt j= grid.ys; j< grid.ys+grid.ym; j++) {
      ierr = this->fromTauc(tauc_a[i][j], &zeta_a[i][j]); CHKERRQ(ierr);
      if(std::isnan(zeta_a[i][j])) {
        PetscPrintf(PETSC_COMM_WORLD,"made a zeta nan tauc=%g zeta=%g\n",tauc_a[i][j],zeta_a[i][j]);
      } 
    }    
  }
  ierr = zeta.end_access(); CHKERRQ(ierr);
  ierr = tauc.end_access(); CHKERRQ(ierr);
  if(communicate) {
    ierr = zeta.update_ghosts(); CHKERRQ(ierr);
  }
  return 0;
}

PetscErrorCode InvTaucParamIdent::toTauc( PetscReal p, PetscReal *value, 
                                          PetscReal *derivative)
{
  if(value != NULL) *value = m_tauc_scale*p;
  if(derivative != NULL) *derivative = m_tauc_scale;
  return 0;
}

PetscErrorCode InvTaucParamIdent::fromTauc( PetscReal tauc, PetscReal *OUTPUT)
{
  *OUTPUT = tauc/m_tauc_scale;
  return 0;
}


PetscErrorCode InvTaucParamSquare::toTauc( PetscReal p, PetscReal *value, 
                                           PetscReal *derivative)
{
  if(value != NULL) *value = m_tauc_scale*p*p;
  if(derivative != NULL) *derivative = m_tauc_scale*2*p;
  return 0;
}

PetscErrorCode InvTaucParamSquare::fromTauc( PetscReal tauc, PetscReal *OUTPUT)
{
  if(tauc<0) {
    tauc = 0;
  }
  *OUTPUT = sqrt(tauc/m_tauc_scale);
  return 0;
}

PetscErrorCode InvTaucParamExp::init( const NCConfigVariable &config ) { 
  PetscErrorCode ierr;
  ierr = InvTaucParameterization::init(config); CHKERRQ(ierr);
  m_tauc_eps = config.get("tauc_param_tauc_eps");
  return 0;
}

PetscErrorCode InvTaucParamExp::toTauc( PetscReal p, PetscReal *value, 
                                           PetscReal *derivative)
{
  if(value != NULL) {
    *value = m_tauc_scale*exp(p);
  } 
  
  if(derivative != NULL) *derivative = m_tauc_scale*exp(p);    
  return 0;
}

PetscErrorCode InvTaucParamExp::fromTauc( PetscReal tauc, PetscReal *OUTPUT)
{
  if(tauc < m_tauc_eps)
  {
    tauc= m_tauc_eps;
  }
  *OUTPUT=log(tauc/m_tauc_scale);

  return 0;
}

PetscErrorCode InvTaucParamTruncatedIdent::init( const NCConfigVariable &config ) {
  PetscErrorCode ierr;
  ierr = InvTaucParameterization::init(config); CHKERRQ(ierr);
  
  PetscReal tauc0 = config.get("tauc_param_trunc_tauc0"); 
  m_tauc0_sq = tauc0*tauc0/(m_tauc_scale*m_tauc_scale);

  m_tauc_eps = config.get("tauc_param_tauc_eps");

  return 0;
}

PetscErrorCode InvTaucParamTruncatedIdent::toTauc( PetscReal p, 
  PetscReal *value, PetscReal *derivative)
{
  PetscReal alpha = sqrt(p*p+4*m_tauc0_sq);
  if(value != NULL) *value = m_tauc_scale*(p+alpha)*0.5;
  if(derivative != NULL) *derivative = m_tauc_scale*(1+p/alpha)*0.5;
  return 0;
}

PetscErrorCode InvTaucParamTruncatedIdent::fromTauc( PetscReal tauc, PetscReal *OUTPUT)
{
  if(tauc < m_tauc_eps)
  {
    tauc= m_tauc_eps;
  }
  PetscReal tauc_dimensionless = tauc/m_tauc_scale;
  *OUTPUT=tauc_dimensionless-m_tauc0_sq/tauc_dimensionless;
  return 0;
}
