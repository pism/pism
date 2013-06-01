// Copyright (C) 2011, 2012, 2013 David Maxwell
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

#include "IPDesignVariableParameterization.hh"
#include "pism_options.hh"
#include <cmath>

PetscErrorCode IPDesignVariableParameterization::set_scales( const NCConfigVariable & config, const char *design_var_name ) { 
  std::string key("design_param_");
  key += design_var_name;
  key += "_scale";
  m_d_scale = config.get(key);
  return 0;
}

PetscErrorCode IPDesignVariableParameterization::convertToDesignVariable( IceModelVec2S &zeta, IceModelVec2S &d, bool communicate ) {
  PetscReal **zeta_a, **d_a;
  PetscErrorCode ierr;
  
  ierr = zeta.get_array(zeta_a); CHKERRQ(ierr);
  ierr = d.get_array(d_a); CHKERRQ(ierr);
  IceGrid &grid = *zeta.get_grid();
  for(PetscInt i= grid.xs; i< grid.xs+grid.xm; i++) {
    for(PetscInt j= grid.ys; j< grid.ys+grid.ym; j++) {
      ierr = this->toDesignVariable(zeta_a[i][j], &d_a[i][j], NULL); CHKERRQ(ierr);
      if(std::isnan(d_a[i][j])) {
        PetscPrintf(PETSC_COMM_WORLD,"made a d nan zeta=%g d=%g\n",zeta_a[i][j],d_a[i][j]);
      } 
    }
  }
  ierr = zeta.end_access(); CHKERRQ(ierr);
  ierr = d.end_access(); CHKERRQ(ierr);
  if(communicate) {
    ierr = d.update_ghosts(); CHKERRQ(ierr);
  }
  return 0;
}

PetscErrorCode IPDesignVariableParameterization::convertFromDesignVariable( IceModelVec2S &d, IceModelVec2S &zeta, bool communicate ) {
  PetscReal **zeta_a, **d_a;
  PetscErrorCode ierr;
  ierr = zeta.get_array(zeta_a); CHKERRQ(ierr);
  ierr = d.get_array(d_a); CHKERRQ(ierr);
  IceGrid &grid = *zeta.get_grid();
  for(PetscInt i= grid.xs; i< grid.xs+grid.xm; i++) {
    for(PetscInt j= grid.ys; j< grid.ys+grid.ym; j++) {
      ierr = this->fromDesignVariable(d_a[i][j], &zeta_a[i][j]); CHKERRQ(ierr);
      if(std::isnan(zeta_a[i][j])) {
        PetscPrintf(PETSC_COMM_WORLD,"made a zeta nan d=%g zeta=%g\n",d_a[i][j],zeta_a[i][j]);
      } 
    }    
  }
  ierr = zeta.end_access(); CHKERRQ(ierr);
  ierr = d.end_access(); CHKERRQ(ierr);
  if(communicate) {
    ierr = zeta.update_ghosts(); CHKERRQ(ierr);
  }
  return 0;
}

PetscErrorCode IPDesignVariableParamIdent::toDesignVariable( PetscReal p, PetscReal *value, 
                                          PetscReal *derivative)
{
  if(value != NULL) *value = m_d_scale*p;
  if(derivative != NULL) *derivative = m_d_scale;
  return 0;
}

PetscErrorCode IPDesignVariableParamIdent::fromDesignVariable( PetscReal d, PetscReal *OUTPUT)
{
  *OUTPUT = d/m_d_scale;
  return 0;
}


PetscErrorCode IPDesignVariableParamSquare::toDesignVariable( PetscReal p, PetscReal *value, 
                                           PetscReal *derivative)
{
  if(value != NULL) *value = m_d_scale*p*p;
  if(derivative != NULL) *derivative = m_d_scale*2*p;
  return 0;
}

PetscErrorCode IPDesignVariableParamSquare::fromDesignVariable( PetscReal d, PetscReal *OUTPUT)
{
  if(d<0) {
    d = 0;
  }
  *OUTPUT = sqrt(d/m_d_scale);
  return 0;
}

PetscErrorCode IPDesignVariableParamExp::set_scales( const NCConfigVariable &config, const char *design_var_name ) { 
  PetscErrorCode ierr;
  ierr = IPDesignVariableParameterization::set_scales(config, design_var_name); CHKERRQ(ierr);

  std::string key("design_param_");
  key += design_var_name;
  key += "_eps";
  m_d_eps = config.get(key);
  return 0;
}

PetscErrorCode IPDesignVariableParamExp::toDesignVariable( PetscReal p, PetscReal *value, 
                                           PetscReal *derivative)
{
  if(value != NULL) {
    *value = m_d_scale*exp(p);
  } 
  
  if(derivative != NULL) *derivative = m_d_scale*exp(p);    
  return 0;
}

PetscErrorCode IPDesignVariableParamExp::fromDesignVariable( PetscReal d, PetscReal *OUTPUT)
{
  if(d < m_d_eps)
  {
    d= m_d_eps;
  }
  *OUTPUT=log(d/m_d_scale);

  return 0;
}


PetscErrorCode IPDesignVariableParamTruncatedIdent::set_scales( const NCConfigVariable &config, const char *design_var_name ) {
  PetscErrorCode ierr;
  ierr = IPDesignVariableParameterization::set_scales(config,design_var_name); CHKERRQ(ierr);

  std::string key("design_param_trunc_");
  key += design_var_name;
  key += "0";
  
  PetscReal d0 = config.get(key); 
  m_d0_sq = d0*d0/(m_d_scale*m_d_scale);


  key = "design_param_";
  key += design_var_name;
  key += "_eps";
  m_d_eps = config.get(key);

  return 0;
}

PetscErrorCode IPDesignVariableParamTruncatedIdent::toDesignVariable( PetscReal p, 
  PetscReal *value, PetscReal *derivative)
{
  PetscReal alpha = sqrt(p*p+4*m_d0_sq);
  if(value != NULL) *value = m_d_scale*(p+alpha)*0.5;
  if(derivative != NULL) *derivative = m_d_scale*(1+p/alpha)*0.5;
  return 0;
}

PetscErrorCode IPDesignVariableParamTruncatedIdent::fromDesignVariable( PetscReal d, PetscReal *OUTPUT)
{
  if(d < m_d_eps)
  {
    d= m_d_eps;
  }
  PetscReal d_dimensionless = d/m_d_scale;
  *OUTPUT=d_dimensionless-m_d0_sq/d_dimensionless;
  return 0;
}
