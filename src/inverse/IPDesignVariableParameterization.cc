// Copyright (C) 2011, 2012, 2013, 2014 David Maxwell
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
#include "PISMConfig.hh"

namespace pism {

//! Initializes the scale parameters of the parameterization.
/*! Every IPDesignVariableParameterization has an associated scale for the design variable
\f$d_{\rm scale}\f$ that equals 1 in internal units.  The scale for a design variable named \a foo
is stored in an PISMConfig file as design_param_foo_scale.  Subclasses may have additional
parameters that are follow the naming convention \a design_param_foo_*.

\param config          The config file to read the scale parameters from.
\param design_var_name The associated name of the design variable, e.g. 'tauc' or 'hardav'
*/
PetscErrorCode IPDesignVariableParameterization::set_scales(const PISMConfig & config, const char *design_var_name ) {
  std::string key("design_param_");
  key += design_var_name;
  key += "_scale";
  m_d_scale = config.get(key);
  return 0;
}

//! Transforms a vector of \f$\zeta\f$ values to a vector of \f$d\f$ values.
PetscErrorCode IPDesignVariableParameterization::convertToDesignVariable(IceModelVec2S &zeta,
                                                                         IceModelVec2S &d,
                                                                         bool communicate) {
  PetscErrorCode ierr;

  ierr = zeta.begin_access(); CHKERRQ(ierr);
  ierr = d.begin_access(); CHKERRQ(ierr);
  IceGrid &grid = *zeta.get_grid();
  for (int i = grid.xs; i < grid.xs + grid.xm; i++) {
    for (int j = grid.ys; j < grid.ys + grid.ym; j++) {
      ierr = this->toDesignVariable(zeta(i, j), &d(i, j), NULL); CHKERRQ(ierr);
      if (std::isnan(d(i, j))) {
        PetscPrintf(PETSC_COMM_WORLD, "made a d nan zeta = %g d = %g\n", zeta(i, j), d(i, j));
      }
    }
  }
  ierr = zeta.end_access(); CHKERRQ(ierr);
  ierr = d.end_access(); CHKERRQ(ierr);
  if (communicate) {
    ierr = d.update_ghosts(); CHKERRQ(ierr);
  }
  return 0;
}

  //! Transforms a vector of \f$d\f$ values to a vector of \f$\zeta\f$ values.
PetscErrorCode IPDesignVariableParameterization::convertFromDesignVariable(IceModelVec2S &d,
                                                                            IceModelVec2S &zeta,
                                                                            bool communicate) {
  PetscErrorCode ierr;
  ierr = zeta.begin_access(); CHKERRQ(ierr);
  ierr = d.begin_access(); CHKERRQ(ierr);
  IceGrid &grid = *zeta.get_grid();
  for (int i = grid.xs; i < grid.xs + grid.xm; i++) {
    for (int j = grid.ys; j < grid.ys + grid.ym; j++) {
      ierr = this->fromDesignVariable(d(i, j), &zeta(i, j)); CHKERRQ(ierr);
      if (std::isnan(zeta(i, j))) {
        PetscPrintf(PETSC_COMM_WORLD, "made a zeta nan d = %g zeta = %g\n", d(i, j), zeta(i, j));
      }
    }
  }
  ierr = zeta.end_access(); CHKERRQ(ierr);
  ierr = d.end_access(); CHKERRQ(ierr);
  if (communicate) {
    ierr = zeta.update_ghosts(); CHKERRQ(ierr);
  }
  return 0;
}

PetscErrorCode IPDesignVariableParamIdent::toDesignVariable(double p, double *value,
                                          double *derivative)
{
  if (value != NULL) *value = m_d_scale*p;
  if (derivative != NULL) *derivative = m_d_scale;
  return 0;
}

PetscErrorCode IPDesignVariableParamIdent::fromDesignVariable(double d, double *OUTPUT)
{
  *OUTPUT = d / m_d_scale;
  return 0;
}


PetscErrorCode IPDesignVariableParamSquare::toDesignVariable(double p, double *value,
                                           double *derivative)
{
  if (value != NULL) *value = m_d_scale*p*p;
  if (derivative != NULL) *derivative = m_d_scale*2*p;
  return 0;
}

PetscErrorCode IPDesignVariableParamSquare::fromDesignVariable(double d, double *OUTPUT)
{
  if (d < 0) {
    d = 0;
  }
  *OUTPUT = sqrt(d / m_d_scale);
  return 0;
}

PetscErrorCode IPDesignVariableParamExp::set_scales(const PISMConfig &config, const char *design_var_name ) {
  PetscErrorCode ierr;
  ierr = IPDesignVariableParameterization::set_scales(config, design_var_name); CHKERRQ(ierr);

  std::string key("design_param_");
  key += design_var_name;
  key += "_eps";
  m_d_eps = config.get(key);
  return 0;
}

PetscErrorCode IPDesignVariableParamExp::toDesignVariable(double p, double *value,
                                           double *derivative)
{
  if (value != NULL) {
    *value = m_d_scale*exp(p);
  }

  if (derivative != NULL) *derivative = m_d_scale*exp(p);
  return 0;
}

PetscErrorCode IPDesignVariableParamExp::fromDesignVariable(double d, double *OUTPUT)
{
  if (d < m_d_eps)
  {
    d = m_d_eps;
  }
  *OUTPUT = log(d / m_d_scale);

  return 0;
}


PetscErrorCode IPDesignVariableParamTruncatedIdent::set_scales(const PISMConfig &config,
                                                               const char *design_var_name ) {
  PetscErrorCode ierr;
  ierr = IPDesignVariableParameterization::set_scales(config, design_var_name); CHKERRQ(ierr);

  std::string key("design_param_trunc_");
  key += design_var_name;
  key += "0";

  double d0 = config.get(key);
  m_d0_sq = d0*d0 / (m_d_scale*m_d_scale);


  key = "design_param_";
  key += design_var_name;
  key += "_eps";
  m_d_eps = config.get(key);

  return 0;
}

PetscErrorCode IPDesignVariableParamTruncatedIdent::toDesignVariable(double p,
                                                                     double *value,
                                                                     double *derivative)
{
  double alpha = sqrt(p*p + 4*m_d0_sq);
  if (value != NULL) *value = m_d_scale*(p + alpha)*0.5;
  if (derivative != NULL) *derivative = m_d_scale*(1 + p / alpha)*0.5;
  return 0;
}

PetscErrorCode IPDesignVariableParamTruncatedIdent::fromDesignVariable(double d, double *OUTPUT)
{
  if (d < m_d_eps)
  {
    d = m_d_eps;
  }
  double d_dimensionless = d / m_d_scale;
  *OUTPUT = d_dimensionless - m_d0_sq / d_dimensionless;
  return 0;
}

} // end of namespace pism
