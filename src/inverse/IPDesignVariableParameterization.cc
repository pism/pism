// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017 David Maxwell
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

#include <cmath>

#include "pism/util/iceModelVec.hh"
#include "IPDesignVariableParameterization.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/error_handling.hh"

namespace pism {
namespace inverse {

//! Initializes the scale parameters of the parameterization.
/*! Every IPDesignVariableParameterization has an associated scale for the design variable
\f$d_{\rm scale}\f$ that equals 1 in internal units.  The scale for a design variable named \a foo
is stored in an Config file as design_param_foo_scale.  Subclasses may have additional
parameters that are follow the naming convention \a design_param_foo_*.

\param config          The config file to read the scale parameters from.
\param design_var_name The associated name of the design variable, e.g. 'tauc' or 'hardav'
*/
void IPDesignVariableParameterization::set_scales(const Config & config,
                                                  const std::string &design_var_name) {
  std::string key("inverse.design.param_");
  key += design_var_name;
  key += "_scale";
  m_d_scale = config.get_double(key);
}

//! Transforms a vector of \f$\zeta\f$ values to a vector of \f$d\f$ values.
void IPDesignVariableParameterization::convertToDesignVariable(IceModelVec2S &zeta,
                                                               IceModelVec2S &d,
                                                               bool communicate) {
  PetscErrorCode ierr;

  IceModelVec::AccessList list{&zeta, &d};

  const IceGrid &grid = *zeta.grid();

  ParallelSection loop(grid.com);
  try {
    for (Points p(grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      this->toDesignVariable(zeta(i, j), &d(i, j), NULL);
      if (std::isnan(d(i, j))) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,
                           "made a d nan zeta = %g d = %g\n",
                           zeta(i, j), d(i, j));
        PISM_CHK(ierr, "PetscPrintf");
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  if (communicate) {
    d.update_ghosts();
  }
}

  //! Transforms a vector of \f$d\f$ values to a vector of \f$\zeta\f$ values.
void IPDesignVariableParameterization::convertFromDesignVariable(IceModelVec2S &d,
                                                                 IceModelVec2S &zeta,
                                                                 bool communicate) {
  PetscErrorCode ierr;
  IceModelVec::AccessList list{&zeta, &d};

  const IceGrid &grid = *zeta.grid();

  ParallelSection loop(grid.com);
  try {
    for (Points p(grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      this->fromDesignVariable(d(i, j), &zeta(i, j));
      if (std::isnan(zeta(i, j))) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,
                           "made a zeta nan d = %g zeta = %g\n",
                           d(i, j), zeta(i, j));
        PISM_CHK(ierr, "PetscPrintf");
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  if (communicate) {
    zeta.update_ghosts();
  }
}

void IPDesignVariableParamIdent::toDesignVariable(double p, double *value,
                                                  double *derivative) {
  if (value != NULL) {
    *value = m_d_scale*p;
  }
  if (derivative != NULL) {
    *derivative = m_d_scale;
  }
}

void IPDesignVariableParamIdent::fromDesignVariable(double d, double *OUTPUT) {
  *OUTPUT = d / m_d_scale;
}


void IPDesignVariableParamSquare::toDesignVariable(double p, double *value,
                                                   double *derivative) {
  if (value != NULL) {
    *value = m_d_scale*p*p;
  }
  if (derivative != NULL) {
    *derivative = m_d_scale*2*p;
  }
}

void IPDesignVariableParamSquare::fromDesignVariable(double d, double *OUTPUT) {
  if (d < 0) {
    d = 0;
  }
  *OUTPUT = sqrt(d / m_d_scale);
}

void IPDesignVariableParamExp::set_scales(const Config &config, const std::string &design_var_name) {
  IPDesignVariableParameterization::set_scales(config, design_var_name);

  std::string key("inverse.design.param_");
  key += design_var_name;
  key += "_eps";
  m_d_eps = config.get_double(key);
}

void IPDesignVariableParamExp::toDesignVariable(double p, double *value,
                                           double *derivative) {
  if (value != NULL) {
    *value = m_d_scale*exp(p);
  }

  if (derivative != NULL) {
    *derivative = m_d_scale*exp(p);
  }
}

void IPDesignVariableParamExp::fromDesignVariable(double d, double *OUTPUT) {
  if (d < m_d_eps) {
    d = m_d_eps;
  }
  *OUTPUT = log(d / m_d_scale);
}


void IPDesignVariableParamTruncatedIdent::set_scales(const Config &config,
                                                     const std::string &design_var_name) {
  IPDesignVariableParameterization::set_scales(config, design_var_name);

  std::string key("inverse.design.param_trunc_");
  key += design_var_name;
  key += "0";

  double d0 = config.get_double(key);
  m_d0_sq = d0*d0 / (m_d_scale*m_d_scale);


  key = "inverse.design.param_";
  key += design_var_name;
  key += "_eps";
  m_d_eps = config.get_double(key);
}

void IPDesignVariableParamTruncatedIdent::toDesignVariable(double p,
                                                           double *value,
                                                           double *derivative) {
  double alpha = sqrt(p*p + 4*m_d0_sq);
  if (value != NULL) {
    *value = m_d_scale*(p + alpha)*0.5;
  }
  if (derivative != NULL) {
    *derivative = m_d_scale*(1 + p / alpha)*0.5;
  }
}

void IPDesignVariableParamTruncatedIdent::fromDesignVariable(double d, double *OUTPUT) {
  if (d < m_d_eps) {
    d = m_d_eps;
  }

  double d_dimensionless = d / m_d_scale;
  *OUTPUT = d_dimensionless - m_d0_sq / d_dimensionless;
}

} // end of namespace inverse
} // end of namespace pism
