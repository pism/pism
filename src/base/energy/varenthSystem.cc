// Copyright (C) 2011, 2012 The PISM Authors
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

#include <gsl/gsl_math.h>

#include "enthalpyConverter.hh"
#include "varcEnthalpyConverter.hh"
#include "varenthSystem.hh"
#include "NCVariable.hh"

varenthSystemCtx::varenthSystemCtx(const NCConfigVariable &config,
                                     IceModelVec3 &my_Enth3, int my_Mz,
                                     string my_prefix, EnthalpyConverter *my_EC)
  : enthSystemCtx(config, my_Enth3, my_Mz, my_prefix), EC(my_EC),
    ice_thickness(0) {

  if (config.get_flag("use_temperature_dependent_thermal_conductivity"))
    k_depends_on_T = true;
  else
    k_depends_on_T = false;
        
}


/*!  Equation (4.37) in \ref GreveBlatter2009 is
  \f[ k(T) = 9.828 e^{−0.0057 T} \f]
  where \f$T\f$ is in Kelvin and the resulting conductivity is in units W m−1 K−1.
 */
PetscScalar varenthSystemCtx::k_from_T(PetscScalar T) {
  return 9.828 * exp(-0.0057 * T);
}


PetscErrorCode varenthSystemCtx::initThisColumn(bool my_ismarginal,
                                                 PetscScalar my_lambda,
                                                 PetscReal my_ice_thickness) {
  PetscErrorCode ierr;
  ierr = enthSystemCtx::initThisColumn(my_ismarginal, my_lambda,
                                       my_ice_thickness); CHKERRQ(ierr);

  ice_thickness = my_ice_thickness;
  return 0;
}
  

//! \brief Assemble the R array for the current column, in the case with temperature-dependent properties.
/*!
In this implementation the conductivity \f$k=k(T)\f$ and the specific heat capacity
\f$c=c(T)\f$ are allowed to depend on temperature.  Thus also \f$R=R(T)\f$
depends on the temperature.  We get the temperature from the enthalpy.
 */
PetscErrorCode varenthSystemCtx::assemble_R() {
  PetscErrorCode ierr;

  const PetscScalar Rfactor = dtTemp / (PetscSqr(dzEQ) * ice_rho);

  for (PetscInt k = 0; k <= ks; k++) {
    if (Enth[k] < Enth_s[k]) {
      // cold case
      const PetscScalar depth = ice_thickness - k * dzEQ;
      PetscScalar T;
      ierr = EC->getAbsTemp(Enth[k], EC->getPressureFromDepth(depth), // FIXME: issue #15
                            T); CHKERRQ(ierr);

      R[k] = ((k_depends_on_T ? k_from_T(T) : ice_k) / EC->c_from_T(T)) * Rfactor;
    } else {
      // temperate case
      R[k] = iceRtemp;
    }
  }

  // R[k] for k > ks are never used
#if (PISM_DEBUG==1)
  for (int k = ks + 1; k < Mz; ++k)
    R[k] = GSL_NAN;
#endif

  return 0;
}

