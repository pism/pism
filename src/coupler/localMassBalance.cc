// Copyright (C) 2009 Ed Bueler
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

#include <petsc.h>
#include <ctime>  // for time(), used to initialize random number gen
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>       // for erfc() in CalovGreveIntegrand()
#include "../base/pism_const.hh"
#include "localMassBalance.hh"



PetscScalar LocalMassBalance::getMassFluxFromTemperatureTimeSeries(
             PetscScalar t, PetscScalar dt, PetscScalar *T, PetscInt N,
             PetscScalar precip) {
  PetscPrintf(PETSC_COMM_WORLD,"LocalMassBalance is a virtual class.  ENDING ...\n");
  PetscEnd();
  return 0.0;
}


PDDMassBalance::PDDMassBalance() {
  pddStdDev       = 5.0;  // K; std dev of daily temp variation; EISMINT-Greenland value
  pddFactorSnow = 0.003;  // m K^-1 day^-1; EISMINT-Greenland value; = 3 mm / (pos degree day)
  pddFactorIce  = 0.008;  // m K^-1 day^-1; EISMINT-Greenland value; = 8 mm / (pos degree day)
  pddRefreezeFrac = 0.6;  // [pure fraction]; EISMINT-Greenland value
}


PetscScalar PDDMassBalance::getMassFluxFromTemperatureTimeSeries(
             PetscScalar t, PetscScalar dt, PetscScalar *T, PetscInt N,
             PetscScalar precip) {
  return 0.0;
}


PetscScalar PDDMassBalance::CalovGreveIntegrand(
             PetscScalar sigma, PetscScalar Tac) {
  const PetscScalar TacC = Tac - 273.15;
  return (sigma / sqrt(2.0 * pi)) * exp(- TacC * TacC / (2.0 * sigma * sigma)) 
           + (TacC / 2.0) * gsl_sf_erfc(- TacC / (sqrt(2.0) * sigma));
}


PetscScalar PDDMassBalance::getPDDSumFromTemperatureTimeSeries(
             PetscScalar t, PetscScalar dt, PetscScalar *T, PetscInt N) {
  // N needs to be ODD for this to work
  if ((N % 2) == 0) {
    PetscPrintf(PETSC_COMM_WORLD,
        "PDDMassBalance::getPDDSum..() needs ODD N for Simpson's rule\n");
    PetscEnd();
  }
  PetscScalar  pdd_sum = 0.0;  // return value has units:  K day
  // FIXME: this is just the default mechanism
  const PetscScalar sperd = 8.64e4, // exact; from UDUNITS
                    h_days = dt / sperd;
  // Simpson's is: (h/3) * sum([1 4 2 4 2 4 ... 4 1] .* [f(x0) f(x1) ... f(xN)])
  for (PetscInt m = 0; m < N; ++m) {
    PetscScalar  coeff = ((m % 2) == 1) ? 4.0 : 2.0;
    if ( (m == 0) || (m == (N-1)) )  coeff = 1.0;
    pdd_sum += coeff * CalovGreveIntegrand(T[m],pddStdDev);  // pass in temp in K
  }
  pdd_sum = (h_days / 3.0) * pdd_sum;
  return pdd_sum;
}

