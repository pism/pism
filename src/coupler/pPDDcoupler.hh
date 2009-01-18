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

#ifndef __pPDDcoupler_hh
#define __pPDDcoupler_hh

#include <petsc.h>
#include "../base/grid.hh"
#include "pccoupler.hh"

//! A derived class of PISMConstAtmosCoupler which provides a PDD to PISM.
/*!
The PDD here is the one already implemented in PISM.  That is, it is the one
from EISMINT-Greenland.  Thus it has various constants parameterizing the 
melt and refreeze processes.
 */
class PISMPDDCoupler : public PISMAtmosphereCoupler {

public:
  PISMPDDCoupler();
  virtual ~PISMPDDCoupler();  // destroys PDD

  PetscErrorCode userOptionsChoosePDD(PetscTruth &userWantsPDD);
  virtual PetscErrorCode initFromOptions(IceGrid* g);
  virtual PetscErrorCode writeCouplingFieldsToFile(const char *filename);


  virtual PetscErrorCode updateSurfTempAndProvide(
             const PetscScalar t_years, const PetscScalar dt_years,
             IceModelVec2 mask, IceModelVec2 surface_elev,
             IceModelVec2* &vst);


protected:
  IceModelVec2 vsurfaccum;
  gsl_rng      *pddRandGen;      // usually NULL; default is expectation integral which
                                 //   does not use actual random numbers
  PetscScalar  pddStdDev,        // K; daily amount of randomness
               pddFactorSnow,    // m day^-1 K^-1; amount of snow melted,
                                 //    as ice equivalent, per positive degree day
               pddFactorIce,     // m day^-1 K^-1; amount of ice melted
                                 //    per positive degree day
               pddRefreezeFrac,  // [pure fraction]; amount of melted snow which refreezes
                                 //    as ice
               pddSummerWarming, // K; amplitude of yearly temperature cycle
               pddSummerPeakDay; // Julian day of summer temperature peak

  PetscScalar getSummerWarming(
             const PetscScalar elevation, const PetscScalar latitude,
             const PetscScalar Tma);
  PetscScalar getTemperatureFromYearlyCycle(
       const PetscScalar summer_warming, const PetscScalar Tma, const PetscScalar day);

  PetscTruth   usingMonthlyTemps;
  IceModelVec2 vmonthlysurftemp[12]; // usually these are not created; if user supplies monthly
                                 //   temperature maps then we will allocate 12 IceModelVec2's 
  PetscErrorCode readMonthlyTemps(const char *filename);
  PetscErrorCode getMonthIndicesFromDay(const PetscScalar day, PetscInt &curr, PetscInt &next);
  PetscScalar getTemperatureFromMonthlyData(
       PetscScalar **currMonthSurfTemps, PetscScalar **nextMonthSurfTemps,
       const PetscInt i, const PetscInt j, const PetscScalar day);

  PetscScalar getSurfaceBalanceFromSnowAndPDD(const PetscScalar snowrate,
                    const PetscScalar dt, const PetscScalar pddsum);
  double CalovGreveIntegrand(const double Tac);  // double because handed to gsl quadrature routine
};


#endif

