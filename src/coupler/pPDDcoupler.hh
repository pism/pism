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
#include <gsl/gsl_rng.h>
#include "../base/grid.hh"
#include "../base/iceModelVec.hh"
#include "pccoupler.hh"

//! A derived class of PISMAtmosphereCoupler which provides a PDD to PISM.
/*!
The PDD here is the one already implemented in PISM.  That is, it is the one
from EISMINT-Greenland.

If <tt>-pdd_monthly_temps</tt> is not used, it reads 'artm' from input file
and interprets it as the (location dependent) mean annual surface temperature 
vannmeansurftemp (K).  Then the surface temperature at a time, vsurftemp, 
comes from a standard yearly cycle without randomness.  Note 
vannmeansurftemp is written to output files as 'annavartm'.

If <tt>-pdd_monthly_temps</tt> is used then it reads 12 monthly temperature
data sets, 'temp_mon0', ..., 'temp_mon11', from input file and linearly
interpolates these to get vsurftemp.

It reads 'acab' from input file and interprets it as ice-equivalent snow
accumulation rate, vsurfaccum.  The PDD is used to convert to ice-equivalent 
net surface mass flux, vsurfmassflux.  Note vsurfaccum is writen to output
files as 'accum'.

It has various constants parameterizing the melt and refreeze processes.  See REFERENCE
 */
class PISMPDDCoupler : public PISMConstAtmosCoupler {

public:
  PISMPDDCoupler();
  virtual ~PISMPDDCoupler();  // destroys PDD

  PetscErrorCode userOptionsChoosePDD(PetscTruth &userWantsPDD);
  
  virtual PetscErrorCode initFromOptions(IceGrid* g);
  virtual PetscErrorCode writeCouplingFieldsToFile(const char *filename);

  virtual PetscErrorCode updateSurfMassFluxAndProvide(
             const PetscScalar t_years, const PetscScalar dt_years, 
             void *iceInfoNeeded, // will be interpreted as type IceInfoNeededByAtmosphereCoupler*
             IceModelVec2* &pvsmf);  // vsmf = pointer to vsurfmassflux

  virtual PetscErrorCode updateSurfTempAndProvide(
             const PetscScalar t_years, const PetscScalar dt_years,
             void *iceInfoNeeded, // will be interpreted as type iceInfoNeededByAtmosphereCoupler*
             IceModelVec2* &pvst);  // vst = pointer to vsurftemp


protected:
  IceModelVec2 vannmeansurftemp,
               vsurfaccum;

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
                    const PetscScalar dt_secs, const PetscScalar pddsum);
  double CalovGreveIntegrand(const double Tac);  // double because handed to gsl quadrature routine

};


#endif

