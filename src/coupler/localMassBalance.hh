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

#ifndef __localMassBalance_hh
#define __localMassBalance_hh

#include <petsc.h>
#include "../base/NCVariable.hh"

//! Base class for a model which computes surface mass flux rate (ice thickness per time) from a precipitation (scalar) and a time series for temperature.
/*!
This is a process model.  It uses a 1D array, with a time dimension, for snow 
temperatures.  This process model does not know its location on the ice sheet, but
simply computes the surface mass balance from three quantities
  - the time interval \f$[t,t+(N-1)\Delta t]\f$,
  - time series of \f$N\f$ values of temperature in the snow at equally-spaced times
    \f$t,t+\Delta t,\dots,t+(N-1)\Delta t]\f$, and
  - a scalar precipation rate which is taken to apply in the whole time interval.

It may be that this base class should be more general.  For instance, to allow as
input a time series for precipation rate.
 */
class LocalMassBalance {

public:
  LocalMassBalance();

  virtual PetscErrorCode init();

  /*! T[0],...,T[N-1] are temperatures (K) at times t, t+dt, ..., t+(N-1)dt 
      Input t,dt in seconds.  Input precip and return value are in 
      ice-equivalent thickness per time (m s-1).  Input precip is amount of snow.
      Rain is not modeled.  If input precip is negative then it is treated 
      directly as ablation and positive degree days are ignored.  */
  virtual PetscScalar getMassFluxFromTemperatureTimeSeries(
             PetscScalar t, PetscScalar dt, PetscScalar *T, PetscInt N,
             PetscScalar precip);

protected:
  NCConfigVariable config;
};


//! A PDD implementation which computes the local mass balance based on an expectation integral to get the number of PDDs.
/*!
The expected number of positive degree days is computed by an integral in \ref CalovGreve05 .
 */
class PDDMassBalance : public LocalMassBalance {

public:
  PDDMassBalance();

  virtual PetscErrorCode init();

  virtual PetscScalar getMassFluxFromTemperatureTimeSeries(
             PetscScalar t, PetscScalar dt, PetscScalar *T, PetscInt N,
             PetscScalar precip);

protected:
  /*! Return value is number of positive degree days (units: K day)  */
  virtual PetscScalar getPDDSumFromTemperatureTimeSeries(
                 PetscScalar t, PetscScalar dt, PetscScalar *T, PetscInt N);
  PetscScalar CalovGreveIntegrand(PetscScalar sigma, PetscScalar Tac);

  PetscScalar  pddStdDev,        // K; daily amount of randomness
               pddFactorSnow,    // m day^-1 K^-1; amount of snow melted,
                                 //    as ice equivalent, per positive degree day
               pddFactorIce,     // m day^-1 K^-1; amount of ice melted
                                 //    per positive degree day
               pddRefreezeFrac;  // [pure fraction]; amount of melted snow which refreezes
                                 //    as ice
};


//! An alternative PDD implementation which computes the local mass balance based on simulating a random process to get the number of PDDs.
/*!
Uses a GSL random number generator.  Significantly slower because new random numbers are
generated for each grid point.

The way the number of positive degree-days are used to produce a surface mass balance
is identical to the more basic class PDDMassBalance.
 */
class PDDrandMassBalance : public PDDMassBalance {

public:
  PDDrandMassBalance(bool repeatable); //! repeatable==true to seed with zero every time.
  virtual ~PDDrandMassBalance();

protected:
  virtual PetscScalar getPDDSumFromTemperatureTimeSeries(
                 PetscScalar t, PetscScalar dt, PetscScalar *T, PetscInt N);
  gsl_rng     *pddRandGen;
};

#endif
