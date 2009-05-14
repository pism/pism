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

//! Base class for a model which computes surface mass flux rate (ice thickness per time) from a precipitation (scalar) and a time series for temperature.
/*!
This is a process model.  It uses a 1D array with snow temperatures.  The fundamental
idea is that this mass balance model does not know the location on the ice sheet, but
simply computes the surface mass balance from three quantities
  - the time interval,
  - time series for temperature in the snow, and
  - a scalar precipation rate which is taken to apply in the whole interval.

It may be that this base class should be more general.  For instance, to allow as
input a time series for precipation rate.

Provides ???

Example:  See PISM???AtmosCoupler for example use  (FIXME: under development).
 */
class LocalMassBalance {

public:
  LocalMassBalance() {}

  /*! T[0],...,T[N-1] are temperatures (K) at times t, t+dt, ..., t+(N-1)dt 
      Input t,dt in seconds.  Input precip and return value are in 
      ice-equivalent thickness per time (m s-1).  */
  virtual PetscScalar getMassFluxFromTemperatureTimeSeries(
             PetscScalar t, PetscScalar dt, PetscScalar *T, PetscInt N,
             PetscScalar precip);
};


//! A PDD implementation which computes the local mass balance.
/*!
Based on the expectation integral in \ref CalovGreve05 .
 */
class PDDMassBalance : public LocalMassBalance {

public:
  PDDMassBalance();

  virtual PetscScalar getMassFluxFromTemperatureTimeSeries(
             PetscScalar t, PetscScalar dt, PetscScalar *T, PetscInt N,
             PetscScalar precip);

protected:
  /*! Return value is number of positive degree days (units: K day)  */
  PetscScalar getPDDSumFromTemperatureTimeSeries(
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

#endif
