// Copyright (C) 2009, 2010, 2011, 2012 Ed Bueler and Constantine Khroulev
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
#include <gsl/gsl_rng.h>
#include "NCVariable.hh"
#include "iceModelVec.hh"  // only needed for FaustoGrevePDDObject


//! A struct which holds degree day factors.
/*!
Degree day factors convert positive degree days (=PDDs) into amount of melt.
 */
struct DegreeDayFactors {
  PetscScalar  snow, //!< m day^-1 K^-1; ice-equivalent amount of snow melted, per PDD
               ice,  //!< m day^-1 K^-1; ice-equivalent amount of ice melted, per PDD
               refreezeFrac;  //!< fraction of melted snow which refreezes as ice
};


//! \brief Base class for a model which computes surface mass flux rate (ice
//! thickness per time) from precipitation and temperature.
/*!
This is a process model.  At each spatial location, it uses a 1D array, with a
time dimension, for the temperature used in melting snow or ice.  At each spatial
location it assumes the precipitation is time-independent.

This process model does not know its location on the ice sheet, but
simply computes the surface mass balance from three quantities:
  - the time interval \f$[t,t+\Delta t]\f$,
  - time series of values of surface temperature at \f$N\f$ equally-spaced
    times in the time interval
  - a scalar precipation rate which is taken to apply in the time interval.

This model also uses degree day factors passed-in in DegreeDayFactors \c ddf,
and the standard deviation \c pddStdDev.  The latter is the standard deviation of the
modeled temperature away from the input temperature time series which contains
the part of location-dependent temperature cycle on the time interval.

\note 
  \li Please avoid using config.get("...") and config.get_flag("...") calls
inside those methods of this class which are called inside loops over 
spatial grids.  Doing otherwise increases computational costs.
  \li This base class should be more general.  For instance, it could allow as
input a time series for precipation rate.
 */
class LocalMassBalance {
public:
  LocalMassBalance(const NCConfigVariable &myconfig)
    : config(myconfig) {}
  virtual ~LocalMassBalance() {}
  virtual PetscErrorCode init() { return 0; };

  /*! Call before getMassFluxFromTemperatureTimeSeries() so that mass balance method can
      decide how to cut up the time interval.  Most implementations will ignore
      t and just use dt.  Input t,dt in seconds.  */
  virtual PetscErrorCode getNForTemperatureSeries(
                PetscScalar t, PetscScalar dt, PetscInt &N) = 0;

  //! Count positive degree days (PDDs).  Returned value in units of K day.
  /*! Inputs T[0],...,T[N-1] are temperatures (K) at times t, t+dt_series, ..., t+(N-1)dt_series.
      Inputs \c t, \c dt_series are in seconds.  */
  virtual PetscScalar getPDDSumFromTemperatureTimeSeries(
                 PetscScalar pddStdDev, PetscScalar pddThresholdTemp,
                 PetscScalar t, PetscScalar dt_series, PetscScalar *T, PetscInt N) = 0;

  /*! Remove rain from precipitation.  Returned value is amount of snow in ice-equivalent m. */
  /*! Inputs \c precip_rate is in ice-equivalent m s-1.  Note
      <tt>dt = N * dt_series</tt> is the full time-step.  */
  virtual PetscScalar getSnowFromPrecipAndTemperatureTimeSeries(
                 PetscScalar precip_rate,
                 PetscScalar t, PetscScalar dt_series, PetscScalar *T, PetscInt N) = 0;

  /*! Input \c dt is in seconds.  Input \c pddsum is in K day.  
      Input \c snow is in ice-equivalent m. 
      Returned mass fluxes, including \c accumulation_rate, \c melt_rate,
      \c runoff_rate, and \c smb (= surface mass balance), are in 
      ice-equivalent thickness per time (m s-1).  */
  virtual PetscErrorCode getMassFluxesFromPDDs(const DegreeDayFactors &ddf,
                                               PetscScalar dt, PetscScalar pddsum,
                                               PetscScalar snow,
                                               PetscScalar &accumulation_rate,
                                               PetscScalar &melt_rate,
                                               PetscScalar &runoff_rate,
                                               PetscScalar &smb_rate) = 0;

protected:
  const NCConfigVariable& config;
};


//! A PDD implementation which computes the local mass balance based on an expectation integral.
/*!
The expected number of positive degree days is computed by an integral in \ref CalovGreve05.

Uses degree day factors which are location-independent.
 */
class PDDMassBalance : public LocalMassBalance {

public:
  PDDMassBalance(const NCConfigVariable& myconfig);
  virtual ~PDDMassBalance() {}

  virtual PetscErrorCode getNForTemperatureSeries(
             PetscScalar t, PetscScalar dt, PetscInt &N);

  virtual PetscScalar getPDDSumFromTemperatureTimeSeries(
                 PetscScalar pddStdDev, PetscScalar pddThresholdTemp,
                 PetscScalar t, PetscScalar dt_series, PetscScalar *T, PetscInt N);

  virtual PetscScalar getSnowFromPrecipAndTemperatureTimeSeries(
                 PetscScalar precip_rate,
                 PetscScalar t, PetscScalar dt_series, PetscScalar *T, PetscInt N);

  virtual PetscErrorCode getMassFluxesFromPDDs(const DegreeDayFactors &ddf,
                                               PetscScalar dt, PetscScalar pddsum,
                                               PetscScalar snow,
                                               PetscScalar &accumulation_rate,
                                               PetscScalar &melt_rate,
                                               PetscScalar &runoff_rate,
                                               PetscScalar &smb_rate);

protected:
  PetscScalar CalovGreveIntegrand(PetscScalar sigma, PetscScalar TacC);

  bool precip_as_snow;          //!< interpret all the precipitation as snow (no rain)
  PetscScalar Tmin,             //!< the temperature below which all precipitation is snow
              Tmax;             //!< the temperature above which all precipitation is rain
};


//! An alternative PDD implementation which simulates a random process to get the number of PDDs.
/*!
Uses a GSL random number generator.  Significantly slower because new random numbers are
generated for each grid point.

The way the number of positive degree-days are used to produce a surface mass balance
is identical to the base class PDDMassBalance.

\note
  \li A more realistic pattern for the variability of surface melting might have correlation 
with appropriate spatial and temporal ranges.
 */
class PDDrandMassBalance : public PDDMassBalance {

public:
  PDDrandMassBalance(const NCConfigVariable& myconfig, bool repeatable); //! repeatable==true to seed with zero every time.
  virtual ~PDDrandMassBalance();

  virtual PetscErrorCode getNForTemperatureSeries(
                PetscScalar t, PetscScalar dt, PetscInt &N);

  virtual PetscScalar getPDDSumFromTemperatureTimeSeries(
               PetscScalar pddStdDev, PetscScalar pddThresholdTemp,
               PetscScalar t, PetscScalar dt_series, PetscScalar *T, PetscInt N);

protected:
  gsl_rng *pddRandGen;
};


/*!
The PDD scheme described by Formula (6) in [\ref Faustoetal2009] requires 
special knowledge of latitude and mean July temp to set degree day factors 
for Greenland.

These formulas are inherited by [\ref Faustoetal2009] from [\ref Greve2005geothermal].
There was, apparently, tuning in [\ref Greve2005geothermal] which mixed ice
dynamical ideas and surface process ideas.  That is, these formulas and parameter
choices arise from looking at margin shape.  This may not be a good source of
PDD parameters.

This may become a derived class of a LocationDependentPDDObject, if the idea
is needed more in the future.
 */
class FaustoGrevePDDObject {

public:
  FaustoGrevePDDObject(IceGrid &g, const NCConfigVariable &myconfig);
  virtual ~FaustoGrevePDDObject() {}

  virtual PetscErrorCode update_temp_mj(IceModelVec2S *surfelev, IceModelVec2S *lat, IceModelVec2S *lon);

  /*! If this method is called, it is assumed that i,j is in the ownership range
      for IceModelVec2S temp_mj. */
  virtual PetscErrorCode setDegreeDayFactors(
              PetscInt i, PetscInt j,
              PetscScalar /* usurf */, PetscScalar lat, PetscScalar /* lon */,
              DegreeDayFactors &ddf);

protected:
  IceGrid &grid;
  const NCConfigVariable &config;
  PetscScalar beta_ice_w, beta_snow_w, T_c, T_w, beta_ice_c, beta_snow_c,
              fresh_water_density, ice_density, pdd_fausto_latitude_beta_w;
  IceModelVec2S temp_mj;
};


#endif
