// Copyright (C) 2009, 2010, 2011, 2012, 2013, 2014 Ed Bueler and Constantine Khroulev
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
  double  snow, //!< m day^-1 K^-1; ice-equivalent amount of snow melted, per PDD
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

  This model also uses degree day factors passed-in in DegreeDayFactors `ddf`,
  and the standard deviation `pddStdDev`.  The latter is the standard deviation of the
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
  LocalMassBalance(const PISMConfig &myconfig);
  virtual ~LocalMassBalance() {}

  virtual unsigned int get_timeseries_length(double dt) = 0;

  //! Count positive degree days (PDDs).  Returned value in units of K day.
  /*! Inputs T[0],...,T[N-1] are temperatures (K) at times t, t+dt_series, ..., t+(N-1)dt_series.
    Inputs `t`, `dt_series` are in seconds.  */
  virtual void get_PDDs(double pddStdDev, double dt_series,
                        double *T, unsigned int N, double *PDDs) = 0;

  /*! Remove rain from precipitation. */
  virtual void get_snow_accumulation(double *precip_rate, double *T,
                                     unsigned int N) = 0;

  /** 
   * Take a step of the PDD model.
   *
   * @param[in] ddf degree day factors
   * @param[in] PDDs number of positive degree days during the time step [K day]
   * @param[in] accumulation total solid (snow) accumulation during the time-step [ice equivalent meters]
   * @param[in,out] snow_depth snow depth [ice equivalent meters]
   * @param[in,out] cumulative_melt [ice equivalent meters]
   * @param[in,out] cumulative_runoff [ice equivalent meters]
   * @param[in,out] cumulative_smb [ice equivalent meters]
   */
  virtual void step(const DegreeDayFactors &ddf,
                    double PDDs,
                    double accumulation,
                    double &snow_depth,
                    double &cumulative_melt,
                    double &cumulative_runoff,
                    double &cumulative_smb) = 0;

protected:
  const PISMConfig& config;
  PISMUnitSystem m_unit_system;
  const double seconds_per_day;
};


//! A PDD implementation which computes the local mass balance based on an expectation integral.
/*!
  The expected number of positive degree days is computed by an integral in \ref CalovGreve05.

  Uses degree day factors which are location-independent.
*/
class PDDMassBalance : public LocalMassBalance {

public:
  PDDMassBalance(const PISMConfig& myconfig);
  virtual ~PDDMassBalance() {}

  virtual unsigned int get_timeseries_length(double dt);
  virtual void get_PDDs(double pddStdDev, double dt_series,
                        double *T, unsigned int N, double *PDDs);

  virtual void get_snow_accumulation(double *precip_rate, double *T,
                                     unsigned int N);

  virtual void step(const DegreeDayFactors &ddf,
                    double PDDs,
                    double accumulation,
                    double &snow_depth,
                    double &cumulative_melt,
                    double &cumulative_runoff,
                    double &cumulative_smb);

protected:
  double CalovGreveIntegrand(double sigma, double TacC);

  bool precip_as_snow,          //!< interpret all the precipitation as snow (no rain)
    refreeze_ice_melt;          //!< refreeze melted ice
  double Tmin,             //!< the temperature below which all precipitation is snow
    Tmax;             //!< the temperature above which all precipitation is rain
  double pdd_threshold_temp; //!< threshold temperature for the PDD computation
};


//! An alternative PDD implementation which simulates a random process to get the number of PDDs.
/*!
  Uses a GSL random number generator.  Significantly slower because new random numbers are
  generated for each grid point.

  The way the number of positive degree-days are used to produce a surface mass balance
  is identical to the base class PDDMassBalance.

  \note
  A more realistic pattern for the variability of surface melting might have correlation 
  with appropriate spatial and temporal ranges.
*/
class PDDrandMassBalance : public PDDMassBalance {

public:
  PDDrandMassBalance(const PISMConfig& myconfig,
                     bool repeatable); //! repeatable==true to seed with zero every time.
  virtual ~PDDrandMassBalance();

  virtual unsigned int get_timeseries_length(double dt);

  virtual void get_PDDs(double pddStdDev, double dt_series,
                        double *T, unsigned int N, double *PDDs);
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
  FaustoGrevePDDObject(IceGrid &g, const PISMConfig &myconfig);
  virtual ~FaustoGrevePDDObject() {}

  virtual PetscErrorCode update_temp_mj(IceModelVec2S *surfelev, IceModelVec2S *lat, IceModelVec2S *lon);

  /*! If this method is called, it is assumed that i,j is in the ownership range
    for IceModelVec2S temp_mj. */
  virtual PetscErrorCode setDegreeDayFactors(int i, int j,
                                             double /* usurf */, double lat, double /* lon */,
                                             DegreeDayFactors &ddf);

protected:
  IceGrid &grid;
  const PISMConfig &config;
  double beta_ice_w, beta_snow_w, T_c, T_w, beta_ice_c, beta_snow_c,
    fresh_water_density, ice_density, pdd_fausto_latitude_beta_w;
  IceModelVec2S temp_mj;
};


#endif
