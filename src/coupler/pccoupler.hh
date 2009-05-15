// Copyright (C) 2008-2009 Ed Bueler and Ricarda Winkelmann
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

#ifndef __pccoupler_hh
#define __pccoupler_hh

#include <petsc.h>
#include <gsl/gsl_rng.h>
#include "../base/pism_const.hh"
#include "../base/grid.hh"
#include "../base/LocalInterpCtx.hh"
#include "../base/iceModelVec.hh"
#include "forcing.hh"
#include "localMassBalance.hh"
#include "monthlyDataMaps.hh"


//! A coupler might need to know these things about IceModel to update.  DEPRECATED MECHANISM.
struct IceInfoNeededByCoupler {
  // "might need to know" for these reasons:
  IceModelVec2 *lat,      // location dependence
               *lon,      // location dependence
               *mask,     // ice surface type dependence (potentially; e.g. ice shelf vs interior)
               *thk,      // thickness dependence (relatively unlikely)
               *surfelev, // surface elevation dependence (surface of ice and of exposed bedrock)
               *topg;     // bed elevation dependence (relatively unlikely)
};


//! A virtual base class for coupling PISM to other climate components.
/*!
  Methods and members here are common to all possible implementations and 
  derived classes.

  <p> All the update... methods should interpret the \c t_years and \c dt_years
  arguments as specifying the time interval (t_years, t_years + dt_years). A
  coupler should provide an estimate of a climate field for this interval, which
  may be, but does not have to be, an average over this interval.</p>

  <p> The \c dt_years argument <b> is not </b> "the time since the last call",
  and <b> should not </b> be used to incrementally update climate fields. </p>
*/
class PISMClimateCoupler {

public:
  PISMClimateCoupler();
  virtual ~PISMClimateCoupler();

  // since climate fields may be in the same file as the one used by
  //   PISM for input, get the info needed to read them; this is normally
  //   a helper routine for derived classes; filename needs to be pre-allocated
  // normally only used by derived classes, but also in src/pcctest.cc
  virtual PetscErrorCode findPISMInputFile(char* filename, LocalInterpCtx* &lic);

  virtual PetscErrorCode initFromOptions(IceGrid* g);
  
  virtual PetscErrorCode updateClimateFields(
             const PetscScalar t_years, const PetscScalar dt_years, 
             IceInfoNeededByCoupler* info);

  // the implementations of this in the base class just terminates; to use,
  //   re-implement in the derived class
  virtual PetscErrorCode writeCouplingFieldsToFile(const PetscScalar t_years, const char *filename);

protected:
  IceGrid* grid;
  PetscErrorCode printIfDebug(const char *message);
  NCConfigVariable config;

private:
  bool PCCDEBUG;  // set in constructor; controls printIfDebug() above
};


/*******************  ATMOSPHERE:  PISMAtmosphereCoupler ********************/

//! A basic derived class of PISMClimateCoupler for coupling PISM to an atmosphere model.
/*!
No files are read to initialize this class, so it is essentially a virtual base class 
for atmosphere model coupling.  No process modeling occurs.  It is expected that only
a derived class of this will actually be used.

Space for time-dependent surface mass flux (vsurfmassflux) and time-dependent surface
temperature (vsurftemp) is allocated.  A pointer to these is provided by appropriate
methods, and that is what IceModel needs.  IceModel has a pointer to an instance of this 
class in IceModel::atmosPCC;

The IceModelVec2 members vsurfmassflux and vsurftemp are exactly the surface fields needed 
as boundary conditions for IceModel.  In particular, vsurfmassflux is the instantaneous net
surface mass balance in ice-thickness-equivalent m s-1.  That is, vsurfmassflux is a term 
in the mass continuity equation.

The IceModelVec2 vsurftemp is the temperature (K) at the ice surface but below firn
processes.  E.g. the "10 m" temperature.  It is the upper surface boundary condition 
for the conservation of energy partial differential equation within the ice.

There is no temperature field present in this basic class which should be interpreted
as the snow temperature.  Derived classes may have that, and may model snow/firn processes.

A serious atmosphere model (CAM in CCSM, POTSDAM-2 in CLIMBER3alpha, etc) could run
non-trivially during these calls.  There might be a snow model in a PISMAtmosphereCoupler
component, or (for example) an existing land model might be used for snow processes.
 */
class PISMAtmosphereCoupler : public PISMClimateCoupler {

public:
  PISMAtmosphereCoupler();
  virtual ~PISMAtmosphereCoupler(); // destroys IceModelVec2 below

  // next three redefine PISMClimateCoupler versions
  virtual PetscErrorCode initFromOptions(IceGrid* g);
  virtual PetscErrorCode writeCouplingFieldsToFile(const PetscScalar t_years, const char *filename);
  virtual PetscErrorCode updateClimateFields(
             const PetscScalar t_years, const PetscScalar dt_years, IceInfoNeededByCoupler* info);

  virtual PetscErrorCode updateSurfMassFluxAndProvide(
             const PetscScalar t_years, const PetscScalar dt_years, IceInfoNeededByCoupler* info,
             IceModelVec2* &pvsmf);  // pvsmf = pointer to vsurfmassflux
  virtual PetscErrorCode updateSurfTempAndProvide(
             const PetscScalar t_years, const PetscScalar dt_years, IceInfoNeededByCoupler* info,
             IceModelVec2* &pvst);  // pvst = pointer to vsurftemp

protected:
  IceModelVec2 vsurfmassflux, vsurftemp; // access these through update...() above
  PetscReal           TsOffset;
  PISMClimateForcing* dTforcing; 
};


//! A derived class of PISMAtmosphereCoupler which reads a constant-in-time surface climate from a NetCDF file.
/*!
There are no redefinitions of updateSurfMassFluxAndProvide() and updateSurfTempAndProvide().
Those procedures function the same way as in PISMAtmosphereCoupler, but here the effect is to provide
access to IceModelVec2 s which are read once, at the beginning of the run, and are not modified.

Some users of this class will read the fields vsurfmassflux and vsurftemp from files at the
beginning of the run and then not change them.  The default is to do this.

Other users of this class will set the vsurfmassflux and/or vsurftemp fields from formulas. 
Examples of this are verification cases and certain simplified geometry experiments (e.g.
EISMINT II).  In this case, set initializeFromFile=false before calling initFromOptions().
 */
class PISMConstAtmosCoupler : public PISMAtmosphereCoupler {

public:
  PISMConstAtmosCoupler();

  virtual PetscErrorCode initFromOptions(IceGrid* g);

  bool initializeFromFile;  // default is true
};


//! A derived class of PISMAtmosphereCoupler which has a snow temperature parameterization and a choice of PDD models for surface mass balance.
/*!
Lots of reading of user options, including choice of monthly temperatures or a parameterization
(the default), choice of \ref CalovGreve05 (default) or random PDD computation, and choice of
parameters in PDD.

The default temperature parameterization will be from \ref Faustoetal2009.
 */
class PISMSnowModelAtmosCoupler : public PISMAtmosphereCoupler {

public:
  PISMSnowModelAtmosCoupler();
  virtual ~PISMSnowModelAtmosCoupler();

  //! Check if user wants a snow temperature and mass balance model.
  /*! Determines if any of a list of related options are set.
      If so, returns true, otherwise false. */
  virtual bool optionsChooseSnowModel();

  virtual PetscErrorCode initFromOptions(IceGrid* g);

  virtual PetscErrorCode writeCouplingFieldsToFile(const PetscScalar t_years, const char *filename);

  virtual PetscErrorCode updateSurfMassFluxAndProvide(
             const PetscScalar t_years, const PetscScalar dt_years, IceInfoNeededByCoupler* info,
             IceModelVec2* &pvsmf);

protected:
  //! Defaults to the \ref Faustoetal2009 scheme.  Called when no monthly temperature maps are available.
  /*!
    Computes the mean annual temperature as function of latitude, longitude, surface elevation, and so on.
    Depends on the current state of IceModel fields, through pointer info.
   */
  virtual PetscErrorCode parameterizedUpdateSnowTemp(IceInfoNeededByCoupler* info);

  //! Instead of temperature parameterization we can used monthly temperature maps read from a file.
  MonthlyDataMaps  *monthlysnowtemps;

  LocalMassBalance *mbscheme;

  //! The snow precipitation rate in ice-equivalent meters per second.
  /*! vsurfmassflux is computed by LocalMassBalance scheme from this rate.  */
  IceModelVec2 vsnowprecip;

  //! The mean annual snow temperature used in the mass balance scheme.
  /*! 
    Note vsurftemp is ice temperature below completion of firn processes, and is used as
    upper boundary condition for the conservation of energy scheme.  This temperature 
    is related is the snow processes model (LocalMassBalance) which converts the snow precipitation rate
    vsnowprecip into surface mass flux vsurfmassflux.  It is frequently the +2 m temperature, above
    the snow surface, from an atmospheric flow/energy model.  This temperature may come from a
    temperature parameterization (\ref Faustoetal2009, for example) or from stored monthly temperature
    (MonthlyDataMaps).

    It is a diagnostic extra output of PISMSnowModelAtmosCoupler.  IceModel never gets a pointer to it.
   */
  IceModelVec2 vsnowtemp_ma,
               vsnowtemp_mj;  //!< The mean July (julian day = 196; July 15) snow temperature used in the mass balance scheme.
};


//! DEPRECATED: A derived class of PISMAtmosphereCoupler which reads monthly surface temperatures from a NetCDF file.
/*!
Stored temperatures must have names \c monsnowtemp1, ...,\c monsnowtemp12 and be in units of K.
These monthly surface temperatures are typically used in a PDD; see PISMPDDCoupler.  These
monthly snow temperatures are for the layers above many firn processes, so they are not directly
meaningful as the upper surface boundary condition for the ice; compare PISMAtmosphereCoupler::vsurftemp.
Instead they are typically the output of an atmospheric model.  They correspond to the
temperature which controls melting and refreezing processes within the snow.

There are no redefinitions of updateSurfMassFluxAndProvide() and updateSurfTempAndProvide().
Those procedures function the same way as in PISMAtmosphereCoupler, so the fields vsurfmassflux
and vsurftemp are not completely initialized.

Furthermore, the read monthly snow temperatures have no effect here, because there is nothing 
like a PDD in which snow temperatures are used to determine mass flux to to model firn processes
which affect the temperature below the firn layer.  Thus this derived class should not be used
directly.  Instead, use PISMPDDCoupler or a more sophisticated derived class.
 */
class PISMMonthlyTempsAtmosCoupler : public PISMAtmosphereCoupler {

public:
  PISMMonthlyTempsAtmosCoupler();
  virtual ~PISMMonthlyTempsAtmosCoupler();

  virtual PetscErrorCode initFromOptions(IceGrid* g);
  virtual PetscErrorCode writeCouplingFieldsToFile(const PetscScalar t_years, const char *filename);

  // default for this boolean is false;  if set to false this derived class completely reverts
  //   to PISMAtmosphereCoupler; set this boolean *before* call to initFromOptions()
  bool readMonthlyTempsFromFile;

  PetscErrorCode setMonthlyTempsFilename(const char* filename); // call *before* initFromOptions()

protected:
  PetscErrorCode readMonthlyTemps();
  PetscErrorCode getMonthIndicesFromDay(PetscScalar day, 
                                        PetscInt &curr, PetscInt &next, PetscScalar &lambda);
  PetscScalar getTemperatureFromMonthlyData(
       PetscScalar **currMonthTemps, PetscScalar **nextMonthTemps, PetscScalar lambda,
       PetscInt i, PetscInt j);

  IceModelVec2 vmonthlytemp[12]; // if readMonthlyTempsFromFile then allocated by initFromOptions()

private:
  char monthlyTempsFile[PETSC_MAX_PATH_LEN];
};


//! DEPRECATED: A derived class, a descendant of PISMAtmosphereCoupler and PISMMonthlyTempsAtmosCoupler,  which provides a positive degree day model to PISM.
/*!
The PDD here defaults to the one used for EISMINT-Greenland.

If <tt>-pdd_monthly_temps</tt> is used then it reads 12 monthly temperature
data sets, 'monsnowtemp1', ..., 'monsnowtemp12', from input file.  These are used in
the PDD to determine surface mass flux, but do not affect the reported surface
temperature, which is supposed to be the temperature below the firn.

It reads 'snowaccum' from input file and interprets it as ice-equivalent snow
accumulation rate, vsnowaccum.  If 'snowaccum' is \e not present, reads 'acab' from input file 
and \e interprets it as ice-equivalent snow accumulation rate, vsnowaccum; the user
is warned in this case. Note vsnowaccum is always writen to output files as 'snowaccum'.

The PDD is used to convert vsnowaccum to ice-equivalent net surface mass flux,
vsurfmassflux.  Monthly snow temperatures are used in this process, if present.  They
are read from a file as described for PISMMonthlyTempsAtmosCoupler.

It has various constants parameterizing the melt and refreeze processes.  See
\ref RitzEISMINT .

This derived class is implemented in src/coupler/pPDDcoupler.cc.
 */
class PISMPDDCoupler : public PISMMonthlyTempsAtmosCoupler {

public:
  PISMPDDCoupler();
  virtual ~PISMPDDCoupler();  // destroys PDD

  PetscErrorCode userOptionsChoosePDD(PetscTruth &userWantsPDD);
  
  virtual PetscErrorCode initFromOptions(IceGrid* g); // reads both vsurftemp and vsnowaccum
  
  virtual PetscErrorCode writeCouplingFieldsToFile(const PetscScalar t_years, const char *filename);

  virtual PetscErrorCode updateSurfMassFluxAndProvide(
             const PetscScalar t_years, const PetscScalar dt_years, IceInfoNeededByCoupler* info,
             IceModelVec2* &pvsmf);  // vsmf = pointer to vsurfmassflux

  // updateSurfTempAndProvide() is *NOT* reimplemented; vsurftemp in PISMAtmosphereCoupler
  //   has the same interpretation here

  PetscScalar  pddStdDev,        // K; daily amount of randomness
               pddFactorSnow,    // m day^-1 K^-1; amount of snow melted,
                                 //    as ice equivalent, per positive degree day
               pddFactorIce,     // m day^-1 K^-1; amount of ice melted
                                 //    per positive degree day
               pddRefreezeFrac,  // [pure fraction]; amount of melted snow which refreezes
                                 //    as ice
               pddSummerWarming, // K; amplitude of yearly temperature cycle
               pddSummerPeakDay; // Julian day of summer temperature peak

  bool         initialize_vsurftemp_FromFile,  // default is true
               initialize_vsnowaccum_FromFile; // default is true

protected:
  IceModelVec2 vsnowaccum,   // vsurfmassflux is computed by PDD from this amount of snow (ice-equivalent units)
               vsurftempPDD; // vsurftemp is temperature below firn; this temperature is the time-dependent snow
                             //    temperature computed from getTemperatureFromYearlyCycle() or from stored
                             //    monthly temperatures; in any case this one is a diagnostic extra output of
                             //    the PDD; IceModel never gets a pointer to it; updated by
                             //    updateSurfMassFluxAndProvide()

protected:
  gsl_rng      *pddRandGen;  // usually NULL; default is expectation integral which
                             //   does not use actual random numbers

  virtual PetscScalar getSummerWarming(
             const PetscScalar elevation, const PetscScalar latitude,
             const PetscScalar Tma);
  virtual PetscScalar getTemperatureFromYearlyCycle(
             const PetscScalar summer_warming, const PetscScalar Tma, const PetscScalar day);
  virtual PetscScalar getSurfaceBalanceFromSnowAndPDD(
             const PetscScalar snowrate, const PetscScalar dt_secs, const PetscScalar pddsum);
  virtual double CalovGreveIntegrand(const double Tac);  // handed to gsl quadrature routine

};


/*******************  OCEAN:  PISMOceanCoupler and derived ********************/

//! A basic derived class of PISMClimateCoupler for coupling PISM to an ocean model.  Essentially virtual.
/*!
It is expected that a derived class of this will actually be used.
 */
class PISMOceanCoupler : public PISMClimateCoupler {

public:
  PISMOceanCoupler();

  ~PISMOceanCoupler();

  virtual PetscErrorCode initFromOptions(IceGrid* g);

  virtual PetscErrorCode writeCouplingFieldsToFile(const PetscScalar t_years, const char *filename);

  // an ocean model could run non-trivially during these calls
  virtual PetscErrorCode updateShelfBaseMassFluxAndProvide(
             const PetscScalar t_years, const PetscScalar dt_years, IceInfoNeededByCoupler* info,
             IceModelVec2* &pvsbmf);  // pvsbmf = pointer to vshelfbasemassflux

  virtual PetscErrorCode updateShelfBaseTempAndProvide(
             const PetscScalar t_years, const PetscScalar dt_years, IceInfoNeededByCoupler* info,
             IceModelVec2* &pvsbt);  // pvsbt = pointer to vshelfbasetemp

  virtual PetscErrorCode updateClimateFields(
             const PetscScalar t_years, const PetscScalar dt_years, IceInfoNeededByCoupler* info);

  virtual PetscErrorCode updateSeaLevelElevation(PetscReal t_years, PetscReal dt_years,
						 PetscReal *new_sea_level);

  bool      reportInitializationToStdOut;  // can turn off report on initialization if
                                           // there can be no floating ice (for example)

protected:
  IceModelVec2 vshelfbasetemp, vshelfbasemassflux;
  PISMClimateForcing*  dSLforcing; // possibly contains sea level offset data, e.g. from sea bed core data
  PetscReal    seaLevel; /*!<  IceModel will read PISMOceanCoupler to determine surface elevation
                               of floating ice, and thus the grounding line.  seaLevel here can
                               be defined as the surface elevation of zero thickness floating ice. */
};


//! A derived class of PISMOceanCoupler for coupling PISM to an ocean model.  Essentially virtual.
class PISMConstOceanCoupler : public PISMOceanCoupler {

public:
  PISMConstOceanCoupler();

  virtual PetscErrorCode initFromOptions(IceGrid* g);

  virtual PetscErrorCode writeCouplingFieldsToFile(const PetscScalar t_years, const char *filename);

  virtual PetscErrorCode updateShelfBaseMassFluxAndProvide(
             const PetscScalar t_years, const PetscScalar dt_years, IceInfoNeededByCoupler* info,
             IceModelVec2* &pvsbmf);  // pvsbmf = pointer to vshelfbasemassflux

  virtual PetscErrorCode updateShelfBaseTempAndProvide(
             const PetscScalar t_years, const PetscScalar dt_years, IceInfoNeededByCoupler* info,
             IceModelVec2* &pvsbt);  // pvsbt = pointer to vshelfbasetemp

  PetscReal constOceanHeatFlux;  // in W m-2; directly converted to constant mass flux
                                 //   by updateShelfBaseMassFluxAndProvide()
};

#endif
