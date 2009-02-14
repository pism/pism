// Copyright (C) 2009 Ed Bueler and Ricarda Winkelmann
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


//! A virtual base class for coupling PISM to other climate components.
/*!
Methods and members here are common to all possible implementations and 
derived classes.
 */
class PISMClimateCoupler {

public:
  PISMClimateCoupler();
  virtual ~PISMClimateCoupler();

  virtual PetscErrorCode initFromOptions(IceGrid* g);

  // since climate fields may be in the same file as the one used by
  //   PISM for input, get the info needed to read them; this is normally
  //   a helper routine for derived classes; filename needs to be pre-allocated
  virtual PetscErrorCode findPISMInputFile(char* filename, LocalInterpCtx* &lic);
  
  virtual PetscErrorCode updateClimateFields(
             const PetscScalar t_years, const PetscScalar dt_years, 
             void *iceInfoNeeded);

  // the implementations of this in the base class just terminates; to use,
  //   re-implement in the derived class
  virtual PetscErrorCode writeCouplingFieldsToFile(const char *filename);

protected:
  IceGrid* grid;
};


/*******************  ATMOSPHERE:  PISMAtmosphereCoupler ********************/

//! An atmosphere model might need to know these things about IceModel to update through PISMAtmosphereCoupler.
struct IceInfoNeededByAtmosphereCoupler {
  // "might need to know" for these reasons:
  IceModelVec2 *lat,      // location dependence
               *lon,      // location dependence
               *mask,     // ice surface type dependence (potentially; e.g. ice shelf vs interior)
               *surfelev; // surface elevation dependence
};


//! A basic derived class of PISMClimateCoupler for coupling PISM to an atmosphere model.
/*!
No files are read to initialize this PISMClimateCoupler, so it is essentially virtual.  Space
for time-dependent surface mass flux (vsurfmassflux) and time-dependent surface temperature 
(vsurftemp) is allocated, however.  A pointer to these is provided by appropriate methods.
It is expected that a derived class of this will actually be used.

The IceModelVec2 members vsurfmassflux and vsurftemp are exactly the surface fields needed 
as boundary conditions for IceModel.  In particular, vsurfmassflux is the instantaneous net
surface mass balance in ice-thickness-equivalent m s-1.  That is vsurfmassflux is a term 
in the mass continuity equation.

And vsurftemp is the temperature (K) at the ice surface but below firn processes.  It is the
upper surface boundary condition for the conservation of energy partial differential equation
within the ice.
 */
class PISMAtmosphereCoupler : public PISMClimateCoupler {

public:
  PISMAtmosphereCoupler();
  virtual ~PISMAtmosphereCoupler(); // destroys IceModelVec2 below

  // next three redefine PISMClimateCoupler versions
  virtual PetscErrorCode initFromOptions(IceGrid* g);
  virtual PetscErrorCode writeCouplingFieldsToFile(const char *filename);
  virtual PetscErrorCode updateClimateFields(
             const PetscScalar t_years, const PetscScalar dt_years, 
             void *iceInfoNeeded);

  // an atmosphere model could run non-trivially during these calls
  virtual PetscErrorCode updateSurfMassFluxAndProvide(
             const PetscScalar t_years, const PetscScalar dt_years, 
             void *iceInfoNeeded, // will be interpreted as type IceInfoNeededByAtmosphereCoupler*
             IceModelVec2* &pvsmf);  // vsmf = pointer to vsurfmassflux
  virtual PetscErrorCode updateSurfTempAndProvide(
             const PetscScalar t_years, const PetscScalar dt_years,
             void *iceInfoNeeded, // will be interpreted as type IceInfoNeededByAtmosphereCoupler*
             IceModelVec2* &pvst);  // vst = pointer to vsurftemp

//FIXME: should be protected:
  IceModelVec2 vsurfmassflux, vsurftemp; // access these through update...() above
};


//! A derived class of PISMAtmosphereCoupler which reads a constant-in-time surface climate from a NetCDF file.
class PISMConstAtmosCoupler : public PISMAtmosphereCoupler {

public:
  PISMConstAtmosCoupler();

  virtual PetscErrorCode initFromOptions(IceGrid* g);

  bool initializeFromFile;  // default is true
};


//! A derived class of PISMAtmosphereCoupler which reads monthly surface temperatures from a NetCDF file.
/*!
Stored temperatures must have names \c temp_mon0, ...,\c temp_mon11 and be in units of K.
These monthly surface temperatures are typically used in a PDD; see PISMPDDCoupler.  These
monthly surface temperatures are not "deep", unlike PISMAtmosphereCoupler::vsurftemp.
Instead they are typically the output of an atmospheric model.  These temperatures are used in
determining snowfall and firn processes.
 */
class PISMMonthlyTempsAtmosCoupler : public PISMAtmosphereCoupler {

public:
  PISMMonthlyTempsAtmosCoupler();
  virtual ~PISMMonthlyTempsAtmosCoupler();

  virtual PetscErrorCode initFromOptions(IceGrid* g);
  virtual PetscErrorCode writeCouplingFieldsToFile(const char *filename);

  bool readMonthlyTempsFromFile;  // default is false
                                  // if set to false this derived class reverts to PISMAtmosphereCoupler
                                  // set this bool *before* call to initFromOptions()
  PetscErrorCode setMonthlyTempsFilename(const char* filename); // call *before* initFromOptions()

protected:
  PetscErrorCode readMonthlyTemps();
/*
  PetscErrorCode getMonthIndicesFromDay(const PetscScalar day, PetscInt &curr, PetscInt &next);
  PetscScalar getTemperatureFromMonthlyData(
       PetscScalar **currMonthSurfTemps, PetscScalar **nextMonthSurfTemps,
       const PetscInt i, const PetscInt j, const PetscScalar day);
*/

private:
  char monthlyTempsFile[PETSC_MAX_PATH_LEN];
  IceModelVec2 vmonthlysurftemp[12]; // if readMonthlyTempsFromFile then allocated by initFromOptions()
};


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
class PISMPDDCoupler : public PISMAtmosphereCoupler {  // implemented in pPDDcoupler.cc

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

  PetscScalar  pddStdDev,        // K; daily amount of randomness
               pddFactorSnow,    // m day^-1 K^-1; amount of snow melted,
                                 //    as ice equivalent, per positive degree day
               pddFactorIce,     // m day^-1 K^-1; amount of ice melted
                                 //    per positive degree day
               pddRefreezeFrac,  // [pure fraction]; amount of melted snow which refreezes
                                 //    as ice
               pddSummerWarming, // K; amplitude of yearly temperature cycle
               pddSummerPeakDay; // Julian day of summer temperature peak

  bool         initialize_vannmeansurftemp_FromFile,  // default is true
               initialize_vsurfaccum_FromFile;        // default is true

// FIXME: should be protected:
  IceModelVec2 vannmeansurftemp, vsurfaccum;

protected:
  gsl_rng      *pddRandGen;      // usually NULL; default is expectation integral which
                                 //   does not use actual random numbers

  PetscTruth   usingMonthlyTemps;
  IceModelVec2 vmonthlysurftemp[12]; // usually these are not created; if user supplies monthly
                                 //   temperature maps then we will allocate 12 IceModelVec2's 
  PetscErrorCode readMonthlyTemps(const char *filename);
  PetscErrorCode getMonthIndicesFromDay(const PetscScalar day, PetscInt &curr, PetscInt &next);
  PetscScalar getTemperatureFromMonthlyData(
       PetscScalar **currMonthSurfTemps, PetscScalar **nextMonthSurfTemps,
       const PetscInt i, const PetscInt j, const PetscScalar day);

  virtual PetscScalar getSummerWarming(
             const PetscScalar elevation, const PetscScalar latitude,
             const PetscScalar Tma);
  virtual PetscScalar getTemperatureFromYearlyCycle(
       const PetscScalar summer_warming, const PetscScalar Tma, const PetscScalar day);

  virtual PetscScalar getSurfaceBalanceFromSnowAndPDD(const PetscScalar snowrate,
                    const PetscScalar dt_secs, const PetscScalar pddsum);
  virtual double CalovGreveIntegrand(const double Tac);  // double because handed to gsl quadrature routine

};


/*******************  OCEAN:  PISMOceanCoupler and derived ********************/

//! An ocean model might need to know these things about IceModel to update through PISMOceanCoupler.
struct IceInfoNeededByOceanCoupler {
  // "might need to know" for these reasons:
  IceModelVec2 *lat,      // location dependence
               *lon,      // location dependence
               *mask,     // ice model type dependence (floating ice shelf vs grounded)
               *thk;      // shelf base elevation dependence; ice thickness gives base elevation
                          //   through floatation criterion; potentially need separate
                          //   base elevation field
};


//! A basic derived class of PISMClimateCoupler for coupling PISM to an ocean model.  Essentially virtual.
/*!
It is expected that a derived class of this will actually be used.
 */
class PISMOceanCoupler : public PISMClimateCoupler {

public:
  PISMOceanCoupler();

  ~PISMOceanCoupler();

  virtual PetscErrorCode initFromOptions(IceGrid* g);

  virtual PetscErrorCode writeCouplingFieldsToFile(const char *filename);

  // an ocean model could run non-trivially during these calls
  virtual PetscErrorCode updateShelfBaseMassFluxAndProvide(
             const PetscScalar t_years, const PetscScalar dt_years, 
             void *iceInfoNeeded, // will be interpreted as type IceInfoNeededByOceanCoupler*
             IceModelVec2* &pvsbmf);  // pvsbmf = pointer to vshelfbasemassflux

  virtual PetscErrorCode updateShelfBaseTempAndProvide(
             const PetscScalar t_years, const PetscScalar dt_years,
             void *iceInfoNeeded, // will be interpreted as type IceInfoNeededByOceanCoupler*
             IceModelVec2* &pvsbt);  // pvsbt = pointer to vshelfbasetemp

  virtual PetscErrorCode updateClimateFields(
             const PetscScalar t_years, const PetscScalar dt_years, 
             void *iceInfoNeeded); // will be interpreted as type IceInfoNeededByOceanCoupler*

  bool      reportInitializationToStdOut;  // can turn off report on initialization if
                                           // there can be no floating ice (for example)

protected:
  IceModelVec2 vshelfbasetemp, vshelfbasemassflux;
};


//! A derived class of PISMOceanCoupler for coupling PISM to an ocean model.  Essentially virtual.
class PISMConstOceanCoupler : public PISMOceanCoupler {

public:
  PISMConstOceanCoupler();

  virtual PetscErrorCode initFromOptions(IceGrid* g);

  virtual PetscErrorCode updateShelfBaseMassFluxAndProvide(
             const PetscScalar t_years, const PetscScalar dt_years, 
             void *iceInfoNeeded, // will be interpreted as type IceInfoNeededByOceanCoupler*
             IceModelVec2* &pvsbmf);  // pvsbmf = pointer to vshelfbasemassflux

  virtual PetscErrorCode updateShelfBaseTempAndProvide(
             const PetscScalar t_years, const PetscScalar dt_years,
             void *iceInfoNeeded, // will be interpreted as type IceInfoNeededByOceanCoupler*
             IceModelVec2* &pvsbt);  // pvsbt = pointer to vshelfbasetemp

  PetscReal constOceanHeatFlux;  // in W m-2; directly converted to constant mass flux
                                 //   by updateShelfBaseMassFluxAndProvide()
};

#endif
