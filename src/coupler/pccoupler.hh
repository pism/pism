// Copyright (C) 2008-2009 Ed Bueler, Constantine Khroulev, and Ricarda Winkelmann
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
#include "../base/pism_const.hh"
#include "../base/grid.hh"
#include "../base/LocalInterpCtx.hh"
#include "../base/iceModelVec.hh"
#include "../base/Timeseries.hh"
#include "../base/PISMVars.hh"

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

  virtual PetscErrorCode initFromOptions(IceGrid* g, const PISMVars &variables);
  
  virtual PetscErrorCode findPISMInputFile(char* filename, LocalInterpCtx* &lic,
					   bool &regrid, int &start);

  //! A virtual method which writes fields associated to the derived class.
  /*!
    The base-class implementation merely terminates with an error.
    Re-implement in the derived class.
  */
  virtual PetscErrorCode writeCouplingFieldsToFile(
             PetscScalar t_years, const char *filename);

  virtual PetscErrorCode updateClimateFields(
             PetscScalar t_years, PetscScalar dt_years);

  virtual PetscErrorCode max_timestep(PetscScalar t_years, PetscScalar &dt_years);

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
as the snow temperature or the air temperature just above the surface.  Derived
classes may have that, and may model snow/firn processes.

A serious atmosphere model (CAM in CCSM, POTSDAM-2 in CLIMBER3alpha, etc) could run
non-trivially during these calls.  There might be a snow model in a PISMAtmosphereCoupler
component, or (for example) an existing land model might be used for snow processes.
 */
class PISMAtmosphereCoupler : public PISMClimateCoupler {

public:
  PISMAtmosphereCoupler();
  virtual ~PISMAtmosphereCoupler(); // destroys IceModelVec2 below

  using PISMClimateCoupler::initFromOptions; 
  virtual PetscErrorCode initFromOptions(IceGrid* g, const PISMVars &variables);

  using PISMClimateCoupler::writeCouplingFieldsToFile;
  virtual PetscErrorCode writeCouplingFieldsToFile(
             PetscScalar t_years, const char *filename);

  using PISMClimateCoupler::updateClimateFields;
  virtual PetscErrorCode updateClimateFields(
             PetscScalar t_years, PetscScalar dt_years);

  virtual PetscErrorCode updateSurfMassFluxAndProvide(
             PetscScalar t_years, PetscScalar dt_years,
             IceModelVec2* &pvsmf);  // pvsmf = pointer to vsurfmassflux

  virtual PetscErrorCode updateSurfTempAndProvide(
             PetscScalar t_years, PetscScalar dt_years,
             IceModelVec2* &pvst);  // pvst = pointer to vsurftemp

  using PISMClimateCoupler::max_timestep;
  virtual PetscErrorCode max_timestep(PetscScalar t_years, PetscScalar &dt_years);

protected:
  IceModelVec2    vsurfmassflux, vsurftemp; // access these through update...()

  Timeseries*     dTforcing; 
  PetscReal       TsOffset;

  PetscTruth      doForceToThickness;
  IceModelVec2    vthktarget;
  IceModelVec2    *ftt_thk;  //!< pointer back to IceModel::vH, needed only if doForceToThickness
  PetscReal       ftt_alphadecay;
};


//! A derived class of PISMAtmosphereCoupler which reads a constant-in-time surface climate from a NetCDF file.
/*!
There are no redefinitions of updateSurfMassFluxAndProvide() and 
updateSurfTempAndProvide().  Those procedures function the same way as in 
PISMAtmosphereCoupler, but here the effect is to provide access to 
IceModelVec2 s which are read once at the beginning of the run.
In particular, the default is to read fields vsurfmassflux and vsurftemp 
from a file at the beginning of the run and then not change them.

Other applications using this class will set the vsurfmassflux and/or 
vsurftemp fields from formulas.  Examples of this are verification cases
and certain simplified geometry experiments (e.g. EISMINT II).  In this 
case, set initializeFromFile=false before calling initFromOptions().
 */
class PISMConstAtmosCoupler : public PISMAtmosphereCoupler {

public:
  PISMConstAtmosCoupler();

  virtual PetscErrorCode initFromOptions(IceGrid* g, const PISMVars &variables);

  bool initializeFromFile;  // default is true
};



/*******************  OCEAN:  PISMOceanCoupler and derived ********************/

//! A basic derived class of PISMClimateCoupler for coupling PISM to an ocean model.  Essentially virtual.
/*!
It is expected that a derived class of this will actually be used.
 */
class PISMOceanCoupler : public PISMClimateCoupler {

public:
  PISMOceanCoupler();

  virtual ~PISMOceanCoupler();

  using PISMClimateCoupler::initFromOptions;
  virtual PetscErrorCode initFromOptions(IceGrid* g, const PISMVars &variables);

  using PISMClimateCoupler::writeCouplingFieldsToFile;
  virtual PetscErrorCode writeCouplingFieldsToFile(PetscScalar t_years, const char *filename);

  using PISMClimateCoupler::updateClimateFields;
  virtual PetscErrorCode updateClimateFields(
             PetscScalar t_years, PetscScalar dt_years);

  // an ocean model could run non-trivially during these calls
  virtual PetscErrorCode updateShelfBaseMassFluxAndProvide(
             PetscScalar t_years, PetscScalar dt_years,
             IceModelVec2* &pvsbmf);  // pvsbmf = pointer to vshelfbasemassflux

  virtual PetscErrorCode updateShelfBaseTempAndProvide(
             PetscScalar t_years, PetscScalar dt_years,
             IceModelVec2* &pvsbt);  // pvsbt = pointer to vshelfbasetemp

  virtual PetscErrorCode updateSeaLevelElevation(PetscReal t_years, PetscReal dt_years,
						 PetscReal *new_sea_level);

  bool      reportInitializationToStdOut;  // can turn off report on initialization if
                                           // there can be no floating ice (for example)

protected:
  IceModelVec2 vshelfbasetemp, vshelfbasemassflux;
  Timeseries*  dSLforcing; // possibly contains sea level offset data, e.g. from sea bed core data
  PetscReal    seaLevel; /*!<  IceModel will read PISMOceanCoupler to determine surface elevation
                               of floating ice, and thus the grounding line.  seaLevel here can
                               be defined as the surface elevation of zero thickness floating ice. */
};


//! A derived class of PISMOceanCoupler for coupling PISM to an ocean model.  Essentially virtual.
class PISMConstOceanCoupler : public PISMOceanCoupler {

public:
  PISMConstOceanCoupler();

  using PISMOceanCoupler::initFromOptions;
  virtual PetscErrorCode initFromOptions(IceGrid* g, const PISMVars &variables);

  using PISMOceanCoupler::writeCouplingFieldsToFile;
  virtual PetscErrorCode writeCouplingFieldsToFile(PetscScalar t_years, const char *filename);

  using PISMOceanCoupler::updateShelfBaseMassFluxAndProvide;
  virtual PetscErrorCode updateShelfBaseMassFluxAndProvide(
             PetscScalar t_years, PetscScalar dt_years,
             IceModelVec2* &pvsbmf);  // pvsbmf = pointer to vshelfbasemassflux

  using PISMOceanCoupler::updateShelfBaseTempAndProvide;
  virtual PetscErrorCode updateShelfBaseTempAndProvide(
             PetscScalar t_years, PetscScalar dt_years,
             IceModelVec2* &pvsbt);  // pvsbt = pointer to vshelfbasetemp

  PetscReal constOceanHeatFlux;  // in W m-2; directly converted to constant mass flux
                                 //   by updateShelfBaseMassFluxAndProvide()
protected:
  IceModelVec2 *thk;  //!< pointer back to IceModel::vH
};

#endif

