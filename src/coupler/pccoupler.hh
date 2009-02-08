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


//! An atmosphere model might need to know these things about IceModel to update through PISMAtmosphereCoupler.
struct IceInfoNeededByAtmosphereCoupler {
  // "might need to know" for these reasons:
  IceModelVec2 *lat,      // location dependence
               *lon,      // location dependence
               *mask,     // ice surface type dependence (potentially; e.g. ice shelf vs interior)
               *surfelev; // surface elevation dependence
};


//! A basic derived class of PISMClimateCoupler for coupling PISM to an atmosphere model.  Essentially virtual.
/*!
It is expected that a derived class of this will actually be used.
 */
class PISMAtmosphereCoupler : public PISMClimateCoupler {

public:
  PISMAtmosphereCoupler();
  virtual ~PISMAtmosphereCoupler(); // destroys IceModelVec2 below

  // next three override PISMClimateCoupler versions
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

protected:
  IceModelVec2 vsurfmassflux, vsurftemp; // access these through update...() above
};


//! A derived class of PISMAtmosphereCoupler which provides a constant-in-time surface climate.
/*!
Reads surface temperature and mass balance from a NetCDF file.
 */
class PISMConstAtmosCoupler : public PISMAtmosphereCoupler {

public:
  PISMConstAtmosCoupler();

  virtual PetscErrorCode initFromOptions(IceGrid* g);

  // because climate is constant, no update occurs in these; they just provide a pointer
  virtual PetscErrorCode updateSurfMassFluxAndProvide(
             const PetscScalar t_years, const PetscScalar dt_years, 
             void *iceInfoNeeded, // will be interpreted as type IceInfoNeededByAtmosphereCoupler*
             IceModelVec2* &pvsmf);  // vsmf = pointer to vsurfmassflux
  virtual PetscErrorCode updateSurfTempAndProvide(
             const PetscScalar t_years, const PetscScalar dt_years,
             void *iceInfoNeeded, // will be interpreted as type IceInfoNeededByAtmosphereCoupler*
             IceModelVec2* &pvst);  // vst = pointer to vsurftemp

  bool initializeFromFile;
};


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

protected:
  IceModelVec2 vshelfbasetemp, vshelfbasemassflux;
};


//! A derived class of PISMOceanCoupler for coupling PISM to an ocean model.  Essentially virtual.
class PISMConstOceanCoupler : public PISMOceanCoupler {

public:
  PISMConstOceanCoupler();

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

