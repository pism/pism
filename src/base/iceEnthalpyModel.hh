// Copyright (C) 2009 Andreas Aschwanden and Ed Bueler
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

#ifndef __iceEnthalpyModel_hh
#define __iceEnthalpyModel_hh

#include <petsc.h>
#include "iceModelVec.hh"
#include "NCVariable.hh"
#include "materials.hh"
#include "iceModel.hh"


//! Glen (1955) and Paterson-Budd (1982) flow law with additional water fraction factor from Lliboutry & Duval (1985).
/*!
See \ref AschwandenBlatter2009.  The basic references are \ref Glen and \ref PatersonBudd 
and \ref LliboutryDuval1985.
 */
class PolyThermalGPBLDIce : public ThermoGlenIce {
public:
  PolyThermalGPBLDIce(MPI_Comm c,const char pre[]);
  virtual PetscErrorCode setFromOptions();
  virtual PetscErrorCode view(PetscViewer viewer) const;

  /* these are not literal reimplementations, but new routines.
  to see where they are needed, do in src/base/:
     $ grep ice->flow *.cc
     $ grep ice->effectiveViscosity *.cc
     $ grep ice->softnessParameter *.cc
     etc
  but arrow (>) must be escaped with backslash
  */
  virtual PetscScalar softnessParameterFromEnth(PetscScalar enthalpy, PetscScalar pressure) const;
  virtual PetscScalar hardnessParameterFromEnth(PetscScalar enthalpy, PetscScalar pressure) const;

  virtual PetscScalar flowFromEnth(PetscScalar stress, PetscScalar enthalpy, PetscScalar pressure,
                                       PetscScalar gs) const; // grainsize arg gs not used

  virtual PetscScalar effectiveViscosityColumnFromEnth(
                PetscScalar thickness,  PetscInt kbelowH, const PetscScalar *zlevels,
                PetscScalar u_x,  PetscScalar u_y, PetscScalar v_x,  PetscScalar v_y,
                const PetscScalar *enthalpy1, const PetscScalar *enthalpy2) const;

  /* these are used in src/base/ssaJed/ stuff only, so not addressed for now:
    integratedStoreSize(), integratedStore(), integratedViscosity()
  */

  NCConfigVariable *config;

protected:
  PetscReal water_frac_coeff;
};



//! Temporary class for development of enthalpy-based polythermal PISM.
/*!
Based for now on Bueler's reading of \ref AschwandenBlatter2009.
 */
class IceEnthalpyModel : public IceModel {

public:
  IceEnthalpyModel(IceGrid &g);

  using IceModel::initFromFile;
  virtual PetscErrorCode initFromFile(const char *);

  using IceModel::write_extra_fields;
  virtual PetscErrorCode write_extra_fields(const char filename[]);

  bool doColdIceMethods;

protected:
  using IceModel::createVecs;
  virtual PetscErrorCode createVecs();
  
  using IceModel::init_physics;
  virtual PetscErrorCode init_physics();

  // PetscErrorCode setEnth3toCTSValue(); ???  desired?
  
  virtual PetscErrorCode setEnth3FromT3_ColdIce();
  
  virtual PetscErrorCode setT3FromEnth3();

  virtual PetscErrorCode setUserLiquidFracFromEnthalpy(IceModelVec3 &useForLiquidFrac);

  using IceModel::velocitySIAStaggered;
  virtual PetscErrorCode velocitySIAStaggered();

  using IceModel::computeEffectiveViscosity;
  virtual PetscErrorCode computeEffectiveViscosity(IceModelVec2 vNuH[2], PetscReal epsilon);

  using IceModel::temperatureStep;
  virtual PetscErrorCode temperatureStep(PetscScalar* vertSacrCount, PetscScalar* bulgeCount);
  
  using IceModel::temperatureAgeStep;
  virtual PetscErrorCode temperatureAgeStep();

  virtual PetscErrorCode enthalpyAndDrainageStep(PetscScalar* vertSacrCount, PetscScalar* bulgeCount);

  virtual PetscErrorCode drainageToBaseModelEnth(
                const PetscScalar thickness, const PetscScalar z, const PetscScalar dz,
                PetscScalar &enthalpy, PetscScalar &Hmelt);

protected: // new data members
  IceModelVec3  Enth3, EnthNew3;
};

#endif

