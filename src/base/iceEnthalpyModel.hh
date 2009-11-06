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
#include "iceModel.hh"
#include "enthalpyConverter.hh"


//! Temporary class for development of enthalpy-based polythermal PISM.
/*!
Based on Bueler's reading of \ref AschwandenBlatter.
 */
class IceEnthalpyModel : public IceModel {

public:
  IceEnthalpyModel(IceGrid &g);

  using IceModel::setFromOptions;
  PetscErrorCode setFromOptions();

  using IceModel::initFromFile;
  virtual PetscErrorCode initFromFile(const char *);

  using IceModel::bootstrapFromFile;
  virtual PetscErrorCode bootstrapFromFile(const char *filename);

  using IceModel::write_extra_fields;
  virtual PetscErrorCode write_extra_fields(const char* filename);

protected:
  using IceModel::createVecs;
  virtual PetscErrorCode createVecs();
  
  using IceModel::init_physics;
  virtual PetscErrorCode init_physics();

  virtual PetscErrorCode setEnth3FromT3_ColdIce();
  
  virtual PetscErrorCode setTnew3FromEnth3();

  virtual PetscErrorCode setLiquidFracFromEnthalpy(IceModelVec3 &useForLiquidFrac);

  virtual PetscErrorCode setCTSFromEnthalpy(IceModelVec3 &useForCTS);

  virtual PetscErrorCode setPATempFromEnthalpy(IceModelVec3 &useForPATemp);

  using IceModel::energyAgeStats;
  virtual PetscErrorCode energyAgeStats(
                    PetscScalar ivol, PetscScalar iarea, bool useHomoTemp, 
                    PetscScalar &gmeltfrac, PetscScalar &gtemp0, PetscScalar &gorigfrac);

  using IceModel::temperatureStep;
  virtual PetscErrorCode temperatureStep(PetscScalar* vertSacrCount, PetscScalar* bulgeCount);
  
  virtual PetscErrorCode enthalpyAndDrainageStep(PetscScalar* vertSacrCount, PetscScalar* bulgeCount);

  virtual PetscErrorCode drainageToBaseModelEnth(EnthalpyConverter &EC,
                PetscScalar L, PetscScalar omega_max,
                PetscScalar thickness, PetscScalar z, PetscScalar dz,
                PetscScalar &enthalpy, PetscScalar &Hmelt);

  using IceModel::updateYieldStressFromHmelt;
  virtual PetscErrorCode updateYieldStressFromHmelt();

  using IceModel::getEffectivePressureOnTill;  // but does not override it; one more arg
  virtual PetscScalar getEffectivePressureOnTill(PetscScalar thk, PetscScalar bwat, PetscScalar bmr,
						 PetscScalar till_pw_frac, PetscScalar max_hmelt) const;

protected: // new data members
  IceModelVec3  EnthNew3;  // NOTE:  Enth3 is an IceModel member, uninitialized and unused within IceModel
  
  PetscTruth    bmr_in_pore_pressure, thk_affects_pore_pressure;
  PetscScalar   bmr_enhance_scale, margin_pore_pressure_reduced,
                margin_pore_pressure_H_high, margin_pore_pressure_H_low;
};

#endif

