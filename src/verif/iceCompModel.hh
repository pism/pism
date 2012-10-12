// Copyright (C) 2004-2012 Jed Brown, Ed Bueler and Constantine Khroulev
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

#ifndef __iceCompModel_hh
#define __iceCompModel_hh

#include "iceModel.hh"
#include "bedrockThermalUnit.hh"

class BTU_Verification : public PISMBedThermalUnit
{
public:
  BTU_Verification(IceGrid &g, const NCConfigVariable &conf, int test, int bii)
    : PISMBedThermalUnit(g, conf) { testname = test; bedrock_is_ice = bii; }
  virtual ~BTU_Verification() {}

  virtual PetscErrorCode get_temp(IceModelVec3BTU* &result);
protected:
  virtual PetscErrorCode bootstrap();
  int testname, bedrock_is_ice;
};

class IceCompModel : public IceModel {

public:
  IceCompModel(IceGrid &g, NCConfigVariable &config, NCConfigVariable &overrides, int mytest);
  virtual ~IceCompModel() {}
  
  // re-defined steps of init() sequence:
  virtual PetscErrorCode set_grid_defaults();     // called by IceModel::grid_setup()
  virtual PetscErrorCode setFromOptions();
  virtual PetscErrorCode createVecs();
  virtual PetscErrorCode allocate_flowlaw();
  virtual PetscErrorCode allocate_stressbalance();
  virtual PetscErrorCode allocate_bedrock_thermal_unit();
  virtual PetscErrorCode allocate_bed_deformation();
  virtual PetscErrorCode allocate_enthalpy_converter();
  virtual PetscErrorCode set_vars_from_options(); // called by IceModel::model_state_setup()

  PetscErrorCode reportErrors();

protected:
  // related to all (or most) tests
  PetscBool   exactOnly;
  int          testname;
  virtual PetscErrorCode additionalAtStartTimestep();
  virtual PetscErrorCode additionalAtEndTimestep();
  PetscErrorCode computeGeometryErrors(    // all tests except K
        PetscScalar &gvolexact, PetscScalar &gareaexact, PetscScalar &gdomeHexact,
        PetscScalar &volerr, PetscScalar &areaerr,
        PetscScalar &gmaxHerr, PetscScalar &gavHerr, PetscScalar &gmaxetaerr,
        PetscScalar &centerHerr);
  virtual PetscErrorCode summary(bool tempAndAge);

  // related to tests A B C D E H
  PetscErrorCode initTestABCDEH();
  PetscErrorCode getCompSourcesTestCDH();  // only for time-dependent compensatory sources
  PetscErrorCode fillSolnTestABCDH();  // only used with exactOnly == PETSC_TRUE
  
  // related to test E
  PetscErrorCode fillSolnTestE();  // only used with exactOnly == PETSC_TRUE
  PetscErrorCode computeBasalVelocityErrors(    // test E only
        PetscScalar &exactmaxspeed,
        PetscScalar &gmaxvecerr, PetscScalar &gavvecerr,
        PetscScalar &gmaxuberr, PetscScalar &gmaxvberr);

  PetscErrorCode reset_thickness_tests_AE();

  // related to test L
  IceModelVec2S   vHexactL;
  PetscErrorCode initTestL();
  PetscErrorCode fillSolnTestL();  // only used with exactOnly == PETSC_TRUE

  // related to tests F G; see iCMthermo.cc
  virtual PetscErrorCode temperatureStep(PetscScalar* vertSacrCount, PetscScalar* bulgeCount);
  PetscErrorCode initTestFG();
  PetscErrorCode getCompSourcesTestFG();
  PetscErrorCode fillSolnTestFG();  // only used with exactOnly == PETSC_TRUE
  PetscErrorCode computeTemperatureErrors(      // tests F and G
                   PetscScalar &gmaxTerr, PetscScalar &gavTerr);
  PetscErrorCode computeBasalTemperatureErrors( // tests F and G
                   PetscScalar &gmaxTerr, PetscScalar &gavTerr, PetscScalar &centerTerr);
  PetscErrorCode computeSigmaErrors(            // tests F and G
                   PetscScalar &gmaxSigmaerr, PetscScalar &gavSigmaerr);
  PetscErrorCode computeSurfaceVelocityErrors(  // tests F and G
                   PetscScalar &gmaxUerr, PetscScalar &gavUerr,  // 2D vector errors
                   PetscScalar &gmaxWerr, PetscScalar &gavWerr); // scalar errors
  
  IceModelVec3   SigmaComp3;

  // related to tests K and O; see iCMthermo.cc
  PetscErrorCode initTestsKO();
  PetscErrorCode fillTemperatureSolnTestsKO();  // used in initialzation
                                                //   and with exactOnly == PETSC_TRUE
  PetscErrorCode fillBasalMeltRateSolnTestO();  // used only with exactOnly == PETSC_TRUE
  PetscErrorCode computeIceBedrockTemperatureErrors( // tests K and O only
                   PetscScalar &gmaxTerr, PetscScalar &gavTerr,
                   PetscScalar &gmaxTberr, PetscScalar &gavTberr);
  PetscErrorCode computeBasalMeltRateErrors( // test O only
                   PetscScalar &gmaxbmelterr, PetscScalar &gminbmelterr);

  // using Van der Veen's exact solution to test CFBC and the part-grid code
  PetscErrorCode test_V_init();

private:
  PetscScalar        f;       // ratio of ice density to bedrock density
  PetscBool         bedrock_is_ice_forK;

  static const PetscScalar ablationRateOutside;

  // see iCMthermo.cc
  static const PetscScalar Ggeo;    // J/m^2 s; geothermal heat flux, assumed constant
  static const PetscScalar ST;      // K m^-1;  surface temperature gradient: T_s = ST * r + Tmin
  static const PetscScalar Tmin;    // K;       minimum temperature (at center)
  static const PetscScalar LforFG;  // m;  exact radius of tests F&G ice sheet
  static const PetscScalar ApforG;  // m;  magnitude A_p of annular perturbation for test G;
  // period t_p is set internally to 2000 years
};

#endif /* __iceCompModel_hh */
