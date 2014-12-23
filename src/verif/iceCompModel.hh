// Copyright (C) 2004-2014 Jed Brown, Ed Bueler and Constantine Khroulev
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

#ifndef __iceCompModel_hh
#define __iceCompModel_hh

#include "iceModel.hh"
#include "bedrockThermalUnit.hh"

namespace pism {

class BTU_Verification : public BedThermalUnit
{
public:
  BTU_Verification(const IceGrid &g, int test, int bii)
    : BedThermalUnit(g) { testname = test; bedrock_is_ice = bii; }
  virtual ~BTU_Verification() {}

  virtual void get_temp(IceModelVec3Custom* &result);
protected:
  virtual void bootstrap();
  int testname, bedrock_is_ice;
};

class IceCompModel : public IceModel {

public:
  IceCompModel(IceGrid &g, Config &config, Config &overrides, int mytest);
  virtual ~IceCompModel() {}
  
  // re-defined steps of init() sequence:
  virtual void set_grid_defaults();     // called by IceModel::grid_setup()
  virtual void setFromOptions();
  virtual void createVecs();
  virtual void allocate_stressbalance();
  virtual void allocate_bedrock_thermal_unit();
  virtual void allocate_bed_deformation();
  virtual void allocate_enthalpy_converter();
  virtual void allocate_couplers();
  virtual void set_vars_from_options(); // called by IceModel::model_state_setup()

  void reportErrors();

protected:
  // related to all (or most) tests
  PetscBool   exactOnly;
  int          testname;
  virtual void additionalAtStartTimestep();
  virtual void additionalAtEndTimestep();
  // all tests except K
  void computeGeometryErrors(double &gvolexact, double &gareaexact, double &gdomeHexact,
                                       double &volerr, double &areaerr,
                                       double &gmaxHerr, double &gavHerr, double &gmaxetaerr,
                                       double &centerHerr);
  virtual void summary(bool tempAndAge);

  // related to tests A B C D E H
  void initTestABCDEH();
  void fillSolnTestABCDH();  // only used with exactOnly == PETSC_TRUE
  
  // related to test E
  void fillSolnTestE();  // only used with exactOnly == PETSC_TRUE

  // test E only
  void computeBasalVelocityErrors(double &exactmaxspeed,
                                            double &gmaxvecerr, double &gavvecerr,
                                            double &gmaxuberr, double &gmaxvberr);

  void reset_thickness_tests_AE();

  // related to test L
  IceModelVec2S   vHexactL;
  void initTestL();
  void fillSolnTestL();  // only used with exactOnly == PETSC_TRUE

  // related to tests F G; see iCMthermo.cc
  virtual void temperatureStep(double* vertSacrCount, double* bulgeCount);
  void initTestFG();
  void getCompSourcesTestFG();
  void fillSolnTestFG();  // only used with exactOnly == PETSC_TRUE
  // tests F and G
  void computeTemperatureErrors(double &gmaxTerr, double &gavTerr);
  // tests F and G
  void computeBasalTemperatureErrors(double &gmaxTerr, double &gavTerr, double &centerTerr);
  // tests F and G
  void compute_strain_heating_errors(double &gmax_strain_heating_err, double &gav_strain_heating_err);

  // tests F and G
  void computeSurfaceVelocityErrors(double &gmaxUerr, double &gavUerr,  // 2D vector errors
                                              double &gmaxWerr, double &gavWerr); // scalar errors
  
  IceModelVec3   strain_heating3_comp;

  // related to tests K and O; see iCMthermo.cc
  void initTestsKO();
  void fillTemperatureSolnTestsKO();  // used in initialzation
  //   and with exactOnly == PETSC_TRUE
  void fillBasalMeltRateSolnTestO();  // used only with exactOnly == PETSC_TRUE
 // tests K and O only
  void computeIceBedrockTemperatureErrors(double &gmaxTerr, double &gavTerr,
                                                    double &gmaxTberr, double &gavTberr);
  // test O only
  void computeBasalMeltRateErrors(double &gmaxbmelterr, double &gminbmelterr);

  // using Van der Veen's exact solution to test CFBC and the part-grid code
  void test_V_init();

  static const double secpera;

private:
  double        f;       // ratio of ice density to bedrock density
  PetscBool         bedrock_is_ice_forK;

  // see iCMthermo.cc
  static const double Ggeo;    // J/m^2 s; geothermal heat flux, assumed constant
  static const double ST;      // K m^-1;  surface temperature gradient: T_s = ST * r + Tmin
  static const double Tmin;    // K;       minimum temperature (at center)
  static const double LforFG;  // m;  exact radius of tests F&G ice sheet
  static const double ApforG;  // m;  magnitude A_p of annular perturbation for test G;
  // period t_p is set internally to 2000 years
};

} // end of namespace pism

#endif /* __iceCompModel_hh */
