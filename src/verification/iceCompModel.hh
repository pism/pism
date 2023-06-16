// Copyright (C) 2004-2017, 2023 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "pism/icemodel/IceModel.hh"

namespace pism {

class IceCompModel : public IceModel {

public:
  IceCompModel(std::shared_ptr<IceGrid> g, std::shared_ptr<Context> ctx, int mytest);
  virtual ~IceCompModel() {}
  
  // re-defined steps of init() sequence:
  virtual void allocate_storage();
  virtual void allocate_bedrock_thermal_unit();
  virtual void allocate_bed_deformation();
  virtual void allocate_couplers();
  virtual void allocate_energy_model();

  // NB! not virtual
  void bootstrap_2d(const File &input_file) __attribute__((noreturn));

  virtual void initialize_2d();

  void reportErrors();

protected:
  // related to all (or most) tests
  int m_testname;

  virtual void post_step_hook();
  // all tests except K
  void computeGeometryErrors(double &gvolexact, double &gareaexact, double &gdomeHexact,
                                       double &volerr, double &areaerr,
                                       double &gmaxHerr, double &gavHerr, double &gmaxetaerr,
                                       double &centerHerr);
  virtual void print_summary(bool tempAndAge);

  // related to tests A B C D H
  void initTestABCDH();

  void reset_thickness_test_A();

  // related to test L
  array::Scalar m_HexactL;
  void initTestL();

  // related to tests F G; see iCMthermo.cc
  virtual void energy_step();
  void initTestFG();
  void getCompSourcesTestFG();

  // tests F and G
  void computeTemperatureErrors(double &gmaxTerr, double &gavTerr);
  // tests F and G
  void computeBasalTemperatureErrors(double &gmaxTerr, double &gavTerr, double &centerTerr);
  // tests F and G
  void compute_strain_heating_errors(double &gmax_strain_heating_err, double &gav_strain_heating_err);

  // tests F and G
  void computeSurfaceVelocityErrors(double &gmaxUerr, double &gavUerr,  // 2D vector errors
                                              double &gmaxWerr, double &gavWerr); // scalar errors
  
  array::Array3D m_strain_heating3_comp;

  // related to tests K and O; see iCMthermo.cc
  void initTestsKO();

 // tests K and O only
  void computeIceBedrockTemperatureErrors(double &gmaxTerr, double &gavTerr,
                                                    double &gmaxTberr, double &gavTberr);
  // test O only
  void computeBasalMeltRateErrors(double &gmaxbmelterr, double &gminbmelterr);

  // using Van der Veen's exact solution to test CFBC and the part-grid code
  void test_V_init();

private:
  double m_f;       // ratio of ice density to bedrock density
  bool m_bedrock_is_ice_forK;

  // see iCMthermo.cc
  static const double m_ST;      // K m^-1;  surface temperature gradient: T_s = ST * r + Tmin
  static const double m_Tmin;    // K;       minimum temperature (at center)
  static const double m_LforFG;  // m;  exact radius of tests F&G ice sheet
  static const double m_ApforG;  // m;  magnitude A_p of annular perturbation for test G;
  // period t_p is set internally to 2000 years
};

} // end of namespace pism

#endif /* __iceCompModel_hh */
