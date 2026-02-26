// Copyright (C) 2020--2026 PISM Authors
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

static char help[] =
  "\nPISM_BLATTER_TEST\n"
  "  Standalone verification test program for PISM's Blatter-Pattyn solver.\n"
  "  Supports tests: xy, xz, cfbc, halfar, vanderveen.\n\n";

#include <petsc.h>

#include "pism/stressbalance/blatter/verification/BlatterTestCase.hh"
#include "pism/stressbalance/blatter/verification/BlatterTestXY.hh"
#include "pism/stressbalance/blatter/verification/BlatterTestXZ.hh"
#include "pism/stressbalance/blatter/verification/BlatterTestCFBC.hh"
#include "pism/stressbalance/blatter/verification/BlatterTestHalfar.hh"
#include "pism/stressbalance/blatter/verification/BlatterTestvanderVeen.hh"
#include "pism/stressbalance/blatter/verification/manufactured_solutions.hh"
#include "pism/util/Context.hh"
#include "pism/util/Grid.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/petscwrappers/PetscInitializer.hh"
#include "pism/util/pism_options.hh"

namespace pism {
namespace stressbalance {

// ============================================================
// Test XY
// ============================================================
class BlatterTestCaseXY : public BlatterTestCase {
public:
  BlatterTestCaseXY(std::shared_ptr<Blatter> solver)
    : BlatterTestCase(solver) {
    m_ice_enthalpy.set(1e5);
  }
protected:
  void initializeGeometry() override {
    m_tauc.set(0.0);
    m_geometry.bed_elevation.set(0.0);
    m_geometry.ice_thickness.set(1.0);
    m_geometry.sea_level_elevation.set(0.0);
    m_geometry.ensure_consistency(0.0);
  }

  Vector2d exactSolution(int /*i*/, int /*j*/, double x, double y,
                         double /*z*/) const override {
    return blatter_xy_exact(x, y);
  }
};

// ============================================================
// Test XZ
// ============================================================
class BlatterTestCaseXZ : public BlatterTestCase {
public:
  BlatterTestCaseXZ(std::shared_ptr<Blatter> solver)
    : BlatterTestCase(solver) {
    m_ice_enthalpy.set(1e5);

    m_s0    = 2000.0;
    m_alpha = 4e-8;
    m_H0    = 1000.0;
    m_A     = m_config->get_number("flow_law.isothermal_Glen.ice_softness");
    m_beta  = units::convert(m_sys, 1.0, "kPa year m-1", "Pa s m-1");
    m_rho   = m_config->get_number("constants.ice.density");
    m_g     = m_config->get_number("constants.standard_gravity");
  }
protected:
  void initializeGeometry() override {
    m_tauc.set(m_beta);

    {
      array::AccessScope list{&m_geometry.bed_elevation};
      for (auto p : m_grid->points()) {
        const int i = p.i(), j = p.j();
        double x = m_grid->x(i);
        m_geometry.bed_elevation(i, j) = m_s0 - m_H0 - m_alpha * x * x;
      }
    }

    m_geometry.ice_thickness.set(m_H0);

    // ensure all ice is grounded
    m_geometry.sea_level_elevation.copy_from(m_geometry.bed_elevation);
    m_geometry.sea_level_elevation.shift(-1.0);

    m_geometry.ensure_consistency(0.0);
  }

  Vector2d exactSolution(int /*i*/, int /*j*/, double x, double /*y*/,
                         double z) const override {
    return blatter_xz_exact(x, z, m_A, m_rho, m_g, m_s0, m_alpha, m_H0, m_beta);
  }

  double m_s0, m_alpha, m_H0, m_A, m_beta, m_rho, m_g;
};

// ============================================================
// Test CFBC
// ============================================================
class BlatterTestCaseCFBC : public BlatterTestCase {
public:
  BlatterTestCaseCFBC(std::shared_ptr<Blatter> solver)
    : BlatterTestCase(solver) {
    m_ice_enthalpy.set(1e5);

    m_H = 1e3;
    m_L = 2.0 * m_grid->Lx();
    double A = m_config->get_number("flow_law.isothermal_Glen.ice_softness");
    m_B = 1.0 / A;  // n == 1
    m_rho_i = m_config->get_number("constants.ice.density");
    m_rho_w = m_config->get_number("constants.sea_water.density");
    m_g     = m_config->get_number("constants.standard_gravity");
  }
protected:
  void initializeGeometry() override {
    m_tauc.set(0.0);

    m_geometry.bed_elevation.set(-m_H);
    m_geometry.ice_thickness.set(m_H);
    m_geometry.ice_surface_elevation.set(0.0);
    m_geometry.cell_type.set(MASK_FLOATING);
    m_geometry.sea_level_elevation.set(0.0);

    // do *not* call ensure_consistency: we want surface elevation at sea level
  }

  Vector2d exactSolution(int /*i*/, int /*j*/, double x, double /*y*/,
                         double z) const override {
    return blatter_xz_cfbc_exact(x, z, m_B, m_L, m_rho_i, m_rho_w, m_g);
  }

  double m_H, m_L, m_B, m_rho_i, m_rho_w, m_g;
};

// ============================================================
// Test Halfar
// ============================================================
class BlatterTestCaseHalfar : public BlatterTestCase {
public:
  BlatterTestCaseHalfar(std::shared_ptr<BlatterTestHalfar> solver)
    : BlatterTestCase(solver), m_halfar(solver) {
    m_ice_enthalpy.set(0.0);
  }
protected:
  void initializeGeometry() override {
    m_tauc.set(0.0);

    {
      array::AccessScope list{&m_geometry.ice_thickness};
      for (auto p : m_grid->points()) {
        const int i = p.i(), j = p.j();
        m_geometry.ice_thickness(i, j) = m_halfar->H_exact(m_grid->x(i));
      }
    }

    m_geometry.bed_elevation.set(0.0);
    m_geometry.sea_level_elevation.set(0.0);
    m_geometry.ensure_consistency(0.0);
  }

  Vector2d exactSolution(int /*i*/, int /*j*/, double x, double /*y*/,
                         double z) const override {
    double u = m_halfar->u_exact(x, z);
    return {u, 0.0};
  }

  std::shared_ptr<BlatterTestHalfar> m_halfar;
};

// ============================================================
// Test van der Veen
// ============================================================
class BlatterTestCaseVanderVeen : public BlatterTestCase {
public:
  BlatterTestCaseVanderVeen(std::shared_ptr<BlatterTestvanderVeen> solver)
    : BlatterTestCase(solver), m_vdv(solver) {
    m_ice_enthalpy.set(0.0);
  }
protected:
  void initializeGeometry() override {
    {
      array::AccessScope list{&m_geometry.ice_thickness, &m_geometry.bed_elevation,
                              &m_tauc};
      for (auto p : m_grid->points()) {
        const int i = p.i(), j = p.j();
        double x = m_grid->x(i);
        m_geometry.ice_thickness(i, j) = m_vdv->H_exact(x);
        m_geometry.bed_elevation(i, j) = m_vdv->b_exact(x);
        m_tauc(i, j) = m_vdv->beta_exact(x);
      }
    }

    // low enough to make it grounded
    m_geometry.sea_level_elevation.set(-100.0);
    m_geometry.ensure_consistency(0.0);
  }

  Vector2d exactSolution(int /*i*/, int /*j*/, double x, double /*y*/,
                         double /*z*/) const override {
    return m_vdv->u_exact(x);
  }

  std::shared_ptr<BlatterTestvanderVeen> m_vdv;
};

} // end of namespace stressbalance
} // end of namespace pism

int main(int argc, char *argv[]) {

  using namespace pism;
  using namespace pism::stressbalance;

  MPI_Comm com = MPI_COMM_WORLD;
  petsc::Initializer petsc(argc, argv, help);

  com = PETSC_COMM_WORLD;

  try {
    std::shared_ptr<Context> ctx = context_from_options(com, "pism_blatter_test");
    auto config = ctx->config();
    auto log    = ctx->log();

    std::string usage =
      "usage of PISM_BLATTER_TEST:\n"
      "  run pism_blatter_test -test <testname> -Mx <number>\n"
      "\n"
      "where <testname> is one of: xy, xz, cfbc, halfar, vanderveen\n";

    bool stop = maybe_show_usage(*log, "pism_blatter_test", usage);

    if (stop) {
      return 0;
    }

    auto testname = options::Keyword("-test",
                                     "Blatter verification test name",
                                     "xy,xz,cfbc,halfar,vanderveen",
                                     "xy");

    unsigned int Mx = config->get_number("grid.Mx");

    auto output_file = config->get_string("output.file");
    bool write_output = config->get_string("output.size") != "none";

    // common config: do not ignore thin ice
    config->set_number("geometry.ice_free_thickness_standard", 0.0);
    config->set_number("stress_balance.ice_free_thickness_standard", 0.0);

    // common config: isothermal Glen flow law
    config->set_string("stress_balance.blatter.flow_law", "isothermal_glen");

    // common config: physical constants
    config->set_number("constants.ice.density", 910.0);
    config->set_number("constants.standard_gravity", 9.81);

    std::string test = testname.value();

    if (test == "xy") {
      // Glen exponent n = 3, ice softness A = 1e-4 Pa^-3 yr^-1
      config->set_number("stress_balance.blatter.Glen_exponent", 3.0);
      config->set_number("flow_law.isothermal_Glen.ice_softness",
                         units::convert(ctx->unit_system(), 1e-4, "Pa-3 year-1", "Pa-3 s-1"));

      // Domain: [0,1] * [0,1] * [0,1]
      unsigned int My = Mx;
      double Lx = 0.5, Ly = 0.5;

      grid::Parameters P(*config, Mx, My, Lx, Ly);
      P.x0 = 0.5;
      P.y0 = 0.5;
      P.z = {0.0, 1.0};
      P.registration = grid::CELL_CORNER;
      P.periodicity  = grid::NOT_PERIODIC;
      P.ownership_ranges_from_options(*config, ctx->size());

      auto grid = std::make_shared<Grid>(ctx, P);

      int blatter_Mz = 2;
      int coarsening_factor = 1;
      auto solver = std::make_shared<BlatterTestXY>(grid, blatter_Mz, coarsening_factor);

      BlatterTestCaseXY testcase(solver);
      testcase.init();
      testcase.run();
      testcase.report("xy");
      if (write_output) {
        testcase.write(output_file);
      }

    } else if (test == "xz") {
      // Glen exponent n = 3, ice softness A = 1e-16 Pa^-3 yr^-1
      config->set_number("stress_balance.blatter.Glen_exponent", 3.0);
      config->set_number("flow_law.isothermal_Glen.ice_softness",
                         units::convert(ctx->unit_system(), 1e-16, "Pa-3 year-1", "Pa-3 s-1"));
      config->set_number("flow_law.Schoof_regularizing_velocity", 1e-5);

      // Set sliding law parameters to make "tauc" equivalent to "beta"
      config->set_flag("basal_resistance.pseudo_plastic.enabled", true);
      config->set_number("basal_resistance.pseudo_plastic.q", 1.0);
      config->set_number("basal_resistance.pseudo_plastic.u_threshold",
                         units::convert(ctx->unit_system(), 1.0, "m / s", "m / year"));

      double Lx = 50e3;
      double H  = 1000.0;

      double dx = (2 * Lx) / (Mx - 1);
      double Ly = dx;

      grid::Parameters P(*config, Mx, 3, Lx, Ly);
      P.x0 = 0.0;
      P.y0 = 0.0;
      P.z = {0.0, H};
      P.registration = grid::CELL_CORNER;
      P.periodicity  = grid::Y_PERIODIC;
      P.ownership_ranges_from_options(*config, ctx->size());

      auto grid = std::make_shared<Grid>(ctx, P);

      int blatter_Mz = Mx;
      int coarsening_factor = 4;
      auto solver = std::make_shared<BlatterTestXZ>(grid, blatter_Mz, coarsening_factor);

      BlatterTestCaseXZ testcase(solver);
      testcase.init();
      testcase.run();
      testcase.report("xz");
      if (write_output) {
        testcase.write(output_file);
      }

    } else if (test == "cfbc") {
      // Glen exponent n = 1 (constant viscosity)
      config->set_number("stress_balance.blatter.Glen_exponent", 1.0);
      config->set_number("flow_law.isothermal_Glen.ice_softness",
                         units::convert(ctx->unit_system(), 1e-3, "Pa-3 year-1", "Pa-3 s-1"));

      double H  = 1e3;
      double L  = 1.0;
      double Lx = 0.5 * L;

      double dx = (2 * Lx) / (Mx - 1);
      double Ly = dx;

      grid::Parameters P(*config, Mx, Mx, Lx, Ly);
      P.x0 = Lx;
      P.y0 = 0.0;
      P.z = {0.0, H};
      P.registration = grid::CELL_CORNER;
      P.periodicity  = grid::Y_PERIODIC;
      P.ownership_ranges_from_options(*config, ctx->size());

      auto grid = std::make_shared<Grid>(ctx, P);

      int blatter_Mz = Mx;
      int coarsening_factor = 1;
      auto solver = std::make_shared<BlatterTestCFBC>(grid, blatter_Mz, coarsening_factor);

      BlatterTestCaseCFBC testcase(solver);
      testcase.init();
      testcase.run();
      testcase.report("cfbc");
      if (write_output) {
        testcase.write(output_file);
      }

    } else if (test == "halfar") {
      config->set_number("stress_balance.blatter.Glen_exponent", 3.0);

      double R0 = 750e3;
      double H0 = 3600.0;

      double padding = 1e3;
      double Lx = R0 / 2 - padding;

      double dx = (2 * Lx) / (Mx - 1);
      double Ly = dx;

      grid::Parameters P(*config, Mx, 3, Lx, Ly);
      P.x0 = Lx;
      P.y0 = 0.0;
      P.z = {0.0, H0};
      P.registration = grid::CELL_CORNER;
      P.periodicity  = grid::Y_PERIODIC;
      P.ownership_ranges_from_options(*config, ctx->size());

      auto grid = std::make_shared<Grid>(ctx, P);

      // Mz from config, default to reasonable value
      int blatter_Mz = config->get_number("stress_balance.blatter.Mz");
      int coarsening_factor = 2;
      auto solver = std::make_shared<BlatterTestHalfar>(grid, blatter_Mz, coarsening_factor);

      BlatterTestCaseHalfar testcase(solver);
      testcase.init();
      testcase.run();
      testcase.report("halfar");
      if (write_output) {
        testcase.write(output_file);
      }

    } else if (test == "vanderveen") {
      config->set_number("stress_balance.blatter.Glen_exponent", 3.0);

      // Set sliding law parameters to make "tauc" equivalent to "beta"
      config->set_flag("basal_resistance.pseudo_plastic.enabled", true);
      config->set_number("basal_resistance.pseudo_plastic.q", 1.0);
      config->set_number("basal_resistance.pseudo_plastic.u_threshold",
                         units::convert(ctx->unit_system(), 1.0, "m / s", "m / year"));

      double Lx = 1e5;
      double dx = (2 * Lx) / (Mx - 1);
      double Ly = dx;

      grid::Parameters P(*config, Mx, 3, Lx, Ly);
      P.x0 = 1.1 * Lx;
      P.y0 = 0.0;
      P.z = {0.0, 1000.0};
      P.registration = grid::CELL_CORNER;
      P.periodicity  = grid::Y_PERIODIC;
      P.ownership_ranges_from_options(*config, ctx->size());

      auto grid = std::make_shared<Grid>(ctx, P);

      int blatter_Mz = 2;
      int coarsening_factor = 1;
      auto solver = std::make_shared<BlatterTestvanderVeen>(grid, blatter_Mz, coarsening_factor);

      BlatterTestCaseVanderVeen testcase(solver);
      testcase.init();
      testcase.run();
      testcase.report("vanderveen");
      if (write_output) {
        testcase.write(output_file);
      }

    } else {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "unknown Blatter test: %s", test.c_str());
    }
  }
  catch (...) {
    handle_fatal_errors(com);
    return 1;
  }

  return 0;
}
