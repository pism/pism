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

#include "pism/stressbalance/blatter/verification/BlatterTestCase.hh"
#include "pism/stressbalance/StressBalance.hh"
#include "pism/util/Context.hh"
#include "pism/util/Logger.hh"
#include "pism/util/Time.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/io/SynchronousOutputWriter.hh"
#include "pism/util/io/IO_Flags.hh"

namespace pism {
namespace stressbalance {

BlatterTestCase::BlatterTestCase(std::shared_ptr<Blatter> solver)
  : m_grid(solver->grid()),
    m_ctx(m_grid->ctx()),
    m_config(m_ctx->config()),
    m_sys(m_ctx->unit_system()),
    m_geometry(m_grid),
    m_tauc(m_grid, "tauc"),
    m_ice_enthalpy(m_grid, "enthalpy", array::WITHOUT_GHOSTS, m_grid->z()),
    m_solver(solver) {

  m_tauc.metadata(0)
      .long_name("yield stress for basal till (plastic or pseudo-plastic model)")
      .units("Pa");

  m_ice_enthalpy.metadata(0)
      .long_name("ice enthalpy (includes sensible heat, latent heat, pressure)")
      .units("J kg^-1");
}

void BlatterTestCase::init() {
  m_solver->init();
  initializeGeometry();
}

void BlatterTestCase::run() {
  m_ctx->log()->message(2, "* Solving the Blatter stress balance ...\n");

  m_geometry.ensure_consistency(m_config->get_number("stress_balance.ice_free_thickness_standard"));

  Inputs inputs;
  inputs.water_column_pressure = nullptr;
  inputs.geometry              = &m_geometry;
  inputs.enthalpy              = &m_ice_enthalpy;
  inputs.basal_yield_stress    = &m_tauc;
  inputs.bc_mask               = nullptr;
  inputs.bc_values             = nullptr;

  m_solver->update(inputs, true);
}

void BlatterTestCase::report(const std::string &testname) {

  auto u_sigma = m_solver->velocity_u_sigma();
  auto v_sigma = m_solver->velocity_v_sigma();

  const auto &Z = u_sigma->levels();
  int Mz = (int)Z.size();

  double maxuerr = 0.0, maxverr = 0.0, maxvecerr = 0.0;
  double avuerr = 0.0, avverr = 0.0, avvecerr = 0.0;
  double exactvelmax = 0.0;

  array::AccessScope list{u_sigma.get(), v_sigma.get(),
                          &m_geometry.bed_elevation, &m_geometry.ice_thickness};

  unsigned int N = 0;

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();

    double x = m_grid->x(i);
    double y = m_grid->y(j);
    double bed = m_geometry.bed_elevation(i, j);
    double H   = m_geometry.ice_thickness(i, j);

    const double *u_col = u_sigma->get_column(i, j);
    const double *v_col = v_sigma->get_column(i, j);

    for (int k = 0; k < Mz; ++k) {
      double z = bed + H * Z[k];

      Vector2d exact = exactSolution(i, j, x, y, z);

      double uerr = fabs(u_col[k] - exact.u);
      double verr = fabs(v_col[k] - exact.v);
      double vecerr = sqrt(uerr * uerr + verr * verr);
      double exactnorm = sqrt(exact.u * exact.u + exact.v * exact.v);

      exactvelmax = std::max(exactvelmax, exactnorm);
      maxuerr     = std::max(maxuerr, uerr);
      maxverr     = std::max(maxverr, verr);
      maxvecerr   = std::max(maxvecerr, vecerr);
      avuerr     += uerr;
      avverr     += verr;
      avvecerr   += vecerr;
      ++N;
    }
  }

  double gexactvelmax = GlobalMax(m_grid->com, exactvelmax);
  double gmaxuerr     = GlobalMax(m_grid->com, maxuerr);
  double gmaxverr     = GlobalMax(m_grid->com, maxverr);
  double gmaxvecerr   = GlobalMax(m_grid->com, maxvecerr);

  unsigned int gN = 0;
  {
    double local_N = (double)N;
    double global_N = GlobalSum(m_grid->com, local_N);
    gN = (unsigned int)global_N;
  }

  double gavuerr   = GlobalSum(m_grid->com, avuerr)   / gN;
  double gavverr   = GlobalSum(m_grid->com, avverr)   / gN;
  double gavvecerr = GlobalSum(m_grid->com, avvecerr) / gN;

  using pism::units::convert;

  m_ctx->log()->message(1,
      "NUMERICAL ERRORS in Blatter velocity (%s) relative to exact solution:\n",
      testname.c_str());
  m_ctx->log()->message(1,
      "velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv\n");
  m_ctx->log()->message(1,
      "           %11.4f%13.5f%10.4f%10.4f%10.4f%10.4f\n",
      convert(m_sys, gmaxvecerr, "m second-1", "m year-1"),
      (gexactvelmax > 0.0 ? (gavvecerr / gexactvelmax) * 100.0 : 0.0),
      convert(m_sys, gmaxuerr, "m second-1", "m year-1"),
      convert(m_sys, gmaxverr, "m second-1", "m year-1"),
      convert(m_sys, gavuerr, "m second-1", "m year-1"),
      convert(m_sys, gavverr, "m second-1", "m year-1"));
  m_ctx->log()->message(1, "NUM ERRORS DONE\n");
}

Vector2d BlatterTestCase::exactSolution(int /*i*/, int /*j*/, double /*x*/, double /*y*/,
                                        double /*z*/) const {
  return {0.0, 0.0};
}

void BlatterTestCase::write(const std::string &filename) {
  auto writer = std::make_shared<SynchronousOutputWriter>(m_grid->com, *m_config);
  writer->initialize({}, true);

  OutputFile file(writer, filename);

  file.define_variable(m_ctx->time()->metadata());
  file.append_time(0.0);

  // Write geometry and 2D averaged velocity (same pattern as SSATestCase::write)
  const array::Array *variables[] = {
    &m_geometry.ice_surface_elevation, &m_geometry.ice_thickness,
    &m_geometry.bed_elevation, &m_tauc, &m_ice_enthalpy,
    &m_solver->velocity()
  };

  // define
  for (const auto *v : variables) {
    for (const auto &m : v->all_metadata()) {
      file.define_variable(m);
    }
  }

  // write
  for (const auto *v : variables) {
    v->write(file);
  }

  // Write solver diagnostics (including 3D sigma velocity fields)
  {
    std::vector<std::shared_ptr<array::Array>> diags;
    for (auto &pair : m_solver->diagnostics()) {
      try {
        diags.push_back(pair.second->compute());
      } catch (RuntimeError &) {
        // ignore: some diagnostics may not be available
      }
    }

    for (auto &v : diags) {
      for (auto &m : v->all_metadata()) {
        file.define_variable(m);
      }
    }
    for (auto &v : diags) {
      v->write(file);
    }
  }

  file.close();
}

} // end of namespace stressbalance
} // end of namespace pism
