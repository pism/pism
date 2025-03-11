// Copyright (C) 2009--2025 Ed Bueler, Constantine Khroulev, and David Maxwell
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

#include "pism/stressbalance/ssa/tests/SSATestCase.hh"
#include "pism/stressbalance/ssa/SSAFD_SNES.hh"
#include "pism/stressbalance/StressBalance.hh"
#include "pism/util/Context.hh"
#include "pism/util/Interpolation1D.hh"
#include "pism/util/io/File.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/pism_utilities.hh"

#include "pism/stressbalance/ssa/SSAFD.hh"
#include "pism/stressbalance/ssa/SSAFEM.hh"

namespace pism {
namespace stressbalance {

SSATestCase::SSATestCase(std::shared_ptr<SSA> ssa)
  : m_grid(ssa->grid()),
    m_ctx(m_grid->ctx()),
    m_config(m_ctx->config()),
    m_sys(m_ctx->unit_system()),
    m_tauc(m_grid, "tauc"),
    m_ice_enthalpy(m_grid, "enthalpy", array::WITH_GHOSTS, m_grid->z(), 1),
    m_bc_values(m_grid, "_bc"), // u_bc and v_bc
    m_bc_mask(m_grid, "bc_mask"),
    m_geometry(m_grid),
    m_ssa(ssa) {

  m_bc_mask.set_interpolation_type(NEAREST);

  // yield stress for basal till (plastic or pseudo-plastic model)
  m_tauc.metadata(0)
      .long_name("yield stress for basal till (plastic or pseudo-plastic model)")
      .units("Pa");

  // enthalpy
  m_ice_enthalpy.metadata(0)
      .long_name("ice enthalpy (includes sensible heat, latent heat, pressure)")
      .units("J kg^-1");

  // dirichlet boundary condition (FIXME: perhaps unused!)
  m_bc_values.metadata(0)
      .long_name("X-component of the SSA velocity boundary conditions")
      .units("m s^-1")
      .output_units("m year^-1");
  m_bc_values.metadata(1)
      .long_name("Y-component of the SSA velocity boundary conditions")
      .units("m s^-1")
      .output_units("m year^-1");

  units::System::Ptr sys  = m_grid->ctx()->unit_system();
  double fill_value =
      units::convert(sys, m_config->get_number("output.fill_value"), "m year^-1", "m second^-1");

  auto large_number = units::convert(m_sys, 1e6, "m year^-1", "m second^-1");

  m_bc_values.metadata(0)["valid_range"] = { -large_number, large_number };
  m_bc_values.metadata(0)["_FillValue"]  = { fill_value };

  m_bc_values.metadata(1)["valid_range"] = { -large_number, large_number };
  m_bc_values.metadata(1)["_FillValue"]  = { fill_value };

  m_bc_values.set(fill_value);

  // Dirichlet B.C. mask
  m_bc_mask.metadata().long_name("grounded_dragging_floating integer mask");

  m_bc_mask.metadata()["flag_values"]   = { 0.0, 1.0 };
  m_bc_mask.metadata()["flag_meanings"] = "no_data ssa.dirichlet_bc_location";
}

std::shared_ptr<SSA> SSATestCase::solver(std::shared_ptr<Grid> grid, const std::string &method) {
  if (method == "fem") {
    return std::make_shared<SSAFEM>(grid);
  }
  if (method == "fd") {
    return std::make_shared<SSAFD>(grid, false);
  }
  return  std::make_shared<SSAFD_SNES>(grid, false);
}

//! Initialize the test case at the start of a run
void SSATestCase::init() {

  m_ssa->init();

  // Allow the subclass to set the coefficients.
  initializeSSACoefficients();
}

//! Solve the SSA
void SSATestCase::run() {
  m_ctx->log()->message(2, "* Solving the SSA stress balance ...\n");

  m_geometry.ensure_consistency(m_config->get_number("stress_balance.ice_free_thickness_standard"));

  Inputs inputs;
  inputs.water_column_pressure = nullptr;
  inputs.geometry              = &m_geometry;
  inputs.enthalpy              = &m_ice_enthalpy;
  inputs.basal_yield_stress    = &m_tauc;
  inputs.bc_mask               = &m_bc_mask;
  inputs.bc_values             = &m_bc_values;

  bool full_update = true;
  m_ssa->update(inputs, full_update);
}

//! Report on the generated solution
void SSATestCase::report(const std::string &testname) {

  m_ctx->log()->message(3, m_ssa->stdout_report());

  double maxvecerr = 0.0, avvecerr = 0.0, avuerr = 0.0, avverr = 0.0, maxuerr = 0.0, maxverr = 0.0;
  double gmaxvecerr = 0.0, gavvecerr = 0.0, gavuerr = 0.0, gavverr = 0.0, gmaxuerr = 0.0,
         gmaxverr = 0.0;

  if (m_config->get_flag("basal_resistance.pseudo_plastic.enabled") &&
      m_config->get_number("basal_resistance.pseudo_plastic.q") != 1.0) {
    m_ctx->log()->message(1, "WARNING: numerical errors not valid for pseudo-plastic till\n");
  }
  m_ctx->log()->message(1, "NUMERICAL ERRORS in velocity relative to exact solution:\n");

  const array::Vector &vel_ssa = m_ssa->velocity();

  array::AccessScope list{ &vel_ssa };

  double exactvelmax = 0, gexactvelmax = 0;
  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    double uexact, vexact;
    double myx = m_grid->x(i), myy = m_grid->y(j);

    exactSolution(i, j, myx, myy, &uexact, &vexact);

    double exactnormsq = sqrt(uexact * uexact + vexact * vexact);
    exactvelmax        = std::max(exactnormsq, exactvelmax);

    // compute maximum errors
    const double uerr   = fabs(vel_ssa(i, j).u - uexact);
    const double verr   = fabs(vel_ssa(i, j).v - vexact);
    avuerr              = avuerr + uerr;
    avverr              = avverr + verr;
    maxuerr             = std::max(maxuerr, uerr);
    maxverr             = std::max(maxverr, verr);
    const double vecerr = sqrt(uerr * uerr + verr * verr);
    maxvecerr           = std::max(maxvecerr, vecerr);
    avvecerr            = avvecerr + vecerr;
  }

  unsigned int N = (m_grid->Mx() * m_grid->My());

  gexactvelmax = GlobalMax(m_grid->com, exactvelmax);
  gmaxuerr     = GlobalMax(m_grid->com, maxuerr);
  gmaxverr     = GlobalMax(m_grid->com, maxverr);
  gavuerr      = GlobalSum(m_grid->com, avuerr);
  gavuerr      = gavuerr / N;
  gavverr      = GlobalSum(m_grid->com, avverr);
  gavverr      = gavverr / N;
  gmaxvecerr   = GlobalMax(m_grid->com, maxvecerr);
  gavvecerr    = GlobalSum(m_grid->com, avvecerr);
  gavvecerr    = gavvecerr / N;

  using pism::units::convert;

  m_ctx->log()->message(
      1, "velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv\n");
  m_ctx->log()->message(1, "           %11.4f%13.5f%10.4f%10.4f%10.4f%10.4f\n",
                        convert(m_sys, gmaxvecerr, "m second-1", "m year-1"),
                        (gavvecerr / gexactvelmax) * 100.0,
                        convert(m_sys, gmaxuerr, "m second-1", "m year-1"),
                        convert(m_sys, gmaxverr, "m second-1", "m year-1"),
                        convert(m_sys, gavuerr, "m second-1", "m year-1"),
                        convert(m_sys, gavverr, "m second-1", "m year-1"));

  m_ctx->log()->message(1, "NUM ERRORS DONE\n");

  report_netcdf(testname, convert(m_sys, gmaxvecerr, "m second-1", "m year-1"),
                (gavvecerr / gexactvelmax) * 100.0,
                convert(m_sys, gmaxuerr, "m second-1", "m year-1"),
                convert(m_sys, gmaxverr, "m second-1", "m year-1"),
                convert(m_sys, gavuerr, "m second-1", "m year-1"),
                convert(m_sys, gavverr, "m second-1", "m year-1"));
}

void SSATestCase::report_netcdf(const std::string &testname, double max_vector, double rel_vector,
                                double max_u, double max_v, double avg_u, double avg_v) {
  auto sys = m_grid->ctx()->unit_system();

  VariableMetadata global_attributes("PISM_GLOBAL", sys);

  options::String filename("-report_file", "NetCDF error report file");

  if (not filename.is_set()) {
    return;
  }

  m_ctx->log()->message(2, "Also writing errors to '%s'...\n", filename->c_str());

  bool append = options::Bool("-append", "Append the NetCDF error report");

  io::Mode mode = io::PISM_READWRITE;
  if (not append) {
    mode = io::PISM_READWRITE_MOVE;
  }

  global_attributes["source"] = std::string("PISM ") + pism::revision;

  // Find the number of records in this file:
  File file(m_grid->com, filename, io::PISM_NETCDF3, mode); // OK to use NetCDF3.
  size_t start = static_cast<size_t>(file.dimension_length("N"));

  io::write_attributes(file, global_attributes, io::PISM_DOUBLE);

  {
    VariableMetadata err("N", sys);
    io::define_timeseries(err, "N", file, io::PISM_DOUBLE);
    io::write_timeseries(file, err, start, { (double)(start + 1) });
  }
  {
    VariableMetadata dx{"dx", sys};
    dx.units("meters");
    io::define_timeseries(dx, "N", file, io::PISM_DOUBLE);
    io::write_timeseries(file, dx, start, { m_grid->dx() });
  }
  {
    VariableMetadata dy{"dy", sys};
    dy.units("meters");
    io::define_timeseries(dy, "N", file, io::PISM_DOUBLE);
    io::write_timeseries(file, dy, start, { m_grid->dy() });
  }
  {
    VariableMetadata test{"test", sys};
    test.units("1");
    io::define_timeseries(test, "N", file, io::PISM_INT);
    io::write_timeseries(file, test, start, { (double)testname[0] });
  }
  {
    VariableMetadata max_velocity{ "max_velocity", sys };
    max_velocity.long_name("maximum ice velocity magnitude error").units("m year^-1");
    io::define_timeseries(max_velocity, "N", file, io::PISM_DOUBLE);
    io::write_timeseries(file, max_velocity, start, { max_vector });
  }
  {
    VariableMetadata rel_velocity{ "relative_velocity", sys };
    rel_velocity.long_name("relative ice velocity magnitude error").units("percent");
    io::define_timeseries(rel_velocity, "N", file, io::PISM_DOUBLE);
    io::write_timeseries(file, rel_velocity, start, { rel_vector });
  }
  {
    VariableMetadata maximum_u{ "maximum_u", sys };
    maximum_u.long_name("maximum error in the X-component of the ice velocity").units("m year^-1");
    io::define_timeseries(maximum_u, "N", file, io::PISM_DOUBLE);
    io::write_timeseries(file, maximum_u, start, { max_u });
  }
  {
    VariableMetadata maximum_v{ "maximum_v", sys };
    maximum_v.long_name("maximum error in the Y-component of the ice velocity").units("m year^-1");
    io::define_timeseries(maximum_v, "N", file, io::PISM_DOUBLE);
    io::write_timeseries(file, maximum_v, start, { max_v });
  }
  {
    VariableMetadata average_u{ "average_u", sys };
    average_u.long_name("average error in the X-component of the ice velocity").units("m year^-1");
    io::define_timeseries(average_u, "N", file, io::PISM_DOUBLE);
    io::write_timeseries(file, average_u, start, { avg_u });
  }
  {
    VariableMetadata average_v{ "average_v", sys };
    average_v.long_name("average error in the Y-component of the ice velocity").units("m year^-1");
    io::define_timeseries(average_v, "N", file, io::PISM_DOUBLE);
    io::write_timeseries(file, average_v, start, { avg_v });
  }
  file.close();
}

void SSATestCase::exactSolution(int /*i*/, int /*j*/, double /*x*/, double /*y*/, double *u,
                                double *v) {
  *u = 0;
  *v = 0;
}

//! Save the computation and data to a file.
void SSATestCase::write(const std::string &filename) {

  // Write results to an output file:
  File file(m_grid->com, filename, io::PISM_NETCDF3, io::PISM_READWRITE_MOVE);
  io::define_time(file, *m_grid->ctx());
  io::append_time(file, *m_config, 0.0);

  m_geometry.ice_surface_elevation.write(file);
  m_geometry.ice_thickness.write(file);
  m_bc_mask.write(file);
  m_tauc.write(file);
  m_geometry.bed_elevation.write(file);
  m_ice_enthalpy.write(file);
  m_bc_values.write(file);

  m_ssa->velocity().write(file);

  // write all diagnostics:
  {
    auto diagnostics = m_ssa->diagnostics();

    for (auto &p : diagnostics) {
      try {
        p.second->compute()->write(file);
      } catch (RuntimeError &e) {
        // ignore errors
      }
    }
  }

  array::Vector tmp(m_grid, "_exact");
  tmp.metadata(0)
      .long_name("X-component of the SSA exact solution")
      .units("m s^-1");
  tmp.metadata(1)
      .long_name("Y-component of the SSA exact solution")
      .units("m s^-1");

  array::AccessScope list(tmp);
  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    exactSolution(i, j, m_grid->x(i), m_grid->y(j),
                  &(tmp(i,j).u), &(tmp(i,j).v));
  }
  tmp.write(file);

  tmp.metadata(0)
    .set_name("u_error")
    .long_name("X-component of the error (exact - computed)")
    .units("m s^-1");
  tmp.metadata(1)
    .set_name("v_error")
    .long_name("Y-component of the error (exact - computed)")
    .units("m s^-1");

  tmp.add(-1.0, m_ssa->velocity());
  tmp.write(file);

  array::Scalar error_mag(m_grid, "error_mag");
  error_mag.metadata(0).long_name("magnitude of the error").units("m s^-1");
  array::compute_magnitude(tmp, error_mag);
  error_mag.write(file);

  file.close();
}

} // end of namespace stressbalance
} // end of namespace pism
