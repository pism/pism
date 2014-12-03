// Copyright (C) 2009--2014 Ed Bueler, Constantine Khroulev, and David Maxwell
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

#include "SSATestCase.hh"
#include "PIO.hh"

#include "SSAFD.hh"
#include "SSAFEM.hh"
#include "pism_options.hh"
#include "Mask.hh"

namespace pism {

//! Initialize the storage for the various coefficients used as input to the SSA
//! (ice elevation, thickness, etc.)
PetscErrorCode SSATestCase::buildSSACoefficients()
{

  const unsigned int WIDE_STENCIL = config.get("grid_max_stencil_width");

  // ice surface elevation
  surface.create(*grid, "usurf", WITH_GHOSTS, WIDE_STENCIL);
  surface.set_attrs("diagnostic", "ice upper surface elevation", "m",
                    "surface_altitude");
  vars.add(surface);

  // land ice thickness
  thickness.create(*grid, "thk", WITH_GHOSTS, WIDE_STENCIL);
  thickness.set_attrs("model_state", "land ice thickness", "m",
                      "land_ice_thickness");
  thickness.metadata().set_double("valid_min", 0.0);
  vars.add(thickness);

  // bedrock surface elevation
  bed.create(*grid, "topg", WITH_GHOSTS, WIDE_STENCIL);
  bed.set_attrs("model_state", "bedrock surface elevation", "m",
                "bedrock_altitude");
  vars.add(bed);

  // yield stress for basal till (plastic or pseudo-plastic model)
  tauc.create(*grid, "tauc", WITH_GHOSTS, WIDE_STENCIL);
  tauc.set_attrs("diagnostic",
                 "yield stress for basal till (plastic or pseudo-plastic model)", "Pa", "");
  vars.add(tauc);

  // enthalpy
  enthalpy.create(*grid, "enthalpy", WITH_GHOSTS, WIDE_STENCIL);
  enthalpy.set_attrs("model_state",
                     "ice enthalpy (includes sensible heat, latent heat, pressure)",
                     "J kg-1", "");
  vars.add(enthalpy);


  // dirichlet boundary condition (FIXME: perhaps unused!)
  vel_bc.create(*grid, "_bc", WITH_GHOSTS, WIDE_STENCIL); // u_bc and v_bc
  vel_bc.set_attrs("intent",
                   "X-component of the SSA velocity boundary conditions",
                   "m s-1", "", 0);
  vel_bc.set_attrs("intent",
                   "Y-component of the SSA velocity boundary conditions",
                   "m s-1", "", 1);
  vel_bc.set_glaciological_units("m year-1");
  vel_bc.metadata(0).set_double("valid_min", grid->convert(-1e6, "m/year", "m/second"));
  vel_bc.metadata(0).set_double("valid_max", grid->convert( 1e6, "m/year", "m/second"));
  vel_bc.metadata(0).set_double("_FillValue", config.get("fill_value", "m/year", "m/s"));
  vel_bc.metadata(1).set_double("valid_min", grid->convert(-1e6, "m/year", "m/second"));
  vel_bc.metadata(1).set_double("valid_max", grid->convert( 1e6, "m/year", "m/second"));
  vel_bc.metadata(1).set_double("_FillValue", config.get("fill_value", "m/year", "m/s"));
  vel_bc.write_in_glaciological_units = true;
  vel_bc.set(config.get("fill_value", "m/year", "m/s"));

  // grounded_dragging_floating integer mask
  ice_mask.create(*grid, "mask", WITH_GHOSTS, WIDE_STENCIL);
  ice_mask.set_attrs("model_state",
                     "grounded_dragging_floating integer mask", "", "");
  std::vector<double> mask_values(4);
  mask_values[0] = MASK_ICE_FREE_BEDROCK;
  mask_values[1] = MASK_GROUNDED;
  mask_values[2] = MASK_FLOATING;
  mask_values[3] = MASK_ICE_FREE_OCEAN;
  ice_mask.metadata().set_doubles("flag_values", mask_values);
  ice_mask.metadata().set_string("flag_meanings",
                                 "ice_free_bedrock grounded_ice floating_ice ice_free_ocean");
  vars.add(ice_mask);

  ice_mask.set(MASK_GROUNDED);

  // Dirichlet B.C. mask
  bc_mask.create(*grid, "bc_mask", WITH_GHOSTS, WIDE_STENCIL);
  bc_mask.set_attrs("model_state",
                    "grounded_dragging_floating integer mask", "", "");
  mask_values.resize(2);
  mask_values[0] = 0;
  mask_values[1] = 1;
  bc_mask.metadata().set_doubles("flag_values", mask_values);
  bc_mask.metadata().set_string("flag_meanings",
                                "no_data ssa_dirichlet_bc_location");
  vars.add(bc_mask);

  melange_back_pressure.create(*grid, "melange_back_pressure_fraction",
                               WITH_GHOSTS, WIDE_STENCIL);
  melange_back_pressure.set_attrs("boundary_condition",
                                  "melange back pressure fraction", "", "");
  melange_back_pressure.set(0.0);

  return 0;
}

SSATestCase::SSATestCase(MPI_Comm com, Config &c)
  : m_com(com), config(c), enthalpyconverter(NULL), ssa(NULL)
{
  // empty
}

SSATestCase::~SSATestCase()
{
  delete enthalpyconverter;
  delete ssa;
}

//! Initialize the test case at the start of a run
PetscErrorCode SSATestCase::init(int Mx, int My, SSAFactory ssafactory)
{
  PetscErrorCode ierr;

  // Set options from command line.
  config.scalar_from_option("ssa_eps",  "epsilon_ssa");
  config.scalar_from_option("ssa_maxi", "max_iterations_ssafd");
  config.scalar_from_option("ssa_rtol", "ssafd_relative_convergence");

  // Subclass builds grid->
  ierr = initializeGrid(Mx,My);

  // Subclass builds ice flow law, basal resistance, etc.
  initializeSSAModel();

  // We setup storage for the coefficients.
  buildSSACoefficients();

  // Allocate the actual SSA solver.
  ssa = ssafactory(*grid, *enthalpyconverter, config);
  ssa->init(vars); // vars was setup preivouisly with buildSSACoefficients

  // Allow the subclass to setup the coefficients.
  initializeSSACoefficients();

  return 0;
}

//! Solve the SSA
PetscErrorCode SSATestCase::run()
{
  // Solve (fast==true means "no update"):
  verbPrintf(2,grid->com,"* Solving the SSA stress balance ...\n");

  bool fast = false;
  ssa->update(fast, melange_back_pressure);

  return 0;
}

//! Report on the generated solution
PetscErrorCode SSATestCase::report(const std::string &testname) {

  std::string ssa_stdout;
  ssa->stdout_report(ssa_stdout);
  verbPrintf(3,grid->com,ssa_stdout.c_str());

  double maxvecerr = 0.0, avvecerr = 0.0,
    avuerr = 0.0, avverr = 0.0, maxuerr = 0.0, maxverr = 0.0;
  double gmaxvecerr = 0.0, gavvecerr = 0.0, gavuerr = 0.0, gavverr = 0.0,
    gmaxuerr = 0.0, gmaxverr = 0.0;

  if (config.get_flag("do_pseudo_plastic_till") &&
      config.get("pseudo_plastic_q") != 1.0) {
    verbPrintf(1,grid->com,
               "WARNING: numerical errors not valid for pseudo-plastic till\n");
  }
  verbPrintf(1,grid->com,
             "NUMERICAL ERRORS in velocity relative to exact solution:\n");


  IceModelVec2V *vel_ssa;
  ssa->get_2D_advective_velocity(vel_ssa);

  IceModelVec::AccessList list;
  list.add(*vel_ssa);

  double exactvelmax = 0, gexactvelmax = 0;
  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double uexact, vexact;
    double myx = grid->x[i], myy = grid->y[j];

    exactSolution(i,j,myx,myy,&uexact,&vexact);

    double exactnormsq=sqrt(uexact*uexact+vexact*vexact);
    exactvelmax = PetscMax(exactnormsq,exactvelmax);

    // compute maximum errors
    const double uerr = PetscAbsReal((*vel_ssa)(i,j).u - uexact);
    const double verr = PetscAbsReal((*vel_ssa)(i,j).v - vexact);
    avuerr = avuerr + uerr;
    avverr = avverr + verr;
    maxuerr = PetscMax(maxuerr,uerr);
    maxverr = PetscMax(maxverr,verr);
    const double vecerr = sqrt(uerr * uerr + verr * verr);
    maxvecerr = PetscMax(maxvecerr,vecerr);
    avvecerr = avvecerr + vecerr;
  }

  unsigned int N = (grid->Mx()*grid->My());
  GlobalMax(grid->com, &exactvelmax,  &gexactvelmax);
  GlobalMax(grid->com, &maxuerr,  &gmaxuerr);
  GlobalMax(grid->com, &maxverr,  &gmaxverr);
  GlobalSum(grid->com, &avuerr,  &gavuerr);
  gavuerr = gavuerr / N;
  GlobalSum(grid->com, &avverr,  &gavverr);
  gavverr = gavverr / N;
  GlobalMax(grid->com, &maxvecerr,  &gmaxvecerr);
  GlobalSum(grid->com, &avvecerr,  &gavvecerr);
  gavvecerr = gavvecerr / N;

  verbPrintf(1,grid->com,
             "velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv\n");
  verbPrintf(1,grid->com,
             "           %11.4f%13.5f%10.4f%10.4f%10.4f%10.4f\n",
             grid->convert(gmaxvecerr, "m/second", "m/year"),
             (gavvecerr/gexactvelmax)*100.0,
             grid->convert(gmaxuerr, "m/second", "m/year"),
             grid->convert(gmaxverr, "m/second", "m/year"),
             grid->convert(gavuerr, "m/second", "m/year"),
             grid->convert(gavverr, "m/second", "m/year"));

  verbPrintf(1,grid->com, "NUM ERRORS DONE\n");

  report_netcdf(testname,
                grid->convert(gmaxvecerr, "m/second", "m/year"),
                (gavvecerr/gexactvelmax)*100.0,
                grid->convert(gmaxuerr, "m/second", "m/year"),
                grid->convert(gmaxverr, "m/second", "m/year"),
                grid->convert(gavuerr, "m/second", "m/year"),
                grid->convert(gavverr, "m/second", "m/year"));

  return 0;
}

PetscErrorCode SSATestCase::report_netcdf(const std::string &testname,
                                          double max_vector,
                                          double rel_vector,
                                          double max_u,
                                          double max_v,
                                          double avg_u,
                                          double avg_v) {
  NCTimeseries err("N", "N", grid->get_unit_system());
  unsigned int start;
  std::string filename;
  bool flag, append;
  NCVariable global_attributes("PISM_GLOBAL", grid->get_unit_system());

  OptionsString("-report_file", "NetCDF error report file",
                filename, flag);

  if (flag == false) {
    return 0;
  }

  err.set_units("1");

  verbPrintf(2, grid->com, "Also writing errors to '%s'...\n", filename.c_str());

  OptionsIsSet("-append", "Append the NetCDF error report",
               append);

  IO_Mode mode = PISM_READWRITE;
  if (append == false) {
    mode = PISM_READWRITE_MOVE;
  }

  global_attributes.set_string("source", std::string("PISM ") + PISM_Revision);

  // Find the number of records in this file:
  PIO nc(*grid, "netcdf3");      // OK to use NetCDF3.
  nc.open(filename, mode);
  start = nc.inq_dimlen("N");

  nc.write_global_attributes(global_attributes);

  // Write the dimension variable:
  nc.write_timeseries(err, (size_t)start, (double)(start + 1), PISM_INT);

  // Always write grid parameters:
  err.set_name("dx");
  err.set_units("meters");
  nc.write_timeseries(err, (size_t)start, grid->dx());
  err.set_name("dy");
  nc.write_timeseries(err, (size_t)start, grid->dy());

  // Always write the test name:
  err.clear_all_strings(); err.clear_all_doubles(); err.set_units("1");
  err.set_name("test");
  nc.write_timeseries(err, (size_t)start, testname[0], PISM_BYTE);

  err.clear_all_strings(); err.clear_all_doubles(); err.set_units("1");
  err.set_name("max_velocity");
  err.set_units("m/year");
  err.set_string("long_name", "maximum ice velocity magnitude error");
  nc.write_timeseries(err, (size_t)start, max_vector);

  err.clear_all_strings(); err.clear_all_doubles(); err.set_units("1");
  err.set_name("relative_velocity");
  err.set_units("percent");
  err.set_string("long_name", "relative ice velocity magnitude error");
  nc.write_timeseries(err, (size_t)start, rel_vector);

  err.clear_all_strings(); err.clear_all_doubles(); err.set_units("1");
  err.set_name("maximum_u");
  err.set_units("m/year");
  err.set_string("long_name", "maximum error in the X-component of the ice velocity");
  nc.write_timeseries(err, (size_t)start, max_u);

  err.clear_all_strings(); err.clear_all_doubles(); err.set_units("1");
  err.set_name("maximum_v");
  err.set_units("m/year");
  err.set_string("long_name", "maximum error in the Y-component of the ice velocity");
  nc.write_timeseries(err, (size_t)start, max_v);

  err.clear_all_strings(); err.clear_all_doubles(); err.set_units("1");
  err.set_name("average_u");
  err.set_units("m/year");
  err.set_string("long_name", "average error in the X-component of the ice velocity");
  nc.write_timeseries(err, (size_t)start, avg_u);

  err.clear_all_strings(); err.clear_all_doubles(); err.set_units("1");
  err.set_name("average_v");
  err.set_units("m/year");
  err.set_string("long_name", "average error in the Y-component of the ice velocity");
  nc.write_timeseries(err, (size_t)start, avg_v);

  nc.close();

  return 0;
}

PetscErrorCode SSATestCase::exactSolution(int /*i*/, int /*j*/,
                                          double /*x*/, double /*y*/,
                                          double *u, double *v)
{
  *u=0; *v=0;
  return 0;
}

//! Save the computation and data to a file.
PetscErrorCode SSATestCase::write(const std::string &filename)
{

  // Write results to an output file:
  PIO pio(*grid, grid->config.get_string("output_format"));
  pio.open(filename, PISM_READWRITE_MOVE);
  pio.def_time(config.get_string("time_dimension_name"),
               grid->time->calendar(),
               grid->time->CF_units_string());
  pio.append_time(config.get_string("time_dimension_name"), 0.0);

  surface.write(pio);
  thickness.write(pio);
  bc_mask.write(pio);
  tauc.write(pio);
  bed.write(pio);
  enthalpy.write(pio);
  vel_bc.write(pio);

  IceModelVec2V *vel_ssa;
  ssa->get_2D_advective_velocity(vel_ssa);
  vel_ssa->write(pio);

  IceModelVec2V exact;
  exact.create(*grid, "_exact", WITHOUT_GHOSTS);
  exact.set_attrs("diagnostic",
                  "X-component of the SSA exact solution",
                  "m s-1", "", 0);
  exact.set_attrs("diagnostic",
                  "Y-component of the SSA exact solution",
                  "m s-1", "", 1);
  exact.set_glaciological_units("m year-1");
  exact.write_in_glaciological_units = true;

  IceModelVec::AccessList list(exact);
  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    exactSolution(i, j, grid->x[i], grid->y[j],
                  &(exact(i,j).u), &(exact(i,j).v));
  }
  exact.write(pio);

  pio.close();
  return 0;
}

} // end of namespace pism
