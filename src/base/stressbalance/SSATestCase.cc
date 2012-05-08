// Copyright (C) 2009--2012 Ed Bueler, Constantine Khroulev, and David Maxwell
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

#include "SSATestCase.hh"
#include "PIO.hh"

#include "SSAFD.hh"
#include "SSAFEM.hh"
#include "pism_options.hh"

//! Initialize the storage for the various coefficients used as input to the SSA
//! (ice elevation, thickness, etc.)  
PetscErrorCode SSATestCase::buildSSACoefficients()
{
  PetscErrorCode ierr;

  const PetscInt WIDE_STENCIL = 2;
  
  // ice surface elevation
  ierr = surface.create(grid, "usurf", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = surface.set_attrs("diagnostic", "ice upper surface elevation", "m", 
                                      "surface_altitude"); CHKERRQ(ierr);
  ierr = vars.add(surface); CHKERRQ(ierr);
  
  // land ice thickness
  ierr = thickness.create(grid, "thk", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = thickness.set_attrs("model_state", "land ice thickness", "m", 
                                  "land_ice_thickness"); CHKERRQ(ierr);
  ierr = thickness.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = vars.add(thickness); CHKERRQ(ierr);

  // bedrock surface elevation
  ierr = bed.create(grid, "topg", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = bed.set_attrs("model_state", "bedrock surface elevation", "m", 
                                          "bedrock_altitude"); CHKERRQ(ierr);
  ierr = vars.add(bed); CHKERRQ(ierr);

  // yield stress for basal till (plastic or pseudo-plastic model)
  ierr = tauc.create(grid, "tauc", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = tauc.set_attrs("diagnostic",  
  "yield stress for basal till (plastic or pseudo-plastic model)", "Pa", "");
      CHKERRQ(ierr);
  ierr = vars.add(tauc); CHKERRQ(ierr);

  // enthalpy
  ierr = enthalpy.create(grid, "enthalpy", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = enthalpy.set_attrs("model_state",
              "ice enthalpy (includes sensible heat, latent heat, pressure)",
              "J kg-1", ""); CHKERRQ(ierr);
  ierr = vars.add(enthalpy); CHKERRQ(ierr);


  // dirichlet boundary condition (FIXME: perhaps unused!)
  ierr = vel_bc.create(grid, "_bc", true,WIDE_STENCIL); CHKERRQ(ierr); // u_bc and v_bc
  ierr = vel_bc.set_attrs("intent", 
            "X-component of the SSA velocity boundary conditions", 
            "m s-1", "", 0); CHKERRQ(ierr);
  ierr = vel_bc.set_attrs("intent", 
            "Y-component of the SSA velocity boundary conditions", 
            "m s-1", "", 1); CHKERRQ(ierr);
  ierr = vel_bc.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  ierr = vel_bc.set_attr("valid_min", convert(-1e6, "m/year", "m/second"), 0); CHKERRQ(ierr);
  ierr = vel_bc.set_attr("valid_max", convert( 1e6, "m/year", "m/second"), 0); CHKERRQ(ierr);
  ierr = vel_bc.set_attr("valid_min", convert(-1e6, "m/year", "m/second"), 1); CHKERRQ(ierr);
  ierr = vel_bc.set_attr("valid_max", convert( 1e6, "m/year", "m/second"), 1); CHKERRQ(ierr);
  ierr = vel_bc.set_attr("_FillValue",convert( 2e6, "m/year", "m/second"), 0); CHKERRQ(ierr);
  ierr = vel_bc.set_attr("_FillValue",convert( 2e6, "m/year", "m/second"), 1); CHKERRQ(ierr);
  vel_bc.write_in_glaciological_units = true;
  ierr = vel_bc.set(convert(2e6, "m/year", "m/second")); CHKERRQ(ierr);
  
  // grounded_dragging_floating integer mask
  ierr = ice_mask.create(grid, "mask", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = ice_mask.set_attrs("model_state", 
          "grounded_dragging_floating integer mask", "", ""); CHKERRQ(ierr);
  vector<double> mask_values(4);
  mask_values[0] = MASK_ICE_FREE_BEDROCK;
  mask_values[1] = MASK_GROUNDED;
  mask_values[2] = MASK_FLOATING;
  mask_values[3] = MASK_ICE_FREE_OCEAN;
  ierr = ice_mask.set_attr("flag_values", mask_values); CHKERRQ(ierr);
  ierr = ice_mask.set_attr("flag_meanings",
                           "ice_free_bedrock grounded_ice floating_ice ice_free_ocean");
  CHKERRQ(ierr);
  ice_mask.output_data_type = PISM_BYTE;
  ierr = vars.add(ice_mask); CHKERRQ(ierr);

  ierr = ice_mask.set(MASK_GROUNDED); CHKERRQ(ierr);

  // Dirichlet B.C. mask
  ierr = bc_mask.create(grid, "bc_mask", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = bc_mask.set_attrs("model_state", 
          "grounded_dragging_floating integer mask", "", ""); CHKERRQ(ierr);
  mask_values.resize(2);
  mask_values[0] = 0;
  mask_values[1] = 1;
  ierr = bc_mask.set_attr("flag_values", mask_values); CHKERRQ(ierr);
  ierr = bc_mask.set_attr("flag_meanings",
                          "no_data ssa_dirichlet_bc_location"); CHKERRQ(ierr);
  bc_mask.output_data_type = PISM_BYTE;
  ierr = vars.add(bc_mask); CHKERRQ(ierr);
  
  return 0;
}


//! Initialize the test case at the start of a run
PetscErrorCode SSATestCase::init(PetscInt Mx, PetscInt My, SSAFactory ssafactory)
{
  PetscErrorCode ierr;

  // Set options from command line.  
  ierr = config.scalar_from_option("ssa_eps",  "epsilon_ssa"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("ssa_maxi", "max_iterations_ssafd"); CHKERRQ(ierr);
  ierr = config.scalar_from_option("ssa_rtol", "ssafd_relative_convergence"); CHKERRQ(ierr);

  // Subclass builds grid.
  ierr = initializeGrid(Mx,My);

  // Subclass builds ice flow law, basal resistance, etc.
  ierr = initializeSSAModel(); CHKERRQ(ierr);

  // We setup storage for the coefficients.
  ierr = buildSSACoefficients(); CHKERRQ(ierr);

  // Allocate the actual SSA solver.
  ssa = ssafactory(grid, *basal, *enthalpyconverter, config);
  ierr = ssa->init(vars); CHKERRQ(ierr); // vars was setup preivouisly with buildSSACoefficients

  // Allow the subclass to setup the coefficients.
  ierr = initializeSSACoefficients(); CHKERRQ(ierr);

  return 0;
}

//! Solve the SSA
PetscErrorCode SSATestCase::run()
{
  PetscErrorCode ierr;
  // Solve (fast==true means "no update"):
  ierr = verbPrintf(2,grid.com,"* Solving the SSA stress balance ...\n"); CHKERRQ(ierr);

  bool fast = false;
  ierr = ssa->update(fast); CHKERRQ(ierr);

  return 0;
}

//! Report on the generated solution
PetscErrorCode SSATestCase::report(string testname) {
  PetscErrorCode  ierr;
    
  string ssa_stdout;
  ierr = ssa->stdout_report(ssa_stdout); CHKERRQ(ierr);
  ierr = verbPrintf(3,grid.com,ssa_stdout.c_str()); CHKERRQ(ierr);
  
  PetscScalar maxvecerr = 0.0, avvecerr = 0.0, 
    avuerr = 0.0, avverr = 0.0, maxuerr = 0.0, maxverr = 0.0;
  PetscScalar gmaxvecerr = 0.0, gavvecerr = 0.0, gavuerr = 0.0, gavverr = 0.0,
    gmaxuerr = 0.0, gmaxverr = 0.0;

  if (config.get_flag("do_pseudo_plastic_till")) {
    ierr = verbPrintf(1,grid.com, 
                    "WARNING: numerical errors not valid for pseudo-plastic till\n"); CHKERRQ(ierr);
  }
  ierr = verbPrintf(1,grid.com, 
                    "NUMERICAL ERRORS in velocity relative to exact solution:\n"); CHKERRQ(ierr);


  IceModelVec2V *vel_ssa;
  ierr = ssa->get_advective_2d_velocity(vel_ssa); CHKERRQ(ierr);
  ierr = vel_ssa->begin_access(); CHKERRQ(ierr);

  PetscScalar exactvelmax = 0, gexactvelmax = 0;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      PetscScalar uexact, vexact;
      PetscScalar myx = grid.x[i], myy = grid.y[j];

      exactSolution(i,j,myx,myy,&uexact,&vexact);

      PetscScalar exactnormsq=sqrt(uexact*uexact+vexact*vexact);
      exactvelmax = PetscMax(exactnormsq,exactvelmax);

      // compute maximum errors
      const PetscScalar uerr = PetscAbsReal((*vel_ssa)(i,j).u - uexact);
      const PetscScalar verr = PetscAbsReal((*vel_ssa)(i,j).v - vexact);
      avuerr = avuerr + uerr;
      avverr = avverr + verr;
      maxuerr = PetscMax(maxuerr,uerr);
      maxverr = PetscMax(maxverr,verr);
      const PetscScalar vecerr = sqrt(uerr * uerr + verr * verr);
      maxvecerr = PetscMax(maxvecerr,vecerr);
      avvecerr = avvecerr + vecerr;
    }
  }
  ierr = vel_ssa->end_access(); CHKERRQ(ierr);


  ierr = PISMGlobalMax(&exactvelmax, &gexactvelmax,grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalMax(&maxuerr, &gmaxuerr, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalMax(&maxverr, &gmaxverr, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&avuerr, &gavuerr, grid.com); CHKERRQ(ierr);
  gavuerr = gavuerr/(grid.Mx*grid.My);
  ierr = PISMGlobalSum(&avverr, &gavverr, grid.com); CHKERRQ(ierr);
  gavverr = gavverr/(grid.Mx*grid.My);
  ierr = PISMGlobalMax(&maxvecerr, &gmaxvecerr, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&avvecerr, &gavvecerr, grid.com); CHKERRQ(ierr);
  gavvecerr = gavvecerr/(grid.Mx*grid.My);

  ierr = verbPrintf(1,grid.com, 
                    "velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv\n");
  CHKERRQ(ierr);
  ierr = verbPrintf(1,grid.com, 
                    "           %11.4f%13.5f%10.4f%10.4f%10.4f%10.4f\n", 
                    gmaxvecerr*report_velocity_scale, (gavvecerr/gexactvelmax)*100.0,
                    gmaxuerr*report_velocity_scale, gmaxverr*report_velocity_scale, gavuerr*report_velocity_scale, 
                    gavverr*report_velocity_scale); CHKERRQ(ierr);

  ierr = verbPrintf(1,grid.com, "NUM ERRORS DONE\n");  CHKERRQ(ierr);

  ierr = report_netcdf(testname,
                       gmaxvecerr*report_velocity_scale,
                       (gavvecerr/gexactvelmax)*100.0,
                       gmaxuerr*report_velocity_scale,
                       gmaxverr*report_velocity_scale,
                       gavuerr*report_velocity_scale,
                       gavverr*report_velocity_scale); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode SSATestCase::report_netcdf(string testname,
                                          double max_vector,
                                          double rel_vector,
                                          double max_u,
                                          double max_v,
                                          double avg_u,
                                          double avg_v) {
  PetscErrorCode ierr;
  NCTimeseries err;
  unsigned int start;
  string filename;
  bool flag, append;
  NCGlobalAttributes global_attributes;

  ierr = PISMOptionsString("-report_file", "NetCDF error report file",
                           filename, flag); CHKERRQ(ierr);

  if (flag == false)
    return 0;

  ierr = verbPrintf(2, grid.com, "Also writing errors to '%s'...\n", filename.c_str());
  CHKERRQ(ierr);

  ierr = PISMOptionsIsSet("-append", "Append the NetCDF error report",
                          append); CHKERRQ(ierr);

  global_attributes.init("global_attributes", grid.com, grid.rank);
  global_attributes.set_string("source", string("PISM ") + PISM_Revision);

  // Find the number of records in this file:
  PIO nc(grid.com, grid.rank, "netcdf3");
  ierr = nc.open(filename, PISM_WRITE, append); CHKERRQ(ierr); // append == true
  ierr = nc.inq_dimlen("N", start); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  ierr = global_attributes.write(filename); CHKERRQ(ierr);

  // Write the dimension variable:
  err.init("N", "N", grid.com, grid.rank);
  ierr = err.write(filename, (size_t)start, (double)(start + 1), PISM_INT); CHKERRQ(ierr);

  // Always write grid parameters:
  err.short_name = "dx";
  ierr = err.set_units("meters"); CHKERRQ(ierr);
  ierr = err.write(filename, (size_t)start, grid.dx); CHKERRQ(ierr);
  err.short_name = "dy";
  ierr = err.write(filename, (size_t)start, grid.dy); CHKERRQ(ierr);

  // Always write the test name:
  err.reset();
  err.short_name = "test";
  ierr = err.write(filename, (size_t)start, testname[0], PISM_BYTE); CHKERRQ(ierr);

  err.reset();
  err.short_name = "max_velocity";
  ierr = err.set_units("m/year"); CHKERRQ(ierr);
  err.set_string("long_name", "maximum ice velocity magnitude error");
  ierr = err.write(filename, (size_t)start, max_vector); CHKERRQ(ierr);

  err.reset();
  err.short_name = "relative_velocity";
  ierr = err.set_units("percent"); CHKERRQ(ierr);
  err.set_string("long_name", "relative ice velocity magnitude error");
  ierr = err.write(filename, (size_t)start, rel_vector); CHKERRQ(ierr);

  err.reset();
  err.short_name = "maximum_u";
  ierr = err.set_units("m/year"); CHKERRQ(ierr);
  err.set_string("long_name", "maximum error in the X-component of the ice velocity");
  ierr = err.write(filename, (size_t)start, max_u); CHKERRQ(ierr);

  err.reset();
  err.short_name = "maximum_v";
  ierr = err.set_units("m/year"); CHKERRQ(ierr);
  err.set_string("long_name", "maximum error in the Y-component of the ice velocity");
  ierr = err.write(filename, (size_t)start, max_v); CHKERRQ(ierr);

  err.reset();
  err.short_name = "average_u";
  ierr = err.set_units("m/year"); CHKERRQ(ierr);
  err.set_string("long_name", "average error in the X-component of the ice velocity");
  ierr = err.write(filename, (size_t)start, avg_u); CHKERRQ(ierr);

  err.reset();
  err.short_name = "average_v";
  ierr = err.set_units("m/year"); CHKERRQ(ierr);
  err.set_string("long_name", "average error in the Y-component of the ice velocity");
  ierr = err.write(filename, (size_t)start, avg_v); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode SSATestCase::exactSolution(PetscInt /*i*/, PetscInt /*j*/, 
                                          PetscReal /*x*/, PetscReal /*y*/,
                                          PetscReal *u, PetscReal *v )
{
  *u=0; *v=0;
  return 0;
}

//! Save the computation and data to a file.
PetscErrorCode SSATestCase::write(const string &filename)
{
  PetscErrorCode ierr;

  // Write results to an output file:
  PIO pio(grid.com, grid.rank, grid.config.get_string("output_format"));
  ierr = pio.open(filename, PISM_WRITE); CHKERRQ(ierr);
  ierr = pio.def_time(config.get_string("time_dimension_name"),
                      config.get_string("calendar"),
                      grid.time->CF_units()); CHKERRQ(ierr);
  ierr = pio.append_time(config.get_string("time_dimension_name"), 0.0); CHKERRQ(ierr);
  ierr = pio.close(); CHKERRQ(ierr);

  ierr = surface.write(filename.c_str()); CHKERRQ(ierr);
  ierr = thickness.write(filename.c_str()); CHKERRQ(ierr);
  ierr = bc_mask.write(filename.c_str()); CHKERRQ(ierr);
  ierr = tauc.write(filename.c_str()); CHKERRQ(ierr);
  ierr = bed.write(filename.c_str()); CHKERRQ(ierr);
  ierr = enthalpy.write(filename.c_str()); CHKERRQ(ierr);
  ierr = vel_bc.write(filename.c_str()); CHKERRQ(ierr);

  IceModelVec2V *vel_ssa;
  ierr = ssa->get_advective_2d_velocity(vel_ssa); CHKERRQ(ierr);
  ierr = vel_ssa->write(filename.c_str()); CHKERRQ(ierr);

  IceModelVec2V exact;
  ierr = exact.create(grid, "_exact", false); CHKERRQ(ierr);
  ierr = exact.set_attrs("diagnostic", 
            "X-component of the SSA exact solution", 
            "m s-1", "", 0); CHKERRQ(ierr);
  ierr = exact.set_attrs("diagnostic", 
            "Y-component of the SSA exact solution", 
            "m s-1", "", 1); CHKERRQ(ierr);
  ierr = exact.set_glaciological_units("m year-1"); CHKERRQ(ierr);

  ierr = exact.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      PetscScalar myx = grid.x[i], myy = grid.y[j];
      exactSolution(i,j,myx,myy,&(exact(i,j).u),&(exact(i,j).v));
    }
  }
  ierr = exact.end_access(); CHKERRQ(ierr);
  ierr = exact.write(filename.c_str()); CHKERRQ(ierr);

  return 0;
}


/*! Initialize a uniform, shallow (3 z-levels), doubly periodic grid with 
half-widths (Lx,Ly) and Mx by My nodes for time-independent computations.*/
PetscErrorCode init_shallow_grid(IceGrid &grid, PetscReal Lx, 
                                 PetscReal Ly, PetscInt Mx, PetscInt My, Periodicity p)
{
  PetscErrorCode ierr;
  
  grid.Lx = Lx;
  grid.Ly = Ly;
  grid.periodicity = p;
  grid.Mx = Mx; grid.My=My; grid.Mz = 3;
  
  grid.compute_nprocs();
  grid.compute_ownership_ranges();
  ierr = grid.compute_vertical_levels(); CHKERRQ(ierr);
  ierr = grid.compute_horizontal_spacing(); CHKERRQ(ierr);
  ierr = grid.createDA(); CHKERRQ(ierr);

  return 0;
}

