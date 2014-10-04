// Copyright (C) 2010, 2011, 2012, 2013, 2014 Ed Bueler and Constantine Khroulev
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
  "\nSIAFD_TEST\n"
  "  Testing program for SIA, time-independent calculations separate from\n"
  "  IceModel. Uses verification test F. Also may be used in a PISM software"
  "(regression) test.\n\n";

#include "pism_const.hh"
#include "pism_options.hh"
#include "iceModelVec.hh"
#include "flowlaws.hh" // IceFlowLaw
#include "PIO.hh"
#include "NCVariable.hh"
#include "PISMStressBalance.hh"
#include "SIAFD.hh"
#include "exactTestsFG.h"
#include "basal_resistance.hh"
#include "enthalpyConverter.hh"
#include "SSB_Modifier.hh"
#include "ShallowStressBalance.hh"
#include "PISMVars.hh"
#include "Mask.hh"

using namespace pism;

PetscErrorCode compute_strain_heating_errors(const Config &config,
                                  IceModelVec3 &strain_heating,
                                  IceModelVec2S &thickness,
                                  IceGrid &grid,
                                  double &gmax_strain_heating_err,
                                  double &gav_strain_heating_err) {

  PetscErrorCode ierr;
  double    max_strain_heating_error = 0.0, av_strain_heating_error = 0.0, avcount = 0.0;
  const int Mz = grid.Mz;

  const double LforFG = 750000; // m

  const double
    ice_rho   = config.get("ice_density"),
    ice_c     = config.get("ice_specific_heat_capacity");

  double   *dummy1, *dummy2, *dummy3, *dummy4, *strain_heating_exact;
  double   junk0, junk1;

  strain_heating_exact = new double[Mz];
  dummy1 = new double[Mz];  dummy2 = new double[Mz];
  dummy3 = new double[Mz];  dummy4 = new double[Mz];

  double *strain_heating_ij;

  IceModelVec::AccessList list;
  list.add(thickness);
  list.add(strain_heating);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double xx = grid.x[i], yy = grid.y[j],
      r = sqrt(PetscSqr(xx) + PetscSqr(yy));
    if ((r >= 1.0) && (r <= LforFG - 1.0)) {  // only evaluate error if inside sheet
      // and not at central singularity
      bothexact(0.0,r,&grid.zlevels[0],Mz,0.0,
                &junk0,&junk1,dummy1,dummy2,dummy3,strain_heating_exact,dummy4);

      for (int k = 0; k < Mz; k++)
        strain_heating_exact[k] *= ice_rho * ice_c; // scale exact strain_heating to J/(s m^3)
      const int ks = grid.kBelowHeight(thickness(i,j));
      ierr = strain_heating.getInternalColumn(i,j,&strain_heating_ij); CHKERRQ(ierr);
      for (int k = 0; k < ks; k++) {  // only eval error if below num surface
        const double _strain_heating_error = PetscAbs(strain_heating_ij[k] - strain_heating_exact[k]);
        max_strain_heating_error = PetscMax(max_strain_heating_error,_strain_heating_error);
        avcount += 1.0;
        av_strain_heating_error += _strain_heating_error;
      }
    }
  }

  delete [] strain_heating_exact;
  delete [] dummy1;  delete [] dummy2;  delete [] dummy3;  delete [] dummy4;

  ierr = GlobalMax(&max_strain_heating_error, &gmax_strain_heating_err, grid.com); CHKERRQ(ierr);
  ierr = GlobalSum(&av_strain_heating_error, &gav_strain_heating_err, grid.com); CHKERRQ(ierr);
  double  gavcount;
  ierr = GlobalSum(&avcount, &gavcount, grid.com); CHKERRQ(ierr);
  gav_strain_heating_err = gav_strain_heating_err/PetscMax(gavcount,1.0);  // avoid div by zero
  return 0;
}


PetscErrorCode computeSurfaceVelocityErrors(IceGrid &grid,
                                            IceModelVec2S &ice_thickness,
                                            IceModelVec3 &u3,
                                            IceModelVec3 &v3,
                                            IceModelVec3 &w3,
                                            double &gmaxUerr,
                                            double &gavUerr,
                                            double &gmaxWerr,
                                            double &gavWerr) {

  PetscErrorCode ierr;
  double    maxUerr = 0.0, maxWerr = 0.0, avUerr = 0.0, avWerr = 0.0;

  const double LforFG = 750000; // m

  IceModelVec::AccessList list;
  list.add(ice_thickness);
  list.add(u3);
  list.add(v3);
  list.add(w3);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double xx = grid.x[i], yy = grid.y[j],
      r = sqrt(PetscSqr(xx) + PetscSqr(yy));
    if ((r >= 1.0) && (r <= LforFG - 1.0)) {  // only evaluate error if inside sheet
      // and not at central singularity
      double radialUex,wex;
      double dummy0,dummy1,dummy2,dummy3,dummy4;
      bothexact(0.0,r,&(ice_thickness(i,j)),1,0.0,
                &dummy0,&dummy1,&dummy2,&radialUex,&wex,&dummy3,&dummy4);

      const double uex = (xx/r) * radialUex;
      const double vex = (yy/r) * radialUex;
      // note that because getValZ does linear interpolation and ice_thickness(i,j) is not exactly at
      // a grid point, this causes nonzero errors even with option -eo
      const double Uerr = sqrt(PetscSqr(u3.getValZ(i,j,ice_thickness(i,j)) - uex)
                               + PetscSqr(v3.getValZ(i,j,ice_thickness(i,j)) - vex));
      maxUerr = PetscMax(maxUerr,Uerr);
      avUerr += Uerr;
      const double Werr = PetscAbs(w3.getValZ(i,j,ice_thickness(i,j)) - wex);
      maxWerr = PetscMax(maxWerr,Werr);
      avWerr += Werr;
    }
  }

  ierr = GlobalMax(&maxUerr, &gmaxUerr, grid.com); CHKERRQ(ierr);
  ierr = GlobalMax(&maxWerr, &gmaxWerr, grid.com); CHKERRQ(ierr);
  ierr = GlobalSum(&avUerr, &gavUerr, grid.com); CHKERRQ(ierr);
  gavUerr = gavUerr/(grid.Mx*grid.My);
  ierr = GlobalSum(&avWerr, &gavWerr, grid.com); CHKERRQ(ierr);
  gavWerr = gavWerr/(grid.Mx*grid.My);
  return 0;
}


PetscErrorCode enthalpy_from_temperature_cold(EnthalpyConverter &EC,
                                              IceGrid &grid,
                                              IceModelVec2S &thickness,
                                              IceModelVec3 &temperature,
                                              IceModelVec3 &enthalpy) {
  PetscErrorCode ierr;

  IceModelVec::AccessList list;
  list.add(temperature);
  list.add(enthalpy);
  list.add(thickness);

  double *T_ij, *E_ij; // columns of these values
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    ierr = temperature.getInternalColumn(i,j,&T_ij); CHKERRQ(ierr);
    ierr = enthalpy.getInternalColumn(i,j,&E_ij); CHKERRQ(ierr);

    for (unsigned int k=0; k<grid.Mz; ++k) {
      double depth = thickness(i,j) - grid.zlevels[k];
      ierr = EC.getEnthPermissive(T_ij[k], 0.0,
                                  EC.getPressureFromDepth(depth),
                                  E_ij[k]); CHKERRQ(ierr);
    }

  }


  ierr = enthalpy.update_ghosts(); CHKERRQ(ierr);
  return 0;
}


//! \brief Set the test F initial state.
PetscErrorCode setInitStateF(IceGrid &grid,
                             EnthalpyConverter &EC,
                             IceModelVec2S *bed,
                             IceModelVec2Int *mask,
                             IceModelVec2S *surface,
                             IceModelVec2S *thickness,
                             IceModelVec3 *enthalpy) {
  PetscErrorCode ierr;
  int        Mz=grid.Mz;
  double     *dummy1, *dummy2, *dummy3, *dummy4, *dummy5;

  double ST = 1.67e-5,
    Tmin = 223.15,  // K
    LforFG = 750000; // m

  dummy1=new double[Mz];  dummy2=new double[Mz];
  dummy3=new double[Mz];  dummy4=new double[Mz];
  dummy5=new double[Mz];

  ierr = bed->set(0); CHKERRQ(ierr);
  ierr = mask->set(MASK_GROUNDED); CHKERRQ(ierr);

  double *T;
  T = new double[grid.Mz];

  IceModelVec::AccessList list;
  list.add(*thickness);
  list.add(*enthalpy);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double r = grid.radius(i, j),
      Ts = Tmin + ST * r;
    if (r > LforFG - 1.0) { // if (essentially) outside of sheet
      (*thickness)(i, j) = 0.0;
      for (int k = 0; k < Mz; k++)
        T[k]=Ts;
    } else {
      r = PetscMax(r, 1.0); // avoid singularity at origin
      bothexact(0.0, r, &grid.zlevels[0], Mz, 0.0,
                &(*thickness)(i, j), dummy5, T, dummy1, dummy2, dummy3, dummy4);
    }
    ierr = enthalpy->setInternalColumn(i, j, T); CHKERRQ(ierr);
  }


  ierr = thickness->update_ghosts(); CHKERRQ(ierr);

  ierr = thickness->copy_to(*surface); CHKERRQ(ierr);

  delete [] dummy1;  delete [] dummy2;  delete [] dummy3;  delete [] dummy4;
  delete [] T; delete [] dummy5;

  ierr = enthalpy_from_temperature_cold(EC, grid, *thickness,
                                        *enthalpy,
                                        *enthalpy); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode reportErrors(const Config &config,
                            IceGrid &grid,
                            IceModelVec2S *thickness,
                            IceModelVec3 *u_sia, IceModelVec3 *v_sia,
                            IceModelVec3 *w_sia,
                            IceModelVec3 *strain_heating) {
  PetscErrorCode ierr;

  // strain_heating errors if appropriate; reported in 10^6 J/(s m^3)
  double max_strain_heating_error, av_strain_heating_error;
  ierr = compute_strain_heating_errors(config, *strain_heating, *thickness,
                            grid,
                            max_strain_heating_error, av_strain_heating_error); CHKERRQ(ierr);

  ierr = verbPrintf(1,grid.com,
                    "Sigma     :      maxSig       avSig\n"); CHKERRQ(ierr);
  ierr = verbPrintf(1,grid.com, "           %12.6f%12.6f\n",
                    max_strain_heating_error*1.0e6, av_strain_heating_error*1.0e6); CHKERRQ(ierr);

  // surface velocity errors if exact values are available; reported in m/year
  double maxUerr, avUerr, maxWerr, avWerr;
  ierr = computeSurfaceVelocityErrors(grid, *thickness,
                                      *u_sia,
                                      *v_sia,
                                      *w_sia,
                                      maxUerr, avUerr,
                                      maxWerr, avWerr); CHKERRQ(ierr);

  ierr = verbPrintf(1,grid.com,
                    "surf vels :     maxUvec      avUvec        maxW         avW\n"); CHKERRQ(ierr);
  ierr = verbPrintf(1,grid.com, "           %12.6f%12.6f%12.6f%12.6f\n",
                    grid.convert(maxUerr, "m/s", "m/year"),
                    grid.convert(avUerr,  "m/s", "m/year"),
                    grid.convert(maxWerr, "m/s", "m/year"),
                    grid.convert(avWerr,  "m/s", "m/year")); CHKERRQ(ierr);

  return 0;
}


int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm    com;  // won't be used except for rank,size

  ierr = PetscInitialize(&argc, &argv, NULL, help); CHKERRQ(ierr);

  com = PETSC_COMM_WORLD;

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  {
    UnitSystem unit_system;
    Config config(com, "pism_config", unit_system),
      overrides(com, "pism_overrides", unit_system);
    ierr = init_config(com, config, overrides); CHKERRQ(ierr);

    config.set_flag("compute_grain_size_using_age", false);

    PetscBool usage_set, help_set;
    ierr = PetscOptionsHasName(NULL, "-usage", &usage_set); CHKERRQ(ierr);
    ierr = PetscOptionsHasName(NULL, "-help", &help_set); CHKERRQ(ierr);
    if ((usage_set==PETSC_TRUE) || (help_set==PETSC_TRUE)) {
      PetscPrintf(com,
                  "\n"
                  "usage of SIAFD_TEST:\n"
                  "  run siafd_test -Mx <number> -My <number> -Mz <number> -o foo.nc\n"
                  "\n");
    }

    IceGrid grid(com, config);

    grid.Lx = grid.Ly = 900e3;
    grid.Lz = 4000;
    grid.Mx = grid.My = 61;
    grid.Mz = 61;

    std::string output_file = "siafd_test_F.nc";
    int tmp = grid.Mz;
    ierr = PetscOptionsBegin(grid.com, "", "SIAFD_TEST options", ""); CHKERRQ(ierr);
    {
      bool flag;
      ierr = OptionsInt("-Mx", "Number of grid points in the X direction",
                            grid.Mx, flag); CHKERRQ(ierr);
      ierr = OptionsInt("-My", "Number of grid points in the X direction",
                            grid.My, flag); CHKERRQ(ierr);
      ierr = OptionsInt("-Mz", "Number of vertical grid levels",
                            tmp, flag); CHKERRQ(ierr);
      ierr = OptionsString("-o", "Set the output file name",
                               output_file, flag); CHKERRQ(ierr);
    }
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    if (tmp > 0) {
      grid.Mz = tmp;
    } else {
      PetscPrintf(grid.com, "PISM ERROR: -Mz %d is invalid (has to be positive).\n",
                  tmp);
      PISMEnd();
    }

    grid.compute_nprocs();
    grid.compute_ownership_ranges();
    ierr = grid.compute_vertical_levels(); CHKERRQ(ierr);
    ierr = grid.compute_horizontal_spacing(); CHKERRQ(ierr);
    ierr = grid.allocate(); CHKERRQ(ierr);

    ierr = setVerbosityLevel(5); CHKERRQ(ierr);
    ierr = grid.printInfo(1); CHKERRQ(ierr);
    //ierr = grid.printVertLevels(1); CHKERRQ(ierr);

    ICMEnthalpyConverter EC(config);
    ThermoGlenArrIce ice(grid.com, "", config, &EC);

    IceModelVec2S ice_surface_elevation, ice_thickness, bed_topography;
    IceModelVec2Int vMask;
    IceModelVec3 enthalpy,
      age;                      // is not used (and need not be allocated)
    const int WIDE_STENCIL = config.get("grid_max_stencil_width");

    Vars vars;

    // ice upper surface elevation
    ierr = ice_surface_elevation.create(grid, "usurf", WITH_GHOSTS, WIDE_STENCIL); CHKERRQ(ierr);
    ierr = ice_surface_elevation.set_attrs("diagnostic", "ice upper surface elevation",
          "m", "surface_altitude"); CHKERRQ(ierr);
    ierr = vars.add(ice_surface_elevation); CHKERRQ(ierr);

    // land ice thickness
    ierr = ice_thickness.create(grid, "thk", WITH_GHOSTS, WIDE_STENCIL); CHKERRQ(ierr);
    ierr = ice_thickness.set_attrs("model_state", "land ice thickness",
          "m", "land_ice_thickness"); CHKERRQ(ierr);
    ice_thickness.metadata().set_double("valid_min", 0.0);
    ierr = vars.add(ice_thickness); CHKERRQ(ierr);

    // bedrock surface elevation
    ierr = bed_topography.create(grid, "topg", WITH_GHOSTS, WIDE_STENCIL); CHKERRQ(ierr);
    ierr = bed_topography.set_attrs("model_state", "bedrock surface elevation",
          "m", "bedrock_altitude"); CHKERRQ(ierr);
    ierr = vars.add(bed_topography); CHKERRQ(ierr);

    // age of the ice; is not used here
    ierr = age.create(grid, "age", WITHOUT_GHOSTS); CHKERRQ(ierr);
    ierr = age.set_attrs("diagnostic", "age of the ice", "s", ""); CHKERRQ(ierr);
    ierr = age.set_glaciological_units("year"); CHKERRQ(ierr);
    age.write_in_glaciological_units = true;
    ierr = vars.add(age); CHKERRQ(ierr);

    // enthalpy in the ice
    ierr = enthalpy.create(grid, "enthalpy", WITH_GHOSTS, WIDE_STENCIL); CHKERRQ(ierr);
    ierr = enthalpy.set_attrs("model_state",
                              "ice enthalpy (includes sensible heat, latent heat, pressure)",
                              "J kg-1", ""); CHKERRQ(ierr);
    ierr = vars.add(enthalpy); CHKERRQ(ierr);

    // grounded_dragging_floating integer mask
    ierr = vMask.create(grid, "mask", WITH_GHOSTS, WIDE_STENCIL); CHKERRQ(ierr);
    ierr = vMask.set_attrs("model_state", "grounded_dragging_floating integer mask",
                           "", ""); CHKERRQ(ierr);
    std::vector<double> mask_values(4);
    mask_values[0] = MASK_ICE_FREE_BEDROCK;
    mask_values[1] = MASK_GROUNDED;
    mask_values[2] = MASK_FLOATING;
    mask_values[3] = MASK_ICE_FREE_OCEAN;
    vMask.metadata().set_doubles("flag_values", mask_values);
    vMask.metadata().set_string("flag_meanings",
                                "ice_free_bedrock grounded_ice floating_ice ice_free_ocean");
    ierr = vars.add(vMask); CHKERRQ(ierr);

    // Create the SIA solver object:

    // We use SIA_Nonsliding and not SIAFD here because we need the z-component
    // of the ice velocity, which is computed using incompressibility of ice in
    // StressBalance::compute_vertical_velocity().
    SIAFD *sia = new SIAFD(grid, EC, config);
    ZeroSliding *no_sliding = new ZeroSliding(grid, EC, config);

    StressBalance stress_balance(grid, no_sliding, sia, config);

    // fill the fields:
    ierr = setInitStateF(grid, EC,
                         &bed_topography, &vMask, &ice_surface_elevation, &ice_thickness,
                         &enthalpy); CHKERRQ(ierr);

    // Allocate the SIA solver:
    ierr = stress_balance.init(vars); CHKERRQ(ierr);

    IceModelVec2S melange_back_pressure;
    ierr = melange_back_pressure.create(grid, "melange_back_pressure", WITHOUT_GHOSTS); CHKERRQ(ierr);
    ierr = melange_back_pressure.set_attrs("boundary_condition",
                                           "melange back pressure fraction", "", ""); CHKERRQ(ierr);
    ierr = melange_back_pressure.set(0.0); CHKERRQ(ierr);

    // Solve (fast==true means "no 3D update and no strain heating computation"):
    bool fast = false;
    ierr = stress_balance.update(fast, 0.0, melange_back_pressure); CHKERRQ(ierr);

    // Report errors relative to the exact solution:
    IceModelVec3 *u_sia, *v_sia, *w_sia, *sigma;
    ierr = stress_balance.get_3d_velocity(u_sia, v_sia, w_sia); CHKERRQ(ierr);

    ierr = stress_balance.get_volumetric_strain_heating(sigma); CHKERRQ(ierr);

    ierr = reportErrors(config, grid,
                        &ice_thickness, u_sia, v_sia, w_sia, sigma); CHKERRQ(ierr);

    // Write results to an output file:
    PIO pio(grid, "guess_mode");

    ierr = pio.open(output_file, PISM_READWRITE_MOVE); CHKERRQ(ierr);
    ierr = pio.def_time(config.get_string("time_dimension_name"),
                        grid.time->calendar(),
                        grid.time->CF_units_string()); CHKERRQ(ierr);
    ierr = pio.append_time(config.get_string("time_dimension_name"), 0.0);
    ierr = pio.close(); CHKERRQ(ierr);

    ierr = ice_surface_elevation.write(output_file); CHKERRQ(ierr);
    ierr = ice_thickness.write(output_file); CHKERRQ(ierr);
    ierr = vMask.write(output_file); CHKERRQ(ierr);
    ierr = bed_topography.write(output_file); CHKERRQ(ierr);
    // ierr = enthalpy.write(output_file); CHKERRQ(ierr);
    ierr = u_sia->write(output_file); CHKERRQ(ierr);
    ierr = v_sia->write(output_file); CHKERRQ(ierr);
    ierr = w_sia->write(output_file); CHKERRQ(ierr);
    ierr = sigma->write(output_file); CHKERRQ(ierr);
  }
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
