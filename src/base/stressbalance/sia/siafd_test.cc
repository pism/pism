// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016 Ed Bueler and Constantine Khroulev
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

#include "SIAFD.hh"
#include "base/basalstrength/basal_resistance.hh"
#include "base/enthalpyConverter.hh"
#include "base/rheology/PatersonBuddCold.hh"
#include "base/stressbalance/PISMStressBalance.hh"
#include "base/stressbalance/SSB_Modifier.hh"
#include "base/stressbalance/ShallowStressBalance.hh"
#include "base/util/IceGrid.hh"
#include "base/util/Mask.hh"
#include "base/util/Context.hh"
#include "base/util/PISMTime.hh"
#include "base/util/PISMVars.hh"
#include "base/util/VariableMetadata.hh"
#include "base/util/error_handling.hh"
#include "base/util/iceModelVec.hh"
#include "base/util/io/PIO.hh"
#include "base/util/petscwrappers/PetscInitializer.hh"
#include "base/util/pism_const.hh"
#include "base/util/pism_options.hh"
#include "verif/tests/exactTestsFG.h"
#include "base/util/io/io_helpers.hh"
#include "base/util/pism_utilities.hh"

namespace pism {

static void compute_strain_heating_errors(const IceModelVec3 &strain_heating,
                                          const IceModelVec2S &thickness,
                                          const IceGrid &grid,
                                          double &gmax_strain_heating_err,
                                          double &gav_strain_heating_err) {
  double    max_strain_heating_error = 0.0, av_strain_heating_error = 0.0, avcount = 0.0;
  const int Mz = grid.Mz();

  const double LforFG = 750000; // m

  const double
    ice_rho   = grid.ctx()->config()->get_double("ice_density"),
    ice_c     = grid.ctx()->config()->get_double("ice_specific_heat_capacity");

  double   junk0, junk1;

  std::vector<double> strain_heating_exact(Mz), dummy1(Mz), dummy2(Mz), dummy3(Mz), dummy4(Mz);

  const double *strain_heating_ij;

  IceModelVec::AccessList list;
  list.add(thickness);
  list.add(strain_heating);

  ParallelSection loop(grid.com);
  try {
    for (Points p(grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double
        xx = grid.x(i),
        yy = grid.y(j),
        r  = sqrt(PetscSqr(xx) + PetscSqr(yy));

      if ((r >= 1.0) && (r <= LforFG - 1.0)) {
        // only evaluate error if inside sheet and not at central
        // singularity
        bothexact(0.0, r, &(grid.z()[0]), Mz, 0.0,
                  &junk0, &junk1, &dummy1[0], &dummy2[0], &dummy3[0],
                  &strain_heating_exact[0], &dummy4[0]);

        for (int k = 0; k < Mz; k++) {
          strain_heating_exact[k] *= ice_rho * ice_c; // scale exact strain_heating to J/(s m^3)
        }
        const int ks = grid.kBelowHeight(thickness(i, j));
        strain_heating_ij = strain_heating.get_column(i, j);
        for (int k = 0; k < ks; k++) {  // only eval error if below num surface
          const double _strain_heating_error = fabs(strain_heating_ij[k] - strain_heating_exact[k]);
          max_strain_heating_error = std::max(max_strain_heating_error, _strain_heating_error);
          avcount += 1.0;
          av_strain_heating_error += _strain_heating_error;
        }
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  gmax_strain_heating_err = GlobalMax(grid.com, max_strain_heating_error);
  gav_strain_heating_err = GlobalSum(grid.com, av_strain_heating_error);
  double gavcount = GlobalSum(grid.com, avcount);
  gav_strain_heating_err = gav_strain_heating_err/std::max(gavcount,1.0);  // avoid div by zero
}


static void computeSurfaceVelocityErrors(const IceGrid &grid,
                                         const IceModelVec2S &ice_thickness,
                                         const IceModelVec3 &u3,
                                         const IceModelVec3 &v3,
                                         const IceModelVec3 &w3,
                                         double &gmaxUerr,
                                         double &gavUerr,
                                         double &gmaxWerr,
                                         double &gavWerr) {
  double    maxUerr = 0.0, maxWerr = 0.0, avUerr = 0.0, avWerr = 0.0;

  const double LforFG = 750000; // m

  IceModelVec::AccessList list;
  list.add(ice_thickness);
  list.add(u3);
  list.add(v3);
  list.add(w3);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double xx = grid.x(i), yy = grid.y(j),
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
      maxUerr = std::max(maxUerr,Uerr);
      avUerr += Uerr;
      const double Werr = fabs(w3.getValZ(i,j,ice_thickness(i,j)) - wex);
      maxWerr = std::max(maxWerr,Werr);
      avWerr += Werr;
    }
  }

  gmaxUerr = GlobalMax(grid.com, maxUerr);
  gmaxWerr = GlobalMax(grid.com, maxWerr);
  gavUerr = GlobalSum(grid.com, avUerr);
  gavUerr = gavUerr/(grid.Mx()*grid.My());
  gavWerr = GlobalSum(grid.com, avWerr);
  gavWerr = gavWerr/(grid.Mx()*grid.My());
}


static void enthalpy_from_temperature_cold(EnthalpyConverter &EC,
                                           IceGrid &grid,
                                           IceModelVec2S &thickness,
                                           IceModelVec3 &temperature,
                                           IceModelVec3 &enthalpy) {

  IceModelVec::AccessList list;
  list.add(temperature);
  list.add(enthalpy);
  list.add(thickness);

  double *T_ij, *E_ij; // columns of these values
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    T_ij = temperature.get_column(i,j);
    E_ij = enthalpy.get_column(i,j);

    for (unsigned int k=0; k<grid.Mz(); ++k) {
      double depth = thickness(i,j) - grid.z(k);
      E_ij[k] = EC.enthalpy_permissive(T_ij[k], 0.0,
                                     EC.pressure(depth));
    }

  }

  enthalpy.update_ghosts();
}


//! \brief Set the test F initial state.
static void setInitStateF(IceGrid &grid,
                          EnthalpyConverter &EC,
                          IceModelVec2S *bed,
                          IceModelVec2Int *mask,
                          IceModelVec2S *surface,
                          IceModelVec2S *thickness,
                          IceModelVec3 *enthalpy) {
  int Mz = grid.Mz();

  double ST = 1.67e-5,
    Tmin = 223.15,  // K
    LforFG = 750000; // m

  std::vector<double> dummy1(Mz), dummy2(Mz), dummy3(Mz), dummy4(Mz), dummy5(Mz);

  bed->set(0);
  mask->set(MASK_GROUNDED);

  std::vector<double> T(Mz);

  IceModelVec::AccessList list;
  list.add(*thickness);
  list.add(*enthalpy);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double r = radius(grid, i, j),
      Ts = Tmin + ST * r;

    if (r > LforFG - 1.0) { // if (essentially) outside of sheet
      (*thickness)(i, j) = 0.0;
      for (int k = 0; k < Mz; k++) {
        T[k]=Ts;
      }
    } else {
      r = std::max(r, 1.0); // avoid singularity at origin
      bothexact(0.0, r, &(grid.z()[0]), Mz, 0.0,
                &(*thickness)(i, j), &dummy5[0], &T[0], &dummy1[0], &dummy2[0], &dummy3[0], &dummy4[0]);
    }
    enthalpy->set_column(i, j, &T[0]);
  }


  thickness->update_ghosts();

  surface->copy_from(*thickness);

  enthalpy_from_temperature_cold(EC, grid, *thickness,
                                 *enthalpy,
                                 *enthalpy);
}

static void reportErrors(const IceGrid &grid,
                         units::System::Ptr unit_system,
                         const IceModelVec2S &thickness,
                         const IceModelVec3 &u_sia,
                         const IceModelVec3 &v_sia,
                         const IceModelVec3 &w_sia,
                         const IceModelVec3 &strain_heating) {

  Logger::ConstPtr log = grid.ctx()->log();

  // strain_heating errors if appropriate; reported in 10^6 J/(s m^3)
  double max_strain_heating_error, av_strain_heating_error;
  compute_strain_heating_errors(strain_heating, thickness,
                                grid,
                                max_strain_heating_error, av_strain_heating_error);

  log->message(1,
               "Sigma     :      maxSig       avSig\n");
  log->message(1, "           %12.6f%12.6f\n",
               max_strain_heating_error*1.0e6, av_strain_heating_error*1.0e6);

  // surface velocity errors if exact values are available; reported in m year-1
  double maxUerr, avUerr, maxWerr, avWerr;
  computeSurfaceVelocityErrors(grid, thickness,
                               u_sia,
                               v_sia,
                               w_sia,
                               maxUerr, avUerr,
                               maxWerr, avWerr);

  log->message(1,
               "surf vels :     maxUvec      avUvec        maxW         avW\n");
  log->message(1, "           %12.6f%12.6f%12.6f%12.6f\n",
               units::convert(unit_system, maxUerr, "m second-1", "m year-1"),
               units::convert(unit_system, avUerr,  "m second-1", "m year-1"),
               units::convert(unit_system, maxWerr, "m second-1", "m year-1"),
               units::convert(unit_system, avWerr,  "m second-1", "m year-1"));
}

} // end of namespace pism

int main(int argc, char *argv[]) {

  using namespace pism;
  using namespace pism::stressbalance;

  MPI_Comm com = MPI_COMM_WORLD;
  petsc::Initializer petsc(argc, argv, help);
  PetscErrorCode ierr;

  com = PETSC_COMM_WORLD;

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  try {
    units::System::Ptr unit_system(new units::System);
    Context::Ptr ctx = context_from_options(com, "siafd_test");
    Config::Ptr config = ctx->config();

    config->set_boolean("compute_grain_size_using_age", false);

    bool
      usage_set = options::Bool("-usage", "print usage info"),
      help_set  = options::Bool("-help", "print help info");
    if (usage_set or help_set) {
      ierr = PetscPrintf(com,
                         "\n"
                         "usage of SIAFD_TEST:\n"
                         "  run siafd_test -Mx <number> -My <number> -Mz <number> -o foo.nc\n"
                         "\n");
      PISM_CHK(ierr, "PetscPrintf");
    }

    options::String output_file("-o", "Set the output file name", "siafd_test_F.nc");

    GridParameters P(config);
    P.Lx = 900e3;
    P.Ly = P.Lx;
    P.horizontal_size_from_options();

    double Lz = 4000.0;
    options::Integer Mz("-Mz", "Number of vertical grid levels", 61);

    P.z = IceGrid::compute_vertical_levels(Lz, Mz, EQUAL);
    P.ownership_ranges_from_options(ctx->size());

    // create grid and set defaults
    IceGrid::Ptr grid(new IceGrid(ctx, P));
    grid->report_parameters();

    setVerbosityLevel(5);

    EnthalpyConverter::Ptr EC(new ColdEnthalpyConverter(*config));
    rheology::PatersonBuddCold ice("sia_", *config, EC);

    IceModelVec2S ice_surface_elevation, ice_thickness, bed_topography;
    IceModelVec2Int cell_type;
    IceModelVec3 enthalpy,
      age;                      // is not used (and need not be allocated)
    const int WIDE_STENCIL = config->get_double("grid_max_stencil_width");

    Vars &vars = grid->variables();

    bed_topography.create(grid, "topg", WITH_GHOSTS, WIDE_STENCIL);
    bed_topography.set_attrs("model_state", "bedrock surface elevation",
                             "m", "bedrock_altitude");
    vars.add(bed_topography);

    // ice upper surface elevation
    ice_surface_elevation.create(grid, "usurf", WITH_GHOSTS, WIDE_STENCIL);
    ice_surface_elevation.set_attrs("diagnostic", "ice upper surface elevation",
                                    "m", "surface_altitude");
    vars.add(ice_surface_elevation);

    // land ice thickness
    ice_thickness.create(grid, "thk", WITH_GHOSTS, WIDE_STENCIL);
    ice_thickness.set_attrs("model_state", "land ice thickness",
                            "m", "land_ice_thickness");
    ice_thickness.metadata().set_double("valid_min", 0.0);
    vars.add(ice_thickness);

    // age of the ice; is not used here
    age.create(grid, "age", WITHOUT_GHOSTS);
    age.set_attrs("diagnostic", "age of the ice", "s", "");
    age.metadata().set_string("glaciological_units", "year");
    age.write_in_glaciological_units = true;
    vars.add(age);

    // enthalpy in the ice
    enthalpy.create(grid, "enthalpy", WITH_GHOSTS, WIDE_STENCIL);
    enthalpy.set_attrs("model_state",
                       "ice enthalpy (includes sensible heat, latent heat, pressure)",
                       "J kg-1", "");
    vars.add(enthalpy);

    // grounded_dragging_floating integer mask
    cell_type.create(grid, "mask", WITH_GHOSTS, WIDE_STENCIL);
    cell_type.set_attrs("model_state", "grounded_dragging_floating integer mask",
                    "", "");
    std::vector<double> mask_values(4);
    mask_values[0] = MASK_ICE_FREE_BEDROCK;
    mask_values[1] = MASK_GROUNDED;
    mask_values[2] = MASK_FLOATING;
    mask_values[3] = MASK_ICE_FREE_OCEAN;
    cell_type.metadata().set_doubles("flag_values", mask_values);
    cell_type.metadata().set_string("flag_meanings",
                                "ice_free_bedrock grounded_ice floating_ice ice_free_ocean");
    vars.add(cell_type);

    // Create the SIA solver object:

    // We use SIA_Nonsliding and not SIAFD here because we need the z-component
    // of the ice velocity, which is computed using incompressibility of ice in
    // StressBalance::compute_vertical_velocity().
    SIAFD *sia = new SIAFD(grid, EC);
    ZeroSliding *no_sliding = new ZeroSliding(grid, EC);

    StressBalance stress_balance(grid, no_sliding, sia);

    // fill the fields:
    setInitStateF(*grid, *EC,
                  &bed_topography, &cell_type, &ice_surface_elevation, &ice_thickness,
                  &enthalpy);

    // Allocate the SIA solver:
    stress_balance.init();

    IceModelVec2S melange_back_pressure;
    melange_back_pressure.create(grid, "melange_back_pressure", WITHOUT_GHOSTS);
    melange_back_pressure.set_attrs("boundary_condition",
                                    "melange back pressure fraction", "", "");
    melange_back_pressure.set(0.0);

    // Solve (fast==true means "no 3D update and no strain heating computation"):
    bool fast = false;
    stress_balance.update(fast, 0.0, melange_back_pressure);

    // Report errors relative to the exact solution:
    const IceModelVec3
      &u3 = stress_balance.velocity_u(),
      &v3 = stress_balance.velocity_v(),
      &w3 = stress_balance.velocity_w();

    const IceModelVec3 &sigma = stress_balance.volumetric_strain_heating();

    reportErrors(*grid, unit_system,
                 ice_thickness, u3, v3, w3, sigma);

    // Write results to an output file:
    PIO pio(grid->com, "netcdf3");

    pio.open(output_file, PISM_READWRITE_MOVE);
    io::define_time(pio, config->get_string("time_dimension_name"),
                    grid->ctx()->time()->calendar(),
                    grid->ctx()->time()->CF_units_string(),
                    unit_system);
    io::append_time(pio, config->get_string("time_dimension_name"), 0.0);
    pio.close();

    ice_surface_elevation.write(output_file);
    ice_thickness.write(output_file);
    cell_type.write(output_file);
    bed_topography.write(output_file);
    
    u3.write(output_file);
    v3.write(output_file);
    w3.write(output_file);
    sigma.write(output_file);
  }
  catch (...) {
    handle_fatal_errors(com);
    return 1;
  }

  return 0;
}
