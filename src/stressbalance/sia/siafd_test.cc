// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 Ed Bueler and Constantine Khroulev
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
#include "pism/basalstrength/basal_resistance.hh"
#include "pism/util/EnthalpyConverter.hh"
#include "pism/rheology/PatersonBuddCold.hh"
#include "pism/stressbalance/StressBalance.hh"
#include "pism/stressbalance/SSB_Modifier.hh"
#include "pism/stressbalance/ShallowStressBalance.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/Mask.hh"
#include "pism/util/Context.hh"
#include "pism/util/Time.hh"
#include "pism/util/VariableMetadata.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/io/PIO.hh"
#include "pism/util/petscwrappers/PetscInitializer.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/pism_options.hh"
#include "pism/verification/tests/exactTestsFG.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/geometry/Geometry.hh"

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
    ice_rho   = grid.ctx()->config()->get_double("constants.ice.density"),
    ice_c     = grid.ctx()->config()->get_double("constants.ice.specific_heat_capacity");

  IceModelVec::AccessList list{&thickness, &strain_heating};

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
        TestFGParameters F = exactFG(0.0, r, grid.z(), 0.0);

        for (int k = 0; k < Mz; k++) {
          F.Sig[k] *= ice_rho * ice_c; // scale exact strain_heating to J/(s m^3)
        }
        const int ks = grid.kBelowHeight(thickness(i, j));
        const double *strain_heating_ij = strain_heating.get_column(i, j);
        for (int k = 0; k < ks; k++) {  // only eval error if below num surface
          const double _strain_heating_error = fabs(strain_heating_ij[k] - F.Sig[k]);
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

  IceModelVec::AccessList list{&ice_thickness, &u3, &v3, &w3};

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double xx = grid.x(i), yy = grid.y(j),
      r = sqrt(PetscSqr(xx) + PetscSqr(yy));
    if ((r >= 1.0) && (r <= LforFG - 1.0)) {  // only evaluate error if inside sheet
      // and not at central singularity
      const double H = ice_thickness(i, j);
      std::vector<double> z(1, H);
      TestFGParameters F = exactFG(0.0, r, z, 0.0);

      const double uex = (xx/r) * F.U[0];
      const double vex = (yy/r) * F.U[0];
      // note that because getValZ does linear interpolation and H is not exactly at
      // a grid point, this causes nonzero errors
      const double Uerr = sqrt(PetscSqr(u3.getValZ(i,j,H) - uex) +
                               PetscSqr(v3.getValZ(i,j,H) - vex));
      maxUerr = std::max(maxUerr,Uerr);
      avUerr += Uerr;
      const double Werr = fabs(w3.getValZ(i,j,H) - F.w[0]);
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
                                           const IceGrid &grid,
                                           const IceModelVec2S &thickness,
                                           const IceModelVec3 &temperature,
                                           IceModelVec3 &enthalpy) {

  IceModelVec::AccessList list{&temperature, &enthalpy, &thickness};

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double *T_ij = temperature.get_column(i,j);
    double *E_ij = enthalpy.get_column(i,j);

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
                          IceModelVec2S &bed,
                          IceModelVec2Int &mask,
                          IceModelVec2S &surface,
                          IceModelVec2S &thickness,
                          IceModelVec3 &enthalpy) {

  double
    ST     = 1.67e-5,
    Tmin   = 223.15,            // K
    LforFG = 750000;            // m

  bed.set(0.0);
  mask.set(MASK_GROUNDED);

  IceModelVec::AccessList list{&thickness, &enthalpy};

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double
      r  = std::max(radius(grid, i, j), 1.0), // avoid singularity at origin
      Ts = Tmin + ST * r;

    if (r > LforFG - 1.0) { // if (essentially) outside of sheet
      thickness(i, j) = 0.0;
      enthalpy.set_column(i, j, Ts);
    } else {
      TestFGParameters F = exactFG(0.0, r, grid.z(), 0.0);

      thickness(i, j) = F.H;
      enthalpy.set_column(i, j, &F.T[0]);
    }
  }


  thickness.update_ghosts();

  surface.copy_from(thickness);

  enthalpy_from_temperature_cold(EC, grid, thickness, enthalpy,
                                 enthalpy);
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

  com = PETSC_COMM_WORLD;

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  try {
    // set default verbosity
    Context::Ptr ctx = context_from_options(com, "siafd_test");
    Config::Ptr config = ctx->config();

    config->set_boolean("stress_balance.sia.grain_size_age_coupling", false);
    config->set_string("stress_balance.sia.flow_law", "arr");

    set_config_from_options(*config);

    std::string usage = "\n"
      "usage of SIAFD_TEST:\n"
      "  run siafd_test -Mx <number> -My <number> -Mz <number> -o foo.nc\n"
      "\n";

    bool stop = show_usage_check_req_opts(*ctx->log(), "siafd_test", {}, usage);

    if (stop) {
      return 0;
    }

    auto output_file = config->get_string("output.file_name");

    GridParameters P(config);
    P.Lx = 900e3;
    P.Ly = P.Lx;
    P.horizontal_size_from_options();

    double Lz = 4000.0;
    unsigned int Mz = config->get_double("grid.Mz");

    P.z = IceGrid::compute_vertical_levels(Lz, Mz, EQUAL);
    P.ownership_ranges_from_options(ctx->size());

    // create grid and set defaults
    IceGrid::Ptr grid(new IceGrid(ctx, P));
    grid->report_parameters();

    EnthalpyConverter::Ptr EC(new ColdEnthalpyConverter(*config));

    const int WIDE_STENCIL = config->get_double("grid.max_stencil_width");

    IceModelVec3
      enthalpy(grid, "enthalpy", WITH_GHOSTS, WIDE_STENCIL),
      age(grid, "age", WITHOUT_GHOSTS);

    Geometry geometry(grid);
    geometry.sea_level_elevation.set(0.0);

    // age of the ice; is not used here
    age.set_attrs("diagnostic", "age of the ice", "s", "");
    age.set(0.0);

    // enthalpy in the ice
    enthalpy.set_attrs("model_state",
                       "ice enthalpy (includes sensible heat, latent heat, pressure)",
                       "J kg-1", "");
    //
    enthalpy.set(EC->enthalpy(263.15, 0.0, EC->pressure(1000.0)));

    // Create the SIA solver object:

    // We use SIA_Nonsliding and not SIAFD here because we need the z-component
    // of the ice velocity, which is computed using incompressibility of ice in
    // StressBalance::compute_vertical_velocity().
    SIAFD *sia = new SIAFD(grid);
    ZeroSliding *no_sliding = new ZeroSliding(grid);

    // stress_balance will de-allocate no_sliding and sia.
    StressBalance stress_balance(grid, no_sliding, sia);

    // fill the fields:
    setInitStateF(*grid, *EC,
                  geometry.bed_elevation,
                  geometry.cell_type,
                  geometry.ice_surface_elevation,
                  geometry.ice_thickness,
                  enthalpy);

    geometry.ensure_consistency(config->get_double("geometry.ice_free_thickness_standard"));

    // Initialize the SIA solver:
    stress_balance.init();

    IceModelVec2S melange_back_pressure;
    melange_back_pressure.create(grid, "melange_back_pressure", WITHOUT_GHOSTS);
    melange_back_pressure.set_attrs("boundary_condition",
                                    "melange back pressure fraction", "", "");
    melange_back_pressure.set(0.0);

    bool full_update = true;

    stressbalance::Inputs inputs;
    inputs.geometry              = &geometry;
    inputs.melange_back_pressure = &melange_back_pressure;
    inputs.enthalpy              = &enthalpy;
    inputs.age                   = &age;

    stress_balance.update(inputs, full_update);

    // Report errors relative to the exact solution:
    const IceModelVec3
      &u3 = stress_balance.velocity_u(),
      &v3 = stress_balance.velocity_v(),
      &w3 = stress_balance.velocity_w();

    const IceModelVec3 &sigma = stress_balance.volumetric_strain_heating();

    reportErrors(*grid, ctx->unit_system(),
                 geometry.ice_thickness, u3, v3, w3, sigma);

    // Write results to an output file:
    PIO file(grid->com, "netcdf3", output_file, PISM_READWRITE_MOVE);
    io::define_time(file, *ctx);
    io::append_time(file, *ctx->config(), ctx->time()->current());

    geometry.ice_surface_elevation.write(file);
    geometry.ice_thickness.write(file);
    geometry.cell_type.write(file);
    geometry.bed_elevation.write(file);

    sia->diffusivity().write(file);
    u3.write(file);
    v3.write(file);
    w3.write(file);
    sigma.write(file);
  }
  catch (...) {
    handle_fatal_errors(com);
    return 1;
  }

  return 0;
}
