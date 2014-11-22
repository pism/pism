// Copyright (C) 2011, 2012, 2013, 2014 Ed Bueler and Constantine Khroulev
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
  "Tests BedThermalUnit using Test K.  Sans IceModel.\n\n";

#include "pism_options.hh"
#include "IceGrid.hh"
#include "PIO.hh"
#include "NCVariable.hh"
#include "bedrockThermalUnit.hh"
#include "PISMTime.hh"
#include "PISMVars.hh"
#include "PISMConfig.hh"

#include "../../verif/tests/exactTestK.h"

#include "PetscInitializer.hh"
#include "error_handling.hh"

using namespace pism;

class BTU_Test : public BedThermalUnit
{
public:
  BTU_Test(IceGrid &g, const Config &conf)
    : BedThermalUnit(g, conf) {}
  virtual ~BTU_Test() {}
  /** Initialize the bedrock temperature field at the beginning of the run. */
  virtual void bootstrap();
};

void BTU_Test::bootstrap() {

  // fill exact bedrock temperature from Test K at time ys
  if (m_Mbz > 1) {
    std::vector<double> zlevels = temp.get_levels();

    IceModelVec::AccessList list(temp);
    double *Tb; // columns of these values
    for (Points p(grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      temp.getInternalColumn(i,j,&Tb);
      for (unsigned int k=0; k < m_Mbz; k++) {
        const double z = zlevels[k];
        double FF; // Test K:  use Tb[k], ignore FF
        exactK(grid.time->start(), z, &Tb[k], &FF, 0);
      }
    }
  }
}


static PetscErrorCode createVecs(IceGrid &grid, Vars &variables) {

  PetscErrorCode ierr;
  IceModelVec2S *bedtoptemp = new IceModelVec2S,
                *ghf        = new IceModelVec2S;

  ghf->create(grid, "bheatflx", WITHOUT_GHOSTS);
  ghf->set_attrs("",
                 "upward geothermal flux at bedrock thermal layer base",
                 "W m-2", "");
  ghf->set_glaciological_units("mW m-2");
  variables.add(*ghf);

  bedtoptemp->create(grid, "bedtoptemp", WITHOUT_GHOSTS);
  bedtoptemp->set_attrs("",
                        "temperature at top of bedrock thermal layer",
                        "K", "");
  variables.add(*bedtoptemp);

  return 0;
}


static PetscErrorCode doneWithIceInfo(Vars &variables) {
  // Get the names of all the variables allocated earlier:
  std::set<std::string> vars = variables.keys();
  std::set<std::string>::iterator i = vars.begin();
  while (i != vars.end()) {
    IceModelVec *var = variables.get(*i);
    delete var;
    i++;
  }
  return 0;
}


int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm com = MPI_COMM_WORLD;
  PetscInitializer petsc(argc, argv, help);

  com = PETSC_COMM_WORLD;

  try {
    UnitSystem unit_system;
    Config config(com, "pism_config", unit_system),
      overrides(com, "pism_overrides", unit_system);

    verbosityLevelFromOptions();
    verbPrintf(2,com, "BTUTEST %s (test program for BedThermalUnit)\n",
               PISM_Revision);
    stop_on_version_option();

    // check required options
    std::vector<std::string> required;
    required.push_back("-Mbz");
    show_usage_check_req_opts(com, "btutest", required,
                              "  btutest -Mbz NN -Lbz 1000.0 [-o OUT.nc -ys A -ye B -dt C -Mz D -Lz E]\n"
                              "where these are required because they are used in BedThermalUnit:\n"
                              "  -Mbz           number of bedrock thermal layer levels to use\n"
                              "  -Lbz 1000.0    depth of bedrock thermal layer (required; Lbz=1000.0 m in Test K)\n"
                              "and these are allowed:\n"
                              "  -o             output file name; NetCDF format\n"
                              "  -ys            start year in using Test K\n"
                              "  -ye            end year in using Test K\n"
                              "  -dt            time step B (= positive float) in years\n"
                              "  -Mz            number of ice levels to use\n"
                              "  -Lz            height of ice/atmospher box\n");

    verbPrintf(2,com,
               "btutest tests BedThermalUnit and IceModelVec3BTU\n");

    // read the config option database:
    init_config(com, config, overrides);
    config.set_string("calendar", "none");

    // when IceGrid constructor is called, these settings are used
    config.set_string("grid_ice_vertical_spacing","equal");
    config.set_string("grid_bed_vertical_spacing","equal");
    config.set_double("start_year", 0.0);
    config.set_double("run_length_years", 1.0);

    // create grid and set defaults
    IceGrid grid(com, config);
    grid.Mz = 41;
    grid.Lz = 4000.0;
    grid.Mx = 3;
    grid.My = 3;
    grid.Lx = 1500e3;
    grid.Ly = grid.Lx;

    // Mbz and Lbz are used by the BedThermalUnit, not by IceGrid
    config.set_double("grid_Mbz", 11); 
    config.set_double("grid_Lbz", 1000); 

    verbPrintf(2,com,
               "  initializing IceGrid from options ...\n");
    bool flag;
    double dt_years = 1.0;
    std::string outname="unnamed_btutest.nc";
    int tmp = grid.Mz;
    ierr = PetscOptionsBegin(grid.com, "", "BTU_TEST options", "");
    PISM_PETSC_CHK(ierr, "PetscOptionsBegin");
    {
      OptionsString("-o", "Output file name", outname, flag);
      OptionsReal("-dt", "Time-step, in years", dt_years, flag);
      OptionsInt("-Mz", "number of vertical layers in ice", tmp, flag);
      OptionsReal("-Lz", "height of ice/atmosphere boxr", grid.Lz, flag);
    }
    ierr = PetscOptionsEnd();
    PISM_PETSC_CHK(ierr, "PetscOptionsEnd");

    if (tmp > 0) {
      grid.Mz = tmp;
    } else {
      throw RuntimeError::formatted("-Mz %d is invalid (has to be positive)", tmp);
    }

    // complete grid initialization based on user options
    grid.compute_nprocs();
    grid.compute_ownership_ranges();
    grid.compute_horizontal_spacing();
    grid.compute_vertical_levels();
    grid.time->init();
    grid.allocate();

    // allocate tools and IceModelVecs
    Vars variables;
    createVecs(grid, variables);

    // these vars are owned by this driver, outside of BedThermalUnit
    IceModelVec2S *bedtoptemp, *ghf;

    // top of bedrock layer temperature; filled from Test K exact values
    bedtoptemp = variables.get_2d_scalar("bedtoptemp");
    // lithosphere (bottom of bedrock layer) heat flux; has constant value
    ghf = variables.get_2d_scalar("bheatflx");

    ghf->set(0.042);  // see Test K

    // initialize BTU object:
    BTU_Test btu(grid, config);

    bool bootstrapping_needed = true; // we know it's true
    btu.init(variables, bootstrapping_needed);
    btu.bootstrap();

    double dt_seconds = unit_system.convert(dt_years, "years", "seconds");

    // worry about time step
    int  N = (int)ceil((grid.time->end() - grid.time->start()) / dt_seconds);
    dt_seconds = (grid.time->end() - grid.time->start()) / (double)N;
    verbPrintf(2,com,
               "  user set timestep of %.4f years ...\n"
               "  reset to %.4f years to get integer number of steps ... \n",
               dt_years, unit_system.convert(dt_seconds, "seconds", "years"));
    double max_dt;
    bool restrict_dt;
    btu.max_timestep(0.0, max_dt, restrict_dt);
    verbPrintf(2,com,
               "  BedThermalUnit reports max timestep of %.4f years ...\n",
               unit_system.convert(max_dt, "seconds", "years"));


    // actually do the time-stepping
    verbPrintf(2,com,"  running ...\n");
    for (int n = 0; n < N; n++) {
      const double time = grid.time->start() + dt_seconds * (double)n;  // time at start of time-step

      // compute exact ice temperature at z=0 at time y
      IceModelVec::AccessList list(*bedtoptemp);
      for (Points p(grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        double TT, FF; // Test K:  use TT, ignore FF
        exactK(time, 0.0, &TT, &FF, 0);
        (*bedtoptemp)(i,j) = TT;
      }
      // we are not communicating anything, which is fine

      // update the temperature inside the thermal layer using bedtoptemp
      btu.update(time, dt_seconds);
      verbPrintf(2,com,".");
    }

    verbPrintf(2,com,"\n  done ...\n");

    // compute final output heat flux G_0 at z=0; reuse ghf for this purpose
    ghf->set_name("bheatflx0");
    ghf->set_attrs("",
                   "upward geothermal flux at ice/bedrock interface",
                   "W m-2", "");
    btu.get_upward_geothermal_flux(*ghf);

    // get, and tell stdout, the correct answer from Test K
    double TT, FF; // Test K:  use FF, ignore TT
    exactK(grid.time->end(), 0.0, &TT, &FF, 0);
    verbPrintf(2,com,
               "  exact Test K reports upward heat flux at z=0, at end time %s, as G_0 = %.7f W m-2;\n",
               grid.time->end_date().c_str(), FF);

    // compute numerical error
    double maxghferr, avghferr;
    ghf->shift(-FF);
    ghf->norm(NORM_INFINITY,maxghferr);
    ghf->norm(NORM_1,avghferr);
    ghf->shift(+FF); // shift it back for writing
    avghferr /= (grid.Mx * grid.My);
    verbPrintf(2,grid.com, 
               "case dt = %.5f:\n", dt_years);
    verbPrintf(1,grid.com, 
               "NUMERICAL ERRORS in upward heat flux at z=0 relative to exact solution:\n");
    verbPrintf(1,grid.com, 
               "bheatflx0  :       max    prcntmax          av\n");
    verbPrintf(1,grid.com, 
               "           %11.7f  %11.7f  %11.7f\n", 
               maxghferr,100.0*maxghferr/FF,avghferr);
    verbPrintf(1,grid.com, "NUM ERRORS DONE\n");

    std::set<std::string> vars;
    btu.add_vars_to_output("big", vars); // "write everything you can"

    PIO pio(grid, grid.config.get_string("output_format"));

    std::string time_name = config.get_string("time_dimension_name");
    pio.open(outname, PISM_READWRITE_MOVE);
    pio.def_time(time_name, grid.time->calendar(),
                 grid.time->CF_units_string());
    pio.append_time(time_name, grid.time->end());

    btu.define_variables(vars, pio, PISM_DOUBLE);
    btu.write_variables(vars, pio);

    bedtoptemp->write(pio);
    ghf->write(pio);

    pio.close();

    doneWithIceInfo(variables);
    verbPrintf(2,com, "done.\n");
  }
  catch (...) {
    handle_fatal_errors(com);
  }

  return 0;
}

