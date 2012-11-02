// Copyright (C) 2011, 2012 Ed Bueler and Constantine Khroulev
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

static char help[] =
  "Tests PISMBedThermalUnit using Test K.  Sans IceModel.\n\n";

#include "pism_options.hh"
#include "IceGrid.hh"
#include "PIO.hh"
#include "NCVariable.hh"
#include "bedrockThermalUnit.hh"
#include "PISMTime.hh"
#include "PISMVars.hh"

#include "../../verif/tests/exactTestK.h"

class BTU_Test : public PISMBedThermalUnit
{
public:
  BTU_Test(IceGrid &g, const NCConfigVariable &conf)
    : PISMBedThermalUnit(g, conf) {}
  virtual ~BTU_Test() {}
protected:
  virtual PetscErrorCode bootstrap();
};

PetscErrorCode BTU_Test::bootstrap() {
  PetscErrorCode ierr;

  // fill exact bedrock temperature from Test K at time ys
  if (Mbz > 1) {
    vector<double> zlevels = temp.get_levels();

    ierr = temp.begin_access(); CHKERRQ(ierr);
    PetscScalar *Tb; // columns of these values
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        ierr = temp.getInternalColumn(i,j,&Tb); CHKERRQ(ierr);
        for (PetscInt k=0; k < Mbz; k++) {
          const PetscReal z = zlevels[k];
          PetscReal FF; // Test K:  use Tb[k], ignore FF
          ierr = exactK(grid.time->start(), z, &Tb[k], &FF, 0); CHKERRQ(ierr);
        }
      }
    }
    ierr = temp.end_access(); CHKERRQ(ierr);
  }

  return 0;
}


static PetscErrorCode createVecs(IceGrid &grid, PISMVars &variables) {

  PetscErrorCode ierr;
  IceModelVec2S *bedtoptemp = new IceModelVec2S,
                *ghf        = new IceModelVec2S;

  ierr = ghf->create(grid, "bheatflx", false); CHKERRQ(ierr);
  ierr = ghf->set_attrs("",
                       "upward geothermal flux at bedrock thermal layer base",
		       "W m-2", ""); CHKERRQ(ierr);
  ierr = ghf->set_glaciological_units("mW m-2");
  ierr = variables.add(*ghf); CHKERRQ(ierr);

  ierr = bedtoptemp->create(grid, "bedtoptemp", false); CHKERRQ(ierr);
  ierr = bedtoptemp->set_attrs("",
                       "temperature at top of bedrock thermal layer",
		       "K", ""); CHKERRQ(ierr);
  ierr = variables.add(*bedtoptemp); CHKERRQ(ierr);

  return 0;
}


static PetscErrorCode doneWithIceInfo(PISMVars &variables) {
  // Get the names of all the variables allocated earlier:
  set<string> vars = variables.keys();
  set<string>::iterator i = vars.begin();
  while (i != vars.end()) {
    IceModelVec *var = variables.get(*i);
    delete var;
    i++;
  }
  return 0;
}


int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm    com;
  PetscMPIInt rank, size;
  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);
  com = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  {
    NCConfigVariable config, overrides;

    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);
    ierr = verbPrintf(2,com, "BTUTEST %s (test program for PISMBedThermalUnit)\n",
		      PISM_Revision); CHKERRQ(ierr);
    ierr = stop_on_version_option(); CHKERRQ(ierr);

    // check required options
    vector<string> required;
    required.push_back("-Mbz");
    ierr = show_usage_check_req_opts(com, "btutest", required,
      "  btutest -Mbz NN -Lbz 1000.0 [-o OUT.nc -ys A -ye B -dt C -Mz D -Lz E]\n"
      "where these are required because they are used in PISMBedThermalUnit:\n"
      "  -Mbz           number of bedrock thermal layer levels to use\n"
      "  -Lbz 1000.0    depth of bedrock thermal layer (required; Lbz=1000.0 m in Test K)\n"
      "and these are allowed:\n"
      "  -o             output file name; NetCDF format\n"
      "  -ys            start year in using Test K\n"
      "  -ye            end year in using Test K\n"
      "  -dt            time step B (= positive float) in years\n"
      "  -Mz            number of ice levels to use\n"
      "  -Lz            height of ice/atmospher box\n"
      ); CHKERRQ(ierr);

    ierr = verbPrintf(2,com,
        "btutest tests PISMBedThermalUnit and IceModelVec3BTU\n"); CHKERRQ(ierr);

    // read the config option database:
    ierr = init_config(com, rank, config, overrides); CHKERRQ(ierr);

    // when IceGrid constructor is called, these settings are used
    config.set_string("grid_ice_vertical_spacing","equal");
    config.set_string("grid_bed_vertical_spacing","equal");
    config.set("start_year", 0.0);
    config.set("run_length_years", 1.0);

    // create grid and set defaults
    IceGrid grid(com, rank, size, config);
    grid.Mz = 41;
    grid.Lz = 4000.0;
    grid.Mx = 3;
    grid.My = 3;
    grid.Lx = 1500e3;
    grid.Ly = grid.Lx;

    // Mbz and Lbz are used by the PISMBedThermalUnit, not by IceGrid
    config.set("grid_Mbz", 11); 
    config.set("grid_Lbz", 1000); 

    ierr = verbPrintf(2,com,
	"  initializing IceGrid from options ...\n"); CHKERRQ(ierr);
    bool flag;
    PetscReal dt_years = 1.0;
    string outname="unnamed_btutest.nc";
    ierr = PetscOptionsBegin(grid.com, "", "BTU_TEST options", ""); CHKERRQ(ierr);
    {
      ierr = PISMOptionsString("-o", "Output file name", outname, flag); CHKERRQ(ierr);
      ierr = PISMOptionsReal("-dt", "Time-step, in years", dt_years, flag); CHKERRQ(ierr);
      ierr = PISMOptionsInt("-Mz", "number of vertical layers in ice", grid.Mz, flag); CHKERRQ(ierr);
      ierr = PISMOptionsReal("-Lz", "height of ice/atmosphere boxr", grid.Lz, flag); CHKERRQ(ierr);
    }
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    // complete grid initialization based on user options
    grid.compute_nprocs();
    grid.compute_ownership_ranges();
    ierr = grid.compute_horizontal_spacing(); CHKERRQ(ierr);
    ierr = grid.compute_vertical_levels(); CHKERRQ(ierr);
    ierr = grid.time->init(); CHKERRQ(ierr);
    ierr = grid.createDA(); CHKERRQ(ierr);

    // allocate tools and IceModelVecs
    PISMVars variables;
    ierr = createVecs(grid, variables); CHKERRQ(ierr);

    // these vars are owned by this driver, outside of PISMBedThermalUnit
    IceModelVec2S *bedtoptemp, *ghf;

    // top of bedrock layer temperature; filled from Test K exact values
    bedtoptemp = dynamic_cast<IceModelVec2S*>(variables.get("bedtoptemp"));
    if (bedtoptemp == NULL) SETERRQ(com, 1, "bedtoptemp is not available");

    // lithosphere (bottom of bedrock layer) heat flux; has constant value
    ghf = dynamic_cast<IceModelVec2S*>(variables.get("bheatflx"));
    if (ghf == NULL) SETERRQ(com, 2, "bheatflx is not available");
    ierr = ghf->set(0.042); CHKERRQ(ierr);  // see Test K

    // initialize BTU object:
    BTU_Test btu(grid, config);

    ierr = btu.init(variables); CHKERRQ(ierr);  // FIXME: this is bootstrapping, really

    PetscReal dt_seconds = convert(dt_years, "years", "seconds");

    // worry about time step
    int  N = (int)ceil((grid.time->end() - grid.time->start()) / dt_seconds);
    dt_seconds = (grid.time->end() - grid.time->start()) / (double)N;
    ierr = verbPrintf(2,com,
                      "  user set timestep of %.4f years ...\n"
                      "  reset to %.4f years to get integer number of steps ... \n",
                      dt_years, convert(dt_seconds, "seconds", "years")); CHKERRQ(ierr);
    PetscReal max_dt;
    bool restrict_dt;
    ierr = btu.max_timestep(0.0, max_dt, restrict_dt); CHKERRQ(ierr);
    ierr = verbPrintf(2,com,
        "  PISMBedThermalUnit reports max timestep of %.4f years ...\n",
                      convert(max_dt, "seconds", "years")); CHKERRQ(ierr);


    // actually do the time-stepping
    ierr = verbPrintf(2,com,"  running ...\n  "); CHKERRQ(ierr);
    for (PetscInt n = 0; n < N; n++) {
      const PetscReal time = grid.time->start() + dt_seconds * (double)n;  // time at start of time-step

      // compute exact ice temperature at z=0 at time y
      ierr = bedtoptemp->begin_access(); CHKERRQ(ierr);
      for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
        for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
          PetscReal TT, FF; // Test K:  use TT, ignore FF
          ierr = exactK(time, 0.0, &TT, &FF, 0); CHKERRQ(ierr);
          (*bedtoptemp)(i,j) = TT;
        }
      }
      ierr = bedtoptemp->end_access(); CHKERRQ(ierr);
      // we are not communicating anything, which is fine

      // update the temperature inside the thermal layer using bedtoptemp
      ierr = btu.update(time, dt_seconds); CHKERRQ(ierr);
      ierr = verbPrintf(2,com,"."); CHKERRQ(ierr);
    }

    ierr = verbPrintf(2,com,"\n  done ...\n"); CHKERRQ(ierr);

    // compute final output heat flux G_0 at z=0; reuse ghf for this purpose
    ierr = ghf->set_name("bheatflx0"); CHKERRQ(ierr);
    ierr = ghf->set_attrs("",
                       "upward geothermal flux at ice/bedrock interface",
		       "W m-2", ""); CHKERRQ(ierr);
    ierr = btu.get_upward_geothermal_flux(*ghf); CHKERRQ(ierr);

    // get, and tell stdout, the correct answer from Test K
    PetscReal TT, FF; // Test K:  use FF, ignore TT
    ierr = exactK(grid.time->end(), 0.0, &TT, &FF, 0); CHKERRQ(ierr);
    ierr = verbPrintf(2,com,
        "  exact Test K reports upward heat flux at z=0, at end time %s, as G_0 = %.7f W m-2;\n",
                      grid.time->end_date().c_str(), FF); CHKERRQ(ierr);

    // compute numerical error
    PetscReal maxghferr, avghferr;
    ierr = ghf->shift(-FF); CHKERRQ(ierr);
    ierr = ghf->norm(NORM_INFINITY,maxghferr); CHKERRQ(ierr);
    ierr = ghf->norm(NORM_1,avghferr); CHKERRQ(ierr);
    ierr = ghf->shift(+FF); CHKERRQ(ierr); // shift it back for writing
    avghferr /= (grid.Mx * grid.My);
    ierr = verbPrintf(2,grid.com, 
                      "case dt = %.5f:\n", dt_years); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, 
                      "NUMERICAL ERRORS in upward heat flux at z=0 relative to exact solution:\n");
                      CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, 
                      "bheatflx0  :       max    prcntmax          av\n"); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, 
                      "           %11.7f  %11.7f  %11.7f\n", 
                      maxghferr,100.0*maxghferr/FF,avghferr); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, "NUM ERRORS DONE\n");  CHKERRQ(ierr);

    map<string, NCSpatialVariable> list;
    set<string> vars;
    btu.add_vars_to_output("big", list); // "write everything you can"

    map<string, NCSpatialVariable>::iterator j = list.begin();
    while (j != list.end()) {
      vars.insert(j->first);
      ++j;
    }

    PIO pio(grid.com, grid.rank, grid.config.get_string("output_format"));

    string time_name = config.get_string("time_dimension_name");
    ierr = pio.open(outname, PISM_WRITE); CHKERRQ(ierr);
    ierr = pio.def_time(time_name, config.get_string("calendar"),
                        grid.time->CF_units()); CHKERRQ(ierr);
    ierr = pio.append_time(time_name, grid.time->end()); CHKERRQ(ierr);

    ierr = btu.define_variables(vars, pio, PISM_DOUBLE); CHKERRQ(ierr);
    ierr = btu.write_variables(vars, pio); CHKERRQ(ierr);

    ierr = bedtoptemp->write(pio); CHKERRQ(ierr);
    ierr = ghf->write(pio); CHKERRQ(ierr);

    ierr = pio.close(); CHKERRQ(ierr);

    ierr = doneWithIceInfo(variables); CHKERRQ(ierr);
    ierr = verbPrintf(2,com, "done.\n"); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}

