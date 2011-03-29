// Copyright (C) 2011 Ed Bueler and Constantine Khroulev
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

#include "pism_const.hh"
#include "grid.hh"
#include "PISMIO.hh"
#include "NCVariable.hh"
#include "bedrockThermalUnit.hh"

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
          ierr = exactK(grid.start_year * secpera, z, &Tb[k], &FF, 0); CHKERRQ(ierr);
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

    // create grid and set defaults
    IceGrid grid(com, rank, size, config);
    grid.Mbz = 11;
    grid.Lbz = 1000.0;
    grid.Mz = 41;
    grid.Lz = 4000.0;
    grid.Mx = 3;
    grid.My = 3;
    grid.Lx = 1500e3;
    grid.Ly = grid.Lx;
    grid.start_year = 0.0;
    grid.end_year = 1.0;

    ierr = verbPrintf(2,com,
	"  initializing IceGrid from options ...\n"); CHKERRQ(ierr);
    bool flag;
    PetscReal dt_years = 1.0;
    string outname="unnamed_btutest.nc";
    ierr = PetscOptionsBegin(grid.com, "", "BTU_TEST options", ""); CHKERRQ(ierr);
    {
      ierr = PISMOptionsString("-o", "Output file name", outname, flag); CHKERRQ(ierr);
      ierr = PISMOptionsReal("-ys", "starting time in years", grid.start_year, flag); CHKERRQ(ierr);
      ierr = PISMOptionsReal("-ye", "starting time in years", grid.end_year, flag); CHKERRQ(ierr);
      ierr = PISMOptionsReal("-dt", "Time-step, in years", dt_years, flag); CHKERRQ(ierr);
      ierr = PISMOptionsInt("-Mz", "number of vertical layers in ice", grid.Mz, flag); CHKERRQ(ierr);
      ierr = PISMOptionsReal("-Lz", "height of ice/atmosphere boxr", grid.Lz, flag); CHKERRQ(ierr);
      // -Mbz is ALSO checked by IceModelVec3BTU
      ierr = PISMOptionsInt("-Mbz", "number of vertical layers in bedrock thermal layer", grid.Mbz, flag); CHKERRQ(ierr);
      // -Lbz is ALSO checked by IceModelVec3BTU
      ierr = PISMOptionsReal("-Lbz", "thickness of bedrock thermal layer", grid.Lbz, flag); CHKERRQ(ierr);
    }
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    // complete grid initialization based on user options
    grid.compute_nprocs();
    grid.compute_ownership_ranges();
    ierr = grid.compute_horizontal_spacing(); CHKERRQ(ierr);
    if (grid.Mbz < 1) grid.Mbz = 1;  // FIXME: just to protect current IceGrid init requirement; removable ...
    ierr = grid.compute_vertical_levels(); CHKERRQ(ierr);
    ierr = grid.createDA(); CHKERRQ(ierr);

    // allocate tools and IceModelVecs
    PISMVars variables;
    ierr = createVecs(grid, variables); CHKERRQ(ierr);

    // these vars are owned by this driver, outside of PISMBedThermalUnit
    IceModelVec2S *bedtoptemp, *ghf;

    // top of bedrock layer temperature; filled from Test K exact values
    bedtoptemp = dynamic_cast<IceModelVec2S*>(variables.get("bedtoptemp"));
    if (bedtoptemp == NULL) SETERRQ(1, "bedtoptemp is not available");

    // lithosphere (bottom of bedrock layer) heat flux; has constant value
    ghf = dynamic_cast<IceModelVec2S*>(variables.get("bheatflx"));
    if (ghf == NULL) SETERRQ(2, "bheatflx is not available");
    ierr = ghf->set(0.042); CHKERRQ(ierr);  // see Test K

    // initialize BTU object:
    BTU_Test btu(grid, config);

    ierr = btu.init(variables); CHKERRQ(ierr);  // FIXME: this is bootstrapping, really

    // worry about time step
    int  N = (int)ceil((grid.end_year - grid.start_year) / dt_years);
    PetscReal old_dt_years = dt_years;
    dt_years = (grid.end_year - grid.start_year) / (double)N;
    ierr = verbPrintf(2,com,
        "  user set timestep of %.4f years ...\n"
        "  reset to %.4f years to get integer number of steps ... \n",
	old_dt_years,dt_years); CHKERRQ(ierr);
    PetscReal max_dt_years;
    ierr = btu.max_timestep(0.0, max_dt_years); CHKERRQ(ierr);
    ierr = verbPrintf(2,com,
        "  PISMBedThermalUnit reports max timestep of %.4f years ...\n",
	max_dt_years); CHKERRQ(ierr);


    // actually do the time-stepping
    ierr = verbPrintf(2,com,"  running ...\n  "); CHKERRQ(ierr);
    for (PetscInt n = 0; n < N; n++) {
      const PetscReal y = grid.start_year + dt_years * (double)n;  // time at start of time-step

      // compute exact ice temperature at z=0 at time y
      ierr = bedtoptemp->begin_access(); CHKERRQ(ierr);
      for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
        for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
          PetscReal TT, FF; // Test K:  use TT, ignor FF
          ierr = exactK(y * secpera, 0.0, &TT, &FF, 0); CHKERRQ(ierr);
          (*bedtoptemp)(i,j) = TT;
        }
      }
      ierr = bedtoptemp->end_access(); CHKERRQ(ierr);
      // we are not communicating anything, which is fine

      // update the temperature inside the thermal layer using bedtoptemp
      ierr = btu.update(y, dt_years); CHKERRQ(ierr);
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
    PetscReal TT, FF; // Test K:  use FF, ignor TT
    ierr = exactK(grid.end_year * secpera, 0.0, &TT, &FF, 0); CHKERRQ(ierr);
    ierr = verbPrintf(2,com,
        "  exact Test K reports upward heat flux at z=0, at end year %.2f, as G_0 = %.7f W m-2;\n",
	grid.end_year,FF); CHKERRQ(ierr);

    // compute numerical error
    PetscReal maxghferr, avghferr;
    ierr = ghf->shift(-FF); CHKERRQ(ierr);
    ierr = ghf->norm(NORM_INFINITY,maxghferr); CHKERRQ(ierr);
    ierr = ghf->norm(NORM_1,avghferr); CHKERRQ(ierr);
    ierr = ghf->shift(+FF); CHKERRQ(ierr); // shift it back for writing
    avghferr /= (grid.Mx * grid.My);
    ierr = verbPrintf(2,grid.com, 
                      "case Mbz = %d, dt = %.5f:\n", grid.Mbz, dt_years); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, 
                      "NUMERICAL ERRORS in upward heat flux at z=0 relative to exact solution:\n");
                      CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, 
                      "bheatflx0  :       max    prcntmax          av\n"); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, 
                      "           %11.7f  %11.7f  %11.7f\n", 
                      maxghferr,100.0*maxghferr/FF,avghferr); CHKERRQ(ierr);
    ierr = verbPrintf(1,grid.com, "NUM ERRORS DONE\n");  CHKERRQ(ierr);

    set<string> vars;
    btu.add_vars_to_output("big", vars); // "write everything you can"

    PISMIO pio(&grid);

    ierr = pio.open_for_writing(outname, false, true); CHKERRQ(ierr);
    // append == true and check_dims == true
    ierr = pio.append_time(grid.end_year); CHKERRQ(ierr);
    ierr = btu.define_variables(vars, pio, NC_DOUBLE); CHKERRQ(ierr);
    ierr = pio.close(); CHKERRQ(ierr);

    ierr = btu.write_variables(vars, outname); CHKERRQ(ierr);
    ierr = bedtoptemp->write(outname.c_str()); CHKERRQ(ierr);
    ierr = ghf->write(outname.c_str()); CHKERRQ(ierr);

    ierr = doneWithIceInfo(variables); CHKERRQ(ierr);
    ierr = verbPrintf(2,com, "done.\n"); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}

