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
  "Tests PISMBedThermalUnit using Test K.  Sans IceModel.\n"
  "Requires an input file only to get 2D grid information, and otherwise ignors\n"
  "3D and thermodynamic information in the file.  Example with one time step:\n"
  "    pisms -Mx 3 -My 3 -Mbz 11 -Lbz 1000 -zb_spacing equal -o foo.nc\n"
  "    btutest -i foo.nc -o bar.nc -ys 0.0 -ye 1.0 -dt 1.0 -Mz 31 -Mbz 11 -Lbz 1000\n"
  "And continuing on refinement path:\n"
  "    pisms -Mx 3 -My 3 -Mbz 21 -Lbz 1000 -zb_spacing equal -o foo.nc\n"
  "    btutest -i foo.nc -o bar.nc -ys 0.0 -ye 1.0 -dt 1.0 -Mz 31 -Mbz 21 -Lbz 1000\n"
  "    pisms -Mx 3 -My 3 -Mbz 41 -Lbz 1000 -zb_spacing equal -o foo.nc\n"
  "    btutest -i foo.nc -o bar.nc -ys 0.0 -ye 1.0 -dt 1.0 -Mz 31 -Mbz 41 -Lbz 1000\n\n";

#include "pism_const.hh"
#include "grid.hh"
#include "PISMIO.hh"
#include "NCVariable.hh"
#include "bedrockThermalUnit.hh"

#include "../../verif/tests/exactTestK.h"


static PetscErrorCode setupIceGridFromFile(string filename, IceGrid &grid) {
  PetscErrorCode ierr;

  PISMIO nc(&grid);
  ierr = nc.get_grid(filename.c_str()); CHKERRQ(ierr);
  grid.compute_nprocs();
  grid.compute_ownership_ranges();
  ierr = grid.createDA(); CHKERRQ(ierr);
  return 0;
}


static PetscErrorCode createVecs(IceGrid &grid, PISMVars &variables) {

  PetscErrorCode ierr;
  IceModelVec2Mask *mask = new IceModelVec2Mask;
  IceModelVec2S *thk = new IceModelVec2S,
                *ghf = new IceModelVec2S;
  IceModelVec3  *enthalpy = new IceModelVec3;

  ierr = mask->create(grid, "mask", true); CHKERRQ(ierr);
  ierr = mask->set_attrs("", "grounded_dragging_floating integer mask",
			      "", ""); CHKERRQ(ierr);
  ierr = variables.add(*mask); CHKERRQ(ierr);

  ierr = thk->create(grid, "thk", true); CHKERRQ(ierr);
  ierr = thk->set_attrs("", "ice thickness",
			      "m", "land_ice_thickness"); CHKERRQ(ierr);
  ierr = variables.add(*thk); CHKERRQ(ierr);

  ierr = enthalpy->create(grid, "enthalpy", false); CHKERRQ(ierr);
  ierr = enthalpy->set_attrs(
     "",
     "ice enthalpy (includes sensible heat, latent heat, pressure)",
     "J kg-1", ""); CHKERRQ(ierr);
  ierr = variables.add(*enthalpy); CHKERRQ(ierr);

  ierr = ghf->create(grid, "bheatflx", false); CHKERRQ(ierr);
  ierr = ghf->set_attrs("",
                       "upward geothermal flux at bedrock thermal layer base",
		       "W m-2", ""); CHKERRQ(ierr);
  ierr = ghf->set_glaciological_units("mW m-2");
  ierr = variables.add(*ghf); CHKERRQ(ierr);

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
    required.push_back("-i");
    required.push_back("-o");
    required.push_back("-ys");
    required.push_back("-ye");
    required.push_back("-dt");
    required.push_back("-Mz");
    required.push_back("-Mbz");
    required.push_back("-Lbz");
    ierr = show_usage_check_req_opts(com, "btutest", required,
      "  btutest -i IN.nc -o OUT.nc -ys A -ye B -dt C -Mz MM -Mbz NN -Lbz 1000.0\n"
      "where:\n"
      "  -i             input file in NetCDF format\n"
      "  -o             output file in NetCDF format\n"
      "  -ys            start year in using Test K\n"
      "  -ye            end year in using Test K\n"
      "  -dt            time step B (= positive float) in years\n"
      "  -Mz            number of ice levels to use; note -Lz option not needed here\n"
      "  -Mbz           number of bedrock thermal layer levels to use\n"
      "  -Lbz 1000.0    depth of bedrock thermal layer (required; Lbz=1000.0 m in Test K)\n"
      ); CHKERRQ(ierr);

    // read the config option database:
    ierr = init_config(com, rank, config, overrides); CHKERRQ(ierr);

    IceGrid grid(com, rank, size, config);

    bool flag;
    PetscReal dt_years = 0.0, ys, ye;
    PetscInt  Mz;
    string inname, outname;
    ierr = PetscOptionsBegin(grid.com, "", "BTU_TEST options", ""); CHKERRQ(ierr);
    {
      ierr = PISMOptionsString("-i", "Input file name", inname, flag); CHKERRQ(ierr);
      ierr = PISMOptionsString("-o", "Output file name", outname, flag); CHKERRQ(ierr);
      ierr = PISMOptionsReal("-ys", "starting time in years", ys, flag); CHKERRQ(ierr);
      ierr = PISMOptionsReal("-ye", "starting time in years", ye, flag); CHKERRQ(ierr);
      ierr = PISMOptionsReal("-dt", "Time-step, in years", dt_years, flag); CHKERRQ(ierr);
      ierr = PISMOptionsInt("-Mz", "number of ice layers", Mz, flag); CHKERRQ(ierr);
      // -Mbz is checked by IceModelVec3BTU (FIXME)
      // -Lbz is checked by IceModelVec3BTU (FIXME)
    }
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    // initialize the computational grid:
    // FIXME:  in fact the axis dimension and variable zb is kept from the input file,
    //         so in order for the output file to make sense the input file has to share
    //         the vertical bedrock grid
    ierr = verbPrintf(2,com,
	"  initializing grid (2D) from NetCDF file %s...\n", inname.c_str()); CHKERRQ(ierr);
    ierr = setupIceGridFromFile(inname,grid); CHKERRQ(ierr);
    grid.start_year = ys;
    grid.end_year   = ye;

    // allocate tools and IceModelVecs
    EnthalpyConverter EC(config);
    PISMVars variables;
    ierr = createVecs(grid, variables); CHKERRQ(ierr);

    // thickness has constant value
    IceModelVec2S *thk;
    thk = dynamic_cast<IceModelVec2S*>(variables.get("thk"));
    if (thk == NULL) SETERRQ(3, "thk is not available");
    ierr = thk->set(3000.0); CHKERRQ(ierr);  // see Test K

    // lithosphere (bottom of bedrock layer) heat flux has constant value
    IceModelVec2S *ghf;
    ghf = dynamic_cast<IceModelVec2S*>(variables.get("bheatflx"));
    if (ghf == NULL) SETERRQ(4, "bheatflx is not available");
    ierr = ghf->set(0.042); CHKERRQ(ierr);  // see Test K

    // we will be filling enthalpy from temperature used in Test K
    IceModelVec3 *enthalpy;
    enthalpy = dynamic_cast<IceModelVec3*>(variables.get("enthalpy"));
    if (enthalpy == NULL) SETERRQ(5, "enthalpy is not available");

    // Initialize BTU object:
    PISMBedThermalUnit btu(grid, EC, config);

    ierr = btu.init(variables); CHKERRQ(ierr);  // FIXME: this is bootstrapping, really

    // worry about time step
    int  N = (int)ceil((ye - ys) / dt_years);
    PetscReal old_dt_years = dt_years;
    dt_years = (ye - ys) / (double)N;
    ierr = verbPrintf(2,com,
        "  user set timestep of %.4f years ...\n"
        "  reset to %.4f years to get integer number of steps ... \n",
	old_dt_years,dt_years); CHKERRQ(ierr);
    PetscReal max_dt_years;
    ierr = btu.max_timestep(0.0, max_dt_years); CHKERRQ(ierr);
    ierr = verbPrintf(2,com,
        "  PISMBedThermalUnit reports max timestep of %.4f years ...\n",
	max_dt_years); CHKERRQ(ierr);

    // fill exact bedrock temperature from Test K at time ys
    PetscReal dzb, Lbz;
    PetscInt  Mbz;
    btu.temp.get_levels(Mbz);
    btu.temp.get_spacing(dzb);
    btu.temp.get_layer_depth(Lbz);
    ierr = btu.temp.begin_access(); CHKERRQ(ierr);
    PetscScalar *Tb; // columns of these values
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        ierr = btu.temp.getInternalColumn(i,j,&Tb); CHKERRQ(ierr);
        for (PetscInt k=0; k < Mbz; k++) {
          const PetscReal z = - Lbz + (double)k * dzb;
          PetscReal FF; // Test K:  use Tb[k], ignor FF
          ierr = exactK(ys * secpera, z, &Tb[k], &FF, 0); CHKERRQ(ierr);
        }
      }
    }
    ierr = btu.temp.end_access(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com,"  running: "); CHKERRQ(ierr);
    // actually ask btu to do its time-step
    for (PetscInt n = 0; n < N; n++) {
      const PetscReal y = ys + dt_years * (double)n;  // time at start of time-step

      // compute exact ice temperature, thus enthalpy, at time y
      ierr = enthalpy->begin_access(); CHKERRQ(ierr);
      PetscScalar *Enthij; // columns of these values
      for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
        for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
          ierr = enthalpy->getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
          for (PetscInt k=0; k<grid.Mz; ++k) {
            const PetscReal z     = grid.zlevels[k],
                            depth = 3000.0 - z; // FIXME task #7297
            PetscReal TT, FF; // Test K:  use TT, ignor FF
            ierr = exactK(y * secpera, z, &TT, &FF, 0); CHKERRQ(ierr);
            ierr = EC.getEnth(TT,0.0,EC.getPressureFromDepth(depth),Enthij[k]); CHKERRQ(ierr);
          }
        }
      }
      ierr = enthalpy->end_access(); CHKERRQ(ierr);
      // we are not communicating anything, which is fine

      // update the temperature inside the thermal layer
      ierr = btu.update(y, dt_years); CHKERRQ(ierr);
      ierr = verbPrintf(2,com,"."); CHKERRQ(ierr);
    }
    ierr = verbPrintf(2,com,"\n  done ...\n"); CHKERRQ(ierr);

    // reuse ghf for final output heat flux G_0 at z=0
    ierr = ghf->set_name("bheatflx0"); CHKERRQ(ierr);
    ierr = ghf->set_attrs("",
                       "upward geothermal flux at ice/bedrock interface",
		       "W m-2", ""); CHKERRQ(ierr);
    ierr = btu.get_upward_geothermal_flux(*ghf); CHKERRQ(ierr);

    // report result to std out, at center of grid
    ierr = ghf->begin_access(); CHKERRQ(ierr);
    const PetscInt ci = grid.Mx/2, cj = grid.My/2;
    if ((ci >= grid.xs) && (ci < grid.xs+grid.xm) && (cj >= grid.ys) && (cj < grid.ys+grid.ym)) {
      ierr = verbPrintf(2,com,
        "  numerical upward heat flux at z=0, at end year %.2f, at center of grid is\n"
        "    G_0 = %.7f W m-2;\n",
	ye,(*ghf)(ci,cj)); CHKERRQ(ierr);
    }
    ierr = ghf->end_access(); CHKERRQ(ierr);

    // compare to the right answer from Test K
    PetscReal TT, FF; // Test K:  use FF, ignor TT
    ierr = exactK(ye * secpera, 0.0, &TT, &FF, 0); CHKERRQ(ierr);
    ierr = verbPrintf(2,com,
        "  exact Test K reports upward heat flux at z=0, at end year %.2f, as G_0 = %.7f W m-2;\n",
	ye,FF); CHKERRQ(ierr);

    set<string> vars;
    btu.add_vars_to_output("big", vars); // "write everything you can"

    PISMIO pio(&grid);

    ierr = pio.open_for_writing(outname, false, true); CHKERRQ(ierr);
    // append == true and check_dims == true
    ierr = pio.append_time(grid.end_year); CHKERRQ(ierr);
    ierr = btu.define_variables(vars, pio, NC_DOUBLE); CHKERRQ(ierr);
    ierr = pio.close(); CHKERRQ(ierr);

    ierr = btu.write_variables(vars, outname); CHKERRQ(ierr);
    ierr = ghf->write(outname.c_str()); CHKERRQ(ierr);

    ierr = doneWithIceInfo(variables); CHKERRQ(ierr);
    ierr = verbPrintf(2,com, "done.\n"); CHKERRQ(ierr);

  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}

