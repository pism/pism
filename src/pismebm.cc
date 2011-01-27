// Copyright (C) 2004-2011 Ed Bueler and Constantine Khroulev
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
  "Ice sheet driver for PISM ice sheet simulations, initialized from data.\n"
  "Uses external serial code providing top-surface boundary conditions..\n";

#include <petsc.h>
#include "grid.hh"
#include "iceModel.hh"

#include "coupler/PCFactory.hh"
#include "coupler/PISMOcean.hh"
#include "PSExternal.hh"

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm    inter_comm, intra_comm;
  PetscMPIInt rank, size;
  bool ebm_side, pism_side;

  // Initialize MPI
  ierr = MPI_Init(&argc, &argv); CHKERRQ(ierr);

  // Get the processor rank within MPI_COMM_WORLD (used to split MPI_COMM_WORLD
  // into "rank == 0" and "everyone else").
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &size); CHKERRQ(ierr);

  if (size == 1) {
    fprintf(stderr, "ERROR: pismebm cannot run on one processor!\n");
    ierr = MPI_Finalize(); CHKERRQ(ierr); 
    exit(1);
  }

  // Split.
  ebm_side  = (rank == 0);
  pism_side = (rank != 0);
  ierr = MPI_Comm_split(MPI_COMM_WORLD,
                        rank != 0, // color
                        rank,      // rank
                        &intra_comm); CHKERRQ(ierr);

  // Create the inter-communicator
  ierr = MPI_Intercomm_create(intra_comm, // local communicator
                              0,          // local leader (rank)
                              MPI_COMM_WORLD, // peer communicator
                              rank == 0 ? 1 : 0, // remote leader (rank in MPI_COMM_WORLD)
                              123, // tag, has to be a "unique" integer
                              &inter_comm); CHKERRQ(ierr); 

  // Now this inter_comm can be used to send data from the processor running
  // EBM to other processors (running PISM).

  // We could also use MPI_COMM_WORLD, but it would be more confusing (it would
  // be harder to keep track of processor ranks).

  if (ebm_side) {
    EBM_driver ebm(inter_comm);
    ierr = ebm.run(); CHKERRQ(ierr);
  }

  if (pism_side) {
    // Initialize PETSC
    PETSC_COMM_WORLD = intra_comm;
    ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);

    ierr = MPI_Comm_rank(intra_comm, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(intra_comm, &size); CHKERRQ(ierr);

    // explicit scoping
    {
      ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);

      ierr = verbPrintf(2,intra_comm, "PISMEBM %s (basic evolution run mode)\n",
                        PISM_Revision); CHKERRQ(ierr);
      ierr = stop_on_version_option(); CHKERRQ(ierr);

      ierr = check_old_option_and_stop(intra_comm, "-boot_from", "-boot_file"); CHKERRQ(ierr); 

      bool iset, bfset;
      ierr = PISMOptionsIsSet("-i", iset); CHKERRQ(ierr);
      ierr = PISMOptionsIsSet("-boot_file", bfset); CHKERRQ(ierr);
      string usage =
        "  pismebm {-i IN.nc|-boot_file IN.nc} [OTHER PISM & PETSc OPTIONS]\n"
        "where:\n"
        "  -i          IN.nc is input file in NetCDF format: contains PISM-written model state\n"
        "  -boot_file  IN.nc is input file in NetCDF format: contains a few fields, from which\n"
        "              heuristics will build initial model state\n"
        "notes:\n"
        "  * one of -i or -boot_file is required\n"
        "  * if -boot_file is used then also '-Mx A -My B -Mz C -Lz D' are required\n";
      if ((iset == PETSC_FALSE) && (bfset == PETSC_FALSE)) {
        ierr = PetscPrintf(intra_comm,
                           "\nPISM ERROR: one of options -i,-boot_file is required\n\n"); CHKERRQ(ierr);
        ierr = show_usage_and_quit(intra_comm, "pismebm", usage.c_str()); CHKERRQ(ierr);
      } else {
        vector<string> required;  required.clear();
        ierr = show_usage_check_req_opts(intra_comm, "pismebm", required, usage.c_str()); CHKERRQ(ierr);
      }

      NCConfigVariable config, overrides;
      ierr = init_config(intra_comm, rank, config, overrides); CHKERRQ(ierr);

      IceGrid g(intra_comm, rank, size, config);
      IceModel m(g, config, overrides);

      // Initialize boundary models:

      PISMSurfaceModel *surface;
      bool flag;
      ierr = PISMOptionsIsSet("-artm_lapse_rate",
                              "Atmospheric lapse rate for the temperature at the top of the ice",
                              flag); CHKERRQ(ierr);
      if (flag) {
        surface = new PSExternal_ALR(g, config, inter_comm);
      } else {
        // The special surface model (does *not* use an atmosphere model)
        surface = new PSExternal(g, config, inter_comm);
      }
      

      // An ocean model
      POFactory po(g, config);
      PISMOceanModel *ocean;

      ierr = PetscOptionsBegin(intra_comm, "", "Options choosing PISM boundary models", ""); CHKERRQ(ierr);
      po.create(ocean);
      ierr = PetscOptionsEnd(); CHKERRQ(ierr);

      m.attach_ocean_model(ocean);
      m.attach_surface_model(surface);
      ierr = m.setExecName("pismebm"); CHKERRQ(ierr);

      ierr = m.init(); CHKERRQ(ierr);

      ierr = m.run(); CHKERRQ(ierr);

      ierr = verbPrintf(2,intra_comm, "... done with run\n"); CHKERRQ(ierr);
      // provide a default output file name if no -o option is given.
      ierr = m.writeFiles("unnamed.nc"); CHKERRQ(ierr);
    }

    // PETSc clean-up
    ierr = PetscFinalize(); CHKERRQ(ierr);
  } // end of "if (pism_side)"

  // De-allocate communicators
  ierr = MPI_Comm_free(&inter_comm); CHKERRQ(ierr);
  ierr = MPI_Comm_free(&intra_comm); CHKERRQ(ierr);

  // MPI clean-up
  ierr = MPI_Finalize(); CHKERRQ(ierr); 
                        
  return 0;
}
