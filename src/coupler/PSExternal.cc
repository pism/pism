// Copyright (C) 2010, 2011 Constantine Khroulev
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

#include "PSExternal.hh"
#include "PISMIO.hh"
#include <sys/file.h>

PSExternal::~PSExternal() {
  int done = 1;

  // tell the EBM driver to stop:
  if (grid.rank == 0) {
    MPI_Send(&done, 1, MPI_INT, 0, TAG_EBM_STOP, inter_comm);
  }
}

//! Initialize the PSExternal model.
PetscErrorCode PSExternal::init(PISMVars &vars) {
  PetscErrorCode ierr;
  string pism_input;
  LocalInterpCtx *lic;
  bool regrid;
  int start;

  usurf = dynamic_cast<IceModelVec2S*>(vars.get("surface_altitude"));
  if (!usurf) { SETERRQ(1, "ERROR: Surface elevation is not available"); }

  topg = dynamic_cast<IceModelVec2S*>(vars.get("bedrock_altitude"));
  if (!topg) { SETERRQ(1, "ERROR: Bed elevation is not available"); }

  ierr = acab.create(grid, "acab", false); CHKERRQ(ierr);
  ierr = acab.set_attrs(
            "climate_from_PISMSurfaceModel",  // FIXME: can we do better?
            "ice-equivalent surface mass balance (accumulation/ablation) rate",
	    "m s-1",  // m *ice-equivalent* per second
	    "land_ice_surface_specific_mass_balance");  // CF standard_name
	    CHKERRQ(ierr);
  ierr = acab.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  acab.write_in_glaciological_units = true;
  acab.set_attr("comment", "positive values correspond to ice gain");

  // annual mean air temperature at "ice surface", at level below all firn
  //   processes (e.g. "10 m" or ice temperatures)
  ierr = artm.create(grid, "artm", false); CHKERRQ(ierr);
  ierr = artm.set_attrs(
            "climate_from_PISMSurfaceModel",  // FIXME: can we do better?
            "annual average ice surface temperature, below firn processes",
            "K", 
            "");  // PROPOSED CF standard_name = land_ice_surface_temperature_below_firn
  CHKERRQ(ierr);

  // artm_0 is the initial condition; artm_0 = artm(t_0) + gamma*usurf(t_0)
  ierr = artm_0.create(grid, "usurf", false); CHKERRQ(ierr);
  ierr = artm_0.set_attrs("internal", "ice upper surface elevation",
                           "m", "surface_altitude"); CHKERRQ(ierr);

  ierr = find_pism_input(pism_input, lic, regrid, start); CHKERRQ(ierr); 

  if (regrid) {
    ierr = artm_0.regrid(pism_input.c_str(), *lic, true); CHKERRQ(ierr);
    ierr =   artm.regrid(pism_input.c_str(), *lic, true); CHKERRQ(ierr);
  } else {
    ierr = artm_0.read(pism_input.c_str(), start); CHKERRQ(ierr);
    ierr =   artm.read(pism_input.c_str(), start); CHKERRQ(ierr);
  }

  delete lic;

  ierr = PetscOptionsBegin(grid.com, "", "PSExternal model options", ""); CHKERRQ(ierr);
  {
    bool flag;
    ierr = PISMOptionsReal("-lapse_rate", "Air temperature lapse rate, degrees K per kilometer",
			   gamma, flag); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-update_interval", "Energy balance model update interval, years",
			   update_interval, flag); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-min_lag", "Minimal lag (should be less than -update_interval but greater than 0)",
                           min_lag, flag); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  gamma = gamma / 1000;         // convert to K/meter

  // Use gamma to compute the initial condition:
  ierr = artm_0.scale(gamma); CHKERRQ(ierr);
  ierr = artm_0.add(1.0, artm); CHKERRQ(ierr);

  // Initialize the EBM driver:

  // Send parameters read from command-line options to the EBM driver:
  ebm_input  = "ebm_input.nc";
  ebm_output = "ebm_output.nc";

  char command[PETSC_MAX_PATH_LEN];
  snprintf(command, PETSC_MAX_PATH_LEN,
           "ebm -i %s -o %s", ebm_input.c_str(), ebm_output.c_str());

  if (grid.rank == 0) {
    MPI_Send(command, PETSC_MAX_PATH_LEN, MPI_CHAR, 0, TAG_EBM_COMMAND, inter_comm);
  }

  return 0;
}

PetscErrorCode PSExternal::ice_surface_mass_flux(PetscReal t_years, PetscReal dt_years,
                                            IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = update(t_years, dt_years); CHKERRQ(ierr);

  ierr = acab.copy_to(result); CHKERRQ(ierr); 

  return 0;
}

PetscErrorCode PSExternal::ice_surface_temperature(PetscReal t_years, PetscReal dt_years,
                                              IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = update(t_years, dt_years); CHKERRQ(ierr);

  ierr = artm.copy_to(result); CHKERRQ(ierr); 

  return 0;
}

//! \brief Update the surface mass balance field by reading from a file created
//! by an EBM. Also, write ice surface elevation and bed topography for an EBM to read.
PetscErrorCode PSExternal::update(PetscReal t_years, PetscReal dt_years) {
  PetscErrorCode ierr;

  if ((fabs(t_years - t) < 1e-12) &&
      (fabs(dt_years - dt) < 1e-12))
    return 0;

  if (t_years + dt_years < t + update_interval)
    return 0;

  // FIXME: insert a check for whether we need to run EBM to pre-compute acab

  t  = t_years;
  dt = dt_years;

  // The actual update:

  // update PISM's b.c.:
  ierr = update_artm(); CHKERRQ(ierr);
  ierr = update_acab(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSExternal::update_acab() {
  PetscErrorCode ierr;
  PISMIO nc(&grid);
  int ebm_status;

  ierr = write_coupling_fields(); CHKERRQ(ierr);

  if (grid.rank == 0) {
    MPI_Send(&ebm_status, 1, MPI_INT, 0, TAG_EBM_RUN, inter_comm);
  }

  if (grid.rank == 0) {
    int sleep_interval = 1,       // seconds
      threshold = 10;

    int wait_counter = 0;

    MPI_Status status;
    while (wait_counter * sleep_interval < threshold) {
      int flag;
      MPI_Iprobe(0, TAG_EBM_STATUS, inter_comm, &flag, &status);

      if (flag)
        break;

      fprintf(stderr, "PISM: Waiting for a message from the EBM driver...\n");
      sleep(sleep_interval);
      wait_counter++;
    }

    if (wait_counter >= threshold) {
      // exited the loop above because of a timeout
      fprintf(stderr, "ERROR: spent %1.1f minutes waiting for the EBM driver... Giving up...\n",
              threshold / 60.0);
      PISMEnd();
    }

    fprintf(stderr, "PISM: Got a status message from EBM\n");

    // receive the EBM status
    MPI_Recv(&ebm_status, 1, MPI_INT,
             0, TAG_EBM_STATUS, inter_comm, NULL);
  }

  // Broadcast status:
  MPI_Bcast(&ebm_status, 1, MPI_INT, 0, grid.com);

  if (ebm_status == EBM_STATUS_FAILED) {
    PetscPrintf(grid.com, "ERROR: EBM failure. Exiting...\n");
    PISMEnd();
  }

  grid_info gi;
  ierr = nc.open_for_reading(ebm_output.c_str()); CHKERRQ(ierr);
  ierr = nc.get_grid_info_2d(gi); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);
  LocalInterpCtx lic(gi, NULL, NULL, grid); // 2D only

  ierr = acab.regrid(ebm_output.c_str(), lic, true); CHKERRQ(ierr);

  return 0;
}

//! Update artm using an atmospheric lapse rate.
PetscErrorCode PSExternal::update_artm() {
  PetscErrorCode ierr;

  ierr = usurf->begin_access(); CHKERRQ(ierr);
  ierr = artm.begin_access(); CHKERRQ(ierr);
  ierr = artm_0.begin_access(); CHKERRQ(ierr); 
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      artm(i,j) = artm_0(i,j) - gamma * (*usurf)(i,j);
    }
  }
  ierr = usurf->end_access(); CHKERRQ(ierr);
  ierr = artm.end_access(); CHKERRQ(ierr);
  ierr = artm_0.end_access(); CHKERRQ(ierr); 

  return 0;
}

//! Write fields that a model PISM is coupled to needs. Currently: usurf and topg.
PetscErrorCode PSExternal::write_coupling_fields() {
  PetscErrorCode ierr;
  PISMIO nc(&grid);

  ierr = nc.open_for_writing(ebm_input.c_str(),
                             false, true); CHKERRQ(ierr);
  // "append" (i.e. do not move the file aside) and check dimensions. Note that
  // we only append a record (by calling append_time()) if the file is empty.

  ierr = nc.append_time(grid.year); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  // write the fields an EBM needs:
  ierr = usurf->write(ebm_input.c_str()); CHKERRQ(ierr);
  ierr = topg->write(ebm_input.c_str()); CHKERRQ(ierr);

  return 0;
}


