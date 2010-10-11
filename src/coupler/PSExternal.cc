#include "PSExternal.hh"
#include <sys/file.h>

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

  // do ghosts and use the wide stencil so that we can do arithmetic (with
  // usurf) easily:
  ierr = usurf_0.create(grid, "usurf", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = usurf_0.set_attrs("", "ice upper surface elevation",
                           "m", "surface_altitude"); CHKERRQ(ierr);

  ierr = find_pism_input(pism_input, lic, regrid, start); CHKERRQ(ierr); 

  if (regrid) {
    ierr = usurf_0.regrid(pism_input.c_str(), *lic, true); CHKERRQ(ierr);
    ierr =    artm.regrid(pism_input.c_str(), *lic, true); CHKERRQ(ierr);
  } else {
    ierr = usurf_0.read(pism_input.c_str(), start); CHKERRQ(ierr);
    ierr =    artm.read(pism_input.c_str(), start); CHKERRQ(ierr);
  }

  delete lic;

  ierr = PetscOptionsBegin(grid.com, "", "Air temp. lapse rate model options", ""); CHKERRQ(ierr);
  {
    bool flag;
    ierr = PISMOptionsReal("-lapse_rate", "Air temperature lapse rate, degrees K per kilometer",
			   gamma, flag); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-ebm_update_interval", "Energy balance model update interval, years",
			   update_interval, flag); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

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

PetscErrorCode PSExternal::lock(string name, int op, int &fd) {
  PetscErrorCode ierr;
  ierr = verbPrintf(3, grid.com, "  Trying to access %s...",
                    name.c_str()); CHKERRQ(ierr);

  if (grid.rank == 0) {
    fd = open(name.c_str(), O_RDONLY);
    if (ierr == -1) {
      printf("ERROR: can't open '%s'.\n", name.c_str());
      PetscEnd();
    }

    ierr = flock(fd, op);    // a "shared" lock
    if (ierr == -1) {
      printf("ERROR: obtaining a lock on '%s' failed. Error code: %d\n",
             name.c_str(), errno);
      PetscEnd();
    }
  } else {
    fd = 0;
  }

  // make other (non-zero) rank processors wait
  ierr = MPI_Barrier(grid.com); CHKERRQ(ierr); 
  ierr = verbPrintf(3, grid.com, " done.\n");
  return 0;
}

PetscErrorCode PSExternal::unlock(int fd) {
  ierr = verbPrintf(3, grid.com, "  Removing a lock from %s...",
                    in_filename.c_str()); CHKERRQ(ierr);

  if (grid.rank == 0) {
    ierr = flock(fd, LOCK_UN);
    if (ierr == -1) {
      printf("ERROR: removing a lock from '%s' failed. Error code: %d\n",
             in_filename.c_str(), errno);
      PetscEnd();
    }

    ierr = close(fd);
    if (ierr == -1) {
      printf("ERROR: closing '%s' failed. Error code: %d\n",
             in_filename.c_str(), errno);
      PetscEnd();
    }
  }

  // make other (non-zero) rank processors wait
  ierr = MPI_Barrier(grid.com); CHKERRQ(ierr); 
  ierr = verbPrintf(3, grid.com, " done.\n");
  return 0;
}

//! \brief Update the surface mass balance field by reading from a file created
//! by an EBM. Also, write ice surface elevation and bed topography for an EBM to read.
/*!
  \note The coupling NetCDF file in_filename has to exist before calling this method.
  
  \note This method will try to open the file and read its grid information; if
  the last time record does not correspond to a time within the current update
  interval, it will wait for an EBM to update the file if possible. If an EBM
  is ahead of PISM, PISM cannot continue and stops.
 
  \note Uses the UNIX file locking mechanism (flock()) to ensure that coupling
  files are completely written before being read. For this to work, an EBM
  needs to apply a shared lock (LOCK_SH) to a file it is currently reading and an
  exclusive lock (LOCK_EX) to a file that it writes the SMB to.

  \note It should be possible to switch to a different locking mechanism by
  re-implementing private methods lock() and unlock().

  \note It might be easiest to implement file locking in an EBM written in
  Fortran by writing wrappers calling C functions open() and flock().
 */
petscerrorcode PSExternal::update(PetscReal t_years, PetscReal dt_years) {
  PetscErrorCode ierr;
  PISMIO nc(grid);
  LocalInterpCtx *lic;
  grid_info g;
  bool file_is_ready = false;
  int fd;

  if ((fabs(t_years - t) < 1e-12) &&
      (fabs(dt_years - dt) < 1e-12))
    return 0;

  if (t_years + dt_years < t + update_interval)
    return 0;

  t  = t_years;
  dt = dt_years;

  // The actual update:

  while ( ! file_is_ready ) {
    double time_found;

    ierr = lock(in_filename.c_str(), LOCK_SH, fd); CHKERRQ(ierr); 
    {
      ierr = nc.open_for_reading(in_filename.c_str()); CHKERRQ(ierr); 
      ierr = nc.get_grid_info_2d(g); CHKERRQ(ierr);
      ierr = nc.close(); CHKERRQ(ierr);
    }
    ierr = unlock(fd); CHKERRQ(ierr); 

    time_found = g.time / secpera; // convert to years

    if (time_found < t_years) {
      ierr = verbPrintf(2,grid.com,
                        "  NOTE: file %s is not ready (time expected: %3.3f years, time found: %3.3f years)\n"
                        "  Will wait 5 seconds and try again...\n",
                        in_filename.c_str(), t_years, time_found); CHKERRQ(ierr);
    } else if (time_found <= t_years + dt_years) {
      file_is_ready = true;
    } else {
      ierr = PetscPrintf(grid.com,
                         "  ERROR: according to contents of %s (time expected: %3.3f years, time found: %3.3f years)\n"
                         "  EBM somehow got ahead of PISM. PISM is confused and cannot continue.\n",
                         in_filename.c_str(), t_years, time_found); CHKERRQ(ierr);
      PetscEnd();
    }

    if ( ! file_is_ready ) {
      ierr = PetscSleep(5); CHKERRQ(ierr);
    }
  } // end of while (!file_is_ready)

  // read an accumulation/ablation rate field provided by an EBM:
  ierr = lock(in_filename.c_str(), LOCK_SH, fd); CHKERRQ(ierr); 
  {
    lic = new LocalInterpCtx(g, NULL, NULL, grid);
    ierr = acab.regrid(in_filename.c_str(), *lic, true); CHKERRQ(ierr);
    delete lic;
  }
  ierr = unlock(fd); CHKERRQ(ierr); 

  // write data for an EBM to read

  // apply an "exclusive" lock
  ierr = lock(out_filename.c_str(), LOCK_EX, fd); CHKERRQ(ierr);
  {
    ierr = nc.open_for_writing(out_filename.c_str(),
                               true, true); CHKERRQ(ierr);
    // "append" (i.e. do not move the file aside) and check dimensions. Note that
    // we only append a record (by calling append_time()) if the file is empty.

    int t_len;
    ierr = nc.get_dim_length("t", &t_len); CHKERRQ(ierr);
    if (t_len == 0) {
      ierr = nc.append_time(grid.year); CHKERRQ(ierr); 
    } else {
      int varid;
      bool exists;
      ierr = nc.find_variable("t", &varid, exists); CHKERRQ(ierr); 
      ierr = nc.put_dimension(varid, 1, &grid.year); CHKERRQ(ierr);
    }
    ierr = nc.close(); CHKERRQ(ierr);

    // write the fields an EBM needs:
    ierr = usurf->write(out_filename.c_str()); CHKERRQ(ierr);
    ierr = topg->write(out_filename.c_str()); CHKERRQ(ierr);
  }
  ierr = unlock(fd); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSExternal::write_model_state(PetscReal t_years, PetscReal dt_years,
                                        string filename) {
  PetscErrorCode ierr;

  // write current usurf and current artm (corrected for this usurf)
  return 0;
}
