// Copyright (C) 2012-2013 PISM Authors
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

#include "PISMVars.hh"
#include "pism_options.hh"
#include "Mask.hh"
#include "PISMHydrology.hh"
#include "hydrology_diagnostics.hh"


PISMHydrology::PISMHydrology(IceGrid &g, const NCConfigVariable &conf)
  : PISMComponent_TS(g, conf)
{
  thk   = NULL;
  bed   = NULL;
  cellarea = NULL;
  bmelt = NULL;
  mask  = NULL;
  inputtobed = NULL;
  variables = NULL;

  PetscErrorCode ierr1, ierr2, ierr3;
  ierr1 = total_input.create(grid, "total_input_hydro", false);
  ierr2 = total_input.set_attrs("internal",
                         "workspace for total input rate into subglacial water layer",
                         "m s-1", "");
  if ((ierr1 != 0) || (ierr2 != 0)) {
      PetscPrintf(grid.com,
        "PISM ERROR: memory allocation failed in PISMHydrology constructor (total_input).\n");
      PISMEnd();
  }

  // *all* PISMHydrology classes have layer of water stored in till
  ierr1 = Wtil.create(grid, "tillwat", false);
  ierr2 = Wtil.set_attrs("model_state",
                     "effective thickness of subglacial water stored in till",
                     "m", "");
  ierr3 = Wtil.set_attr("valid_min", 0.0);
  if ((ierr1 != 0) || (ierr2 != 0) || (ierr3 != 0)) {
      PetscPrintf(grid.com,
        "PISM ERROR: memory allocation failed in PISMHydrology constructor (Wtil).\n");
      PISMEnd();
  }
}


PetscErrorCode PISMHydrology::init(PISMVars &vars) {
  PetscErrorCode ierr;
  string itbfilename;  // itb = input_to_bed
  bool itbfile_set, itbperiod_set, itbreference_set;
  bool i_set, bootstrap;
  PetscReal itbperiod_years = 0.0, itbreference_year = 0.0;

  ierr = verbPrintf(4, grid.com,
    "entering PISMHydrology::init() ...\n"); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, "",
            "Options controlling the base class PISMHydrology", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsString("-input_to_bed_file",
      "A time- and space-dependent file with amount of water (depth per time) which should be put at the ice sheet bed at the given location at the given time",
      itbfilename, itbfile_set); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-input_to_bed_period",
      "The period (i.e. duration before repeat), in years, of -input_to_bed_file data",
      itbperiod_years, itbperiod_set); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-input_to_bed_reference_year",
      "The reference year for periodizing the -input_to_bed_file data",
      itbreference_year, itbreference_set); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-i", "PISM input file", i_set); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-boot_file", "PISM bootstrapping file",
                            bootstrap); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  variables = &vars;

  // the following are IceModelVec pointers into IceModel generally and are read by code in the
  // update() method at the current PISMHydrology time

  thk = dynamic_cast<IceModelVec2S*>(vars.get("thk"));
  if (thk == NULL) SETERRQ(grid.com, 1, "thk is not available to PISMHydrology");

  bed = dynamic_cast<IceModelVec2S*>(vars.get("topg"));
  if (bed == NULL) SETERRQ(grid.com, 1, "topg is not available to PISMHydrology");

  bmelt = dynamic_cast<IceModelVec2S*>(vars.get("bmelt"));
  if (bmelt == NULL) SETERRQ(grid.com, 1, "bmelt is not available to PISMHydrology");

  cellarea = dynamic_cast<IceModelVec2S*>(vars.get("cell_area"));
  if (cellarea == NULL) SETERRQ(grid.com, 1, "cell_area is not available to PISMHydrology");

  mask = dynamic_cast<IceModelVec2Int*>(vars.get("mask"));
  if (mask == NULL) SETERRQ(grid.com, 1, "mask is not available to PISMHydrology");

  // the following inputtobed is not related to IceModel; we must read it ourselves

  if (itbfile_set) {
    inputtobed_period = (itbperiod_set) ? itbperiod_years : 0.0;
    inputtobed_reference_time = (itbreference_set) ? grid.convert(itbreference_year, "years", "seconds") : 0.0;

    unsigned int buffer_size = (unsigned int) config.get("climate_forcing_buffer_size"),
                 n_records = 1;

    PIO nc(grid.com, grid.rank, "netcdf3", grid.get_unit_system());
    ierr = nc.open(itbfilename, PISM_NOWRITE); CHKERRQ(ierr);
    ierr = nc.inq_nrecords("inputtobed", "", n_records); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);

    // if -..._period is not set, make n_records the minimum of the
    // buffer size and the number of available records. Otherwise try
    // to keep all available records in memory.
    if (itbperiod_set == false) {
      n_records = PetscMin(n_records, buffer_size);
    }

    if (n_records == 0) {
      PetscPrintf(grid.com, "PISM ERROR: can't find 'inputtobed' in -input_to_bed file with name '%s'.\n",
                  itbfilename.c_str());
      PISMEnd();
    }

    ierr = verbPrintf(2,grid.com,
      "    option -input_to_bed_file seen ... creating 'inputtobed' variable ...\n"); CHKERRQ(ierr);
    ierr = verbPrintf(2,grid.com,
      "    allocating buffer space for n = %d 'inputtobed' records ...\n", n_records); CHKERRQ(ierr);
    inputtobed = new IceModelVec2T;
    inputtobed->set_n_records(n_records);
    ierr = inputtobed->create(grid, "inputtobed", false); CHKERRQ(ierr);
    ierr = inputtobed->set_attrs("climate_forcing",
                                 "amount of water (depth per time like bmelt) which should be put at the ice sheet bed",
                                 "m s-1", ""); CHKERRQ(ierr);
    ierr = verbPrintf(2,grid.com,
      "    reading 'inputtobed' variable from file '%s' ...\n",itbfilename.c_str()); CHKERRQ(ierr);
    ierr = inputtobed->init(itbfilename, inputtobed_period, inputtobed_reference_time); CHKERRQ(ierr);
  }

  // initialize till water layer thickness from the context if present,
  //   otherwise from -i or -boot_file, otherwise with constant value
  IceModelVec2S *Wtil_input = dynamic_cast<IceModelVec2S*>(vars.get("tillwat"));
  if (Wtil_input != NULL) { // a variable called "tillwat" is already in context
    ierr = Wtil.copy_from(*Wtil_input); CHKERRQ(ierr);
  } else if (i_set || bootstrap) {
    string filename;
    int start;
    ierr = find_pism_input(filename, bootstrap, start); CHKERRQ(ierr);
    if (i_set) {
      ierr = Wtil.read(filename, start); CHKERRQ(ierr);
    } else {
      ierr = Wtil.regrid(filename,
                         config.get("bootstrapping_tillwat_value_no_var")); CHKERRQ(ierr);
    }
  } else {
    ierr = Wtil.set(config.get("bootstrapping_tillwat_value_no_var")); CHKERRQ(ierr);
  }

  // whether or not we could initialize from file, we could be asked to regrid from file
  ierr = regrid(Wtil); CHKERRQ(ierr);
  return 0;
}


void PISMHydrology::get_diagnostics(map<string, PISMDiagnostic*> &dict,
                                    map<string, PISMTSDiagnostic*> &/*ts_dict*/) {
  dict["bwat"] = new PISMHydrology_bwat(this, grid, *variables);
  dict["bwp"] = new PISMHydrology_bwp(this, grid, *variables);
  dict["bwprel"] = new PISMHydrology_bwprel(this, grid, *variables);
  dict["effbwp"] = new PISMHydrology_effbwp(this, grid, *variables);
  dict["hydroinput"] = new PISMHydrology_hydroinput(this, grid, *variables);
  dict["wallmelt"] = new PISMHydrology_wallmelt(this, grid, *variables);
}


void PISMHydrology::add_vars_to_output(string /*keyword*/, set<string> &result) {
  result.insert("tillwat");
}


PetscErrorCode PISMHydrology::define_variables(set<string> vars, const PIO &nc,
                                               PISM_IO_Type nctype) {
  PetscErrorCode ierr;
  if (set_contains(vars, "tillwat")) {
    ierr = Wtil.define(nc, nctype); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode PISMHydrology::write_variables(set<string> vars, const PIO &nc) {
  PetscErrorCode ierr;
  if (set_contains(vars, "tillwat")) {
    ierr = Wtil.write(nc); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode PISMHydrology::regrid(IceModelVec2S &myvar) {
  PetscErrorCode ierr;
  bool file_set, vars_set;
  string file;
  set<string> vars;

  ierr = PetscOptionsBegin(grid.com, "", "PISMHydrology regridding options", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsString("-regrid_file", "regridding file name",file, file_set); CHKERRQ(ierr);
    ierr = PISMOptionsStringSet("-regrid_vars", "comma-separated list of regridding variables",
                                "", vars, vars_set); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (file_set && vars_set && set_contains(vars, myvar.string_attr("short_name"))) {
    ierr = verbPrintf(2, grid.com, "  regridding '%s' from file '%s' ...\n",
                      myvar.string_attr("short_name").c_str(), file.c_str()); CHKERRQ(ierr);
    ierr = myvar.regrid(file, true); CHKERRQ(ierr);
  }
  return 0;
}


//! Update the overburden pressure from ice thickness.
/*!
Uses the standard hydrostatic (shallow) approximation of overburden pressure,
  \f[ P_0 = \rho_i g H \f]
Accesses H=thk from PISMVars, which points into IceModel.
 */
PetscErrorCode PISMHydrology::overburden_pressure(IceModelVec2S &result) {
  PetscErrorCode ierr;
  // FIXME issue #15
  ierr = result.copy_from(*thk); CHKERRQ(ierr);  // copies into ghosts if result has them
  ierr = result.scale(config.get("ice_density") * config.get("standard_gravity")); CHKERRQ(ierr);
  return 0;
}


//! Return the effective thickness of the water stored in till.
PetscErrorCode PISMHydrology::till_water_thickness(IceModelVec2S &result) {
  PetscErrorCode ierr = Wtil.copy_to(result); CHKERRQ(ierr);
  return 0;
}


//! Set the wall melt rate to zero.  (The most basic subglacial hydrologies have no lateral flux or potential gradient.)
PetscErrorCode PISMHydrology::wall_melt(IceModelVec2S &result) {
  PetscErrorCode ierr = result.set(0.0); CHKERRQ(ierr);
  return 0;
}


/*!
Checks \f$0 \le W_{til} \le W_{til}^{max} =\f$hydrology_tillwat_max.
 */
PetscErrorCode PISMHydrology::check_Wtil_bounds() {
  PetscErrorCode ierr;
  PetscReal tillwatmax = config.get("hydrology_tillwat_max");
  ierr = Wtil.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (Wtil(i,j) < 0.0) {
        PetscPrintf(grid.com,
           "PISMHydrology ERROR: negative till water effective layer thickness Wtil(i,j) = %.6f m\n"
           "            at (i,j)=(%d,%d)\n"
           "ENDING ... \n\n", Wtil(i,j), i, j);
        PISMEnd();
      }
      if (Wtil(i,j) > tillwatmax) {
        PetscPrintf(grid.com,
           "PISMHydrology ERROR: till water effective layer thickness Wtil(i,j) = %.6f m exceeds\n"
           "            hydrology_tillwat_max = %.6f at (i,j)=(%d,%d)\n"
           "ENDING ... \n\n", Wtil(i,j), tillwatmax, i, j);
        PISMEnd();
      }
    }
  }
  ierr = Wtil.end_access(); CHKERRQ(ierr);
  return 0;
}


//! Compute the total water input rate into the basal hydrology layer in the ice-covered region, allowing time-varying input from a file.
/*!
The user can specify the total of en- and supra-glacial drainage contributions
to subglacial hydrology in a time-dependent input file using option -input_to_bed.
This method includes that possible input along with `bmelt` to get the total water
input into the subglacial hydrology.

This method crops the input rate to the ice-covered region.  It
also uses hydrology_const_bmelt if that is requested.

Call this method using the current \e hydrology time step.  This method
may be called many times per IceModel time step.  See update() method
in derived classes of PISMHydrology.
 */
PetscErrorCode PISMHydrology::get_input_rate(
                  PetscReal hydro_t, PetscReal hydro_dt, IceModelVec2S &result) {
  PetscErrorCode ierr;
  bool      use_const   = config.get_flag("hydrology_use_const_bmelt");
  PetscReal const_bmelt = config.get("hydrology_const_bmelt");

  if (inputtobed != NULL) {
    ierr = inputtobed->update(hydro_t, hydro_dt); CHKERRQ(ierr);
    ierr = inputtobed->interp(hydro_t + hydro_dt/2.0); CHKERRQ(ierr);
    ierr = inputtobed->begin_access(); CHKERRQ(ierr);
  }
  ierr = bmelt->begin_access(); CHKERRQ(ierr);
  ierr = mask->begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  MaskQuery m(*mask);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (m.icy(i, j)) {
        result(i,j) = (use_const) ? const_bmelt : (*bmelt)(i,j);
        if (inputtobed != NULL) {
          result(i,j) += (*inputtobed)(i,j);
        }
      } else
        result(i,j) = 0.0;
    }
  }
  ierr = bmelt->end_access(); CHKERRQ(ierr);
  ierr = mask->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  if (inputtobed != NULL) {
    ierr = inputtobed->end_access(); CHKERRQ(ierr);
  }
  return 0;
}

