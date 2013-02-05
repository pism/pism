// Copyright (C) 2012-2013 PISM Authors
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

#include "PISMHydrology.hh"
#include "PISMVars.hh"
#include "pism_options.hh"
#include "Mask.hh"


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

  PetscErrorCode ierr1, ierr2;
  ierr1 = total_input.create(grid, "total_input_hydro", false);
  ierr2 = total_input.set_attrs("internal",
                         "workspace for total input rate into subglacial water layer",
                         "m s-1", "");
  if ((ierr1 != 0) || (ierr2 != 0)) {
      PetscPrintf(grid.com,
        "PISM ERROR: memory allocation failed in PISMHydrology constructor.\n");
      PISMEnd();
  }
}


PetscErrorCode PISMHydrology::init(PISMVars &vars) {
  PetscErrorCode ierr;
  string itbfilename;  // itb = input_to_bed
  bool itbfile_set, itbperiod_set, itbreference_set;
  PetscReal itbperiod_years = 0.0, itbreference_year = 0.0;

  ierr = verbPrintf(4, grid.com,
    "entering PISMHydrology::init() ...\n"); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, "",
            "Options controlling the base class PISMHydrology", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsIsSet("-report_mass_accounting",
      "Report to stdout on mass accounting in hydrology models", report_mass_accounting); CHKERRQ(ierr);
    ierr = PISMOptionsString("-input_to_bed_file",
      "A time- and space-dependent file with amount of water (depth per time) which should be put at the ice sheet bed at the given location at the given time",
      itbfilename, itbfile_set); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-input_to_bed_period",
      "The period (i.e. duration before repeat), in years, of -input_to_bed_file data",
      itbperiod_years, itbperiod_set); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-input_to_bed_reference_year",
      "The reference year for periodizing the -input_to_bed_file data",
      itbreference_year, itbreference_set); CHKERRQ(ierr);
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
    inputtobed_period = (itbperiod_set) ? grid.time->years_to_seconds(itbperiod_years) : 0.0;
    inputtobed_reference_time = (itbreference_set) ? grid.time->years_to_seconds(itbreference_year) : 0.0;

    unsigned int buffer_size = (unsigned int) config.get("climate_forcing_buffer_size"),
                 n_records = 1;

    PIO nc(grid.com, grid.rank, "netcdf3");
    ierr = nc.open(itbfilename, PISM_NOWRITE); CHKERRQ(ierr);
    ierr = nc.inq_nrecords("inputtobed", "", n_records); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);
    n_records = PetscMin(n_records, buffer_size);
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
    ierr = inputtobed->init(itbfilename); CHKERRQ(ierr);
  }
  return 0;
}


void PISMHydrology::get_diagnostics(map<string, PISMDiagnostic*> &dict) {
  dict["enwat"] = new PISMHydrology_enwat(this, grid, *variables);
  dict["bwp"] = new PISMHydrology_bwp(this, grid, *variables);
  dict["bwprel"] = new PISMHydrology_bwprel(this, grid, *variables);
  dict["effbwp"] = new PISMHydrology_effbwp(this, grid, *variables);
  dict["hydroinput"] = new PISMHydrology_hydroinput(this, grid, *variables);
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


//! Set the englacial storage to zero.  (The most basic subglacial hydrologies have no englacial storage.)
PetscErrorCode PISMHydrology::englacial_water_thickness(IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = result.set(0.0); CHKERRQ(ierr);
  return 0;
}


//! The only reason to restrict the time step taken by the calling model (i.e. IceModel) is if there is a time-dependent input file IceModelVec2T *inputtobed.
PetscErrorCode PISMHydrology::max_timestep(PetscReal my_t, PetscReal &my_dt, bool &restrict_dt) {
  if (inputtobed == NULL) {
    my_dt = -1;
    restrict_dt = false;
  } else {
    // "periodize" the forcing data
    my_t = grid.time->mod(my_t - inputtobed_reference_time, inputtobed_period);
    my_dt = inputtobed->max_timestep(my_t);
    // if the user specifies periodized forcing, limit time-steps so that PISM
    // never tries to average data over an interval that begins in one period and
    // ends in the next one.
    if (inputtobed_period > 1e-6)
      my_dt = PetscMin(my_dt, inputtobed_period - my_t);
    restrict_dt = (my_dt > 0);
  }
  return 0;
}


//! Compute the total water input rate into the basal hydrology layer in the ice-covered region, allowing time-varying input from a file.
/*!
The user can specify the total of en- and supra-glacial drainage contributions
to subglacial hydrology in a time-dependent input file using option -input_to_bed.
This method includes that possible input along with \c bmelt to get the total water
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
    ierr = inputtobed->at_time(hydro_t + hydro_dt/2.0); CHKERRQ(ierr);
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


//! Update the water thickness based on boundary requirements.  Do mass accounting.
/*!
At ice free locations and ocean locations we require that the water thickness
is zero at the end of each time step.  Also we require that any negative water
thicknesses be set to zero (i.e. projection to enforce \f$W\ge 0\f$).

This method takes care of these requirements by altering Wnew appropriately.
And we account for the mass changes that these alterations represent.
 */
PetscErrorCode PISMHydrology::boundary_mass_changes(IceModelVec2S &Wnew,
            PetscReal &icefreelost, PetscReal &oceanlost, PetscReal &negativegain) {
  PetscErrorCode ierr;
  PetscReal fresh_water_density = config.get("fresh_water_density");
  PetscReal my_icefreelost = 0.0, my_oceanlost = 0.0, my_negativegain = 0.0;
  MaskQuery M(*mask);
  ierr = Wnew.begin_access(); CHKERRQ(ierr);
  ierr = mask->begin_access(); CHKERRQ(ierr);
  ierr = cellarea->begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscReal dmassdz = (*cellarea)(i,j) * fresh_water_density; // kg m-1
      if (Wnew(i,j) < 0.0) {
        my_negativegain += -Wnew(i,j) * dmassdz;
        Wnew(i,j) = 0.0;
      }
      if (M.ice_free_land(i,j) && (Wnew(i,j) > 0.0)) {
        my_icefreelost += Wnew(i,j) * dmassdz;
        Wnew(i,j) = 0.0;
      }
      if (M.ocean(i,j) && (Wnew(i,j) > 0.0)) {
        my_oceanlost += Wnew(i,j) * dmassdz;
        Wnew(i,j) = 0.0;
      }
    }
  }
  ierr = Wnew.end_access(); CHKERRQ(ierr);
  ierr = mask->end_access(); CHKERRQ(ierr);
  ierr = cellarea->end_access(); CHKERRQ(ierr);

  // make global over all proc domains (i.e. whole glacier/ice sheet)
  ierr = PISMGlobalSum(&my_icefreelost, &icefreelost, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&my_oceanlost, &oceanlost, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&my_negativegain, &negativegain, grid.com); CHKERRQ(ierr);

  // this reporting is redundant for the simpler models but shows short time step
  // reporting for nontrivially-distributed (possibly adaptive) hydrology models
  ierr = verbPrintf(4, grid.com,
    "  mass losses in hydrology time step:\n"
    "     land margin loss = %.3e kg, ocean margin loss = %.3e kg, (W<0) gain = %.3e kg\n",
    icefreelost, oceanlost, negativegain); CHKERRQ(ierr);
  return 0;
}


PISMTillCanHydrology::PISMTillCanHydrology(IceGrid &g, const NCConfigVariable &conf,
                                           bool Whasghosts)
    : PISMHydrology(g, conf)
{
  if (allocate(Whasghosts) != 0) {
    PetscPrintf(grid.com,
      "PISM ERROR: allocation failed in PISMTillCanHydrology constructor.\n");
    PISMEnd();
  }
}


PetscErrorCode PISMTillCanHydrology::allocate(bool Whasghosts) {
  PetscErrorCode ierr;
  if (Whasghosts) {
    ierr = W.create(grid, "bwat", true, 1); CHKERRQ(ierr);
  } else {
    ierr = W.create(grid, "bwat", false); CHKERRQ(ierr);
  }
  ierr = W.set_attrs("model_state",
                     "thickness of subglacial water layer",
                     "m", ""); CHKERRQ(ierr);
  ierr = W.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PISMTillCanHydrology::init(PISMVars &vars) {
  PetscErrorCode ierr;
  ierr = verbPrintf(2, grid.com,
    "* Initializing the 'tillcan' subglacial hydrology model...\n"); CHKERRQ(ierr);
  ierr = PISMHydrology::init(vars); CHKERRQ(ierr);

  // initialize water layer thickness from the context if present,
  //   otherwise from -i or -boot_file, otherwise with constant value
  bool i_set, bootstrap;
  ierr = PetscOptionsBegin(grid.com, "",
            "Options controlling the 'tillcan' subglacial hydrology model", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsIsSet("-i", "PISM input file", i_set); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-boot_file", "PISM bootstrapping file",
                            bootstrap); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);
  IceModelVec2S *W_input = dynamic_cast<IceModelVec2S*>(vars.get("bwat"));
  if (W_input != NULL) { // a variable called "bwat" is already in context
    ierr = W.copy_from(*W_input); CHKERRQ(ierr);
  } else if (i_set || bootstrap) {
    string filename;
    int start;
    ierr = find_pism_input(filename, bootstrap, start); CHKERRQ(ierr);
    if (i_set) {
      ierr = W.read(filename, start); CHKERRQ(ierr);
    } else {
      ierr = W.regrid(filename,
                      config.get("bootstrapping_bwat_value_no_var")); CHKERRQ(ierr);
    }
  } else {
    ierr = W.set(config.get("bootstrapping_bwat_value_no_var")); CHKERRQ(ierr);
  }

  // whether or not we could initialize from file, we could be asked to regrid from file
  ierr = regrid(W); CHKERRQ(ierr);

  // add bwat to the variables in the context if it is not already there
  if (vars.get("bwat") == NULL) {
    ierr = vars.add(W); CHKERRQ(ierr);
  }
  return 0;
}


void PISMTillCanHydrology::add_vars_to_output(string /*keyword*/, map<string,NCSpatialVariable> &result) {
  result["bwat"] = W.get_metadata();
}


PetscErrorCode PISMTillCanHydrology::define_variables(set<string> vars, const PIO &nc,
                                                 PISM_IO_Type nctype) {
  PetscErrorCode ierr;
  if (set_contains(vars, "bwat")) {
    ierr = W.define(nc, nctype); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode PISMTillCanHydrology::write_variables(set<string> vars, const PIO &nc) {
  PetscErrorCode ierr;
  if (set_contains(vars, "bwat")) {
    ierr = W.write(nc); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode PISMTillCanHydrology::subglacial_water_thickness(IceModelVec2S &result) {
  PetscErrorCode ierr = W.copy_to(result); CHKERRQ(ierr);
  return 0;
}


//! Computes pressure diagnostically.
/*!
  \f[ P = \lambda P_o \max\{1,W / W_{crit}\} \f]
where \f$\lambda\f$=till_pw_fraction, \f$P_o = \rho_i g H\f$, \f$W_{crit}\f$=hydrology_bwat_max.
 */
PetscErrorCode PISMTillCanHydrology::subglacial_water_pressure(IceModelVec2S &result) {
  PetscErrorCode ierr;

#if (PISM_DEBUG==1)
  ierr = check_W_bounds(); CHKERRQ(ierr); // check:  W \le bwat_max = Wcrit
#endif

  ierr = overburden_pressure(result); CHKERRQ(ierr);

  double bwat_max = config.get("hydrology_bwat_max"),
    till_pw_fraction = config.get("till_pw_fraction");

  ierr = W.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      // P = lambda (W/W_0) P_o
      result(i,j) = till_pw_fraction * (W(i,j) / bwat_max) * result(i,j);
    }
  }
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = W.end_access(); CHKERRQ(ierr);
  return 0;
}


/*!
Checks \f$0 \le W \le W_{crit} =\f$hydrology_bwat_max.
 */
PetscErrorCode PISMTillCanHydrology::check_W_bounds() {
  PetscErrorCode ierr;
  PetscReal bwat_max = config.get("hydrology_bwat_max");
  ierr = W.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (W(i,j) < 0.0) {
        PetscPrintf(grid.com,
           "PISMTillCanHydrology ERROR: negative subglacial water layer thickness W(i,j) = %.6f m\n"
           "            at (i,j)=(%d,%d)\n"
           "ENDING ... \n\n", W(i,j),i,j);
        PISMEnd();
      }
      if (W(i,j) > bwat_max) {
        PetscPrintf(grid.com,
           "PISMTillCanHydrology ERROR: subglacial water layer thickness W(i,j) = %.6f m exceeds\n"
           "            hydrology_bwat_max = %.6f at (i,j)=(%d,%d)\n"
           "ENDING ... \n\n", W(i,j),bwat_max,i,j);
        PISMEnd();
      }
    }
  }
  ierr = W.end_access(); CHKERRQ(ierr);
  return 0;
}


//! Update the water thickness.
/*!
Does an explicit (forward Euler) step of the integration
  \f[ \frac{dW}{dt} = \text{bmelt} - C \f]
but subject to the inequalities
  \f[ 0 \le W \le W_0 \f]
where \f$C=\f$hydrology_bwat_decay_rate and \f$W_0\f$=hydrology_bwat_max.
 */
PetscErrorCode PISMTillCanHydrology::update(PetscReal icet, PetscReal icedt) {
  // if asked for the identical time interval as last time, then do nothing
  if ((fabs(icet - t) < 1e-6) && (fabs(icedt - dt) < 1e-6))
    return 0;
  t = icet;
  dt = icedt;

  PetscErrorCode ierr;

  ierr = get_input_rate(icet,icedt,total_input); CHKERRQ(ierr);

  PetscReal bwat_max        = config.get("hydrology_bwat_max"),
            bwat_decay_rate = config.get("hydrology_bwat_decay_rate");
  PetscReal icefreelost = 0, oceanlost = 0, negativegain = 0;

  ierr = W.begin_access(); CHKERRQ(ierr);
  ierr = total_input.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      W(i,j) = pointwise_update(W(i,j), total_input(i,j) * icedt, bwat_decay_rate * icedt, bwat_max);
    }
  }
  ierr = W.end_access(); CHKERRQ(ierr);
  ierr = total_input.end_access(); CHKERRQ(ierr);

  // following should *not* alter W, and it should report all zeros by design;
  // this hydrology is *not* distributed
  ierr = boundary_mass_changes(W,icefreelost,oceanlost,negativegain); CHKERRQ(ierr);
  if (report_mass_accounting) {
    ierr = verbPrintf(2, grid.com,
      " 'tillcan' hydrology mass losses:\n"
      "     ice free land lost = %.3e kg, ocean lost = %.3e kg, negative bmelt gain = %.3e kg\n",
      icefreelost, oceanlost, negativegain); CHKERRQ(ierr);
  }
  return 0;
}


PISMDiffuseOnlyHydrology::PISMDiffuseOnlyHydrology(IceGrid &g, const NCConfigVariable &conf)
    : PISMTillCanHydrology(g, conf, true)
{
  if (allocateWnew() != 0) {
    PetscPrintf(grid.com,
      "PISM ERROR: allocation of Wnew failed in PISMDiffuseOnlyHydrology constructor.\n");
    PISMEnd();
  }
}


PetscErrorCode PISMDiffuseOnlyHydrology::allocateWnew() {
  PetscErrorCode ierr;
  // also need temporary space during update
  ierr = Wnew.create(grid, "Wnew-internal", false); CHKERRQ(ierr);
  ierr = Wnew.set_attrs("internal",
                     "new thickness of subglacial water layer during update",
                     "m", ""); CHKERRQ(ierr);
  ierr = Wnew.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PISMDiffuseOnlyHydrology::init(PISMVars &vars) {
  PetscErrorCode ierr;
  ierr = PISMTillCanHydrology::init(vars); CHKERRQ(ierr);
  ierr = verbPrintf(2, grid.com,
    "  using the diffusive water layer variant ...\n"); CHKERRQ(ierr);
  return 0;
}


//! Explicit time step for diffusion of subglacial water layer bwat.
/*!
This model adds a contrived lateral diffusion to the PISMTillCanHydrology
model.  See equation (11) in \ref BBssasliding , namely
  \f[W_t = K \nabla^2 W.\f]
The diffusion constant \f$K\f$ is chosen so that the fundamental solution (Green's
function) of this equation has standard deviation \f$\sigma=L\f$ at time t=\c diffusion_time.
Note that \f$2 \sigma^2 = 4 K t\f$.

The time step restriction for the explicit method for this equation is believed
to be so rare that if it is triggered there is a stdout warning.
 */
PetscErrorCode PISMDiffuseOnlyHydrology::update(PetscReal icet, PetscReal icedt) {
  // if asked for the identical time interval as last time, then do nothing
  if ((fabs(icet - t) < 1e-6) && (fabs(icedt - dt) < 1e-6))
    return 0;
  t = icet;
  dt = icedt;

  PetscErrorCode ierr;
  const PetscReal L = config.get("hydrology_bwat_diffusion_distance");
  if (L <= 0.0)  {
    ierr = PISMTillCanHydrology::update(icet,icedt); CHKERRQ(ierr);
    return 0;
  }

  const PetscReal
    diffusion_time  = config.get("hydrology_bwat_diffusion_time", "years", "seconds"), // convert to seconds
    bwat_max        = config.get("hydrology_bwat_max"),
    bwat_decay_rate = config.get("hydrology_bwat_decay_rate"),
    K               = L * L / (2.0 * diffusion_time);

  PetscReal hdt;
  PetscInt NN;
  hdt = (1.0 / (grid.dx*grid.dx)) + (1.0 / (grid.dy*grid.dy));
  hdt = 1.0 / (2.0 * K * hdt);
  NN = int(ceil(dt / hdt));
  hdt = dt / NN;
  if (NN > 1) {
    verbPrintf(2,grid.com,
      "PISMDiffuseOnlyHydrology WARNING: more than one time step per ice dynamics time step\n"
      "   ... NN = %d > 1 ... THIS IS BELIEVED TO BE RARE\n",NN);
  }

  PetscReal icefreelost = 0, oceanlost = 0, negativegain = 0;

  PetscReal  Rx = K * hdt / (grid.dx * grid.dx),
             Ry = K * hdt / (grid.dy * grid.dy),
             oneM4R = 1.0 - 2.0 * Rx - 2.0 * Ry;
  for (PetscInt n=0; n<NN; ++n) {
    if ((inputtobed != NULL) || (n==0)) {
      ierr = get_input_rate(icet + n * hdt, hdt, total_input); CHKERRQ(ierr);
    }

    // time-splitting: first, Euler step on source terms
    ierr = W.begin_access(); CHKERRQ(ierr);
    ierr = total_input.begin_access(); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        W(i,j) = pointwise_update(W(i,j), total_input(i,j) * icedt, bwat_decay_rate * icedt, bwat_max);
      }
    }
    ierr = W.end_access(); CHKERRQ(ierr);
    ierr = total_input.end_access(); CHKERRQ(ierr);

    // valid ghosts for diffusion below
    ierr = W.update_ghosts(); CHKERRQ(ierr);

    // time-splitting: second, diffusion by first-order explicit
    ierr = W.begin_access(); CHKERRQ(ierr);
    ierr = Wnew.begin_access(); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        Wnew(i,j) = oneM4R * W(i,j) + Rx * (W(i+1,j  ) + W(i-1,j  ))
                                    + Ry * (W(i  ,j+1) + W(i  ,j-1));
        // no check of upper bound here because maximum principle applies to step
      }
    }
    ierr = W.end_access(); CHKERRQ(ierr);
    ierr = Wnew.end_access(); CHKERRQ(ierr);

    ierr = boundary_mass_changes(Wnew,icefreelost,oceanlost,negativegain); CHKERRQ(ierr);
    if (report_mass_accounting) {
      ierr = verbPrintf(2, grid.com,
        " 'diffuseonly' hydrology mass losses:\n"
        "     ice free land lost = %.3e kg, ocean lost = %.3e kg, negative bmelt gain = %.3e kg\n",
        icefreelost, oceanlost, negativegain); CHKERRQ(ierr);
    }

    ierr = Wnew.update_ghosts(W); CHKERRQ(ierr);
  }
  return 0;
}


PISMHydrology_enwat::PISMHydrology_enwat(PISMHydrology *m, IceGrid &g, PISMVars &my_vars)
    : PISMDiag<PISMHydrology>(m, g, my_vars) {
  vars[0].init_2d("enwat", grid);
  set_attrs("effective thickness of englacial water", "", "m", "m", 0);
}


PetscErrorCode PISMHydrology_enwat::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "enwat", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);
  ierr = model->englacial_water_thickness(*result); CHKERRQ(ierr);
  output = result;
  return 0;
}


PISMHydrology_bwp::PISMHydrology_bwp(PISMHydrology *m, IceGrid &g, PISMVars &my_vars)
    : PISMDiag<PISMHydrology>(m, g, my_vars) {
  vars[0].init_2d("bwp", grid);
  set_attrs("pressure of water in subglacial layer", "", "Pa", "Pa", 0);
}


PetscErrorCode PISMHydrology_bwp::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "bwp", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);
  result->write_in_glaciological_units = true;
  ierr = model->subglacial_water_pressure(*result); CHKERRQ(ierr);
  output = result;
  return 0;
}


PISMHydrology_bwprel::PISMHydrology_bwprel(PISMHydrology *m, IceGrid &g, PISMVars &my_vars)
    : PISMDiag<PISMHydrology>(m, g, my_vars) {
  vars[0].init_2d("bwprel", grid);
  set_attrs("pressure of water in subglacial layer as fraction of the overburden pressure", "",
            "", "", 0);
  vars[0].set("_FillValue", grid.config.get("fill_value"));
}


PetscErrorCode PISMHydrology_bwprel::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  PetscReal fill = grid.config.get("fill_value");
  IceModelVec2S *Po     = new IceModelVec2S,
                *result = new IceModelVec2S;
  ierr = result->create(grid, "bwprel", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);
  ierr = Po->create(grid, "Po_temporary", false); CHKERRQ(ierr);
  ierr = Po->set_metadata(vars[0], 0); CHKERRQ(ierr);

  ierr = model->subglacial_water_pressure(*result); CHKERRQ(ierr);
  ierr = model->overburden_pressure(*Po); CHKERRQ(ierr);

  ierr = result->begin_access(); CHKERRQ(ierr);
  ierr = Po->begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if ((*Po)(i,j) > 0.0)
        (*result)(i,j) /= (*Po)(i,j);
      else
        (*result)(i,j) = fill;
    }
  }
  ierr = result->end_access(); CHKERRQ(ierr);
  ierr = Po->end_access(); CHKERRQ(ierr);

  output = result;
  return 0;
}


PISMHydrology_effbwp::PISMHydrology_effbwp(PISMHydrology *m, IceGrid &g, PISMVars &my_vars)
    : PISMDiag<PISMHydrology>(m, g, my_vars) {
  vars[0].init_2d("effbwp", grid);
  set_attrs("effective pressure of water in subglacial layer (overburden pressure minus water pressure)",
            "", "Pa", "Pa", 0);
}


PetscErrorCode PISMHydrology_effbwp::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec2S *P      = new IceModelVec2S,
                *result = new IceModelVec2S;
  ierr = result->create(grid, "effbwp", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);
  ierr = P->create(grid, "P_temporary", false); CHKERRQ(ierr);
  ierr = P->set_metadata(vars[0], 0); CHKERRQ(ierr);

  ierr = model->subglacial_water_pressure(*P); CHKERRQ(ierr);
  ierr = model->overburden_pressure(*result); CHKERRQ(ierr);
  ierr = result->add(-1.0,*P); CHKERRQ(ierr);  // result <-- result + (-1.0) P = Po - P

  output = result;
  return 0;
}


PISMHydrology_hydroinput::PISMHydrology_hydroinput(PISMHydrology *m, IceGrid &g, PISMVars &my_vars)
    : PISMDiag<PISMHydrology>(m, g, my_vars) {
  vars[0].init_2d("hydroinput", grid);
  set_attrs("total water input into subglacial hydrology layer",
            "", "m s-1", "m a-1", 0);
}


PetscErrorCode PISMHydrology_hydroinput::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "hydroinput", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);
  result->write_in_glaciological_units = true;
  // the value reported diagnostically is merely the last value filled
  ierr = (model->total_input).copy_to(*result); CHKERRQ(ierr);
  output = result;
  return 0;
}

