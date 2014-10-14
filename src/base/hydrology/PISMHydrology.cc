// Copyright (C) 2012-2014 PISM Authors
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
#include "error_handling.hh"

namespace pism {

Hydrology::Hydrology(IceGrid &g, const Config &conf)
  : Component_TS(g, conf)
{
  thk   = NULL;
  bed   = NULL;
  cellarea = NULL;
  bmelt = NULL;
  mask  = NULL;
  inputtobed = NULL;
  variables = NULL;
  hold_bmelt = false;

  PetscErrorCode ierr1, ierr2;
  ierr1 = total_input.create(grid, "total_input", WITHOUT_GHOSTS);
  ierr2 = total_input.set_attrs("internal",
                         "hydrology model workspace for total input rate into subglacial water layer",
                         "m s-1", "");
  if ((ierr1 != 0) || (ierr2 != 0)) {
      PetscPrintf(grid.com,
        "PISM ERROR: memory allocation failed in Hydrology constructor (total_input).\n");
      PISMEnd();
  }

  ierr1 = bmelt_local.create(grid, "bmelt", WITHOUT_GHOSTS);
  ierr2 = bmelt_local.set_attrs("internal",
                         "hydrology model workspace for bmelt",
                         "m s-1", "");
  if ((ierr1 != 0) || (ierr2 != 0)) {
      PetscPrintf(grid.com,
        "PISM ERROR: memory allocation failed in Hydrology constructor (bmelt_local).\n");
      PISMEnd();
  }

  // *all* Hydrology classes have layer of water stored in till as a state variable
  ierr1 = Wtil.create(grid, "tillwat", WITHOUT_GHOSTS);
  ierr2 = Wtil.set_attrs("model_state",
                     "effective thickness of subglacial water stored in till",
                     "m", "");
  Wtil.metadata().set_double("valid_min", 0.0);
  if ((ierr1 != 0) || (ierr2 != 0)) {
      PetscPrintf(grid.com,
        "PISM ERROR: memory allocation failed in Hydrology constructor (Wtil).\n");
      PISMEnd();
  }
}


Hydrology::~Hydrology() {
  // empty
}


PetscErrorCode Hydrology::init(Vars &vars) {
  PetscErrorCode ierr;
  std::string itbfilename,  // itb = input_to_bed
              bmeltfilename;
  bool bmeltfile_set, itbfile_set, itbperiod_set, itbreference_set;
  bool i_set, bootstrap;
  double itbperiod_years = 0.0, itbreference_year = 0.0;

  ierr = verbPrintf(4, grid.com,
    "entering Hydrology::init() ...\n"); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, "",
            "Options controlling the base class Hydrology", ""); CHKERRQ(ierr);
  {
    ierr = OptionsString("-hydrology_bmelt_file",
      "Read time-independent values for bmelt from a file; replaces bmelt computed through conservation of energy",
      bmeltfilename, bmeltfile_set); CHKERRQ(ierr);
    ierr = OptionsString("-hydrology_input_to_bed_file",
      "A time- and space-dependent file with amount of water (depth per time) which should be added to the amount of water at the ice sheet bed at the given location at the given time; adds to bmelt",
      itbfilename, itbfile_set); CHKERRQ(ierr);
    ierr = OptionsReal("-hydrology_input_to_bed_period",
      "The period (i.e. duration before repeat), in years, of -hydrology_input_to_bed_file data",
      itbperiod_years, itbperiod_set); CHKERRQ(ierr);
    ierr = OptionsReal("-hydrology_input_to_bed_reference_year",
      "The reference year for periodizing the -hydrology_input_to_bed_file data",
      itbreference_year, itbreference_set); CHKERRQ(ierr);
    ierr = OptionsIsSet("-i", "PISM input file", i_set); CHKERRQ(ierr);
    ierr = OptionsIsSet("-boot_file", "PISM bootstrapping file",
                            bootstrap); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  variables = &vars;

  // the following are IceModelVec pointers into IceModel generally and are read by code in the
  // update() method at the current Hydrology time

  thk      = vars.get_2d_scalar("thk");
  bed      = vars.get_2d_scalar("topg");
  bmelt    = vars.get_2d_scalar("bmelt");
  cellarea = vars.get_2d_scalar("cell_area");
  mask     = vars.get_2d_mask("mask");


  if (bmeltfile_set) {
    ierr = verbPrintf(2, grid.com,
       "  option -hydrology_bmelt_file seen; reading bmelt from '%s'.\n", bmeltfilename.c_str()); CHKERRQ(ierr);
    ierr = bmelt_local.regrid(bmeltfilename, CRITICAL); CHKERRQ(ierr);
    hold_bmelt = true;
  }


  if (itbfile_set) {
    inputtobed_period = (itbperiod_set) ? itbperiod_years : 0.0;
    inputtobed_reference_time = (itbreference_set) ? grid.convert(itbreference_year, "years", "seconds") : 0.0;

    unsigned int buffer_size = (unsigned int) config.get("climate_forcing_buffer_size"),
                 n_records = 1;

    PIO nc(grid.com, "netcdf3", grid.get_unit_system());
    ierr = nc.open(itbfilename, PISM_READONLY); CHKERRQ(ierr);
    ierr = nc.inq_nrecords("inputtobed", "", n_records); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);

    // if -..._period is not set, make n_records the minimum of the
    // buffer size and the number of available records. Otherwise try
    // to keep all available records in memory.
    if (itbperiod_set == false) {
      n_records = PetscMin(n_records, buffer_size);
    }

    if (n_records == 0) {
      PetscPrintf(grid.com, "PISM ERROR: can't find 'inputtobed' in -hydrology_input_to_bed file with name '%s'.\n",
                  itbfilename.c_str());
      PISMEnd();
    }

    ierr = verbPrintf(2,grid.com,
      "    option -hydrology_input_to_bed_file seen ... creating 'inputtobed' variable ...\n"); CHKERRQ(ierr);
    ierr = verbPrintf(2,grid.com,
      "    allocating buffer space for n = %d 'inputtobed' records ...\n", n_records); CHKERRQ(ierr);
    inputtobed = new IceModelVec2T;
    inputtobed->set_n_records(n_records);
    ierr = inputtobed->create(grid, "inputtobed", WITHOUT_GHOSTS); CHKERRQ(ierr);
    ierr = inputtobed->set_attrs("climate_forcing",
                                 "amount of water (depth per time like bmelt) which should be put at the ice sheet bed",
                                 "m s-1", ""); CHKERRQ(ierr);
    ierr = verbPrintf(2,grid.com,
      "    reading 'inputtobed' variable from file '%s' ...\n",itbfilename.c_str()); CHKERRQ(ierr);
    ierr = inputtobed->init(itbfilename, inputtobed_period, inputtobed_reference_time); CHKERRQ(ierr);
  }

  // initialize till water layer thickness from the context if present,
  //   otherwise from -i or -boot_file, otherwise with constant value
  IceModelVec2S *Wtil_input = NULL;

  try {
    // FIXME: this is not an exceptional situation...
    Wtil_input = vars.get_2d_scalar("tillwat");
    ierr = Wtil.copy_from(*Wtil_input); CHKERRQ(ierr);
  } catch (RuntimeError) {
    if (i_set || bootstrap) {
      std::string filename;
      int start;
      ierr = find_pism_input(filename, bootstrap, start); CHKERRQ(ierr);
      if (i_set) {
        ierr = Wtil.read(filename, start); CHKERRQ(ierr);
      } else {
        ierr = Wtil.regrid(filename, OPTIONAL,
                           config.get("bootstrapping_tillwat_value_no_var")); CHKERRQ(ierr);
      }
    } else {
      ierr = Wtil.set(config.get("bootstrapping_tillwat_value_no_var")); CHKERRQ(ierr);
    }
  }

  // whether or not we could initialize from file, we could be asked to regrid from file
  ierr = regrid("Hydrology", &Wtil); CHKERRQ(ierr);
  return 0;
}


void Hydrology::get_diagnostics(std::map<std::string, Diagnostic*> &dict,
                                    std::map<std::string, TSDiagnostic*> &/*ts_dict*/) {
  dict["bwat"] = new Hydrology_bwat(this, grid, *variables);
  dict["bwp"] = new Hydrology_bwp(this, grid, *variables);
  dict["bwprel"] = new Hydrology_bwprel(this, grid, *variables);
  dict["effbwp"] = new Hydrology_effbwp(this, grid, *variables);
  dict["hydrobmelt"] = new Hydrology_hydrobmelt(this, grid, *variables);
  dict["hydroinput"] = new Hydrology_hydroinput(this, grid, *variables);
  dict["wallmelt"] = new Hydrology_wallmelt(this, grid, *variables);
}


void Hydrology::add_vars_to_output(const std::string &/*keyword*/, std::set<std::string> &result) {
  result.insert("tillwat");
}


PetscErrorCode Hydrology::define_variables(const std::set<std::string> &vars, const PIO &nc,
                                               IO_Type nctype) {
  PetscErrorCode ierr;
  if (set_contains(vars, "tillwat")) {
    ierr = Wtil.define(nc, nctype); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode Hydrology::write_variables(const std::set<std::string> &vars, const PIO &nc) {
  PetscErrorCode ierr;
  if (set_contains(vars, "tillwat")) {
    ierr = Wtil.write(nc); CHKERRQ(ierr);
  }
  return 0;
}


//! Update the overburden pressure from ice thickness.
/*!
Uses the standard hydrostatic (shallow) approximation of overburden pressure,
  \f[ P_0 = \rho_i g H \f]
Accesses H=thk from Vars, which points into IceModel.
 */
PetscErrorCode Hydrology::overburden_pressure(IceModelVec2S &result) {
  PetscErrorCode ierr;
  // FIXME issue #15
  ierr = result.copy_from(*thk); CHKERRQ(ierr);  // copies into ghosts if result has them
  ierr = result.scale(config.get("ice_density") * config.get("standard_gravity")); CHKERRQ(ierr);
  return 0;
}


//! Return the effective thickness of the water stored in till.
PetscErrorCode Hydrology::till_water_thickness(IceModelVec2S &result) {
  PetscErrorCode ierr = Wtil.copy_to(result); CHKERRQ(ierr);
  return 0;
}


//! Set the wall melt rate to zero.  (The most basic subglacial hydrologies have no lateral flux or potential gradient.)
PetscErrorCode Hydrology::wall_melt(IceModelVec2S &result) {
  PetscErrorCode ierr = result.set(0.0); CHKERRQ(ierr);
  return 0;
}


/*!
Checks \f$0 \le W_{til} \le W_{til}^{max} =\f$hydrology_tillwat_max.
 */
PetscErrorCode Hydrology::check_Wtil_bounds() {
  double tillwat_max = config.get("hydrology_tillwat_max");

  IceModelVec::AccessList list(Wtil);
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (Wtil(i,j) < 0.0) {
      PetscPrintf(grid.com,
                  "Hydrology ERROR: negative till water effective layer thickness Wtil(i,j) = %.6f m\n"
                  "            at (i,j)=(%d,%d)\n"
                  "ENDING ... \n\n", Wtil(i,j), i, j);
      PISMEnd();
    }
    if (Wtil(i,j) > tillwat_max) {
      PetscPrintf(grid.com,
                  "Hydrology ERROR: till water effective layer thickness Wtil(i,j) = %.6f m exceeds\n"
                  "            hydrology_tillwat_max = %.6f at (i,j)=(%d,%d)\n"
                  "ENDING ... \n\n", Wtil(i,j), tillwat_max, i, j);
      PISMEnd();
    }
  }
  return 0;
}


//! Compute the total water input rate into the basal hydrology layer in the ice-covered region, allowing time-varying input from a file.
/*!
The user can specify the total of en- and supra-glacial drainage contributions
to subglacial hydrology in a time-dependent input file using option -hydrology_input_to_bed.
This method includes that possible input along with `bmelt` to get the total water
input into the subglacial hydrology.

This method crops the input rate to the ice-covered region.  It
also uses hydrology_const_bmelt if that is requested.

Call this method using the current \e hydrology time step.  This method
may be called many times per IceModel time step.  See update() method
in derived classes of Hydrology.
 */
PetscErrorCode Hydrology::get_input_rate(double hydro_t, double hydro_dt,
                                         IceModelVec2S &result) {
  PetscErrorCode ierr;
  bool   use_const   = config.get_flag("hydrology_use_const_bmelt");
  double const_bmelt = config.get("hydrology_const_bmelt");

  IceModelVec::AccessList list;
  if (inputtobed != NULL) {
    ierr = inputtobed->update(hydro_t, hydro_dt); CHKERRQ(ierr);
    ierr = inputtobed->interp(hydro_t + hydro_dt/2.0); CHKERRQ(ierr);
    list.add(*inputtobed);
  }

  if (!hold_bmelt) {
    ierr = bmelt->copy_to(bmelt_local); CHKERRQ(ierr);
  }

  list.add(bmelt_local);
  list.add(*mask);
  list.add(result);
  MaskQuery m(*mask);
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m.icy(i, j)) {
      result(i,j) = (use_const) ? const_bmelt : bmelt_local(i,j);
      if (inputtobed != NULL) {
        result(i,j) += (*inputtobed)(i,j);
      }
    } else
      result(i,j) = 0.0;
  }
  return 0;
}


} // end of namespace pism
