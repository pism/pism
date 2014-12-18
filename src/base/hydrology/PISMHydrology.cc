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

Hydrology::Hydrology(IceGrid &g)
  : Component_TS(g)
{
  thk        = NULL;
  bed        = NULL;
  cellarea   = NULL;
  bmelt      = NULL;
  mask       = NULL;
  inputtobed = NULL;
  hold_bmelt = false;

  total_input.create(m_grid, "total_input", WITHOUT_GHOSTS);
  total_input.set_attrs("internal",
                        "hydrology model workspace for total input rate into subglacial water layer",
                        "m s-1", "");

  bmelt_local.create(m_grid, "bmelt", WITHOUT_GHOSTS);
  bmelt_local.set_attrs("internal",
                        "hydrology model workspace for bmelt",
                        "m s-1", "");

  // *all* Hydrology classes have layer of water stored in till as a state variable
  Wtil.create(m_grid, "tillwat", WITHOUT_GHOSTS);
  Wtil.set_attrs("model_state",
                 "effective thickness of subglacial water stored in till",
                 "m", "");
  Wtil.metadata().set_double("valid_min", 0.0);
}


Hydrology::~Hydrology() {
  // empty
}


void Hydrology::init() {
  std::string itbfilename,  // itb = input_to_bed
              bmeltfilename;
  bool bmeltfile_set, itbfile_set, itbperiod_set, itbreference_set;
  bool i_set, bootstrap;
  double itbperiod_years = 0.0, itbreference_year = 0.0;

  verbPrintf(4, m_grid.com,
             "entering Hydrology::init() ...\n");

  {
    OptionsString("-hydrology_bmelt_file",
                  "Read time-independent values for bmelt from a file; replaces bmelt computed through conservation of energy",
                  bmeltfilename, bmeltfile_set);
    OptionsString("-hydrology_input_to_bed_file",
                  "A time- and space-dependent file with amount of water (depth per time) which should be added to the amount of water at the ice sheet bed at the given location at the given time; adds to bmelt",
                  itbfilename, itbfile_set);
    OptionsReal("-hydrology_input_to_bed_period",
                "The period (i.e. duration before repeat), in years, of -hydrology_input_to_bed_file data",
                itbperiod_years, itbperiod_set);
    OptionsReal("-hydrology_input_to_bed_reference_year",
                "The reference year for periodizing the -hydrology_input_to_bed_file data",
                itbreference_year, itbreference_set);
    i_set = OptionsIsSet("-i", "PISM input file");
    bootstrap = OptionsIsSet("-boot_file", "PISM bootstrapping file");
  }

  // the following are IceModelVec pointers into IceModel generally and are read by code in the
  // update() method at the current Hydrology time

  thk      = m_grid.variables().get_2d_scalar("thk");
  bed      = m_grid.variables().get_2d_scalar("topg");
  bmelt    = m_grid.variables().get_2d_scalar("bmelt");
  cellarea = m_grid.variables().get_2d_scalar("cell_area");
  mask     = m_grid.variables().get_2d_mask("mask");


  if (bmeltfile_set) {
    verbPrintf(2, m_grid.com,
               "  option -hydrology_bmelt_file seen; reading bmelt from '%s'.\n", bmeltfilename.c_str());
    bmelt_local.regrid(bmeltfilename, CRITICAL);
    hold_bmelt = true;
  }


  if (itbfile_set) {
    inputtobed_period = (itbperiod_set) ? itbperiod_years : 0.0;
    inputtobed_reference_time = (itbreference_set) ? m_grid.convert(itbreference_year, "years", "seconds") : 0.0;

    unsigned int buffer_size = (unsigned int) m_config.get("climate_forcing_buffer_size");

    PIO nc(m_grid.com, "netcdf3", m_grid.config.get_unit_system());
    nc.open(itbfilename, PISM_READONLY);
    unsigned int n_records = nc.inq_nrecords("inputtobed", "");
    nc.close();

    // if -..._period is not set, make n_records the minimum of the
    // buffer size and the number of available records. Otherwise try
    // to keep all available records in memory.
    if (itbperiod_set == false) {
      n_records = std::min(n_records, buffer_size);
    }

    if (n_records == 0) {
      throw RuntimeError::formatted("can't find 'inputtobed' in -hydrology_input_to_bed file with name '%s'",
                                    itbfilename.c_str());
    }

    verbPrintf(2,m_grid.com,
               "    option -hydrology_input_to_bed_file seen ... creating 'inputtobed' variable ...\n");
    verbPrintf(2,m_grid.com,
               "    allocating buffer space for n = %d 'inputtobed' records ...\n", n_records);
    inputtobed = new IceModelVec2T;
    inputtobed->set_n_records(n_records);
    inputtobed->create(m_grid, "inputtobed", WITHOUT_GHOSTS);
    inputtobed->set_attrs("climate_forcing",
                          "amount of water (depth per time like bmelt) which should be put at the ice sheet bed",
                          "m s-1", "");
    verbPrintf(2,m_grid.com,
               "    reading 'inputtobed' variable from file '%s' ...\n",itbfilename.c_str());
    inputtobed->init(itbfilename, inputtobed_period, inputtobed_reference_time);
  }

  // initialize till water layer thickness from the context if present,
  //   otherwise from -i or -boot_file, otherwise with constant value
  IceModelVec2S *Wtil_input = NULL;

  try {
    // FIXME: this is not an exceptional situation...
    Wtil_input = m_grid.variables().get_2d_scalar("tillwat");
    Wtil.copy_from(*Wtil_input);
  } catch (RuntimeError) {
    if (i_set || bootstrap) {
      std::string filename;
      int start;
      find_pism_input(filename, bootstrap, start);
      if (i_set) {
        Wtil.read(filename, start);
      } else {
        Wtil.regrid(filename, OPTIONAL,
                    m_config.get("bootstrapping_tillwat_value_no_var"));
      }
    } else {
      Wtil.set(m_config.get("bootstrapping_tillwat_value_no_var"));
    }
  }

  // whether or not we could initialize from file, we could be asked to regrid from file
  regrid("Hydrology", &Wtil);
}


void Hydrology::get_diagnostics(std::map<std::string, Diagnostic*> &dict,
                                    std::map<std::string, TSDiagnostic*> &/*ts_dict*/) {
  dict["bwat"] = new Hydrology_bwat(this, m_grid);
  dict["bwp"] = new Hydrology_bwp(this, m_grid);
  dict["bwprel"] = new Hydrology_bwprel(this, m_grid);
  dict["effbwp"] = new Hydrology_effbwp(this, m_grid);
  dict["hydrobmelt"] = new Hydrology_hydrobmelt(this, m_grid);
  dict["hydroinput"] = new Hydrology_hydroinput(this, m_grid);
  dict["wallmelt"] = new Hydrology_wallmelt(this, m_grid);
}


void Hydrology::add_vars_to_output(const std::string &/*keyword*/, std::set<std::string> &result) {
  result.insert("tillwat");
}


void Hydrology::define_variables(const std::set<std::string> &vars, const PIO &nc,
                                               IO_Type nctype) {
  if (set_contains(vars, "tillwat")) {
    Wtil.define(nc, nctype);
  }
}


void Hydrology::write_variables(const std::set<std::string> &vars, const PIO &nc) {
  if (set_contains(vars, "tillwat")) {
    Wtil.write(nc);
  }
}


//! Update the overburden pressure from ice thickness.
/*!
Uses the standard hydrostatic (shallow) approximation of overburden pressure,
  \f[ P_0 = \rho_i g H \f]
Accesses H=thk from Vars, which points into IceModel.
 */
void Hydrology::overburden_pressure(IceModelVec2S &result) {
  // FIXME issue #15
  result.copy_from(*thk);  // copies into ghosts if result has them
  result.scale(m_config.get("ice_density") * m_config.get("standard_gravity"));
}


//! Return the effective thickness of the water stored in till.
void Hydrology::till_water_thickness(IceModelVec2S &result) {
  Wtil.copy_to(result);
}


//! Set the wall melt rate to zero.  (The most basic subglacial hydrologies have no lateral flux or potential gradient.)
void Hydrology::wall_melt(IceModelVec2S &result) {
  result.set(0.0);
}


/*!
Checks \f$0 \le W_{til} \le W_{til}^{max} =\f$hydrology_tillwat_max.
 */
void Hydrology::check_Wtil_bounds() {
  double tillwat_max = m_config.get("hydrology_tillwat_max");

  IceModelVec::AccessList list(Wtil);
  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (Wtil(i,j) < 0.0) {
      throw RuntimeError::formatted("Hydrology: negative till water effective layer thickness Wtil(i,j) = %.6f m\n"
                                    "at (i,j)=(%d,%d)", Wtil(i,j), i, j);
    }

    if (Wtil(i,j) > tillwat_max) {
      throw RuntimeError::formatted("Hydrology: till water effective layer thickness Wtil(i,j) = %.6f m exceeds\n"
                                    "hydrology_tillwat_max = %.6f at (i,j)=(%d,%d)",
                                    Wtil(i,j), tillwat_max, i, j);
    }
  }
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
void Hydrology::get_input_rate(double hydro_t, double hydro_dt,
                                         IceModelVec2S &result) {
  bool   use_const   = m_config.get_flag("hydrology_use_const_bmelt");
  double const_bmelt = m_config.get("hydrology_const_bmelt");

  IceModelVec::AccessList list;
  if (inputtobed != NULL) {
    inputtobed->update(hydro_t, hydro_dt);
    inputtobed->interp(hydro_t + hydro_dt/2.0);
    list.add(*inputtobed);
  }

  if (!hold_bmelt) {
    bmelt->copy_to(bmelt_local);
  }

  list.add(bmelt_local);
  list.add(*mask);
  list.add(result);
  MaskQuery m(*mask);
  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m.icy(i, j)) {
      result(i,j) = (use_const) ? const_bmelt : bmelt_local(i,j);
      if (inputtobed != NULL) {
        result(i,j) += (*inputtobed)(i,j);
      }
    } else {
      result(i,j) = 0.0;
    }
  }
}


} // end of namespace pism
