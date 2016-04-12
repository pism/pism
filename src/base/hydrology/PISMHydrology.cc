// Copyright (C) 2012-2016 PISM Authors
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

#include "PISMHydrology.hh"
#include "base/util/Mask.hh"
#include "base/util/PISMVars.hh"
#include "base/util/error_handling.hh"
#include "base/util/iceModelVec2T.hh"
#include "base/util/io/PIO.hh"
#include "base/util/pism_options.hh"
#include "hydrology_diagnostics.hh"
#include "base/util/pism_utilities.hh"
#include "base/util/IceModelVec2CellType.hh"

namespace pism {
namespace hydrology {

Hydrology::Hydrology(IceGrid::ConstPtr g)
  : Component_TS(g) {
  m_inputtobed = NULL;
  m_hold_bmelt = false;

  m_total_input.create(m_grid, "total_input", WITHOUT_GHOSTS);
  m_total_input.set_attrs("internal",
                        "hydrology model workspace for total input rate into subglacial water layer",
                        "m s-1", "");

  m_bmelt_local.create(m_grid, "bmelt", WITHOUT_GHOSTS);
  m_bmelt_local.set_attrs("internal",
                        "hydrology model workspace for bmelt",
                        "m s-1", "");

  // *all* Hydrology classes have layer of water stored in till as a state variable
  m_Wtil.create(m_grid, "tillwat", WITHOUT_GHOSTS);
  m_Wtil.set_attrs("model_state",
                 "effective thickness of subglacial water stored in till",
                 "m", "");
  m_Wtil.metadata().set_double("valid_min", 0.0);
}


Hydrology::~Hydrology() {
  // empty
}


void Hydrology::init() {

  m_log->message(4,
             "entering Hydrology::init() ...\n");

  options::String bmelt_file("-hydrology_bmelt_file",
                             "Read time-independent values for bmelt from a file;"
                             " replaces bmelt computed through conservation of energy");
  // itb = input_to_bed
  options::String itb_file("-hydrology_input_to_bed_file",
                           "A time- and space-dependent file with amount of water"
                           " (depth per time) which should be added to the amount of water"
                           " at the ice sheet bed at the given location at the given time;"
                           " adds to bmelt");

  options::Real itb_period_years("-hydrology_input_to_bed_period",
                                 "The period (i.e. duration before repeat), in years,"
                                 " of -hydrology_input_to_bed_file data", 0.0);

  options::Real itb_reference_year("-hydrology_input_to_bed_reference_year",
                                   "The reference year for periodizing the"
                                   " -hydrology_input_to_bed_file data", 0.0);

  // the following are IceModelVec pointers into IceModel generally and are read by code in the
  // update() method at the current Hydrology time

  if (bmelt_file.is_set()) {
    m_log->message(2,
               "  option -hydrology_bmelt_file seen; reading bmelt from '%s'.\n", bmelt_file->c_str());
    m_bmelt_local.regrid(bmelt_file, CRITICAL);
    m_hold_bmelt = true;
  }


  if (itb_file.is_set()) {
    m_inputtobed_period = itb_period_years;
    m_inputtobed_reference_time = units::convert(m_sys, itb_reference_year, "years", "seconds");

    unsigned int buffer_size = (unsigned int) m_config->get_double("climate_forcing_buffer_size");

    PIO nc(m_grid->com, "netcdf3");
    nc.open(itb_file, PISM_READONLY);
    unsigned int n_records = nc.inq_nrecords("inputtobed", "", m_sys);
    nc.close();

    // if -..._period is not set, make n_records the minimum of the
    // buffer size and the number of available records. Otherwise try
    // to keep all available records in memory.
    if (not itb_period_years.is_set()) {
      n_records = std::min(n_records, buffer_size);
    }

    if (n_records == 0) {
      throw RuntimeError::formatted("can't find 'inputtobed' in -hydrology_input_to_bed"
                                    " file with name '%s'",
                                    itb_file->c_str());
    }

    m_log->message(2,
               "    option -hydrology_input_to_bed_file seen ... creating 'inputtobed' variable ...\n");
    m_log->message(2,
               "    allocating buffer space for n = %d 'inputtobed' records ...\n", n_records);
    m_inputtobed = new IceModelVec2T;
    m_inputtobed->set_n_records(n_records);
    m_inputtobed->create(m_grid, "inputtobed");
    m_inputtobed->set_attrs("climate_forcing",
                            "amount of water (depth per time like bmelt)"
                            " which should be put at the ice sheet bed",
                            "m s-1", "");
    m_log->message(2,
               "    reading 'inputtobed' variable from file '%s' ...\n",
               itb_file->c_str());
    m_inputtobed->init(itb_file, m_inputtobed_period, m_inputtobed_reference_time);
  }

  bool bootstrap = false;
  int start = 0;
  std::string filename;
  bool use_input_file = find_pism_input(filename, bootstrap, start);

  if (use_input_file) {
    if (bootstrap) {
      m_Wtil.regrid(filename, OPTIONAL,
                    m_config->get_double("bootstrapping_tillwat_value_no_var"));
    } else {
      m_Wtil.read(filename, start);
    }
  } else {
    m_Wtil.set(m_config->get_double("bootstrapping_tillwat_value_no_var"));
  }

  // whether or not we could initialize from file, we could be asked to regrid from file
  regrid("Hydrology", m_Wtil);
}


void Hydrology::get_diagnostics_impl(std::map<std::string, Diagnostic::Ptr> &dict,
                                     std::map<std::string, TSDiagnostic::Ptr> &/*ts_dict*/) {
  dict["bwat"]       = Diagnostic::Ptr(new Hydrology_bwat(this));
  dict["bwp"]        = Diagnostic::Ptr(new Hydrology_bwp(this));
  dict["bwprel"]     = Diagnostic::Ptr(new Hydrology_bwprel(this));
  dict["effbwp"]     = Diagnostic::Ptr(new Hydrology_effbwp(this));
  dict["hydrobmelt"] = Diagnostic::Ptr(new Hydrology_hydrobmelt(this));
  dict["hydroinput"] = Diagnostic::Ptr(new Hydrology_hydroinput(this));
  dict["wallmelt"]   = Diagnostic::Ptr(new Hydrology_wallmelt(this));
}


void Hydrology::add_vars_to_output_impl(const std::string &/*keyword*/, std::set<std::string> &result) {
  result.insert("tillwat");
}


void Hydrology::define_variables_impl(const std::set<std::string> &vars, const PIO &nc,
                                      IO_Type nctype) {
  if (set_contains(vars, "tillwat")) {
    m_Wtil.define(nc, nctype);
  }
}


void Hydrology::write_variables_impl(const std::set<std::string> &vars, const PIO &nc) {
  if (set_contains(vars, "tillwat")) {
    m_Wtil.write(nc);
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
  const IceModelVec2S *thk = m_grid->variables().get_2d_scalar("thk");

  result.copy_from(*thk);  // copies into ghosts if result has them
  result.scale(m_config->get_double("ice_density") * m_config->get_double("standard_gravity"));
}


//! Return the effective thickness of the water stored in till.
void Hydrology::till_water_thickness(IceModelVec2S &result) {
  result.copy_from(m_Wtil);
}


//! Set the wall melt rate to zero.  (The most basic subglacial hydrologies have no lateral flux or potential gradient.)
void Hydrology::wall_melt(IceModelVec2S &result) {
  result.set(0.0);
}


/*!
Checks \f$0 \le W_{til} \le W_{til}^{max} =\f$hydrology_tillwat_max.
 */
void Hydrology::check_Wtil_bounds() {
  double tillwat_max = m_config->get_double("hydrology_tillwat_max");

  IceModelVec::AccessList list(m_Wtil);
  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (m_Wtil(i,j) < 0.0) {
        throw RuntimeError::formatted("Hydrology: negative till water effective layer thickness Wtil(i,j) = %.6f m\n"
                                      "at (i,j)=(%d,%d)", m_Wtil(i,j), i, j);
      }

      if (m_Wtil(i,j) > tillwat_max) {
        throw RuntimeError::formatted("Hydrology: till water effective layer thickness Wtil(i,j) = %.6f m exceeds\n"
                                      "hydrology_tillwat_max = %.6f at (i,j)=(%d,%d)",
                                      m_Wtil(i,j), tillwat_max, i, j);
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

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
  bool   use_const   = m_config->get_boolean("hydrology_use_const_bmelt");
  double const_bmelt = m_config->get_double("hydrology_const_bmelt");

  IceModelVec::AccessList list;
  if (m_inputtobed != NULL) {
    m_inputtobed->update(hydro_t, hydro_dt);
    m_inputtobed->interp(hydro_t + hydro_dt/2.0);
    list.add(*m_inputtobed);
  }

  const IceModelVec2S        &bmelt = *m_grid->variables().get_2d_scalar("bmelt");
  const IceModelVec2CellType &mask  = *m_grid->variables().get_2d_cell_type("mask");

  if (not m_hold_bmelt) {
    m_bmelt_local.copy_from(bmelt);
  }

  list.add(m_bmelt_local);
  list.add(mask);
  list.add(result);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.icy(i, j)) {
      result(i,j) = (use_const) ? const_bmelt : m_bmelt_local(i,j);
      if (m_inputtobed != NULL) {
        result(i,j) += (*m_inputtobed)(i,j);
      }
    } else {
      result(i,j) = 0.0;
    }
  }
}


} // end of namespace hydrology
} // end of namespace pism
