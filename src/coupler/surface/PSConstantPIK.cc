// Copyright (C) 2008-2014 PISM Authors
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

#include "PSConstantPIK.hh"
#include "PIO.hh"
#include "PISMVars.hh"
#include "IceGrid.hh"
#include <stdexcept>

namespace pism {

///// Constant-in-time surface model for accumulation,
///// ice surface temperature parameterized as in PISM-PIK dependent on latitude and surface elevation


PSConstantPIK::PSConstantPIK(IceGrid &g, const Config &conf)
  : SurfaceModel(g, conf)
{
  PetscErrorCode ierr = allocate_PSConstantPIK(); CHKERRCONTINUE(ierr);
  if (ierr != 0) {
    throw std::runtime_error("PSConstantPIK allocation failed");
  }
}

void PSConstantPIK::attach_atmosphere_model(AtmosphereModel *input)
{
  delete input;
}

PetscErrorCode PSConstantPIK::allocate_PSConstantPIK() {
  PetscErrorCode ierr;

  climatic_mass_balance.create(grid, "climatic_mass_balance", WITHOUT_GHOSTS);
  climatic_mass_balance.set_attrs("climate_state",
                                  "constant-in-time surface mass balance (accumulation/ablation) rate",
                                  "kg m-2 s-1",
                                  "land_ice_surface_specific_mass_balance_flux");
  climatic_mass_balance.set_glaciological_units("kg m-2 year-1");
  climatic_mass_balance.write_in_glaciological_units = true;

  ice_surface_temp.create(grid, "ice_surface_temp", WITHOUT_GHOSTS);
  ice_surface_temp.set_attrs("climate_state",
                             "constant-in-time ice temperature at the ice surface",
                             "K", "");
  return 0;
}


void PSConstantPIK::init(Vars &vars) {
  bool do_regrid = false;
  int start = -1;

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  verbPrintf(2, grid.com,
             "* Initializing the constant-in-time surface processes model PSConstantPIK.\n"
             "  It reads surface mass balance directly from the file and holds it constant.\n"
             "  Ice upper-surface temperature is parameterized as in Martin et al. 2011, Eqn. 2.0.2.\n"
             "  Any choice of atmosphere coupler (option '-atmosphere') is ignored.\n");

  usurf = vars.get_2d_scalar("surface_altitude");
  lat   = vars.get_2d_scalar("latitude");

  // find PISM input file to read data from:
  find_pism_input(input_file, do_regrid, start);

  // read snow precipitation rate from file
  verbPrintf(2, grid.com,
             "    reading surface mass balance rate 'climatic_mass_balance' from %s ... \n",
             input_file.c_str());
  if (do_regrid) {
    climatic_mass_balance.regrid(input_file, CRITICAL); // fails if not found!
  } else {
    climatic_mass_balance.read(input_file, start); // fails if not found!
  }

  // parameterizing the ice surface temperature 'ice_surface_temp'
  verbPrintf(2, grid.com,
             "    parameterizing the ice surface temperature 'ice_surface_temp' ... \n");
}

void PSConstantPIK::update(double my_t, double my_dt)
{
  if ((fabs(my_t - m_t) < 1e-12) &&
      (fabs(my_dt - m_dt) < 1e-12)) {
    return;
  }

  m_t  = my_t;
  m_dt = my_dt;

  IceModelVec::AccessList list;
  list.add(ice_surface_temp);
  list.add(*usurf);
  list.add(*lat);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    ice_surface_temp(i,j) = 273.15 + 30 - 0.0075 * (*usurf)(i,j) - 0.68775 * (*lat)(i,j)*(-1.0);
  }
}

void PSConstantPIK::get_diagnostics(std::map<std::string, Diagnostic*> &/*dict*/,
                                    std::map<std::string, TSDiagnostic*> &/*ts_dict*/)
{
  // empty (does not have an atmosphere model)
}

void PSConstantPIK::ice_surface_mass_flux(IceModelVec2S &result) {
  climatic_mass_balance.copy_to(result);
}

void PSConstantPIK::ice_surface_temperature(IceModelVec2S &result) {
  ice_surface_temp.copy_to(result);
}

void PSConstantPIK::add_vars_to_output(const std::string &/*keyword*/, std::set<std::string> &result) {
  result.insert("climatic_mass_balance");
  result.insert("ice_surface_temp");
  // does not call atmosphere->add_vars_to_output().
}

void PSConstantPIK::define_variables(const std::set<std::string> &vars, const PIO &nc, IO_Type nctype) {
  SurfaceModel::define_variables(vars, nc, nctype);

  if (set_contains(vars, "ice_surface_temp")) {
    ice_surface_temp.define(nc, nctype);
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    climatic_mass_balance.define(nc, nctype);
  }
}

void PSConstantPIK::write_variables(const std::set<std::string> &vars, const PIO &nc) {
  if (set_contains(vars, "ice_surface_temp")) {
    ice_surface_temp.write(nc);
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    climatic_mass_balance.write(nc);
  }
}

} // end of namespace pism
