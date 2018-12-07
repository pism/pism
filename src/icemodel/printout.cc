// Copyright (C) 2004-2018 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <cstring>
#include <cstdlib>

#include <petscsys.h>

#include "IceModel.hh"

#include "pism/stressbalance/StressBalance.hh"

#include "pism/util/IceGrid.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Time.hh"

#include "pism/util/pism_utilities.hh"

namespace pism {


//! Because of the -skip mechanism it is still possible that we can have CFL violations: count them.
/*! This applies to the horizontal part of the 3D advection problem solved by AgeModel and the
horizontal part of the 3D convection-diffusion problems solved by EnthalpyModel and
TemperatureModel.
*/
unsigned int count_CFL_violations(const IceModelVec3 &u3,
                                  const IceModelVec3 &v3,
                                  const IceModelVec2S &ice_thickness,
                                  double dt) {

  if (dt == 0.0) {
    return 0;
  }

  IceGrid::ConstPtr grid = u3.grid();

  const double
    CFL_x = grid->dx() / dt,
    CFL_y = grid->dy() / dt;

  IceModelVec::AccessList list{&ice_thickness, &u3, &v3};

  unsigned int CFL_violation_count = 0;
  ParallelSection loop(grid->com);
  try {
    for (Points p(*grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const int ks = grid->kBelowHeight(ice_thickness(i,j));

      const double
        *u = u3.get_column(i, j),
        *v = v3.get_column(i, j);

      // check horizontal CFL conditions at each point
      for (int k = 0; k <= ks; k++) {
        if (fabs(u[k]) > CFL_x) {
          CFL_violation_count += 1;
        }
        if (fabs(v[k]) > CFL_y) {
          CFL_violation_count += 1;
        }
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return (unsigned int)GlobalMax(grid->com, CFL_violation_count);
}

void IceModel::print_summary(bool tempAndAge) {

  const IceModelVec3
    &u3 = m_stress_balance->velocity_u(),
    &v3 = m_stress_balance->velocity_v();

  unsigned int n_CFL_violations = count_CFL_violations(u3, v3, m_geometry.ice_thickness,
                                                       tempAndAge ? dt_TempAge : m_dt);

  // report CFL violations
  if (n_CFL_violations > 0.0) {
    const double CFLviolpercent = 100.0 * n_CFL_violations / (m_grid->Mx() * m_grid->My() * m_grid->Mz());
    // at default verbosity level, only report CFL viols if above:
    const double CFLVIOL_REPORT_VERB2_PERCENT = 0.1;
    if (CFLviolpercent > CFLVIOL_REPORT_VERB2_PERCENT ||
        m_log->get_threshold() > 2) {
      char tempstr[90] = "";
      snprintf(tempstr,90,
              "  [!CFL#=%d (=%5.2f%% of 3D grid)] ",
              n_CFL_violations,CFLviolpercent);
      m_stdout_flags = tempstr + m_stdout_flags;
    }
  }

  // get maximum diffusivity
  double max_diffusivity = m_stress_balance->max_diffusivity();
  // get volumes in m^3 and areas in m^2
  double volume = ice_volume(0.0);
  double area = ice_area(0.0);

  double meltfrac = 0.0;
  if (tempAndAge or m_log->get_threshold() >= 3) {
    meltfrac = compute_temperate_base_fraction(area);
  }

  // main report: 'S' line
  print_summary_line(false, tempAndAge, m_dt,
                   volume, area, meltfrac, max_diffusivity);
}


//! Print a line to stdout which summarizes the state of the modeled ice sheet at the end of the time step.
/*!
This method is for casual inspection of model behavior, and to provide the user
with some indication of the state of the run.

Generally, two lines are printed to stdout, the first starting with a space
and the second starting with the character 'S' in the left-most column (column 1).

The first line shows flags for which processes executed, and the length of the
time step (and/or substeps under option -skip).  See IceModel::run()
for meaning of these flags.

If printPrototype is TRUE then the first line does not appear and
the second line has alternate appearance.  Specifically, different column 1
characters are printed:
  - 'P' line gives names of the quantities reported in the 'S' line, the
    "prototype", while
  - 'U' line gives units of these quantities.
This column 1 convention allows automatic tools to read PISM stdout
and produce time-series.  The 'P' and 'U' lines are intended to appear once at
the beginning of the run, while an 'S' line appears at every time step.

These quantities are reported in this base class version:
  - `time` is the current model time
  - `ivol` is the total ice sheet volume
  - `iarea` is the total area occupied by positive thickness ice
  - `max_diffusivity` is the maximum diffusivity
  - `max_hor_vel` is the maximum diffusivity

Configuration parameters `output.runtime.time_unit_name`, `output.runtime.volume_scale_factor_log10`,
and `output.runtime.area_scale_factor_log10` control the appearance and units.

For more description and examples, see the PISM User's Manual.
Derived classes of IceModel may redefine this method and print alternate
information.
 */
void IceModel::print_summary_line(bool printPrototype,  bool tempAndAge,
                                double delta_t,
                                double volume,  double area,
                                double /* meltfrac */,  double max_diffusivity) {
  const bool do_energy = m_config->get_boolean("energy.enabled");
  const int log10scalevol  = static_cast<int>(m_config->get_double("output.runtime.volume_scale_factor_log10")),
            log10scalearea = static_cast<int>(m_config->get_double("output.runtime.area_scale_factor_log10"));
  const std::string time_units = m_config->get_string("output.runtime.time_unit_name");
  const bool use_calendar = m_config->get_boolean("output.runtime.time_use_calendar");

  const double scalevol  = pow(10.0, static_cast<double>(log10scalevol)),
               scalearea = pow(10.0, static_cast<double>(log10scalearea));
  char  volscalestr[10] = "     ", areascalestr[10] = "   "; // blank when 10^0 = 1 scaling
  if (log10scalevol != 0) {
    snprintf(volscalestr, sizeof(volscalestr), "10^%1d_", log10scalevol);
  }
  if (log10scalearea != 0) {
    snprintf(areascalestr, sizeof(areascalestr), "10^%1d_", log10scalearea);
  }

  if (printPrototype == true) {
    m_log->message(2,
               "P         time:       ivol      iarea  max_diffusivity  max_hor_vel\n");
    m_log->message(2,
               "U         %s   %skm^3  %skm^2         m^2 s^-1       m/%s\n",
               time_units.c_str(),volscalestr,areascalestr,time_units.c_str());
    return;
  }

  // this version keeps track of what has been done so as to minimize stdout:
  // FIXME: turn these static variables into class members.
  static std::string stdout_flags_count0;
  static int         mass_cont_sub_counter = 0;
  static double      mass_cont_sub_dtsum   = 0.0;
  if (mass_cont_sub_counter == 0) {
    stdout_flags_count0 = m_stdout_flags;
  }
  if (delta_t > 0.0) {
    mass_cont_sub_counter++;
    mass_cont_sub_dtsum += delta_t;
  }

  if ((tempAndAge == true) || (!do_energy) || (m_log->get_threshold() > 2)) {
    char tempstr[90]    = "";

    const double major_dt = m_time->convert_time_interval(mass_cont_sub_dtsum, time_units);
    if (mass_cont_sub_counter <= 1) {
      snprintf(tempstr,90, " (dt=%.5f)", major_dt);
    } else {
      snprintf(tempstr,90, " (dt=%.5f in %d substeps; av dt_sub_mass_cont=%.5f)",
               major_dt, mass_cont_sub_counter, major_dt / mass_cont_sub_counter);
    }
    stdout_flags_count0 += tempstr;

    if (delta_t > 0.0) { // avoids printing an empty line if we have not done anything
      stdout_flags_count0 += "\n";
      m_log->message(2, stdout_flags_count0);
    }

    if (use_calendar) {
      snprintf(tempstr,90, "%12s", m_time->date().c_str());
    } else {
      snprintf(tempstr,90, "%.3f", m_time->convert_time_interval(m_time->current(), time_units));
    }



    const CFLData cfl = m_stress_balance->max_timestep_cfl_2d();
    std::string velocity_units = "meters / (" + time_units + ")";
    const double maxvel = units::convert(m_sys, std::max(cfl.u_max, cfl.v_max),
                                         "m second-1", velocity_units);

    m_log->message(2,
               "S %s:   %8.5f  %9.5f     %12.5f %12.5f\n",
               tempstr,
               volume/(scalevol*1.0e9), area/(scalearea*1.0e6),
               max_diffusivity, maxvel);

    mass_cont_sub_counter = 0;
    mass_cont_sub_dtsum = 0.0;
  }
}


} // end of namespace pism
