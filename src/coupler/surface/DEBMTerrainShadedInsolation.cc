/* Copyright (C) 2026 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */


/*!
 * Reference:
 *
 * A. Chokshi, “Solshade: Terrain-aware Solar Illumination Modelling using Digital
 * Elevation Models and Orbital Geometry,” Journal of Open Source Software, vol. 11, no.
 * 122, p. 9944, Jun. 2026, doi: 10.21105/joss.09944.
 */

#include <cassert>
#include <cmath>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscsystypes.h>
#include <petscvec.h>
#include <vector>

#include "pism/coupler/surface/Insolation.hh"

#include "pism/coupler/surface/terrain_insolation_kernel.hh"

#include "pism/util/Component.hh"
#include "pism/util/Config.hh"
#include "pism/util/Context.hh"
#include "pism/util/Grid.hh"
#include "pism/util/Logger.hh"
#include "pism/util/Profiling.hh"
#include "pism/util/SunPosition.hh"
#include "pism/util/array/Array3D.hh"
#include "pism/util/array/Scalar.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/petscwrappers/DM.hh"
#include "pism/util/petscwrappers/Vec.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/projection.hh"
#include "pism/util/Time.hh"
#include "pism/util/io/IO_Flags.hh"
#include "pism/util/MaxTimestep.hh"

namespace pism {
namespace surface {

void DEBMTerrainShadedInsolation::init_impl(const InputOptions &options) {
  m_t_last = time().current();
  if (options.type == INIT_RESTART or options.type == INIT_BOOTSTRAP) {
    File input_file(m_grid->com, options.filename, io::PISM_NETCDF3, io::PISM_READONLY);

    if (input_file.variable_exists(m_time_name)) {
      auto t_length = input_file.nrecords(m_time_name, "", m_sys);
      auto t_last   = t_length > 0 ? t_length - 1 : 0;

      auto t =
          io::read_timeseries_variable(input_file, m_time_name, time().units(), m_sys, t_last, 1);
      m_t_last = t[0];
    }
  }

  {
    auto log    = m_grid->ctx()->log();
    auto config = m_grid->ctx()->config();

    switch (options.type) {
    case INIT_RESTART:
      // read bed elevation and uplift rate from file
      log->message(2, "    reading bed topography and uplift from %s ... \n",
                   options.filename.c_str());
      // re-starting
      m_horizon->read(options.filename, options.record); // fails if not found!
      break;
    case INIT_BOOTSTRAP:
    case INIT_OTHER:
    default: {
      // do nothing
    } break;
    }

    // process -regrid_file and -regrid_vars
    io::regrid(*config, "insolation", *m_horizon, true);
  }
}

void DEBMTerrainShadedInsolation::update_impl(double t, double dt,
                                              const array::Scalar1 &surface_elevation) {

  update_shading(surface_elevation);

  update_daily_insolation(t + 0.5 * dt, surface_elevation);
}

void DEBMTerrainShadedInsolation::begin_access() const {
  m_insolation.begin_access();
}

void DEBMTerrainShadedInsolation::end_access() const {
  m_insolation.end_access();
}

void DEBMTerrainShadedInsolation::surface_insolation_energy_series_impl(
    int i, int j, const std::vector<OrbitalParameters> &orbital, double dt,
    std::vector<double> &result) const {
  // the daily-mean insolation rate (W m^-2) at this cell was computed in update_shading()
  double rate = m_insolation(i, j);
  for (size_t k = 0; k < orbital.size(); ++k) {
    result[k] = rate * dt;
  }
}

MaxTimestep DEBMTerrainShadedInsolation::max_timestep_impl(double t,
                                                           const CFLData * /*cfl_data*/) const {

  if (t < m_t_last) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "time %f is less than the previous time %f",
                                  t, m_t_last);
  }

  if (m_update_interval == 0.0) {
    return {};
  }

  // Find the smallest time of the form m_t_last + k * m_update_interval that is greater
  // than t
  double k = std::ceil((t - m_t_last) / m_update_interval);

  double t_next = m_t_last + k * m_update_interval;

  double dt_max = m_update_interval;
  if (t < t_next) {
    dt_max = t_next - t;
  }

  if (dt_max < m_t_eps) {
    dt_max = m_update_interval;
  }

  return MaxTimestep(dt_max, "bed_def");
}

std::set<VariableMetadata> DEBMTerrainShadedInsolation::state_impl() const {
  auto variables = array::metadata({ &horizon_altitude() });

  VariableMetadata T(m_time_name, m_sys);
  T.long_name("time of the last update of the terrain-shaded insolation model")
      .set_time_dependent(true)
      .units(time().units())
      .set_string("calendar", time().calendar());

  variables.insert(T);

  return variables;
}

void DEBMTerrainShadedInsolation::write_state_impl(const OutputFile &output) const {
  horizon_altitude().write(output);

  auto t_length = output.time_dimension_length();
  auto t_start = t_length > 0 ? t_length - 1 : 0;

  output.write_timeseries(m_time_name, { t_start }, { 1 }, { m_t_last });
}

DiagnosticList DEBMTerrainShadedInsolation::spatial_diagnostics_impl() const {
  DiagnosticList result = { { "insolation", Diagnostic::wrap(insolation()) },
                            { "horizon", Diagnostic::wrap(horizon_altitude()) } };

  if (sky_view_enabled()) {
    result["sky_view_factor"] = Diagnostic::wrap(sky_view_factor());
  }

  return result;
}

DEBMTerrainShadedInsolation::DEBMTerrainShadedInsolation(
    std::shared_ptr<const Grid> grid, std::function<double(double)> atmosphere_transmissivity)
    : Insolation(grid),
      m_time_name(time().variable_name() + "_insolation"),
      m_t_last(time().current()),
      m_insolation(grid, "insolation"),
      m_y_azimuth(m_grid, "y_azimuth"),
      m_orbital_parameters(*grid->ctx()),
      m_transmissivity(atmosphere_transmissivity) {

  m_update_interval = m_config->get_number("bed_deformation.update_interval", "seconds");
  m_t_eps           = m_config->get_number("time_stepping.resolution", "seconds");

  m_insolation.metadata(0)
      .long_name("daily mean terrain-shaded surface insolation")
      .units("W m^-2");

  auto config = m_grid->ctx()->config();

  m_n_directions =
      static_cast<int>(config->get_number("terrain_shaded_insolation.horizon.n_directions"));
  m_max_distance     = config->get_number("terrain_shaded_insolation.horizon.max_distance");
  m_step             = config->get_number("terrain_shaded_insolation.horizon.step");
  m_insolation_dt    = config->get_number("terrain_shaded_insolation.dt");
  m_solar_constant   = config->get_number("surface.debm_simple.solar_constant");
  m_diffuse_fraction = config->get_number("terrain_shaded_insolation.diffuse_fraction");

  if (not(m_step > 0.0)) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "terrain_shaded_insolation.horizon.step must be positive");
  }

  if (m_n_directions < 1) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "terrain_shaded_insolation.horizon.n_directions must be positive");
  }

  // azimuth sample directions (radians), clockwise from north
  std::vector<double> azimuth(m_n_directions);
  for (int k = 0; k < m_n_directions; ++k) {
    azimuth[k] = 2.0 * M_PI * k / m_n_directions;
  }

  m_horizon = std::make_shared<array::Array3D>(m_grid, "horizon", array::WITHOUT_GHOSTS, azimuth);
  m_horizon->metadata(0)
      .long_name("terrain horizon elevation angle as a function of azimuth")
      .units("radian");

  if (config->get_flag("terrain_shaded_insolation.use_sky_view_factor")) {
    m_sky_view = std::make_shared<array::Scalar>(m_grid, "sky_view_factor");
    m_sky_view->metadata(0)
        .long_name("sky-view factor (fraction of the diffuse sky hemisphere visible "
                   "from the terrain-shaded, tilted surface)")
        .units("1");
  }

  // Allocate the scatter to all ranks and the vector that will hold local copies of the
  // DEM:
  PetscErrorCode ierr;
  ierr = DMDAGlobalToNaturalAllCreate(*m_grid->get_dm(1, 0), m_scatter.rawptr());
  PISM_CHK(ierr, "DMDAGlobalToNaturalAllCreate");

  ierr = VecCreateSeq(PETSC_COMM_SELF, static_cast<PetscInt>(m_grid->Mx() * m_grid->My()),
                      m_dem_local.rawptr());
  PISM_CHK(ierr, "VecCreateSeq");

  // use projection information to compute the azimuth of the Y axis of the grid:
  compute_y_azimuth(m_y_azimuth);
}

const array::Scalar& DEBMTerrainShadedInsolation::insolation() const {
  return m_insolation;
}

const array::Array3D &DEBMTerrainShadedInsolation::horizon_altitude() const {
  return *m_horizon;
}

const array::Scalar &DEBMTerrainShadedInsolation::sky_view_factor() const {
  if (not sky_view_enabled()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "sky view factor is not available");
  }
  return *m_sky_view;
}

bool DEBMTerrainShadedInsolation::sky_view_enabled() const {
  return m_sky_view != nullptr;
}

/*!
 * Given the azimuth of the Y axis, convert the gradient of a function `f` in the X-Y
 * coordinate system (vector (f_x, f_y)) to its gradient in the east-north-up system (f_e,
 * f_n).
 */
static void xy_to_en(double y_azimuth, double f_x, double f_y, double &f_e, double &f_n) {
  // unit vector pointing north:
  double north_x = std::sin(-y_azimuth);
  double north_y = std::cos(-y_azimuth);
  // unit vector pointing east (90 degrees from north, simplified using trigonometric
  // identities):
  double east_x = north_y;
  double east_y = -north_x;

  // f_e and f_n are directional derivatives of f computed as dot products between vectors
  // pointing east and north and the gradient of f in the x-y coordinate system
  //
  // See https://en.wikipedia.org/wiki/Directional_derivative
  f_e = east_x * f_x + east_y * f_y;
  f_n = north_x * f_x + north_y * f_y;
}

/*!
 * Finite difference approximations of partial derivatives of a gridded field.
 *
 * Uses centered differences in the interior and one-sided differences at domain
 * boundaries.
 */
class FD {
public:
  FD(int Mx_, int My_, double dx_, double dy_) {
    Mx = Mx_;
    My = My_;
    dx = dx_;
    dy = dy_;
  }

  inline double x(const array::Scalar1 &F, int i, int j) const {
    int ip = i < Mx - 1 ? i + 1 : i;
    int im = i > 0 ? i - 1 : i;

    return (F(ip, j) - F(im, j)) / ((ip - im) * dx);
  }

  inline double y(const array::Scalar1 &F, int i, int j) const {
    int jp = j < My - 1 ? j + 1 : j;
    int jm = j > 0 ? j - 1 : j;

    return (F(i, jp) - F(i, jm)) / ((jp - jm) * dy);
  }

private:
  int Mx, My;
  double dx, dy;
};


/*!
 * Update the sky view factor and the map of horizon altitude angles (in radians) at each
 * grid point.
 */
void DEBMTerrainShadedInsolation::update_shading(const array::Scalar1 &surface_elevation) {
  auto log = m_grid->ctx()->log();

  log->message(2, "* Updating the horizon map...\n");
  double start = get_time(m_grid->com);

  const auto &profiling = m_grid->ctx()->profiling();

  const int Mx = static_cast<int>(m_grid->Mx());
  const int My = static_cast<int>(m_grid->My());
  const double dx = m_grid->dx();
  const double dy = m_grid->dy();

  FD diff(Mx, My, dx, dy);

  // Scatter the full DEM to every rank. After this block each rank has a copy of the
  // global surface elevation that lets it compute shading at its owned cells without
  // ghost communication.
  {
    profiling.begin("insolation.scatter_dem");
    PetscErrorCode ierr;
    auto dm = surface_elevation.dm();
    petsc::TemporaryGlobalVec dem_global(dm);

    // make sure surface_elevation is actually "local":
    assert(surface_elevation.stencil_width() > 0);

    // Note: we use DMLocalToGlobal because surface_elevation is ghosted (local)
    ierr = DMLocalToGlobal(*dm, surface_elevation.vec(), INSERT_VALUES, dem_global);
    PISM_CHK(ierr, "DMLocalToGlobal");

    ierr = VecScatterBegin(m_scatter, dem_global, m_dem_local, INSERT_VALUES, SCATTER_FORWARD);
    PISM_CHK(ierr, "VecScatterBegin");
    ierr = VecScatterEnd(m_scatter, dem_global, m_dem_local, INSERT_VALUES, SCATTER_FORWARD);
    PISM_CHK(ierr, "VecScatterEnd");
    profiling.end("insolation.scatter_dem");
  }

  profiling.begin("insolation.horizon_altitude");

  const auto &azimuth = m_horizon->levels();

  petsc::VecArray dem(m_dem_local);
  array::AccessScope scope{ &surface_elevation, m_horizon.get(), &m_y_azimuth };

  if (sky_view_enabled()) {
    scope.add(*m_sky_view);
  }

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();

    double y_azimuth = m_y_azimuth(i, j);

    double *horizon_altitude = m_horizon->get_column(i, j);

    for (int k = 0; k < m_n_directions; ++k) {
      horizon_altitude[k] = terrain::ray_horizon(dem.get(), Mx, My, dx, dy, i, j,
                                                 azimuth[k] - y_azimuth, m_step, m_max_distance);
    }

    if (sky_view_enabled()) {
      // Surface elevation derivatives in the east and north directions:
      double s_e = 0.0, s_n = 0.0;
      xy_to_en(y_azimuth, diff.x(surface_elevation, i, j), diff.y(surface_elevation, i, j), s_e,
               s_n);

      // Compute the upward-pointing normal to the surface in the east-north-up system:
      double Ne = -s_e, Nn = -s_n, Nu = 1.0;
      // Scale it to get the unit normal:
      {
        // Note that norm != 0.0 because Nu == 1
        double norm = std::sqrt(Ne * Ne + Nn * Nn + Nu * Nu);

        Ne /= norm;
        Nn /= norm;
        Nu /= norm;
      }

      double slope  = std::acos(pism::clip(Nu, -1.0, 1.0));
      double aspect = std::atan2(Ne, Nn);
      (*m_sky_view)(i, j) =
          terrain::sky_view_factor(horizon_altitude, azimuth.data(), m_n_directions, slope, aspect);
    }
  } // end of the loop over grid points
  profiling.end("insolation.horizon_altitude");

  double end = get_time(m_grid->com);
  log->message(2, "* Updated the horizon map in %f s.\n", end - start);
}

//! Periodic linear interpolation of a horizon column at the given azimuth (radians).
double DEBMTerrainShadedInsolation::horizon_altitude(int i, int j, double azimuth) const {
  const double *column = m_horizon->get_column(i, j);
  const int n = m_n_directions;
  const double two_pi = 2.0 * M_PI;
  const double da = two_pi / n;

  double a = azimuth - two_pi * std::floor(azimuth / two_pi); // wrap to [0, 2*pi)
  double x = a / da;
  int k = static_cast<int>(std::floor(x));
  double frac = x - k;

  int k0 = k % n;
  int k1 = (k + 1) % n;

  return column[k0] * (1.0 - frac) + column[k1] * frac;
}

// The per-timestep shadow test, Lambertian cosine projection (max(0, normal . sun)), and
// inverse-square distance scaling are adapted from solshade's compute_flux_timeseries
// (solshade/irradiance.py).
//
// Here they are integrated over the diurnal cycle to compute daily energy.
void DEBMTerrainShadedInsolation::update_daily_insolation(double time,
                                                          const array::Scalar1 &surface_elevation) {

  auto parameters = m_orbital_parameters.compute(time);
  double declination = parameters.solar_declination;
  double distance_factor = parameters.distance_factor;

  const double seconds_per_day = 86400.0;

  SunPositionCalculator sun_position(declination);

  // number of sub-daily samples used to integrate the diurnal cycle
  int M = static_cast<int>(std::lround(seconds_per_day / m_insolation_dt));
  if (M < 1) {
    M = 1;
  }

  // Pre-compute cos() and sin() of M hour angles:
  {
    std::vector<double> hour_angles(M);
    for (int m = 0; m < M; ++m) {
      // hour angle sweeps the full day, midpoint rule over [-pi, pi)
      hour_angles[m] = -M_PI + 2.0 * M_PI * (m + 0.5) / M;
    }
    sun_position.set_hour_angles(hour_angles);
  }

  const double dt = seconds_per_day / M;

  const auto &profiling = m_grid->ctx()->profiling();
  profiling.begin("insolation.daily");

  const auto &latitude = m_grid->latitude();

  array::AccessScope scope{ &latitude, &m_insolation, &surface_elevation, m_horizon.get(),
                            &m_y_azimuth };

  bool use_sky_view = sky_view_enabled();
  if (use_sky_view) {
    scope.add(*m_sky_view);
  }

  FD diff((int)m_grid->Mx(), (int)m_grid->My(), m_grid->dx(), m_grid->dy());

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();

    sun_position.set_latitude(latitude(i, j) * (M_PI / 180.0));

    // Surface elevation derivatives in the east and north directions:
    double s_e = 0.0, s_n = 0.0;
    xy_to_en(m_y_azimuth(i, j), diff.x(surface_elevation, i, j), diff.y(surface_elevation, i, j),
             s_e, s_n);

    // Compute the upward-pointing normal to the surface in the east-north-up system:
    double Ne = -s_e, Nn = -s_n, Nu = 1.0;
    // Scale it to get the unit normal:
    {
      // Note that norm != 0.0 because Nu == 1
      double norm = std::sqrt(Ne * Ne + Nn * Nn + Nu * Nu);

      Ne /= norm;
      Nn /= norm;
      Nu /= norm;
    }

    const double f_diff = use_sky_view ? m_diffuse_fraction : 0.0;
    const double svf = use_sky_view ? (*m_sky_view)(i, j) : 0.0;

    double toa_insolation_energy = 0.0;
    // loop over hour angles:
    for (int hour_angle_index = 0; hour_angle_index < M; ++hour_angle_index) {

      double altitude = 0.0, azimuth = 0.0;

      // vector pointing toward the center of the sun in the east-north-up coordinate
      // system:
      double solar_vector[3] = {0.0, 0.0, 0.0};

      sun_position.compute_at_set_hour_angle(hour_angle_index, altitude, azimuth, solar_vector);

      if (altitude <= 0.0) {
        // Sun below the astronomical horizon: no direct and no diffuse contribution
        continue;
      }

      // Note: from equation (13) in Sproul2007 sin(altitude) is equal to the Z ("up")
      // component of the solar vector.
      double sin_altitude = solar_vector[2];
      // top-of-atmosphere horizontal irradiance, the basis for the diffuse component
      double toa_horizontal = m_solar_constant * distance_factor * sin_altitude;

      // diffuse: isotropic sky scaled by the sky-view factor; reaches shadowed cells too
      toa_insolation_energy += f_diff * toa_horizontal * svf * dt;

      // direct beam: only when the Sun clears the local horizon and lights the surface
      if (altitude > horizon_altitude(i, j, azimuth)) {
        double Se = solar_vector[0];
        double Sn = solar_vector[1];
        double Su = solar_vector[2];

        double mu = Ne * Se + Nn * Sn + Nu * Su;
        if (mu > 0.0) {
          toa_insolation_energy += (1.0 - f_diff) * m_solar_constant * distance_factor * mu * dt;
        }
      }
    } // end of the loop over hour angles

    // store the daily-mean insolation rate (W m-2)
    m_insolation(i, j) =
        m_transmissivity(surface_elevation(i, j)) * toa_insolation_energy / seconds_per_day;
  } // end of the loop over grid points

  profiling.end("insolation.daily");
}

}
}
