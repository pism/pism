/* Copyright (C) 2018, 2019, 2020, 2021 PISM Authors
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

#include "ScalarForcing.hh"

#include <algorithm>            // std::min_element(), std::max_element()
#include <cassert>              // assert()
#include <cmath>                // std::floor()

#include "pism/util/ConfigInterface.hh"
#include "pism/util/Context.hh"
#include "pism/util/Time.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/Logger.hh"
#include "pism/util/io/File.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/VariableMetadata.hh"

namespace pism {

//! \brief Report the range of a time-series stored in `data`.
static void report_range(const std::vector<double> &data,
                         const units::System::Ptr unit_system,
                         const VariableMetadata &metadata,
                         const Logger &log) {
  double min, max;

  // min_element and max_element return iterators; "*" is used to get
  // the value corresponding to this iterator
  min = *std::min_element(data.begin(), data.end());
  max = *std::max_element(data.begin(), data.end());

  units::Converter c(unit_system,
                     metadata.get_string("units"),
                     metadata.get_string("glaciological_units"));
  min = c(min);
  max = c(max);

  std::string spacer(metadata.get_name().size(), ' ');

  log.message(2,
              "  FOUND  %s / %-60s\n"
              "         %s \\ min,max = %9.3f,%9.3f (%s)\n",
              metadata.get_name().c_str(),
              metadata.get_string("long_name").c_str(), spacer.c_str(), min, max,
              metadata.get_string("glaciological_units").c_str());
}

void ScalarForcing::initialize(const Context &ctx,
                               const std::string &filename,
                               const std::string &variable_name,
                               const std::string &units,
                               const std::string &glaciological_units,
                               const std::string &long_name,
                               bool periodic) {
  try {
    auto unit_system = ctx.unit_system();

    VariableMetadata variable(variable_name, unit_system);

    variable.set_string("units", units);
    variable.set_string("glaciological_units", glaciological_units);
    variable.set_string("long_name", long_name);

    double forcing_t0{};
    double forcing_t1{};

    // Read data from a NetCDF file
    {
      ctx.log()->message(2,
                         "  reading %s (%s) from file '%s'...\n",
                         long_name.c_str(), variable_name.c_str(), filename.c_str());

      File file(ctx.com(), filename, PISM_NETCDF3, PISM_READONLY);

      // Read forcing data. The read_timeseries() call will ensure that variable_name is a
      // scalar variable.
      std::vector<double> data{};
      io::read_timeseries(file, variable, *ctx.log(), data);

      // The following line relies on the fact that this is a scalar variable.
      std::string time_name = file.dimensions(variable_name)[0];

      std::vector<double> times{};
      std::vector<double> bounds{};
      io::read_time_info(*ctx.log(), unit_system,
                         file, time_name, ctx.time()->units_string(),
                         times, bounds);
      size_t N = times.size();

      // Compute values used to extend data read from file
      //
      // Initialize using values corresponding to constant extrapolation:
      double v0 = data.front();
      double v1 = data.back();
      if (periodic) {
        double b0 = bounds.front();
        double b1 = bounds.back();

        // compute the value at the beginning and end of the period:
        {
          double dt1 = times.front() - b0;
          double dt2 = b1 - times.back();

          // interpolation factor
          double alpha = dt1 + dt2 > 0 ? dt2 / (dt1 + dt2) : 0.0;

          v0 = (1.0 - alpha) * data.back() + alpha * data.front();
          v1 = v0;
        }

        {
          m_period       = b1 - b0;
          m_period_start = b0;

          assert(m_period > 0.0);
        }
      }

      // Note: this should take care of file with one record as well.
      m_times.clear();
      m_values.clear();
      {
        if (bounds.front() < times.front()) {
          m_times.emplace_back(bounds.front());
          m_values.emplace_back(v0);
        }
        for (size_t k = 0; k < N; ++k) {
          m_times.emplace_back(times[k]);
          m_values.emplace_back(data[k]);
        }
        if (bounds.back()  > times.back()) {
          m_times.emplace_back(bounds.back());
          m_values.emplace_back(v1);
        }
      }

      // Save the time interval covered by this forcing.
      forcing_t0 = bounds.front();
      forcing_t1 = bounds.back();
    } // end of the block initializing data

    // validate resulting times
    if (not is_increasing(m_times)) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "m_times have to be strictly increasing (this is a bug)");
    }

    report_range(m_values, unit_system, variable, *ctx.log());

    // Set up interpolation.
    {
      m_acc = gsl_interp_accel_alloc();
      m_spline = gsl_spline_alloc(gsl_interp_linear, m_times.size());
      gsl_spline_init(m_spline, m_times.data(), m_values.data(), m_times.size());
    }

    bool extrapolate = ctx.config()->get_flag("input.forcing.time_extrapolation");

    if (not (extrapolate or periodic)) {
      check_forcing_duration(*ctx.time(), forcing_t0, forcing_t1);
    }

  } catch (RuntimeError &e) {
    e.add_context("reading %s (%s) from '%s'",
                  long_name.c_str(), variable_name.c_str(), filename.c_str());
    throw;
  }
}

ScalarForcing::ScalarForcing(const Context &ctx,
                             const std::string &prefix,
                             const std::string &variable_name,
                             const std::string &units,
                             const std::string &glaciological_units,
                             const std::string &long_name)
  : m_period(0.0),
    m_period_start(0.0),
    m_acc(nullptr),
    m_spline(nullptr) {

  Config::ConstPtr config = ctx.config();

  auto filename = config->get_string(prefix + ".file");
  bool periodic = config->get_flag(prefix + ".periodic");

  if (filename.empty()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "%s.file is required", prefix.c_str());
  }

  initialize(ctx, filename, variable_name, units, glaciological_units,
             long_name, periodic);
}

ScalarForcing::ScalarForcing(const Context &ctx,
                             const std::string &filename,
                             const std::string &variable_name,
                             const std::string &units,
                             const std::string &glaciological_units,
                             const std::string &long_name,
                             bool periodic)
  : m_period(0.0),
    m_period_start(0.0),
    m_acc(nullptr),
    m_spline(nullptr) {

  initialize(ctx, filename, variable_name, units, glaciological_units,
             long_name, periodic);
}

ScalarForcing::~ScalarForcing() {
  if (m_spline != nullptr) {
    gsl_spline_free(m_spline);
  }
  if (m_acc != nullptr) {
    gsl_interp_accel_free(m_acc);
  }
}

double ScalarForcing::value(double t) const {
  double T = t;

  if (m_period > 0.0) {
    // compute time since the period start
    T -= m_period_start;

    // remove an integer number of periods
    T -= std::floor(T / m_period) * m_period;

    // add the period start back
    T += m_period_start;
  }

  if (T > m_times.back()) {
    return m_values.back();
  }

  if (T < m_times.front()) {
    return m_values.front();
  }

  return gsl_spline_eval(m_spline, T, m_acc);
}

// Integrate from a to b, interpreting data as *not* periodic
double ScalarForcing::integral(double a, double b) const {
  assert(b >= a);

  double dt = b - a;

  double t0 = m_times.front();
  double t1 = m_times.back();
  double v0 = m_values[0];
  double v1 = m_values.back();

  // both points are to the left of [t0, t1]
  if (b <= t0) {
    return v0 * (b - a);
  }

  // both points are to the right of [t0, t1]
  if (a >= t1) {
    return v1 * (b - a);
  }

  // a is to the left of [t0, t1]
  if (a < t0) {
    return v0 * (t0 - a) + integral(t0, b);
  }

  // b is to the right of [t0, t1]
  if (b > t1) {
    return integral(a, t1) + v1 * (b - t1);
  }

  // both points are inside [t0, t1]
  size_t ai = gsl_interp_bsearch(m_times.data(), a, 0, m_times.size() - 1);
  size_t bi = gsl_interp_bsearch(m_times.data(), b, 0, m_times.size() - 1);

  // gsl_interp_bsearch() above returns the index i of the array ‘m_times’ such that
  // ‘m_times[i] <= t < m_times[i+1]’.  The index is searched for in the
  // range [0, m_times.size() - 1] (inclusive).

  double v_a = gsl_spline_eval(m_spline, a, m_acc);
  double v_b = gsl_spline_eval(m_spline, b, m_acc);

  if (ai == bi) {
    return 0.5 * (v_a + v_b) * dt;
  }

  double result = 0.0;
  // integrate from a to the data point just to its right
  result += 0.5 * (v_a + m_values[ai + 1]) * (m_times[ai + 1] - a);

  // integrate over (possibly zero) intervals between a and b
  for (size_t k = ai + 1; k < bi; ++k) {
    result += 0.5 * (m_values[k] + m_values[k + 1]) * (m_times[k + 1] - m_times[k]);
  }

  result += 0.5 * (m_values[bi] + v_b) * (b - m_times[bi]);

  return result;
}

double ScalarForcing::average(double t, double dt) const {
  if (dt == 0.0) {
    return value(t);
  }

  if (dt < 0.0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "negative interval length");
  }

  if (m_period > 0.0) {
    double a = t;
    double b = t + dt;
    double t0 = m_period_start;
    double P = m_period;

    double N = std::floor((a - t0) / P);
    double M = std::floor((b - t0) / P);
    double delta = a - t0 - P * N;
    double gamma = b - t0 - P * M;

    double N_periods = M - (N + 1);

    if (N_periods >= 0.0) {
      return (integral(t0 + delta, t0 + P) +
              N_periods * integral(t0, t0 + P) +
              integral(t0, t0 + gamma)) / dt;
    }

    return integral(t0 + delta, t0 + gamma) / dt;

  } else {
    return integral(t, t + dt) / dt;
  }
}


} // end of namespace pism
