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

ScalarForcing::ScalarForcing(const Context &ctx,
                             const std::string &prefix,
                             const std::string &variable_name,
                             const std::string &units,
                             const std::string &glaciological_units,
                             const std::string &long_name)
  : m_period(0.0), m_period_start(0.0), m_acc(nullptr), m_spline(nullptr) {

  Config::ConstPtr config = ctx.config();

  auto unit_system = ctx.unit_system();

  VariableMetadata variable(variable_name, unit_system);

  variable.set_string("units", units);
  variable.set_string("glaciological_units", glaciological_units);
  variable.set_string("long_name", long_name);

  // Read data from a NetCDF file
  {
    auto filename = config->get_string(prefix + ".file");

    if (filename.empty()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "%s.file is required", prefix.c_str());
    }

    ctx.log()->message(2,
                       "  reading %s data from forcing file %s...\n",
                       variable_name.c_str(), filename.c_str());

    try {
      bool periodic = config->get_flag(prefix + ".periodic");
      auto time_units = ctx.time()->units_string();

      File file(ctx.com(), filename, PISM_NETCDF3, PISM_READONLY);

      io::read_timeseries(file, variable, *ctx.log(), m_values);

      std::string time_name{};
      {
        auto dims = file.dimensions(variable_name);
        // io::read_timeseries() already ensured that dims.size() == 1
        time_name = dims[0];
      }

      if (periodic) {
        // require time bounds and override times using time bounds

        std::string time_bounds_name = file.read_text_attribute(time_name, "bounds");

        if (time_bounds_name.empty()) {
          throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                        "please provide time bounds for '%s' in '%s'",
                                        variable_name.c_str(), file.filename().c_str());
        }

        VariableMetadata time_bounds(time_bounds_name, unit_system);
        time_bounds.set_string("units", time_units);

        std::vector<double> bounds{};
        io::read_time_bounds(file, time_bounds, *ctx.log(), bounds);

        auto N = m_values.size();

        if (2 * N != bounds.size()) {
          throw RuntimeError(PISM_ERROR_LOCATION, "expected two bounds per value");
        }

        // Add points to data stored in RAM: at the very beginning and the very end of the
        // period. This will simplify interpolation.
        std::vector<double> t(N + 2);
        std::vector<double> v(N + 2);

        // compute the value at the beginning and end of the period:
        double v0 = 0.0;
        {
          // interval lengths at the beginning and the end of the period
          double dt_f = bounds[2 * 0 + 1] - bounds[0 + 0];
          double dt_l = bounds[2 * (N - 1) + 1] - bounds[2 * (N - 1) + 0];

          // interpolation factor
          double alpha = dt_l + dt_f > 0 ? dt_l / (dt_l + dt_f) : 0.0;

          v0 = (1.0 - alpha) * m_values.back() + alpha * m_values.front();
        }

        t.front() = bounds.front();
        v.front() = v0;

        for (size_t k = 0; k < N; ++k) {
          t[k + 1] = 0.5 * (bounds[2 * k + 0] + bounds[2 * k + 1]);
          v[k + 1] = m_values[k];
        }

        t.back() = bounds.back();
        v.back() = v0;          // note: v.front() == v.back()


        m_times  = t;
        m_values = v;
        {
          m_period       = m_times.back() - m_times.front();
          m_period_start = m_times.front();

          assert(m_period > 0.0);
        }
      } else {
        // ignore time bounds and use times

        VariableMetadata time_dimension(time_name, unit_system);
        time_dimension.set_string("units", time_units);

        io::read_timeseries(file, time_dimension, *ctx.log(), m_times);

        if (not is_increasing(m_times)) {
          throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                        "dimension '%s' has to be strictly increasing (read from '%s').",
                                        time_name.c_str(), file.filename().c_str());

        }

        {
          m_period = 0.0;
          m_period_start = 0.0;
        }
      }
      // LCOV_EXCL_START
    } catch (RuntimeError &e) {
      e.add_context("while reading %s (%s) from '%s'",
                    long_name.c_str(), variable_name.c_str(), filename.c_str());
      throw;
    }
    // LCOV_EXCL_STOP
  } // end of the block reading data

  report_range(m_values, unit_system, variable, *ctx.log());

  // set up interpolation
  if (m_times.size() > 1) {
    m_acc = gsl_interp_accel_alloc();
    m_spline = gsl_spline_alloc(gsl_interp_linear, m_times.size());
    gsl_spline_init(m_spline, m_times.data(), m_values.data(), m_times.size());
  }
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

  if (m_times.size() == 1) {
    return m_values[0];
  }

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
  // range [0, m_times.size() - 1]

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
