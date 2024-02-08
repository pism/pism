// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2021, 2023, 2024 PISM Authors
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

#include "ClimateIndex.hh"
#include "pism/util/Time.hh"
#include "pism/util/Grid.hh"

namespace pism {
namespace atmosphere {

ClimateIndex::ClimateIndex(std::shared_ptr<const Grid> g)
  : AtmosphereModel(g),
    // Reference fields for the mean annual and mean summer near-surface air temperature
    m_air_temp_annual(m_grid, "air_temp_annual_ref"),
    m_air_temp_annual_ref(m_grid, "air_temp_annual_ref"),
    m_air_temp_summer(m_grid, "air_temp_summer_ref"),
    m_air_temp_summer_ref(m_grid, "air_temp_summer_ref"),

    // Anaomaly temperature fields for Climate index 0 (e.g. LGM), interglacial index 1 (e.g. LIG) and interglacial index 1X (e.g. mPWP)
    m_air_temp_anomaly_annual_0(m_grid, "air_temp_anomaly_annual_0"),
    m_air_temp_anomaly_annual_1(m_grid, "air_temp_anomaly_annual_1"),
    m_air_temp_anomaly_annual_1X(m_grid, "air_temp_anomaly_annual_1X"),

    m_air_temp_anomaly_summer_0(m_grid, "air_temp_anomaly_summer_0"),
    m_air_temp_anomaly_summer_1(m_grid, "air_temp_anomaly_summer_1"),
    m_air_temp_anomaly_summer_1X(m_grid, "air_temp_anomaly_summer_1X"),

    //Reference precipitation field
    m_precipitation(m_grid, "precipitation_ref"),
    m_precipitation_ref(m_grid, "precipitation_ref"),

    m_precipitation_anomaly_0(m_grid, "precipitation_anomaly_0"),
    m_precipitation_anomaly_1(m_grid, "precipitation_anomaly_1"),
    m_precipitation_anomaly_1X(m_grid, "precipitation_anomaly_1X"),

    // Spatial precipitation scaling factor
    m_spatial_precip_scaling(m_grid, "precip_scaling_factor") {

    auto climate_index_file = m_config->get_string("climate_index.file");
    if (not climate_index_file.empty()) {
        m_climate_index.reset(
            new ClimateIndexWeights(*g->ctx()));
    } else {
          m_log->message(2,
                  "* 'Climate Index Weights' Index File not given. Please specify under -climate_index_file \n");
    }

    auto scaling_file = m_config->get_string("atmosphere.yearly_cycle.scaling.file");
    if (not scaling_file.empty()) {
        m_A.reset(new ScalarForcing(*g->ctx(),
                                "atmosphere.yearly_cycle.scaling",
                                "amplitude_scaling",
                                "1", "1",
                                "temperature amplitude scaling"));
    }

  m_snow_temp_summer_day = m_config->get_number("atmosphere.fausto_air_temp.summer_peak_day");

  m_air_temp_annual.metadata(0)
        .long_name("mean annual near-surface air temperature (without sub-year time-dependence or forcing)")
        .units("K")
        .set_time_independent(true);
  m_air_temp_annual.metadata()["source"] = m_reference;

  m_air_temp_annual_ref.metadata(0)
        .long_name("mean annual near-surface air temperature (without sub-year time-dependence or forcing)")
        .units("K")
        .set_time_independent(true);
  m_air_temp_annual_ref.metadata()["source"] = m_reference;

  m_air_temp_summer.metadata(0)
        .long_name("mean summer (NH: July/ SH: January) near-surface air temperature (without sub-year time-dependence or forcing)")
        .units("K")
        .set_time_independent(true);
  m_air_temp_summer.metadata()["source"] = m_reference;

  m_air_temp_summer_ref.metadata(0)
        .long_name("mean summer (NH: July/ SH: January) near-surface air temperature (without sub-year time-dependence or forcing)")
        .units("K")
        .set_time_independent(true);
  m_air_temp_summer_ref.metadata()["source"] = m_reference;

  m_precipitation.metadata(0)
      .long_name("precipitation rate")
      .units("kg m-2 second-1")
      .output_units("kg m-2 year-1")
      .set_time_independent(true);

  m_precipitation_ref.metadata(0)
      .long_name("precipitation rate")
      .units("kg m-2 second-1")
      .output_units("kg m-2 year-1")
      .set_time_independent(true);

  m_air_temp_anomaly_annual_0.metadata(0)
        .long_name("mean annual near-surface air temperature (without sub-year time-dependence or forcing) for Climate index 0 (e.g. LGM)")
        .units("K")
        .set_time_independent(true);
  m_air_temp_anomaly_annual_0.metadata()["source"] = m_reference;

  m_air_temp_anomaly_annual_1.metadata(0)
        .long_name("mean annual near-surface air temperature (without sub-year time-dependence or forcing) for interglacial index 1 (e.g. LIG)")
        .units("K")
        .set_time_independent(true);
  m_air_temp_anomaly_annual_1.metadata()["source"] = m_reference;

  m_air_temp_anomaly_annual_1X.metadata(0)
        .long_name("mean PD annual near-surface air temperature (without sub-year time-dependence or forcing) for interglacial index 1X (e.g. mPWP)")
        .units("K")
        .set_time_independent(true);
  m_air_temp_anomaly_annual_1X.metadata()["source"] = m_reference;

  // Paleo time slice temperature summer
  m_air_temp_anomaly_summer_0.metadata(0)
        .long_name("mean summer (NH: July/ SH: January) near-surface air temperature (without sub-year time-dependence or forcing) for Climate index 0 (e.g. LGM)")
        .units("K")
        .set_time_independent(true);
  m_air_temp_anomaly_summer_0.metadata()["source"] = m_reference;

  m_air_temp_anomaly_summer_1.metadata(0)
        .long_name("mean summer (NH: July/ SH: January) near-surface air temperature (without sub-year time-dependence or forcing) for interglacial index 1 (e.g. LIG)")
        .units("K")
        .set_time_independent(true);
  m_air_temp_anomaly_summer_1.metadata()["source"] = m_reference;

  m_air_temp_anomaly_summer_1X.metadata(0)
        .long_name("mean summer (NH: July/ SH: January) near-surface air temperature (without sub-year time-dependence or forcing) for interglacial index 1X (e.g. mPWP)")
        .units("K")
        .set_time_independent(true);
  m_air_temp_anomaly_summer_1X.metadata()["source"] = m_reference;

  m_precipitation_anomaly_0.metadata(0)
      .long_name("precipitation rate")
      .units("kg m-2 second-1")
      .output_units("kg m-2 year-1")
      .set_time_independent(true);
  m_precipitation_anomaly_0.metadata()["source"] = m_reference;

  m_precipitation_anomaly_1.metadata(0)
      .long_name("precipitation rate anomaly")
      .units("kg m-2 second-1")
      .output_units("kg m-2 year-1")
      .set_time_independent(true);
  m_precipitation_anomaly_1.metadata()["source"] = m_reference;

  m_precipitation_anomaly_1X.metadata(0)
      .long_name("precipitation rate")
      .units("kg m-2 second-1")
      .output_units("kg m-2 year-1")
      .set_time_independent(true);
  m_precipitation_anomaly_1X.metadata()["source"] = m_reference;

  use_precip_scaling = m_config->get_flag("atmosphere.climate_index.precip_scaling.use");

  // Spatial precip scaling factor
  m_spatial_precip_scaling.metadata(0)
        .long_name("spatial scaling factor with temperature for precipitation")
        .units("K-1")
        .set_time_independent(true);
  m_spatial_precip_scaling.metadata()["source"] = m_reference;
}

void ClimateIndex::init_impl(const Geometry &geometry) {
  (void) geometry;
  m_log->message(2,
            "**** Initializing the 'Climate Index' atmosphere model...\n");

  auto input_file = m_config->get_string("atmosphere.climate_index.climate_snapshots.file");

  if (input_file.empty()) {
      throw RuntimeError(PISM_ERROR_LOCATION,
                      "Please specify an '-atmosphere.climate_index.climate_snapshots' input file\n"
                      "using -atmosphere_climate_snapshots_file or a command-line option.\n");
  }

  m_log->message(2,
                "  Reading mean annual air temperature, mean July air temperature, and\n"
                "  precipitation fields from '%s'...\n", input_file.c_str());

  // Reference fields
  m_precipitation.regrid(input_file, io::Default::Nil());
  m_precipitation_ref.regrid(input_file, io::Default::Nil());
  m_air_temp_annual.regrid(input_file, io::Default::Nil());
  m_air_temp_annual_ref.regrid(input_file, io::Default::Nil());

  try {
    m_air_temp_summer.regrid(input_file, io::Default::Nil());
    m_air_temp_summer_ref.regrid(input_file, io::Default::Nil());
    use_cos = true;
    m_log->message(2,
                " * mean summer temperatures are provided, thus will use a cosinus function\n"
                " for representing yearly cycles \n");
  } catch (...) {
    use_cos = false;
    m_log->message(2,
                " * no mean summer temperatures provided: uses annual mean temperatures only (no yearly cycle) \n");
  }

  // Annual anomaly for Paleo time slices 0=Glacial, 1=glacial, 1X= Super Interglacial e.g. mPWP
  m_air_temp_anomaly_annual_0.regrid(input_file, io::Default::Nil());
  m_air_temp_anomaly_annual_1.regrid(input_file, io::Default::Nil());
  // Summer anomaly
  if (use_cos) {
    m_air_temp_anomaly_summer_0.regrid(input_file, io::Default::Nil());
    m_air_temp_anomaly_summer_1.regrid(input_file, io::Default::Nil());
  }

  try {
    m_air_temp_anomaly_annual_1X.regrid(input_file, io::Default::Nil());
    if (use_cos) { m_air_temp_anomaly_summer_1X.regrid(input_file, io::Default::Nil()); }
    use_1X = true;
    m_log->message(2,
        "* 1X slice found in input file. Will use it for scaling the atmosphere forcing\n");
  } catch (...) {
    use_1X = false;
  }

  precip_scaling_file = m_config->get_string("atmosphere.climate_index.precip_scaling.spatial_linear_factor.file");
  // If a file is giving for spatial scaling, then it will use it. Otherwise by default it uses the linear, uniform scaling factor from config

  if (not use_precip_scaling) {
    m_log->message(2,
        "*  no scaling method used for precipitation\n"
        "   thus it will use the precipitation snapshots from atmosphere.climate_index.climate_snapshots.file if exist\n");
    m_precipitation_anomaly_0.regrid(input_file, io::Default::Nil());
    m_precipitation_anomaly_1.regrid(input_file, io::Default::Nil());
    try {
      m_precipitation_anomaly_1X.regrid(input_file, io::Default::Nil());
    } catch (...) {}
  } else if (not precip_scaling_file.empty()) {
    m_spatial_precip_scaling.regrid(precip_scaling_file, io::Default::Nil());
    m_log->message(2,
                "*  - scaling file given for precipitation scaling in -atmosphere.climate_index.precip_scaling.spatial_linear_factor\n"
                "    thus Climate Index forcing is using temperature anomalies to calculate\n"
                "    precipitation anomalies using a spatially distributed scaling factor from '%s'...\n", precip_scaling_file.c_str());
  } else {
    m_preciplinfactor = m_config->get_number("atmosphere.climate_index.precip_scaling.uniform_linear_factor");
    m_log->message(2,
                  "*   -climate_index_precip_scaling is set to the default uniform scaling factor,\n"
                  "    thus precipitation anomalies are calculated using linear scaling\n"
                  "    with air temperature anomalies (%.3f percent per degree).\n",
                  m_preciplinfactor * 100.0);
  }
}

void ClimateIndex::update_impl(const Geometry &geometry, double t, double dt) {
  (void) geometry;
  (void) t;
  (void) dt;

  double w0{0.0}, w1{0.0}, w1X{0.0};
  m_climate_index->update_weights(t, dt, w0, w1, w1X);

  m_log->message(3,
             "**** Updated weigths in atmo: m_w0 = '%f', m_w1 = '%f', m_w1X = '%f' ****\n", w0, w1, w1X);

  array::AccessScope scope{ &m_air_temp_annual, &m_air_temp_annual_ref, &m_air_temp_anomaly_annual_0,
    &m_air_temp_anomaly_annual_1, &m_precipitation, &m_precipitation_ref };

  if (use_cos) {
    scope.add({ &m_air_temp_summer, &m_air_temp_summer_ref, &m_air_temp_anomaly_summer_0,
                &m_air_temp_anomaly_summer_1 });
  }

  if (use_1X) {
    scope.add(m_air_temp_anomaly_annual_1X);
    if (use_cos) {
      scope.add(m_air_temp_anomaly_summer_1X);
    }
  }

  if (not use_precip_scaling) {
    scope.add({&m_precipitation_anomaly_0, &m_precipitation_anomaly_1});
    if (use_1X) {
      scope.add(m_precipitation_anomaly_1X);
    }
  } else if (not precip_scaling_file.empty()) {
    scope.add(m_spatial_precip_scaling);
  }

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double annual_anomaly = 0.0;
    double summer_anomaly = 0.0;

    if (use_1X) {
      annual_anomaly = w0 * m_air_temp_anomaly_annual_0(i, j) + w1 * m_air_temp_anomaly_annual_1(i, j) + w1X * (m_air_temp_anomaly_annual_1X(i, j) - m_air_temp_anomaly_annual_1(i, j));
      if (use_cos) {
        summer_anomaly = w0 * m_air_temp_anomaly_summer_0(i, j) + w1 * m_air_temp_anomaly_summer_1(i, j) + w1X * (m_air_temp_anomaly_summer_1X(i, j) - m_air_temp_anomaly_summer_1(i, j));
      }
    } else {
      annual_anomaly = w0 * m_air_temp_anomaly_annual_0(i, j) + w1 * m_air_temp_anomaly_annual_1(i, j);
      if (use_cos) {
        summer_anomaly = w0 * m_air_temp_anomaly_summer_0(i, j) + w1 * m_air_temp_anomaly_summer_1(i, j);
      }
    }
    m_air_temp_annual(i, j) = m_air_temp_annual_ref(i, j)  + annual_anomaly;
    if (use_cos) {
      m_air_temp_summer(i, j) = m_air_temp_summer_ref(i, j)  + summer_anomaly;
    }
    if (not use_precip_scaling) {
      if (use_1X) {
        m_precipitation(i, j) = m_precipitation_ref(i, j) + w0 * m_precipitation_anomaly_0(i, j) + w1 * m_precipitation_anomaly_1(i, j) + w1X * (m_precipitation_anomaly_1X(i, j) - m_precipitation_anomaly_1(i, j));
      } else {
        m_precipitation(i, j) = m_precipitation_ref(i, j) + w0 * m_precipitation_anomaly_0(i, j) + w1 * m_precipitation_anomaly_1(i, j);
      }
    } else if (not precip_scaling_file.empty()) {
      m_precipitation(i, j) = m_precipitation_ref(i, j) * (1 + annual_anomaly*m_spatial_precip_scaling(i, j));
    } else {
      m_precipitation(i, j) = m_precipitation_ref(i, j) * (1 + annual_anomaly*m_preciplinfactor);
    }
  } // end of the loop over grid points
}

void ClimateIndex::define_model_state_impl(const File &output) const {
  m_precipitation.define(output, io::PISM_DOUBLE);
}

void ClimateIndex::write_model_state_impl(const File &output) const {
  m_precipitation.write(output);
}

//! Copies the stored precipitation field into result.
const array::Scalar& ClimateIndex::precipitation_impl() const {
  return m_precipitation;
}

//! Copies the stored mean annual near-surface air temperature field into result.
const array::Scalar& ClimateIndex::air_temperature_impl() const {
  return m_air_temp_annual;
}

//! Copies the stored mean summer near-surface air temperature field into result.
const array::Scalar& ClimateIndex::mean_summer_temp() const {
  return m_air_temp_summer;
}

// same as YearlyCycle
void ClimateIndex::init_timeseries_impl(const std::vector<double> &ts) const {

  const double
    summerday_fraction = time().day_of_the_year_to_year_fraction(m_snow_temp_summer_day);

  size_t N = ts.size();

  m_ts_times.resize(N);
  m_cosine_cycle.resize(N);
  for (unsigned int k = 0; k < m_ts_times.size(); k++) {
    double tk = time().year_fraction(ts[k]) - summerday_fraction;

    m_ts_times[k] = ts[k];
    m_cosine_cycle[k] = cos(2.0 * M_PI * tk);
  }

  if (m_A) {
    for (unsigned int k = 0; k < ts.size(); ++k) {
      m_cosine_cycle[k] *= m_A->value(ts[k]);
    }
  }
}

MaxTimestep ClimateIndex::max_timestep_impl(double t) const {
  (void) t;
  return MaxTimestep("atmosphere cosine_yearly_cycle");
}

void ClimateIndex::precip_time_series_impl(int i, int j, std::vector<double> &result) const {
  result.resize(m_ts_times.size());
  for (unsigned int k = 0; k < m_ts_times.size(); k++) {
    result[k] = m_precipitation(i,j);
  }
}

void ClimateIndex::temp_time_series_impl(int i, int j, std::vector<double> &result) const {
  result.resize(m_ts_times.size());
  for (unsigned int k = 0; k < m_ts_times.size(); ++k) {
    if (use_cos) {
      result[k] = m_air_temp_annual(i, j) + (m_air_temp_summer(i, j) - m_air_temp_annual(i, j)) * m_cosine_cycle[k];
    } else {
      result[k] = m_air_temp_annual(i, j);
    }
  }
}

void ClimateIndex::begin_pointwise_access_impl() const {
  m_air_temp_annual.begin_access();
  if (use_cos) { m_air_temp_summer.begin_access(); }
  m_precipitation.begin_access();
}

void ClimateIndex::end_pointwise_access_impl() const {
  m_air_temp_annual.end_access();
  if (use_cos) { m_air_temp_summer.end_access(); }
  m_precipitation.end_access();
}

namespace diagnostics {

/*! @brief Mean summer near-surface air temperature. */
class MeanSummerTemperature : public Diag<ClimateIndex>
{
public:
  MeanSummerTemperature(const ClimateIndex *m) : Diag<ClimateIndex>(m) {
    m_vars = { { m_sys, "air_temp_summer" } };
    m_vars[0]
        .long_name("mean summer near-surface air temperature used in the cosine yearly cycle")
        .units("Kelvin");
  }
private:
  std::shared_ptr<array::Array> compute_impl() const {
    auto result = allocate<array::Scalar>("air_temp_summer");

    result->copy_from(model->mean_summer_temp());

    return result;
  }
};
} // end of namespace diagnostics

DiagnosticList ClimateIndex::diagnostics_impl() const {
  DiagnosticList result = AtmosphereModel::diagnostics_impl();

  result["air_temp_summer"] = Diagnostic::Ptr(new diagnostics::MeanSummerTemperature(this));

  return result;
}

} // end of namespace atmosphere
} // end of namespace pism
