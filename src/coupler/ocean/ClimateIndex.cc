// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2021, 2023 PISM Authors
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
#include "pism/util/Grid.hh"

namespace pism {
namespace ocean {

ClimateIndex::ClimateIndex(std::shared_ptr<const Grid> grid)
: CompleteOceanModel(grid, std::shared_ptr<OceanModel>()),

    m_theta_ocean_ref(m_grid, "theta_ocean_ref"),
    m_salinity_ocean_ref(m_grid, "salinity_ocean_ref"),

    // Anaomaly temperature fields for glacial index 0 (e.g. LGM), interglacial index 1 (e.g. LIG) and interglacial index 1X (e.g. mPWP)
    m_theta_ocean_anomaly_0(m_grid, "theta_ocean_anomaly_0"),
    m_theta_ocean_anomaly_1(m_grid, "theta_ocean_anomaly_1"),
    m_theta_ocean_anomaly_1X(m_grid, "theta_ocean_anomaly_1X"),

    m_salinity_ocean_anomaly_0(m_grid, "salinity_ocean_anomaly_0"),
    m_salinity_ocean_anomaly_1(m_grid, "salinity_ocean_anomaly_1"),
    m_salinity_ocean_anomaly_1X(m_grid, "salinity_ocean_anomaly_1X") {

    auto climate_index_file = m_config->get_string("climate_index.file");
    if (not climate_index_file.empty()) {
        m_climate_index.reset(
            new ClimateIndexWeights(*grid->ctx()));

    } else {
        m_log->message(2,
                "* 'Climate Index Weight' Index File not given. Please specify under -climate_index_file \n");
    }

    m_theta_ocean_ref.metadata(0)
            .long_name("potential temperature of the adjacent ocean")
            .units("Kelvin")
            .set_time_independent(true);
    m_theta_ocean_ref.metadata()["source"] = m_reference;

    // Paleo time slice temperature data annual
    m_theta_ocean_anomaly_0.metadata(0)
        .long_name("absolute potential temperature anomaly of the adjacent ocean")
        .units("Kelvin")
        .set_time_independent(true);
    m_theta_ocean_anomaly_0.metadata()["source"] = m_reference;

    m_theta_ocean_anomaly_1.metadata(0)
            .long_name("potential temperature of the adjacent ocean")
            .units("Kelvin")
            .set_time_independent(true);
    m_theta_ocean_anomaly_1.metadata()["source"] = m_reference;

    m_theta_ocean_anomaly_1X.metadata(0)
            .long_name("potential temperature of the adjacent ocean")
            .units("Kelvin")
            .set_time_independent(true);
    m_theta_ocean_anomaly_1X.metadata()["source"] = m_reference;

    m_salinity_ocean_ref.metadata(0)
        .long_name("salinity of the adjacent ocean")
        .units("g/kg")
        .set_time_independent(true);
    m_salinity_ocean_ref.metadata()["source"] = m_reference;

    // Paleo time slice temperature data annual
    m_salinity_ocean_anomaly_0.metadata(0)
        .long_name("salinity of the adjacent ocean")
        .units("g/kg")
        .set_time_independent(true);
    m_salinity_ocean_anomaly_0.metadata()["source"] = m_reference;

    m_salinity_ocean_anomaly_1.metadata(0)
        .long_name("salinity of the adjacent ocean")
        .units("g/kg")
        .set_time_independent(true);
    m_salinity_ocean_anomaly_1.metadata()["source"] = m_reference;

    m_salinity_ocean_anomaly_1X.metadata(0)
        .long_name("salinity of the adjacent ocean")
        .units("g/kg")
        .set_time_independent(true);
    m_salinity_ocean_anomaly_1X.metadata()["source"] = m_reference;
}

void ClimateIndex::init_forcing() {
    m_log->message(2,
                "**** Initializing the 'Climate Index' ocean model...\n");

    m_w0 = m_w1 = m_w1X = 0.0; // initialise the weights

    auto input_file = m_config->get_string("ocean.climate_index.climate_snapshots.file");
    if (input_file.empty()) {
        throw RuntimeError(PISM_ERROR_LOCATION,
                        "Please specify an ocean snapshots input file\n"
                        "using -ocean_climate_snapshots_file or a command-line option.");
    }
    m_log->message(2,
                    "  Reading salinity and "
                    "theta fields from '%s'...\n", input_file.c_str());

    // Reference fields
    m_theta_ocean_ref.regrid(input_file, io::Default::Nil());
    m_salinity_ocean_ref.regrid(input_file, io::Default::Nil());

    // Annual anomaly for Paleo time slices 0=Glacial, 1=Interglacial, 1X= Super InterGlacial e.g. mPWP
    m_theta_ocean_anomaly_0.regrid(input_file, io::Default::Nil());
    m_theta_ocean_anomaly_1.regrid(input_file, io::Default::Nil());

    m_salinity_ocean_anomaly_0.regrid(input_file, io::Default::Nil());
    m_salinity_ocean_anomaly_1.regrid(input_file, io::Default::Nil());

    try {
        m_theta_ocean_anomaly_1X.regrid(input_file, io::Default::Nil());
        m_salinity_ocean_anomaly_1X.regrid(input_file, io::Default::Nil());
        use_1X = true;
    } catch (...) {
        use_1X = false;
    }
}

void ClimateIndex::update_forcing(double t, double dt, array::Scalar &theta_ocean, array::Scalar &salinity_ocean) {

    m_climate_index->update_weights(t, dt, m_w0, m_w1, m_w1X);

    m_log->message(2,
             "**** Updated weights in ocean: m_w0 = '%f', m_w1 = '%f', m_w1X = '%f' ****\n", m_w0, m_w1, m_w1X);

    theta_ocean.begin_access();
    salinity_ocean.begin_access();

    m_theta_ocean_ref.begin_access();
    m_theta_ocean_anomaly_0.begin_access();
    m_theta_ocean_anomaly_1.begin_access();

    m_salinity_ocean_ref.begin_access();
    m_salinity_ocean_anomaly_0.begin_access();
    m_salinity_ocean_anomaly_1.begin_access();

    if (use_1X) {
        m_theta_ocean_anomaly_1X.begin_access();
        m_salinity_ocean_anomaly_1X.begin_access();
    }

    for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        if (use_1X) {
            theta_ocean(i, j) = m_theta_ocean_ref(i, j) + m_w0 * m_theta_ocean_anomaly_0(i, j) + m_w1 * m_theta_ocean_anomaly_1(i, j) + m_w1X * (m_theta_ocean_anomaly_1X(i, j) - m_theta_ocean_anomaly_1(i, j));
            salinity_ocean(i, j) = m_salinity_ocean_ref(i, j) + m_w0 * m_salinity_ocean_anomaly_0(i, j) + m_w1 * m_salinity_ocean_anomaly_1(i, j) + m_w1X * (m_salinity_ocean_anomaly_1X(i, j) - m_salinity_ocean_anomaly_1(i, j));
        } else {
            theta_ocean(i, j) = m_theta_ocean_ref(i, j) + m_w0 * m_theta_ocean_anomaly_0(i, j) + m_w1 * m_theta_ocean_anomaly_1(i, j);
            salinity_ocean(i, j) = m_salinity_ocean_ref(i, j) + m_w0 * m_salinity_ocean_anomaly_0(i, j) + m_w1 * m_salinity_ocean_anomaly_1(i, j);
        }
    }

    theta_ocean.end_access();
    salinity_ocean.end_access();

    m_theta_ocean_ref.end_access();
    m_theta_ocean_anomaly_0.end_access();
    m_theta_ocean_anomaly_1.end_access();

    m_salinity_ocean_ref.end_access();
    m_salinity_ocean_anomaly_0.end_access();
    m_salinity_ocean_anomaly_1.end_access();

    if (use_1X) {
        m_theta_ocean_anomaly_1X.end_access();
        m_salinity_ocean_anomaly_1X.end_access();
    }
}

} // end of namespace ocean
} // end of namespace pism
