// Copyright (C) 2012-2018 Constantine Khrulev, Ricarda Winkelmann, Ronja Reese, Torsten
// Albrecht, and Matthias Mengel
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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
//
// Please cite this model as:
// 1.
// Antarctic sub-shelf melt rates via PICO
// R. Reese, T. Albrecht, M. Mengel, X. Asay-Davis and R. Winkelmann
// The Cryosphere, 12, 1969-1985, (2018)
// DOI: 10.5194/tc-12-1969-2018
//
// 2.
// A box model of circulation and melting in ice shelf caverns
// D. Olbers & H. Hellmer
// Ocean Dynamics (2010), Volume 60, Issue 1, pp 141–153
// DOI: 10.1007/s10236-009-0252-z

#include <gsl/gsl_math.h> // GSL_NAN

#include "pism/util/ConfigInterface.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/Mask.hh"
#include "pism/util/Vars.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/Time.hh"
#include "pism/geometry/Geometry.hh"

#include "pism/coupler/util/options.hh"

#include "Pico.hh"
#include "PicoGeometry.hh"
#include "PicoPhysics.hh"

namespace pism {
namespace ocean {

Pico::Pico(IceGrid::ConstPtr g)
    : CompleteOceanModel(g, NULL), m_geometry(new PicoGeometry(g)) {

  ForcingOptions opt(*m_grid->ctx(), "ocean.pico");

  {
    unsigned int buffer_size = m_config->get_double("climate_forcing.buffer_size");
    unsigned int evaluations_per_year = m_config->get_double("climate_forcing.evaluations_per_year");
    bool periodic = opt.period > 0;

    PIO file(m_grid->com, "netcdf3", opt.filename, PISM_READONLY);

    m_theta_ocean = IceModelVec2T::ForcingField(m_grid,
                                                file,
                                                "theta_ocean",
                                                "", // no standard name
                                                buffer_size,
                                                evaluations_per_year,
                                                periodic);

    m_salinity_ocean = IceModelVec2T::ForcingField(m_grid,
                                                   file,
                                                   "salinity_ocean",
                                                   "", // no standard name
                                                   buffer_size,
                                                   evaluations_per_year,
                                                   periodic);
  }

  m_theta_ocean->set_attrs("climate_forcing",
                           "potential temperature of the adjacent ocean",
                           "Kelvin", "");

  m_salinity_ocean->set_attrs("climate_forcing",
                              "salinity of the adjacent ocean",
                              "g/kg", "");

  m_basin_mask.create(m_grid, "basins", WITH_GHOSTS);
  m_basin_mask.set_attrs("climate_forcing", "mask determines basins for PICO", "", "");

  // computed salinity in ocean boxes
  m_Soc.create(m_grid, "pico_salinity", WITHOUT_GHOSTS);
  m_Soc.set_attrs("model_state", "ocean salinity field", "g/kg", "ocean salinity field");
  m_Soc.metadata().set_double("_FillValue", 0.0);

  // salinity input for box 1
  m_Soc_box0.create(m_grid, "pico_salinity_box0", WITHOUT_GHOSTS);
  m_Soc_box0.set_attrs("model_state", "ocean base salinity field", "g/kg", "ocean base salinity field");
  m_Soc_box0.metadata().set_double("_FillValue", 0.0);

  // computed temperature in ocean boxes
  m_Toc.create(m_grid, "pico_temperature", WITHOUT_GHOSTS);
  m_Toc.set_attrs("model_state", "ocean temperature field", "K", "ocean temperature field");
  m_Toc.metadata().set_double("_FillValue", 0.0);

  // temperature input for box 1
  m_Toc_box0.create(m_grid, "pico_temperature_box0", WITHOUT_GHOSTS);
  m_Toc_box0.set_attrs("model_state", "ocean base temperature", "K", "ocean base temperature");
  m_Toc_box0.metadata().set_double("_FillValue", 0.0);

  m_T_star.create(m_grid, "pico_T_star", WITHOUT_GHOSTS);
  m_T_star.set_attrs("model_state", "T_star field", "degree C", "T_star field");
  m_T_star.metadata().set_double("_FillValue", 0.0);

  m_overturning.create(m_grid, "pico_overturning", WITHOUT_GHOSTS);
  m_overturning.set_attrs("model_state", "cavity overturning", "m^3 s-1", "cavity overturning");
  m_overturning.metadata().set_double("_FillValue", 0.0);

  m_basal_melt_rate.create(m_grid, "pico_basal_melt_rate", WITHOUT_GHOSTS);
  m_basal_melt_rate.set_attrs("model_state", "PICO sub-shelf melt rate", "m/s", "PICO sub-shelf melt rate");
  m_basal_melt_rate.metadata().set_string("glaciological_units", "m year-1");
  m_basal_melt_rate.metadata().set_double("_FillValue", 0.0);

  m_shelf_base_temperature->metadata().set_double("_FillValue", 0.0);

  m_n_basins = 0;

  m_n_boxes  = m_config->get_double("ocean.pico.number_of_boxes");
}


Pico::~Pico() {
  // empty
}


void Pico::init_impl(const Geometry &geometry) {
  (void) geometry;

  m_log->message(2, "* Initializing the Potsdam Ice-shelf Cavity mOdel for the ocean ...\n");

  ForcingOptions opt(*m_grid->ctx(), "ocean.pico");

  m_theta_ocean->init(opt.filename, opt.period, opt.reference_time);
  m_salinity_ocean->init(opt.filename, opt.period, opt.reference_time);

  m_basin_mask.regrid(opt.filename, CRITICAL);

  // FIXME: m_n_basins is a misnomer
  m_n_basins = m_basin_mask.max() + 1;

  m_log->message(4, "PICO basin min=%f,max=%f\n", m_basin_mask.min(), m_basin_mask.max());

  PicoPhysics physics(*m_config);

  m_log->message(2, "  -Using %d drainage basins and values: \n"
                    "   gamma_T= %.2e, overturning_coeff = %.2e... \n",
                 m_n_basins, physics.gamma_T(), physics.overturning_coeff());

  m_log->message(2, "  -Depth of continental shelf for computation of temperature and salinity input\n"
                    "   is set for whole domain to continental_shelf_depth=%.0f meter\n",
                 physics.continental_shelf_depth());

  // read time-independent data right away:
  if (m_theta_ocean->n_records() == 1 and m_salinity_ocean->n_records() == 1) {
    m_theta_ocean->update(m_grid->ctx()->time()->current(), 0.0);
    m_salinity_ocean->update(m_grid->ctx()->time()->current(), 0.0);
  }
}

void Pico::define_model_state_impl(const PIO &output) const {

  m_basin_mask.define(output);
  m_Soc_box0.define(output);
  m_Toc_box0.define(output);
  m_overturning.define(output);

  OceanModel::define_model_state_impl(output);
}

void Pico::write_model_state_impl(const PIO &output) const {

  m_basin_mask.write(output);
  m_Soc_box0.write(output);
  m_Toc_box0.write(output);
  m_overturning.write(output);

  OceanModel::define_model_state_impl(output);
}

void Pico::update_impl(const Geometry &geometry, double t, double dt) {

  m_theta_ocean->update(t, dt);
  m_salinity_ocean->update(t, dt);

  m_theta_ocean->average(t, dt);
  m_salinity_ocean->average(t, dt);

  // set values that will be used outside of floating ice areas
  {
    double T_fill_value   = m_config->get_double("constants.fresh_water.melting_point_temperature");
    double Toc_fill_value = m_Toc.metadata().get_double("_FillValue");
    double Soc_fill_value = m_Soc.metadata().get_double("_FillValue");
    double M_fill_value   = m_basal_melt_rate.metadata().get_double("_FillValue");
    double O_fill_value   = m_overturning.metadata().get_double("_FillValue");

    m_shelf_base_temperature->set(T_fill_value);
    m_basal_melt_rate.set(M_fill_value);
    m_Toc.set(Toc_fill_value);
    m_Soc.set(Soc_fill_value);
    m_overturning.set(O_fill_value);
    m_T_star.set(Toc_fill_value);
  }

  PicoPhysics physics(*m_config);

  const IceModelVec2S &ice_thickness    = geometry.ice_thickness;
  const IceModelVec2CellType &cell_type = geometry.cell_type;
  const IceModelVec2S &bed_elevation    = geometry.bed_elevation;

  // Geometric part of PICO
  m_geometry->update(bed_elevation, cell_type);

  // FIXME: m_n_shelves is not really the number of shelves.
  m_n_shelves = m_geometry->ice_shelf_mask().max() + 1;

  // Physical part of PICO
  {

    // prepare ocean input temperature and salinity
    {
      std::vector<double> basin_temperature(m_n_basins);
      std::vector<double> basin_salinity(m_n_basins);

      compute_ocean_input_per_basin(physics, m_basin_mask, m_geometry->continental_shelf_mask(), *m_salinity_ocean,
                                    *m_theta_ocean, basin_temperature, basin_salinity); // per basin

      set_ocean_input_fields(physics, ice_thickness, cell_type, m_basin_mask, m_geometry->ice_shelf_mask(),
                             basin_temperature, basin_salinity, m_Toc_box0, m_Soc_box0); // per shelf
    }

    // Use the Beckmann-Goosse parameterization to set reasonable values throughout the
    // domain.
    beckmann_goosse(physics,
                    ice_thickness,                // input
                    m_geometry->ice_shelf_mask(), // input
                    cell_type,                    // input
                    m_Toc_box0,                   // input
                    m_Soc_box0,                   // input
                    m_basal_melt_rate,
                    *m_shelf_base_temperature,
                    m_Toc,
                    m_Soc);

    // In ice shelves, replace Beckmann-Goosse values using the Olbers and Hellmer model.
    process_box1(physics,
                 ice_thickness,                             // input
                 m_geometry->ice_shelf_mask(),              // input
                 m_geometry->box_mask(),                    // input
                 m_Toc_box0,                                // input
                 m_Soc_box0,                                // input
                 m_basal_melt_rate,
                 *m_shelf_base_temperature,
                 m_T_star,
                 m_Toc,
                 m_Soc,
                 m_overturning);

    process_other_boxes(physics,
                        ice_thickness,                // input
                        m_geometry->ice_shelf_mask(), // input
                        m_geometry->box_mask(),       // input
                        m_basal_melt_rate,
                        *m_shelf_base_temperature,
                        m_T_star,
                        m_Toc,
                        m_Soc);
  }

  m_shelf_base_mass_flux->copy_from(m_basal_melt_rate);
  m_shelf_base_mass_flux->scale(physics.ice_density());

  m_melange_back_pressure_fraction->set(0.0);
}


MaxTimestep Pico::max_timestep_impl(double t) const {
  (void) t;

  return MaxTimestep("ocean pico");
}


//! Compute temperature and salinity input from ocean data by averaging.

//! We average the ocean data over the continental shelf reagion for each basin.
//! We use dummy ocean data if no such average can be calculated.


void Pico::compute_ocean_input_per_basin(const PicoPhysics &physics, const IceModelVec2Int &basin_mask,
                                         const IceModelVec2Int &continental_shelf_mask,
                                         const IceModelVec2S &salinity_ocean, const IceModelVec2S &theta_ocean,
                                         std::vector<double> &temperature, std::vector<double> &salinity) {

  std::vector<int> count(m_n_basins, 0);

  temperature.resize(m_n_basins);
  salinity.resize(m_n_basins);
  for (int basin_id = 0; basin_id < m_n_basins; basin_id++) {
    temperature[basin_id] = 0.0;
    salinity[basin_id]    = 0.0;
  }

  IceModelVec::AccessList list{ &theta_ocean, &salinity_ocean, &basin_mask, &continental_shelf_mask };

  // compute the sum for each basin for region that intersects with the continental shelf
  // area and is not covered by an ice shelf. (continental shelf mask excludes ice shelf
  // areas)
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (continental_shelf_mask.as_int(i, j) == 2) {
      int basin_id = basin_mask.as_int(i, j);

      count[basin_id] += 1;
      salinity[basin_id] += salinity_ocean(i, j);
      temperature[basin_id] += theta_ocean(i, j);
    }
  }

  // Divide by number of grid cells if more than zero cells belong to the basin. if no
  // ocean_contshelf_mask values intersect with the basin, count is zero. In such case,
  // use dummy temperature and salinity. This could happen, for example, if the ice shelf
  // front advances beyond the continental shelf break.
  for (int basin_id = 0; basin_id < m_n_basins; basin_id++) {

    count[basin_id]       = GlobalSum(m_grid->com, count[basin_id]);
    salinity[basin_id]    = GlobalSum(m_grid->com, salinity[basin_id]);
    temperature[basin_id] = GlobalSum(m_grid->com, temperature[basin_id]);

    // if basin is not dummy basin 0 or there are no ocean cells in this basin to take the mean over.
    if (basin_id > 0 && count[basin_id] == 0) {
      m_log->message(2, "PICO ocean WARNING: basin %d contains no cells with ocean data on continental shelf\n"
                        "(no values with ocean_contshelf_mask=2).\n"
                        "No mean salinity or temperature values are computed, instead using\n"
                        "the standard values T_dummy =%.3f, S_dummy=%.3f.\n"
                        "This might bias your basal melt rates, check your input data carefully.\n",
                     basin_id, physics.T_dummy(), physics.S_dummy());

      temperature[basin_id] = physics.T_dummy();
      salinity[basin_id]    = physics.S_dummy();

    } else {

      salinity[basin_id] /= count[basin_id];
      temperature[basin_id] /= count[basin_id];

      m_log->message(5, "  %d: temp =%.3f, salinity=%.3f\n", basin_id, temperature[basin_id], salinity[basin_id]);
    }
  }
}

//! Set ocean ocean input from box 0 as boundary condition for box 1.

//! Set ocean temperature and salinity (Toc_box0, Soc_box0)
//! from box 0 (in front of the ice shelf) as inputs for
//! box 1, which is the ocean box adjacent to the grounding line.
//!
//! We enforce that Toc_box0 is always at least the local pressure melting point.
void Pico::set_ocean_input_fields(const PicoPhysics &physics, const IceModelVec2S &ice_thickness,
                                  const IceModelVec2CellType &mask, const IceModelVec2Int &basin_mask,
                                  const IceModelVec2Int &shelf_mask, const std::vector<double> basin_temperature,
                                  const std::vector<double> basin_salinity, IceModelVec2S &Toc_box0,
                                  IceModelVec2S &Soc_box0) {

  IceModelVec::AccessList list{ &ice_thickness, &basin_mask, &Soc_box0, &Toc_box0, &mask, &shelf_mask };

  std::vector<std::vector<int> > n_shelf_cells_per_basin(m_n_shelves, std::vector<int>(m_n_basins, 0));
  std::vector<int> n_shelf_cells(m_n_shelves, 0);

  // 1) count the number of cells in each shelf
  // 2) count the number of cells in the intersection of each shelf with all the basins
  {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();
      int s = shelf_mask.as_int(i, j);
      int b = basin_mask.as_int(i, j);
      n_shelf_cells_per_basin[s][b]++;
      n_shelf_cells[s]++;
    }

    for (int s = 0; s < m_n_shelves; s++) {
      n_shelf_cells[s] = GlobalSum(m_grid->com, n_shelf_cells[s]);
      for (int b = 0; b < m_n_basins; b++) {
        n_shelf_cells_per_basin[s][b] = GlobalSum(m_grid->com, n_shelf_cells_per_basin[s][b]);
      }
    }
  }

  // now set potential temperature and salinity box 0:

  int low_temperature_counter = 0;
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // make sure all temperatures are zero at the beginning of each time step
    Toc_box0(i, j) = 0.0; // in K
    Soc_box0(i, j) = 0.0; // in psu

    int s = shelf_mask.as_int(i, j);

    if (mask.as_int(i, j) == MASK_FLOATING and s > 0) {
      // note: shelf_mask = 0 in lakes

      // weighted input depending on the number of shelf cells in each basin
      for (int b = 1; b < m_n_basins; b++) { //Note: b=0 yields nan
        Toc_box0(i, j) += basin_temperature[b] * n_shelf_cells_per_basin[s][b] / (double)n_shelf_cells[s];
        Soc_box0(i, j) += basin_salinity[b] * n_shelf_cells_per_basin[s][b] / (double)n_shelf_cells[s];
      }

      double theta_pm = physics.theta_pm(Soc_box0(i, j), physics.pressure(ice_thickness(i, j)));

      // temperature input for grounding line box should not be below pressure melting point
      if (Toc_box0(i, j) < theta_pm) {
        // Setting Toc_box0 a little higher than theta_pm ensures that later equations are well solvable.
        Toc_box0(i, j) = theta_pm + 0.001;
        low_temperature_counter += 1;
      }
    }
  }

  low_temperature_counter = GlobalSum(m_grid->com, low_temperature_counter);
  if (low_temperature_counter > 0) {
    m_log->message(2, "PICO ocean warning: temperature has been below pressure melting temperature in %d cases,\n"
                      "setting it to pressure melting temperature\n",
                   low_temperature_counter);
  }
}

/*!
 * Use the simpler parameterization due to [@ref BeckmannGoosse2003] to set default
 * sub-shelf temperature and melt rate values.
 *
 * At grid points containing floating ice not connected to the ocean, set the basal melt
 * rate to zero and set basal temperature to the pressure melting point.
 */
void Pico::beckmann_goosse(const PicoPhysics &physics,
                           const IceModelVec2S &ice_thickness,
                           const IceModelVec2Int &shelf_mask,
                           const IceModelVec2CellType &cell_type,
                           const IceModelVec2S &Toc_box0,
                           const IceModelVec2S &Soc_box0,
                           IceModelVec2S &basal_melt_rate,
                           IceModelVec2S &basal_temperature,
                           IceModelVec2S &Toc,
                           IceModelVec2S &Soc) {

  const double T0          = m_config->get_double("constants.fresh_water.melting_point_temperature"),
               beta_CC     = m_config->get_double("constants.ice.beta_Clausius_Clapeyron"),
               g           = m_config->get_double("constants.standard_gravity"),
               ice_density = m_config->get_double("constants.ice.density");

  IceModelVec::AccessList list{ &ice_thickness, &cell_type, &shelf_mask,      &Toc_box0,          &Soc_box0,
                                &Toc,           &Soc,       &basal_melt_rate, &basal_temperature };

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (cell_type.floating_ice(i, j)) {
      if (shelf_mask.as_int(i, j) > 0) {
        double pressure = physics.pressure(ice_thickness(i, j));

        basal_melt_rate(i, j) =
            physics.melt_rate_beckmann_goosse(physics.theta_pm(Soc_box0(i, j), pressure), Toc_box0(i, j));
        basal_temperature(i, j) = physics.T_pm(Soc_box0(i, j), pressure);

        // diagnostic outputs
        Toc(i, j) = Toc_box0(i, j); // in Kelvin
        Soc(i, j) = Soc_box0(i, j); // in psu
      } else {
        // Floating ice cells not connected to the ocean.
        const double pressure = ice_density * g * ice_thickness(i, j); // FIXME issue #15

        basal_temperature(i, j) = T0 - beta_CC * pressure;
        basal_melt_rate(i, j)   = 0.0;
      }
    }
  }
}


void Pico::process_box1(const PicoPhysics &physics,
                        const IceModelVec2S &ice_thickness,
                        const IceModelVec2Int &shelf_mask,
                        const IceModelVec2Int &box_mask,
                        const IceModelVec2S &Toc_box0,
                        const IceModelVec2S &Soc_box0,
                        IceModelVec2S &basal_melt_rate,
                        IceModelVec2S &basal_temperature,
                        IceModelVec2S &T_star,
                        IceModelVec2S &Toc,
                        IceModelVec2S &Soc,
                        IceModelVec2S &overturning) {

  std::vector<double> box1_area(m_n_shelves);

  compute_box_area(1, shelf_mask, box_mask, box1_area);

  IceModelVec::AccessList list{ &ice_thickness, &shelf_mask, &box_mask,    &T_star,          &Toc_box0,          &Toc,
                                &Soc_box0,      &Soc,        &overturning, &basal_melt_rate, &basal_temperature };

  int n_Toc_failures = 0;

  // basal melt rate, ambient temperature and salinity and overturning calculation
  // for each box1 grid cell.
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    int shelf_id = shelf_mask.as_int(i, j);

    if (box_mask.as_int(i, j) == 1 and shelf_id > 0) {

      const double pressure = physics.pressure(ice_thickness(i, j));

      T_star(i, j) = physics.T_star(Soc_box0(i, j), Toc_box0(i, j), pressure);

      auto Toc_box1 = physics.Toc_box1(box1_area[shelf_id], T_star(i, j), Soc_box0(i, j), Toc_box0(i, j));

      // This can only happen if T_star > 0.25*p_coeff, in particular T_star > 0
      // which can only happen for values of Toc_box0 close to the local pressure melting point
      if (Toc_box1.failed) {
        m_log->message(5, "PICO ocean WARNING: negative square root argument at %d, %d\n"
                          "probably because of positive T_star=%f \n"
                          "Not aborting, but setting square root to 0... \n",
                       i, j, T_star(i, j));

        n_Toc_failures += 1;
      }

      Toc(i, j) = Toc_box1.value;
      Soc(i, j) = physics.Soc_box1(Toc_box0(i, j), Soc_box0(i, j), Toc(i, j)); // in psu

      overturning(i, j) = physics.overturning(Soc_box0(i, j), Soc(i, j), Toc_box0(i, j), Toc(i, j));

      // main outputs
      basal_melt_rate(i, j)    = physics.melt_rate(physics.theta_pm(Soc(i, j), pressure), Toc(i, j));
      basal_temperature(i, j) = physics.T_pm(Soc(i, j), pressure);
    }
  }

  n_Toc_failures = GlobalSum(m_grid->com, n_Toc_failures);
  if (n_Toc_failures > 0) {
    m_log->message(2, "PICO ocean warning: square-root argument for temperature calculation "
                      "has been negative in %d cases!\n",
                   n_Toc_failures);
  }
}

void Pico::process_other_boxes(const PicoPhysics &physics,
                               const IceModelVec2S &ice_thickness,
                               const IceModelVec2Int &shelf_mask,
                               const IceModelVec2Int &box_mask,
                               IceModelVec2S &basal_melt_rate,
                               IceModelVec2S &basal_temperature,
                               IceModelVec2S &T_star,
                               IceModelVec2S &Toc,
                               IceModelVec2S &Soc) {

  std::vector<double> overturning(m_n_shelves, 0.0);
  std::vector<double> salinity(m_n_shelves, 0.0);
  std::vector<double> temperature(m_n_shelves, 0.0);

  // get average overturning from box 1 that is used as input later
  compute_box_average(1, m_overturning, shelf_mask, box_mask, overturning);

  std::vector<bool> use_beckmann_goosse(m_n_shelves);

  IceModelVec::AccessList list{ &ice_thickness, &shelf_mask,      &box_mask,           &T_star,   &Toc,
                                &Soc,           &basal_melt_rate, &basal_temperature };

  // Iterate over all boxes i for i > 1
  for (int box = 2; box <= m_n_boxes; ++box) {

    compute_box_average(box - 1, Toc, shelf_mask, box_mask, temperature);
    compute_box_average(box - 1, Soc, shelf_mask, box_mask, salinity);

    // find all the shelves where we should fall back to the Beckmann-Goosse
    // parameterization
    for (int s = 1; s < m_n_shelves; ++s) {
      if (salinity[s] == 0.0 or temperature[s] == 0.0 or overturning[s] == 0.0) {
        use_beckmann_goosse[s] = true;
      } else {
        use_beckmann_goosse[s] = false;
      }
    }

    std::vector<double> box_area;
    compute_box_area(box, shelf_mask, box_mask, box_area);

    int n_beckmann_goosse_cells = 0;

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      int shelf_id = shelf_mask.as_int(i, j);

      if (box_mask.as_int(i, j) == box and shelf_id > 0) {

        if (use_beckmann_goosse[shelf_id]) {
          n_beckmann_goosse_cells += 1;
          continue;
        }

        // Get the input from previous box
        const double
          S_previous       = salinity[shelf_id],
          T_previous       = temperature[shelf_id],
          overturning_box1 = overturning[shelf_id];

        {
          double pressure = physics.pressure(ice_thickness(i, j));

          // diagnostic outputs
          T_star(i, j) = physics.T_star(S_previous, T_previous, pressure);
          Toc(i, j)    = physics.Toc(box_area[shelf_id], T_previous, T_star(i, j), overturning_box1, S_previous);
          Soc(i, j)    = physics.Soc(S_previous, T_previous, Toc(i, j));

          // main outputs: basal melt rate and temperature
          basal_melt_rate(i, j)   = physics.melt_rate(physics.theta_pm(Soc(i, j), pressure), Toc(i, j));
          basal_temperature(i, j) = physics.T_pm(Soc(i, j), pressure);
        }
      }
    } // loop over grid points

    n_beckmann_goosse_cells = GlobalSum(m_grid->com, n_beckmann_goosse_cells);
    if (n_beckmann_goosse_cells > 0) {
      m_log->message(2, "PICO ocean WARNING: box %d, no boundary data from previous box in %d case(s)!\n"
                        "switching to Beckmann Goosse (2003) meltrate calculation\n",
                     box, n_beckmann_goosse_cells);
    }

  } // loop over boxes
}


// Write diagnostic variables to extra files if requested
DiagnosticList Pico::diagnostics_impl() const {

  DiagnosticList result = {
    { "basins",                 Diagnostic::wrap(m_basin_mask) },
    { "pico_overturning",       Diagnostic::wrap(m_overturning) },
    { "pico_salinity_box0",     Diagnostic::wrap(m_Soc_box0) },
    { "pico_temperature_box0",  Diagnostic::wrap(m_Toc_box0) },
    { "pico_box_mask",          Diagnostic::wrap(m_geometry->box_mask()) },
    { "pico_shelf_mask",        Diagnostic::wrap(m_geometry->ice_shelf_mask()) },
    { "pico_ice_rise_mask",     Diagnostic::wrap(m_geometry->ice_rise_mask()) },
    { "pico_basal_melt_rate",   Diagnostic::wrap(m_basal_melt_rate) },
    { "pico_contshelf_mask",    Diagnostic::wrap(m_geometry->continental_shelf_mask()) },
    { "pico_salinity",          Diagnostic::wrap(m_Soc) },
    { "pico_temperature",       Diagnostic::wrap(m_Toc) },
    { "pico_T_star",            Diagnostic::wrap(m_T_star) },
    { "pico_basal_temperature", Diagnostic::wrap(*m_shelf_base_temperature) },
  };

  return combine(result, OceanModel::diagnostics_impl());
}

/*!
 * For each shelf, compute average of a given field over the box with id `box_id`.
 *
 * This method is used to get inputs from a previous box for the next one.
 */
void Pico::compute_box_average(int box_id,
                               const IceModelVec2S &field,
                               const IceModelVec2Int &shelf_mask,
                               const IceModelVec2Int &box_mask,
                               std::vector<double> &result) {

  IceModelVec::AccessList list{ &field, &shelf_mask, &box_mask };

  std::vector<int> n_cells_per_box(m_n_shelves, 0);

  // fill results with zeros
  result.resize(m_n_shelves);
  for (int s = 0; s < m_n_shelves; ++s) {
    result[s] = 0.0;
  }

  // compute the sum of field in each shelf's box box_id
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    int shelf_id = shelf_mask.as_int(i, j);

    if (box_mask.as_int(i, j) == box_id) {
      n_cells_per_box[shelf_id] += 1;
      result[shelf_id] += field(i, j);
    }
  }

  // compute the global sum and average
  for (int s = 0; s < m_n_shelves; ++s) {
    auto n_cells = GlobalSum(m_grid->com, n_cells_per_box[s]);

    result[s] = GlobalSum(m_grid->com, result[s]);

    if (n_cells > 0) {
      result[s] /= (double)n_cells;
    }
  }
}

/*!
 * For all shelves compute areas of boxes with id `box_id`.
 *
 * @param[in] box_mask box index mask
 * @param[out] result resulting box areas.
 *
 * Note: shelf and box indexes start from 1.
 */
void Pico::compute_box_area(int box_id,
                            const IceModelVec2Int &shelf_mask,
                            const IceModelVec2Int &box_mask,
                            std::vector<double> &result) {
  result.resize(m_n_shelves);

  IceModelVec::AccessList list{ &shelf_mask, &box_mask };

  auto cell_area = m_grid->cell_area();

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    int shelf_id = shelf_mask.as_int(i, j);

    if (shelf_id > 0 and box_mask.as_int(i, j) == box_id) {
      result[shelf_id] += cell_area;
    }
  }

  // compute global sums
  for (int s = 1; s < m_n_shelves; ++s) {
    result[s] = GlobalSum(m_grid->com, result[s]);
  }
}

} // end of namespace ocean
} // end of namespace pism
