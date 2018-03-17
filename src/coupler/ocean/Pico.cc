// Copyright (C) 2012-2018 Ricarda Winkelmann, Ronja Reese, Torsten Albrecht
// and Matthias Mengel
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
// The Cryosphere Discussions (2017)
// DOI: 10.5194/tc-2017-70
//
// 2.
// A box model of circulation and melting in ice shelf caverns
// D. Olbers & H. Hellmer
// Ocean Dynamics (2010), Volume 60, Issue 1, pp 141–153
// DOI: 10.1007/s10236-009-0252-z

#include <gsl/gsl_math.h>       // GSL_NAN

#include "Pico.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/Mask.hh"
#include "pism/util/Vars.hh"
#include "pism/util/iceModelVec.hh"

namespace pism {
namespace ocean {

Pico::Pico(IceGrid::ConstPtr g) : PGivenClimate<CompleteOceanModel, CompleteOceanModel>(g, NULL) {

  m_option_prefix = "-ocean_pico";

  // will be de-allocated by the parent's destructor
  m_theta_ocean    = new IceModelVec2T;
  m_salinity_ocean = new IceModelVec2T;

  m_fields["theta_ocean"]    = m_theta_ocean;
  m_fields["salinity_ocean"] = m_salinity_ocean;

  process_options();

  m_exicerises_set = m_config->get_boolean("ocean.pico.exclude_icerises");

  std::map<std::string, std::string> standard_names;
  set_vec_parameters(standard_names);

  m_Mx = m_grid->Mx(), m_My = m_grid->My();

  m_theta_ocean->create(m_grid, "theta_ocean");
  m_theta_ocean->set_attrs("climate_forcing", "absolute potential temperature of the adjacent ocean", "Kelvin", "");

  m_salinity_ocean->create(m_grid, "salinity_ocean");
  m_salinity_ocean->set_attrs("climate_forcing", "salinity of the adjacent ocean", "g/kg", "");

  m_cbasins.create(m_grid, "basins", WITH_GHOSTS);
  m_cbasins.set_attrs("climate_forcing", "mask determines basins for PICO", "", "");

  // mask to identify ice shelves
  m_shelf_mask.create(m_grid, "pico_shelf_mask", WITH_GHOSTS);
  m_shelf_mask.set_attrs("model_state", "mask for individual ice shelves", "", "");

  // mask to identify the ocean boxes
  m_ocean_box_mask.create(m_grid, "pico_ocean_box_mask", WITH_GHOSTS);
  m_ocean_box_mask.set_attrs("model_state", "mask displaying ocean box model grid", "", "");

  // mask to identify the ice rises
  m_icerise_mask.create(m_grid, "pico_icerise_mask", WITH_GHOSTS);
  m_icerise_mask.set_attrs("model_state", "mask displaying ice rises", "", "");

  // mask displaying continental shelf - region where mean salinity and ocean temperature is calculated
  m_ocean_contshelf_mask.create(m_grid, "pico_ocean_contshelf_mask", WITH_GHOSTS);
  m_ocean_contshelf_mask.set_attrs("model_state", "mask displaying ocean region for parameter input", "", "");

  // mask displaying open ocean - ice-free regions below sea-level except 'holes' in ice shelves
  m_ocean_mask.create(m_grid, "pico_ocean_mask", WITH_GHOSTS);
  m_ocean_mask.set_attrs("model_state", "mask displaying open ocean", "", "");

  // mask displaying subglacial lakes - floating regions with no connection to the ocean
  m_lake_mask.create(m_grid, "pico_lake_mask", WITH_GHOSTS);
  m_lake_mask.set_attrs("model_state", "mask displaying subglacial lakes", "", "");

  // mask with distance (in boxes) to grounding line
  m_DistGL.create(m_grid, "pico_dist_grounding_line", WITH_GHOSTS);
  m_DistGL.set_attrs("model_state", "mask displaying distance to grounding line", "", "");

  // mask with distance (in boxes) to ice front
  m_DistIF.create(m_grid, "pico_dist_iceshelf_front", WITH_GHOSTS);
  m_DistIF.set_attrs("model_state", "mask displaying distance to ice shelf calving front", "", "");

  // computed salinity in ocean boxes
  m_Soc.create(m_grid, "pico_Soc", WITHOUT_GHOSTS);
  m_Soc.set_attrs("model_state", "ocean salinity field", "g/kg", "ocean salinity field");

  // salinity input for box 1
  m_Soc_box0.create(m_grid, "pico_salinity_box0", WITHOUT_GHOSTS);
  m_Soc_box0.set_attrs("model_state", "ocean base salinity field", "g/kg", "ocean base salinity field");

  // computed temperature in ocean boxes
  m_Toc.create(m_grid, "pico_Toc", WITHOUT_GHOSTS);
  m_Toc.set_attrs("model_state", "ocean temperature field", "K", "ocean temperature field");

  // temperature input for box 1
  m_Toc_box0.create(m_grid, "pico_temperature_box0", WITHOUT_GHOSTS);
  m_Toc_box0.set_attrs("model_state", "ocean base temperature", "K", "ocean base temperature");

  // in ocean box i: T_star = aS_{i-1} + b -c p_i - T_{i-1} with T_{-1} = Toc_box0 and S_{-1}=Soc_box0
  // FIXME convert to internal field
  m_T_star.create(m_grid, "pico_T_star", WITHOUT_GHOSTS);
  m_T_star.set_attrs("model_state", "T_star field", "degree C", "T_star field");

  m_overturning.create(m_grid, "pico_overturning", WITHOUT_GHOSTS);
  m_overturning.set_attrs("model_state", "cavity overturning", "m^3 s-1", "cavity overturning"); // no CF standard_name?

  m_basalmeltrate_shelf.create(m_grid, "pico_bmelt_shelf", WITHOUT_GHOSTS);
  m_basalmeltrate_shelf.set_attrs("model_state", "PICO sub-shelf melt rate", "m/s", "PICO sub-shelf melt rate");
  m_basalmeltrate_shelf.metadata().set_string("glaciological_units", "m year-1");
  //basalmeltrate_shelf.write_in_glaciological_units = true;

  // TODO: this may be initialized to NA, it should only have valid values below ice shelves.
  m_T_pressure_melting.create(m_grid, "pico_T_pressure_melting", WITHOUT_GHOSTS);
  m_T_pressure_melting.set_attrs(
      "model_state", "pressure melting temperature at ice shelf base", "Kelvin",
      "pressure melting temperature at ice shelf base"); // no CF standard_name? // This is the in-situ pressure melting point


  // Initialize this early so that we can check the validity of the "basins" mask read from a file
  // in Pico::init_impl(). This number is hard-wired, so I don't think it matters that it did not
  // come from Pico::Constants.
  m_n_basins = 20;

  // This will be re-set by identify_shelf_mask()
  m_n_shelves = 1;
}


Pico::~Pico() {
  // empty
}


void Pico::init_impl() {

  m_t = m_dt = GSL_NAN; // every re-init restarts the clock

  m_log->message(2, "* Initializing the Potsdam Ice-shelf Cavity mOdel for the ocean ...\n");

  m_theta_ocean->init(m_filename, m_bc_period, m_bc_reference_time);
  m_salinity_ocean->init(m_filename, m_bc_period, m_bc_reference_time);

  m_cbasins.regrid(m_filename, CRITICAL);

  m_log->message(4, "PICO basin min=%f,max=%f\n", m_cbasins.min(), m_cbasins.max());

  BoxModel box_model(*m_config);

  m_n_basins = m_config->get_double("ocean.pico.number_of_basins");
  m_n_boxes  = m_config->get_double("ocean.pico.number_of_boxes");

  m_Toc_box0_vec.resize(m_n_basins);
  m_Soc_box0_vec.resize(m_n_basins);

  m_log->message(2, "  -Using %d drainage basins and values: \n"
                    "   gamma_T= %.2e, overturning_coeff = %.2e... \n",
                 m_n_basins, box_model.gamma_T(), box_model.overturning_coeff());

  m_log->message(2, "  -Depth of continental shelf for computation of temperature and salinity input\n"
                    "   is set for whole domain to continental_shelf_depth=%.0f meter\n",
                 box_model.continental_shelf_depth());

  round_basins(m_cbasins);

  // Range basins_range = cbasins.range();

  // if (basins_range.min < 0 or basins_range.max > numberOfBasins - 1) {
  //   throw RuntimeError::formatted(PISM_ERROR_LOCATION,
  //                                 "Some basin numbers in %s read from %s are invalid:"
  //                                 "allowed range is [0, %d], found [%d, %d]",
  //                                 cbasins.get_name().c_str(), m_filename.c_str(),
  //                                 numberOfBasins - 1,
  //                                 (int)basins_range.min, (int)basins_range.max);
  // }

  // read time-independent data right away:
  if (m_theta_ocean->get_n_records() == 1 && m_salinity_ocean->get_n_records() == 1) {
    update(m_grid->ctx()->time()->current(), 0); // dt is irrelevant
  }
}

void Pico::define_model_state_impl(const PIO &output) const {

  m_cbasins.define(output);
  m_ocean_box_mask.define(output);
  m_shelf_mask.define(output);
  m_Soc_box0.define(output);
  m_Toc_box0.define(output);
  m_overturning.define(output);
  //basalmeltrate_shelf.define(output);

  OceanModel::define_model_state_impl(output);
}

void Pico::write_model_state_impl(const PIO &output) const {

  m_cbasins.write(output);
  m_ocean_box_mask.write(output);
  m_shelf_mask.write(output);
  m_Soc_box0.write(output);
  m_Toc_box0.write(output);
  m_overturning.write(output);
  //basalmeltrate_shelf.write(output);

  OceanModel::define_model_state_impl(output);
}

void Pico::update_impl(double my_t, double my_dt) {

  // Make sure that sea water salinity and sea water potential
  // temperature fields are up to date:
  update_internal(my_t, my_dt);

  m_theta_ocean->average(m_t, m_dt);
  m_salinity_ocean->average(m_t, m_dt);

  BoxModel model(*m_config);

  // Geometric part of PICO
  // define the ocean boxes below the ice shelves
  identifyMASK(m_ocean_contshelf_mask, "ocean_continental_shelf");
  if (m_exicerises_set) {
    identifyMASK(m_icerise_mask, "icerises");
  }
  identifyMASK(m_ocean_mask, "ocean");
  identifyMASK(m_lake_mask, "lakes");
  identify_shelf_mask();
  round_basins(m_cbasins);
  compute_distances();
  identify_ocean_box_mask(model);

  m_ocean_box_mask.update_ghosts();

  // Physical part of PICO

  // prepare ocean input temperature and salinity
  {
    compute_ocean_input_per_basin(model,
                                  m_cbasins, m_ocean_contshelf_mask,
                                  *m_salinity_ocean, *m_theta_ocean); // per basin

    const IceModelVec2S &ice_thickness = *m_grid->variables().get_2d_scalar("land_ice_thickness");
    const IceModelVec2CellType &mask   = *m_grid->variables().get_2d_cell_type("mask");

    set_ocean_input_fields(ice_thickness, mask, m_cbasins, m_shelf_mask, model,
                           m_Toc_box0, m_Soc_box0);        // per shelf


    process_box1(ice_thickness, m_shelf_mask, m_ocean_box_mask, // inputs
                 m_Toc_box0, m_Soc_box0, model,                    // inputs
                 m_T_star, m_Toc, m_Soc, m_basalmeltrate_shelf, // outputs
                 m_overturning, m_T_pressure_melting);          // outputs

    process_other_boxes(ice_thickness, m_shelf_mask, model,        // inputs
                        m_ocean_box_mask, m_T_star, m_Toc,      // outputs
                        m_Soc, m_basalmeltrate_shelf, m_T_pressure_melting); // outputs

    //Assumes that mass flux is proportional to the shelf-base heat flux.
    process_missing_cells(model, m_shelf_mask, m_ocean_box_mask, ice_thickness, // inputs
                          m_Toc_box0, m_Soc_box0, // inputs
                          m_Toc, m_Soc, m_basalmeltrate_shelf, m_T_pressure_melting); // outputs

  }

  m_shelf_base_temperature->copy_from(m_T_pressure_melting); // in-situ freezing point at the ice shelf base
  m_shelf_base_mass_flux->copy_from(m_basalmeltrate_shelf);
  m_shelf_base_mass_flux->scale(model.ice_density());

  m_sea_level_elevation->set(0.0);
  m_melange_back_pressure_fraction->set(0.0);
}


//! Compute temperature and salinity input from ocean data by averaging.

//! We average the ocean data over the continental shelf reagion for each basin.
//! We use dummy ocean data if no such average can be calculated.


void Pico::compute_ocean_input_per_basin(const BoxModel &box_model,
                                         const IceModelVec2Int &basin_mask,
                                         const IceModelVec2Int &continental_shelf_mask,
                                         const IceModelVec2S &salinity_ocean,
                                         const IceModelVec2S &theta_ocean) {

  std::vector<double> count(m_n_basins, 0.0);
  std::vector<double> Tval(m_n_basins, 0.0);
  std::vector<double> Sval(m_n_basins, 0.0);

  IceModelVec::AccessList list{&theta_ocean, &salinity_ocean, &basin_mask, &continental_shelf_mask };

  // compute the sum for each basin for region that intersects with the
  // continental shelf area and is not covered by an ice shelf.
  // (continental shelf mask excludes ice shelf areas)
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (continental_shelf_mask(i, j) == INNER) {
      int basin_id = basin_mask.as_int(i, j);
      count[basin_id] += 1;
      Sval[basin_id] += salinity_ocean(i, j);
      Tval[basin_id] += theta_ocean(i, j);
    }
  }

  // Divide by number of grid cells if more than zero cells belong to the basin.
  // if no ocean_contshelf_mask values intersect with the basin, m_count is zero.
  // In such case, use dummy temperature and salinity. This could happen, for
  // example, if the ice shelf front advances beyond the continental shelf break.
  for (int basin_id = 0; basin_id < m_n_basins; basin_id++) {

    count[basin_id] = GlobalSum(m_grid->com, count[basin_id]);
    Sval[basin_id]  = GlobalSum(m_grid->com, Sval[basin_id]);
    Tval[basin_id]  = GlobalSum(m_grid->com, Tval[basin_id]);

    // if basin is not dummy basin 0 or there are no ocean cells in this basin to take the mean over.
    // FIXME: the following warning occurs once at initialization before input is available.
    // Please ignore this very first warning for now.
    if (basin_id > 0 && count[basin_id] == 0) {
      m_log->message(2, "PICO ocean WARNING: basin %d contains no cells with ocean data on continental shelf\n"
                        "(no values with ocean_contshelf_mask=2).\n"
                        "No mean salinity or temperature values are computed, instead using\n"
                        "the standard values T_dummy =%.3f, S_dummy=%.3f.\n"
                        "This might bias your basal melt rates, check your input data carefully.\n",
                     basin_id, box_model.T_dummy(), box_model.S_dummy());
      m_Toc_box0_vec[basin_id] = box_model.T_dummy();
      m_Soc_box0_vec[basin_id] = box_model.S_dummy();
    } else {
      Sval[basin_id] = Sval[basin_id] / count[basin_id];
      Tval[basin_id] = Tval[basin_id] / count[basin_id];

      m_Toc_box0_vec[basin_id] = Tval[basin_id];
      m_Soc_box0_vec[basin_id] = Sval[basin_id];
      m_log->message(5, "  %d: temp =%.3f, salinity=%.3f\n", basin_id, m_Toc_box0_vec[basin_id],
                     m_Soc_box0_vec[basin_id]);
    }
  }
}

//! Set ocean ocean input from box 0 as boundary condition for box 1.

//! Set ocean temperature and salinity (Toc_box0, Soc_box0)
//! from box 0 (in front of the ice shelf) as boundary condition for
//! box 1, which is the ocean box adjacent to the grounding line.
//! Toc_box0 and Soc_box0 were computed in function compute_ocean_input_per_basin.
//! We enforce that Toc_box0 is always at least the local pressure melting point.
void Pico::set_ocean_input_fields(const IceModelVec2S &ice_thickness,
                                  const IceModelVec2CellType &mask,
                                  const IceModelVec2Int &basin_mask,
                                  const IceModelVec2Int &shelf_mask,
                                  const BoxModel &box_model,
                                  IceModelVec2S &Toc_box0,
                                  IceModelVec2S &Soc_box0) {

  IceModelVec::AccessList list{ &ice_thickness, &basin_mask, &Soc_box0, &Toc_box0, &mask, &shelf_mask };

  std::vector<std::vector<int> > n_shelf_cells_per_basin(m_n_shelves,
                                                         std::vector<int>(m_n_basins, 0));
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

    if (mask.as_int(i, j) == MASK_FLOATING and shelf_mask.as_int(i, j) > 0) {
      // note: shelf_mask = 0 in lakes

      int s = shelf_mask.as_int(i, j);
      // weighted input depending on the number of shelf cells in each basin
      for (int b = 1; b < m_n_basins; b++) { //Note: b=0 yields nan
        Toc_box0(i, j) += m_Toc_box0_vec[b] * n_shelf_cells_per_basin[s][b] / (double)n_shelf_cells[s];
        Soc_box0(i, j) += m_Soc_box0_vec[b] * n_shelf_cells_per_basin[s][b] / (double)n_shelf_cells[s];
      }

      double theta_pm = box_model.theta_pm(Soc_box0(i, j), box_model.pressure(ice_thickness(i, j)));

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
    m_log->message(2,
                   "PICO ocean warning: temperature has been below pressure melting temperature in %d cases,\n"
                   "setting it to pressure melting temperature\n",
                   low_temperature_counter);
  }
}

//! compute the area of a box consisting of N cells dx by dy meters in size
double f_area(double N, double dx, double dy){
  //FIXME this assumes rectangular cell areas, adjust with real areas from projection
  return N * dx * dy;
}


//! Compute the basal melt for each ice shelf cell in box 1

//! Here are the core physical equations of the PICO model (for box1):
//! We here calculate basal melt rate, ambient ocean temperature and salinity
//! and overturning within box1. We calculate the average over the box 1 input for box 2.

void Pico::process_box1(const IceModelVec2S &ice_thickness,
                        const IceModelVec2Int &shelf_mask,
                        const IceModelVec2Int &box_mask,
                        const IceModelVec2S &Toc_box0,
                        const IceModelVec2S &Soc_box0,
                        const BoxModel &box_model,
                        IceModelVec2S &T_star,
                        IceModelVec2S &Toc,
                        IceModelVec2S &Soc,
                        IceModelVec2S &basal_melt_rate,
                        IceModelVec2S &overturning,
                        IceModelVec2S &T_pressure_melting) {

  std::vector<double> box1_area(m_n_shelves);
  const IceModelVec2S *cell_area = m_grid->variables().get_2d_scalar("cell_area");
  compute_box_area(1, shelf_mask, box_mask, *cell_area, box1_area);

  IceModelVec::AccessList list{&ice_thickness, &shelf_mask, &box_mask, &T_star,
      &Toc_box0, &Toc, &Soc_box0, &Soc, &overturning, &basal_melt_rate, &T_pressure_melting};

  int n_Toc_failures = 0;

  // basal melt rate, ambient temperature and salinity and overturning calculation
  // for each box1 grid cell.
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    int shelf_id = shelf_mask.as_int(i, j);

    if (box_mask.as_int(i, j) == 1 and shelf_id > 0.0) {

      // pressure in dbar, 1dbar = 10000 Pa = 1e4 kg m-1 s-2
      const double pressure = box_model.pressure(ice_thickness(i, j));

      T_star(i, j) = box_model.T_star(Soc_box0(i, j), Toc_box0(i, j), pressure);

      auto Toc_box1 = box_model.Toc_box1(box1_area[shelf_id],
                                         T_star(i, j), Soc_box0(i, j), Toc_box0(i, j));

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
      Soc(i, j) = box_model.Soc_box1(Toc_box0(i, j), Soc_box0(i, j), Toc(i, j)); // in psu

      overturning(i, j) = box_model.overturning(Soc_box0(i, j), Soc(i, j),
                                                Toc_box0(i, j), Toc(i, j));

      // main outputs
      basal_melt_rate(i, j) = box_model.melt_rate(box_model.theta_pm(Soc(i, j), pressure),
                                                  Toc(i, j));
      T_pressure_melting(i, j) = box_model.T_pm(Soc(i, j), pressure);

    } else {
      T_star(i, j)      = 0.0;    // in Kelvin
      Toc(i, j)         = 273.15; // in Kelvin
      Soc(i, j)         = 0.0;    // in psu
      overturning(i, j) = 0.0;

      basal_melt_rate(i, j)    = 0.0;
      T_pressure_melting(i, j) = 0.0;
    }
  }

  n_Toc_failures = GlobalSum(m_grid->com, n_Toc_failures);
  if (n_Toc_failures > 0) {
    m_log->message(2, "PICO ocean warning: square-root argument for temperature calculation "
                      "has been negative in %d cases!\n",
                   n_Toc_failures);
  }
}

//! Compute the basal melt for each ice shelf cell in boxes other than box1

//! Here are the core physical equations of the PICO model:
//! We here calculate basal melt rate, ambient ocean temperature and salinity.
//! Overturning is only calculated for box 1.
//! We calculate the average values over box i as input for box i+1.

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
 * @param[in] cell_area cell area, m^2
 * @param[out] result resulting box areas.
 *
 * Note: shelf and box indexes start from 1.
 */
void Pico::compute_box_area(int box_id,
                            const IceModelVec2Int &shelf_mask,
                            const IceModelVec2Int &box_mask,
                            const IceModelVec2S &cell_area,
                            std::vector<double> &result) {
  result.resize(m_n_shelves);

  IceModelVec::AccessList list{&shelf_mask, &box_mask, &cell_area};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    int shelf_id = shelf_mask.as_int(i, j);

    if (shelf_id > 0 and box_mask.as_int(i, j) == box_id) {
      result[shelf_id] += cell_area(i, j);
    }
  }

  // compute global sums
  for (int s = 1; s < m_n_shelves; ++s) {
    result[s] = GlobalSum(m_grid->com, result[s]);
  }
}

void Pico::process_other_boxes(const IceModelVec2S &ice_thickness,
                               const IceModelVec2Int &shelf_mask,
                               const BoxModel &box_model,
                               IceModelVec2Int &box_mask,
                               IceModelVec2S &T_star,
                               IceModelVec2S &Toc,
                               IceModelVec2S &Soc,
                               IceModelVec2S &basal_melt_rate,
                               IceModelVec2S &T_pressure_melting) {

  std::vector<double> overturning(m_n_shelves, 0.0);
  std::vector<double> salinity(m_n_shelves, 0.0);
  std::vector<double> temperature(m_n_shelves, 0.0);

  // get average overturning from box 1 that is used as input later
  compute_box_average(1, m_overturning, shelf_mask, box_mask, overturning);

  const IceModelVec2S &cell_area = *m_grid->variables().get_2d_scalar("cell_area");

  IceModelVec::AccessList list{
    &ice_thickness, &shelf_mask, &box_mask, &T_star, &Toc, &Soc,
      &basal_melt_rate, &T_pressure_melting, &cell_area};

  // Iterate over all boxes i for i > 1
  // box number = numberOfBoxes+1 is used as identifier for Beckmann Goose calculation
  // for cells with missing input and excluded in loop here, i.e. box <= m_n_boxes.
  for (int box = 2; box <= m_n_boxes; ++box) {

    compute_box_average(box - 1, Toc, shelf_mask, box_mask, temperature);
    compute_box_average(box - 1, Soc, shelf_mask, box_mask, salinity);

    std::vector<double> box_area;
    compute_box_area(box, shelf_mask, box_mask, cell_area, box_area);

    int countGl0 = 0;

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      int shelf_id = shelf_mask.as_int(i, j);

      if (box_mask.as_int(i, j) == box and shelf_id > 0) {

        // Get the input from previous box
        //
        // Note: overturning is only computed in box 1
        const double
          S_previous       = salinity[shelf_id],
          T_previous       = temperature[shelf_id],
          overturning_box1 = overturning[shelf_id];

        // if there are no boundary values from the box before
        if (S_previous == 0.0 or overturning_box1 == 0.0 or T_previous == 0.0) {

          // set mask to Beckmann Goose identifier, will be handled in calculate_basal_melt_missing_cells
          box_mask(i, j) = m_n_boxes + 1;

          // flag to print warning later
          countGl0 += 1;

        } else {

          double pressure = box_model.pressure(ice_thickness(i, j));

          // diagnostic outputs
          T_star(i, j) = box_model.T_star(S_previous, T_previous, pressure);
          Toc(i, j)    = box_model.Toc(box_area[shelf_id], T_previous, T_star(i, j), overturning_box1, S_previous);
          Soc(i, j)    = box_model.Soc(S_previous, T_previous, Toc(i, j));

          // main outputs: basal melt rate and temperature
          basal_melt_rate(i, j) = box_model.melt_rate(box_model.theta_pm(Soc(i, j), pressure), Toc(i,j));
          T_pressure_melting(i, j) = box_model.T_pm(Soc(i, j), pressure);
        }
      }
      // no else-case, since calculate_basal_melt_box1() and
      // calculate_basal_melt_missing_cells() cover all other cases and we would overwrite
      // those results here.
    }   // loop over grid points

    countGl0 = GlobalSum(m_grid->com, countGl0);
    if (countGl0 > 0) {
      m_log->message(2, "PICO ocean WARNING: box %d, no boundary data from previous box in %d case(s)!\n"
                     "switching to Beckmann Goose (2003) meltrate calculation\n",
                     box, countGl0);
    }

  } // loop over boxes

  // FIXME: we should not modify the box mask here
  box_mask.update_ghosts();
}


//! Compute the basal melt for ice shelf cells with missing input data

//! This covers cells that could not be related to ocean boxes
//! or where input data is missing.
//! Such boxes are identified with the ocean_box_mask value numberOfBoxes+1.
//! For those boxes use the [@ref BeckmannGoosse2003] meltrate parametrization, which
//! only depends on local ocean inputs.
//! We use the open ocean temperature and salinity as input here.
//! Also set basal melt rate to zero if shelf_id is zero, which is mainly
//! at the computational domain boundary.

void Pico::process_missing_cells(const BoxModel &box_model,
                                 const IceModelVec2Int &shelf_mask,
                                 const IceModelVec2Int &box_mask,
                                 const IceModelVec2S &ice_thickness,
                                 const IceModelVec2S &Toc_box0,
                                 const IceModelVec2S &Soc_box0,
                                 IceModelVec2S &Toc,
                                 IceModelVec2S &Soc,
                                 IceModelVec2S &basal_melt_rate,
                                 IceModelVec2S &T_pressure_melting) {

  IceModelVec::AccessList list{
    &ice_thickness, &shelf_mask, &box_mask, &Toc_box0, &Toc, &Soc_box0, &Soc,
    &basal_melt_rate, &T_pressure_melting
  };

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    int shelf_id = shelf_mask.as_int(i, j);

    // mainly at the boundary of computational domain,
    // or through erroneous basin mask
    if (shelf_id == 0) {
      basal_melt_rate(i, j) = 0.0;
    }

    // cell with missing data identifier numberOfBoxes+1, as set in routines before
    if (shelf_id > 0 and box_mask.as_int(i, j) == (m_n_boxes + 1)) {

      double pressure = box_model.pressure(ice_thickness(i, j));

      basal_melt_rate(i, j) = box_model.melt_rate_beckmann_goose(box_model.theta_pm(Soc_box0(i, j), pressure),
                                                                 Toc_box0(i, j));
      T_pressure_melting(i, j) = box_model.T_pm(Soc_box0(i, j), pressure);

      // diagnostic outputs
      Toc(i, j) = Toc_box0(i, j); // in Kelvin
      Soc(i, j) = Soc_box0(i, j); // in psu
    }
  }
}

// Write diagnostic variables to extra files if requested
std::map<std::string, Diagnostic::Ptr> Pico::diagnostics_impl() const {

  DiagnosticList result = {
    { "basins",                    Diagnostic::wrap(m_cbasins) },
    { "pico_overturning",          Diagnostic::wrap(m_overturning) },
    { "pico_salinity_box0",        Diagnostic::wrap(m_Soc_box0) },
    { "pico_temperature_box0",     Diagnostic::wrap(m_Toc_box0) },
    { "pico_ocean_box_mask",       Diagnostic::wrap(m_ocean_box_mask) },
    { "pico_shelf_mask",           Diagnostic::wrap(m_shelf_mask) },
    { "pico_bmelt_shelf",          Diagnostic::wrap(m_basalmeltrate_shelf) },
    { "pico_icerise_mask",         Diagnostic::wrap(m_icerise_mask) },
    { "pico_ocean_contshelf_mask", Diagnostic::wrap(m_ocean_contshelf_mask) },
    { "pico_ocean_mask",           Diagnostic::wrap(m_ocean_mask) },
    { "pico_lake_mask",            Diagnostic::wrap(m_lake_mask) },
    { "pico_dist_grounding_line",  Diagnostic::wrap(m_DistGL) },
    { "pico_dist_iceshelf_front",  Diagnostic::wrap(m_DistIF) },
    { "pico_salinity",             Diagnostic::wrap(m_Soc) },
    { "pico_temperature",          Diagnostic::wrap(m_Toc) },
    { "pico_T_star",               Diagnostic::wrap(m_T_star) },
    { "pico_T_pressure_melting",   Diagnostic::wrap(m_T_pressure_melting) },
  };

  return combine(result, OceanModel::diagnostics_impl());
}

} // end of namespace ocean
} // end of namespace pism
