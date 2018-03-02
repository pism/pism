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

Pico::Constants::Constants(const Config &config) {

  // standard value for Antarctic basin mask
  default_numberOfBasins = static_cast<int>(config.get_double("ocean.pico.number_of_basins"));

  // maximum number of boxes (applies for big ice shelves)
  default_numberOfBoxes = static_cast<int>(config.get_double("ocean.pico.number_of_boxes"));

  // threshold between deep ocean and continental shelf
  continental_shelf_depth = config.get_double("ocean.pico.continental_shelf_depth");

  // value for ocean temperature around Antarctica if no other data available (cold conditions)
  T_dummy = -1.5 + config.get_double("constants.fresh_water.melting_point_temperature");
  S_dummy = 34.7; // value for ocean salinity around Antarctica if no other data available (cold conditions)

  earth_grav = config.get_double("constants.standard_gravity");
  rhoi       = config.get_double("constants.ice.density");
  rhow       = config.get_double("constants.sea_water.density");
  rho_star   = 1033;        // kg/m^3
  nu         = rhoi / rhow; // no unit

  latentHeat = config.get_double("constants.fresh_water.latent_heat_of_fusion"); //Joule / kg
  c_p_ocean  = 3974.0;                 // J/(K*kg), specific heat capacity of ocean mixed layer
  lambda     = latentHeat / c_p_ocean; // °C, NOTE K vs °C

  // Valus for linearized potential freezing point (from Xylar Asay-Davis, should be in Asay-Davis et al 2016, but not correct in there )
  a = -0.0572;         // K/psu
  b = 0.0788 + 273.15; // K
  c = 7.77e-4;         // K/dbar

  // in-situ pressure melting point from Jenkins et al. 2010 paper
  as = -0.0573;         // K/psu
  bs = 0.0832 + 273.15; // K
  cs = 7.53e-4;         // K/dbar

  // in-situ pressure melting point from Olbers & Hellmer 2010 paper
  // as          = -0.057;       // K/psu
  // bs          = 0.0832 + 273.15;       // K
  // cs          = 7.64e-4;      // K/dbar

  alpha = 7.5e-5; // 1/K
  beta  = 7.7e-4; // 1/psu

  default_gamma_T = config.get_double("ocean.pico.heat_exchange_coefficent"); // m s-1, best fit value in paper
  default_overturning_coeff =
      config.get_double("ocean.pico.overturning_coefficent"); // m6 kg−1 s−1, best fit value in paper

  // for shelf cells where normal box model is not calculated,
  // used in calculate_basal_melt_missing_cells(), compare POConstantPIK
  // m/s, thermal exchange velocity for Beckmann-Goose parameterization
  // this is the same meltFactor as in POConstantPIK
  meltFactor = config.get_double("ocean.pik_melt_factor");
}

Pico::Pico(IceGrid::ConstPtr g) : PGivenClimate<CompleteOceanModel, CompleteOceanModel>(g, NULL) {

  m_option_prefix = "-ocean_pico";

  // will be de-allocated by the parent's destructor
  m_theta_ocean    = new IceModelVec2T;
  m_salinity_ocean = new IceModelVec2T;

  m_fields["theta_ocean"]    = m_theta_ocean;
  m_fields["salinity_ocean"] = m_salinity_ocean;

  process_options();

  m_exicerises_set = options::Bool("-exclude_icerises", "exclude ice rises in PICO");

  std::map<std::string, std::string> standard_names;
  set_vec_parameters(standard_names);

  m_Mx = m_grid->Mx(), m_My = m_grid->My(), m_dx = m_grid->dx(), m_dy = m_grid->dy();

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
  m_Soc.set_attrs("model_state", "ocean salinity field", "", "ocean salinity field"); //NOTE unit=psu

  // salinity input for box 1
  m_Soc_box0.create(m_grid, "pico_salinity_box0", WITHOUT_GHOSTS);
  m_Soc_box0.set_attrs("model_state", "ocean base salinity field", "", "ocean base salinity field"); //NOTE unit=psu

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
  m_numberOfBasins = 20;
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

  Constants cc(*m_config);
  initBasinsOptions(cc);

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


//! initialize PICO model variables, can be user-defined.

//! numberOfBasins: number of drainage basins for PICO model
//!                 FIXME: we should infer that from the read-in basin mask
//! numberOfBoxes: maximum number of ocean boxes for PICO model
//!                for smaller shelves, the model may use less.
//! gamma_T: turbulent heat exchange coefficient for ice-ocean boundary layer
//! overturning_coeff: coefficient that scales strength of overturning circulation
//! continental_shelf_depth: threshold for definition of continental shelf area
//!                          area shallower than threshold is used for ocean input

void Pico::initBasinsOptions(const Constants &cc) {

  m_log->message(5, "starting initBasinOptions\n");

  m_numberOfBasins = cc.default_numberOfBasins;
  m_numberOfBoxes  = cc.default_numberOfBoxes;

  m_Toc_box0_vec.resize(m_numberOfBasins);
  m_Soc_box0_vec.resize(m_numberOfBasins);

  counter_boxes.resize(m_numberOfShelves, std::vector<double>(2, 0));
  // FIXME: the three vectors below will later on be resized to numberOfShelves,
  // is the line here still necessary?
  m_mean_salinity_boundary_vector.resize(m_numberOfBasins);
  m_mean_temperature_boundary_vector.resize(m_numberOfBasins);
  m_mean_overturning_box1_vector.resize(m_numberOfBasins);

  m_gamma_T           = cc.default_gamma_T;
  m_overturning_coeff = cc.default_overturning_coeff;
  m_log->message(2, "  -Using %d drainage basins and values: \n"
                    "   gamma_T= %.2e, overturning_coeff = %.2e... \n",
                 m_numberOfBasins, m_gamma_T, m_overturning_coeff);

  m_continental_shelf_depth = cc.continental_shelf_depth;
  m_log->message(2, "  -Depth of continental shelf for computation of temperature and salinity input\n"
                    "   is set for whole domain to continental_shelf_depth=%.0f meter\n",
                 m_continental_shelf_depth);
}

void Pico::update_impl(double my_t, double my_dt) {

  // Make sure that sea water salinity and sea water potential
  // temperature fields are up to date:
  update_internal(my_t, my_dt);

  m_theta_ocean->average(m_t, m_dt);
  m_salinity_ocean->average(m_t, m_dt);

  Constants cc(*m_config);

  test();

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
  identify_ocean_box_mask(cc);


  // Physical part of PICO

  // prepare ocean input temperature and salinity
  compute_ocean_input_per_basin(cc); // per basin
  set_ocean_input_fields(cc);        // per shelf

  //basal melt rates underneath ice shelves
  calculate_basal_melt_box1(cc);
  calculate_basal_melt_other_boxes(cc);
  calculate_basal_melt_missing_cells(cc); //Assumes that mass flux is proportional to the shelf-base heat flux.

  m_shelf_base_temperature->copy_from(m_T_pressure_melting); // in-situ freezing point at the ice shelf base
  m_shelf_base_mass_flux->copy_from(m_basalmeltrate_shelf);
  m_shelf_base_mass_flux->scale(cc.rhoi);

  m_sea_level_elevation->set(0.0);
  m_melange_back_pressure_fraction->set(0.0);
}


//! Compute temperature and salinity input from ocean data by averaging.

//! We average over ocean_contshelf_mask for each basin.
//! We use dummy ocean data if no such average can be calculated.
//!

void Pico::compute_ocean_input_per_basin(const Constants &cc) {

  m_log->message(5, "starting compute_ocean_input_per_basin routine \n");

  std::vector<double> count(m_numberOfBasins, 0.0);
  std::vector<double> Tval(m_numberOfBasins, 0.0);
  std::vector<double> Sval(m_numberOfBasins, 0.0);

  IceModelVec::AccessList list{ m_theta_ocean, m_salinity_ocean, &m_cbasins, &m_ocean_contshelf_mask };

  // compute the sum for each basin
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m_ocean_contshelf_mask(i, j) == INNER) {
      int shelf_id = m_cbasins(i, j);
      count[shelf_id] += 1;
      Sval[shelf_id] += (*m_salinity_ocean)(i, j);
      Tval[shelf_id] += (*m_theta_ocean)(i, j);
    }
  }

  // Divide by number of grid cells if more than zero cells belong to the basin.
  // if no ocean_contshelf_mask values intersect with the basin, m_count is zero.
  // in such case, use dummy temperature and salinity. This could happen, for
  // example, if the ice shelf front advances beyond the continental shelf break.
  for (int shelf_id = 0; shelf_id < m_numberOfBasins; shelf_id++) {

    count[shelf_id] = GlobalSum(m_grid->com, count[shelf_id]);
    Sval[shelf_id]  = GlobalSum(m_grid->com, Sval[shelf_id]);
    Tval[shelf_id]  = GlobalSum(m_grid->com, Tval[shelf_id]);

    // if basin is not dummy basin 0 or there are no ocean cells in this basin to take the mean over.
    // FIXME: the following warning occurs once at initialization before input is available.
    // Please ignore this very first warning for now.
    if (shelf_id > 0 && count[shelf_id] == 0) {
      m_log->message(2, "PICO ocean WARNING: basin %d contains no cells with ocean data on continental shelf\n"
                        "(no values with ocean_contshelf_mask=2).\n"
                        "No mean salinity or temperature values are computed, instead using\n"
                        "the standard values T_dummy =%.3f, S_dummy=%.3f.\n"
                        "This might bias your basal melt rates, check your input data carefully.\n",
                     shelf_id, cc.T_dummy, cc.S_dummy);
      m_Toc_box0_vec[shelf_id] = cc.T_dummy;
      m_Soc_box0_vec[shelf_id] = cc.S_dummy;
    } else {
      Sval[shelf_id] = Sval[shelf_id] / count[shelf_id];
      Tval[shelf_id] = Tval[shelf_id] / count[shelf_id];

      m_Toc_box0_vec[shelf_id] = Tval[shelf_id];
      m_Soc_box0_vec[shelf_id] = Sval[shelf_id];
      m_log->message(5, "  %d: temp =%.3f, salinity=%.3f\n", shelf_id, m_Toc_box0_vec[shelf_id],
                     m_Soc_box0_vec[shelf_id]);
    }
  }
}


//! Set ocean ocean input from box 0 as boundary condition for box 1.

//! Set ocean temperature and salinity (Toc_box0, Soc_box0)
//! from box 0 (in front of the ice shelf) as boundary condition for
//! box 1, which is the ocean box adjacent to the grounding line.
//! Toc_box0 and Soc_box0 were computed in function compute_ocean_input_per_basin.
//! We enforce that Toc_box0 is always at least the local pressure melting point.

void Pico::set_ocean_input_fields(const Constants &cc) {

  m_log->message(5, "starting set_ocean_input_fields routine\n");

  const IceModelVec2S &ice_thickness = *m_grid->variables().get_2d_scalar("land_ice_thickness");
  const IceModelVec2CellType &mask   = *m_grid->variables().get_2d_cell_type("mask");

  IceModelVec::AccessList list{ &ice_thickness, &m_cbasins, &m_Soc_box0, &m_Toc_box0, &m_Toc, &mask, &m_shelf_mask };

  // compute for each shelf the number of cells  within each basin
  std::vector<std::vector<double> > lcounter_shelf_cells_in_basin(m_numberOfShelves,
                                                                  std::vector<double>(m_numberOfBasins));
  std::vector<std::vector<double> > counter_shelf_cells_in_basin(m_numberOfShelves,
                                                                 std::vector<double>(m_numberOfBasins));

  // compute the number of all shelf cells
  std::vector<double> lcounter_shelf_cells(m_numberOfShelves);
  std::vector<double> counter_shelf_cells(m_numberOfShelves);

  for (int shelf_id = 0; shelf_id < m_numberOfShelves; shelf_id++) {
    lcounter_shelf_cells[shelf_id] = 0;
    for (int basin_id = 0; basin_id < m_numberOfBasins; basin_id++) {
      lcounter_shelf_cells_in_basin[shelf_id][basin_id] = 0;
    }
  }

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    int shelf_id = m_shelf_mask(i, j);
    int basin_id = m_cbasins(i, j);
    lcounter_shelf_cells_in_basin[shelf_id][basin_id]++;
    lcounter_shelf_cells[shelf_id]++;
  }

  for (int shelf_id = 0; shelf_id < m_numberOfShelves; shelf_id++) {
    counter_shelf_cells[shelf_id] = GlobalSum(m_grid->com, lcounter_shelf_cells[shelf_id]);
    for (int basin_id = 0; basin_id < m_numberOfBasins; basin_id++) {
      counter_shelf_cells_in_basin[shelf_id][basin_id] =
          GlobalSum(m_grid->com, lcounter_shelf_cells_in_basin[shelf_id][basin_id]);
    }
  }


  // now set temp and salinity box 0:
  double counterTpmp = 0.0, lcounterTpmp = 0.0;

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // make sure all temperatures are zero at the beginning of each timestep
    m_Toc(i, j)      = 273.15; // in K
    m_Toc_box0(i, j) = 0.0;    // in K
    m_Soc_box0(i, j) = 0.0;    // in psu


    if (mask(i, j) == MASK_FLOATING && m_shelf_mask(i, j) > 0) { // shelf_mask = 0 in lakes

      int shelf_id = m_shelf_mask(i, j);
      // weighted input depending on the number of shelf cells in each basin
      for (int basin_id = 1; basin_id < m_numberOfBasins; basin_id++) { //Note: basin_id=0 yields nan
        m_Toc_box0(i, j) +=
            m_Toc_box0_vec[basin_id] * counter_shelf_cells_in_basin[shelf_id][basin_id] / counter_shelf_cells[shelf_id];
        m_Soc_box0(i, j) +=
            m_Soc_box0_vec[basin_id] * counter_shelf_cells_in_basin[shelf_id][basin_id] / counter_shelf_cells[shelf_id];
      }


      // pressure in dbar, 1dbar = 10000 Pa = 1e4 kg m-1 s-2
      const double pressure = cc.rhoi * cc.earth_grav * ice_thickness(i, j) * 1e-4;
      const double T_pmt = cc.a * m_Soc_box0(i, j) + cc.b - cc.c * pressure; // in Kelvin, here potential freezing point

      // temperature input for grounding line box should not be below pressure melting point
      if (m_Toc_box0(i, j) < T_pmt) {
        // Setting Toc_box0 a little higher than T_pmt ensures that later equations are well solvable.
        m_Toc_box0(i, j) = T_pmt + 0.001;
        lcounterTpmp += 1;
      }

    } // end if herefloating
  }

  counterTpmp = GlobalSum(m_grid->com, lcounterTpmp);
  if (counterTpmp > 0) {
    m_log->message(2, "PICO ocean warning: temperature has been below pressure melting temperature in %.0f cases,\n"
                      "setting it to pressure melting temperature\n",
                   counterTpmp);
  }
}


//! Compute the basal melt for each ice shelf cell in box 1

//! Here are the core physical equations of the PICO model (for box1):
//! We here calculate basal melt rate, ambient ocean temperature and salinity
//! and overturning within box1. We calculate the average over the box 1 input for box 2.

void Pico::calculate_basal_melt_box1(const Constants &cc) {

  m_log->message(5, "starting basal calculate_basal_melt_box1 routine\n");

  std::vector<double> lcounter_edge_of_box1_vector(m_numberOfShelves);
  std::vector<double> lmean_salinity_box1_vector(m_numberOfShelves);
  std::vector<double> lmean_temperature_box1_vector(m_numberOfShelves);
  std::vector<double> lmean_meltrate_box1_vector(m_numberOfShelves);
  std::vector<double> lmean_overturning_box1_vector(m_numberOfShelves);

  for (int shelf_id = 0; shelf_id < m_numberOfShelves; shelf_id++) {
    lcounter_edge_of_box1_vector[shelf_id]  = 0.0;
    lmean_salinity_box1_vector[shelf_id]    = 0.0;
    lmean_temperature_box1_vector[shelf_id] = 0.0;
    lmean_meltrate_box1_vector[shelf_id]    = 0.0;
    lmean_overturning_box1_vector[shelf_id] = 0.0;
  }

  m_mean_salinity_boundary_vector.resize(m_numberOfShelves);
  m_mean_temperature_boundary_vector.resize(m_numberOfShelves);
  m_mean_overturning_box1_vector.resize(m_numberOfShelves);


  const IceModelVec2S &ice_thickness = *m_grid->variables().get_2d_scalar("land_ice_thickness");

  IceModelVec::AccessList list{
    &ice_thickness, &m_shelf_mask, &m_ocean_box_mask, &m_T_star, &m_Toc_box0, &m_Toc, &m_Soc_box0, &m_Soc,
    &m_overturning, &m_basalmeltrate_shelf, &m_T_pressure_melting
  };

  double countHelpterm = 0, lcountHelpterm = 0;

  m_ocean_box_mask.update_ghosts();


  // basal melt rate, ambient temperature and salinity and overturning calculation
  // for each box1 grid cell.
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    int shelf_id = m_shelf_mask(i, j);

    // Make sure everything is at default values at the beginning of each timestep
    m_T_star(i, j) = 0.0;    // in Kelvin
    m_Toc(i, j)    = 273.15; // in Kelvin
    m_Soc(i, j)    = 0.0;    // in psu

    m_basalmeltrate_shelf(i, j) = 0.0;
    m_overturning(i, j)         = 0.0;
    m_T_pressure_melting(i, j)  = 0.0;

    if ((m_ocean_box_mask(i, j) == 1) && (shelf_id > 0.0)) {

      // pressure in dbar, 1dbar = 10000 Pa = 1e4 kg m-1 s-2
      const double pressure = cc.rhoi * cc.earth_grav * ice_thickness(i, j) * 1e-4;
      m_T_star(i, j) = cc.a * m_Soc_box0(i, j) + cc.b - cc.c * pressure - m_Toc_box0(i, j); // in Kelvin

      //FIXME this assumes rectangular cell areas, adjust with real areas from projection
      double area_box1 = (counter_boxes[shelf_id][1] * m_dx * m_dy);

      double g1 = area_box1 * m_gamma_T;
      double s1 = m_Soc_box0(i, j) / (cc.nu * cc.lambda);

      // These are the coefficients for solving the quadratic temperature equation
      // trough the p-q formula.
      double p_coeff = g1 / (m_overturning_coeff * cc.rho_star * (cc.beta * s1 - cc.alpha)); // in 1 / (1/K) = K
      double q_coeff =
          (g1 * m_T_star(i, j)) / (m_overturning_coeff * cc.rho_star * (cc.beta * s1 - cc.alpha)); // in K / (1/K) = K^2

      // This can only happen if T_star > 0.25*p_coeff, in particular T_star > 0
      // which can only happen for values of Toc_box0 close to the local pressure melting point
      if ((0.25 * PetscSqr(p_coeff) - q_coeff) < 0.0) {

        m_log->message(5, "PICO ocean WARNING: negative square root argument at %d, %d\n"
                          "probably because of positive T_star=%f \n"
                          "Not aborting, but setting square root to 0... \n",
                       i, j, m_T_star(i, j));

        q_coeff = 0.25 * PetscSqr(p_coeff);
        lcountHelpterm += 1;
      }

      // temperature for box 1, p-q formula
      m_Toc(i, j) = m_Toc_box0(i, j) - (-0.5 * p_coeff + sqrt(0.25 * PetscSqr(p_coeff) - q_coeff)); // in Kelvin
      // salinity for box 1
      m_Soc(i, j) =
          m_Soc_box0(i, j) - (m_Soc_box0(i, j) / (cc.nu * cc.lambda)) * (m_Toc_box0(i, j) - m_Toc(i, j)); // in psu

      // potential pressure melting pointneeded to calculate thermal driving
      // using coefficients for potential temperature
      double potential_pressure_melting_point = cc.a * m_Soc(i, j) + cc.b - cc.c * pressure;

      // basal melt rate for box 1
      m_basalmeltrate_shelf(i, j) =
          (-m_gamma_T / (cc.nu * cc.lambda)) * (potential_pressure_melting_point - m_Toc(i, j)); // in m/s

      m_overturning(i, j) =
          m_overturning_coeff * cc.rho_star *
          (cc.beta * (m_Soc_box0(i, j) - m_Soc(i, j)) - cc.alpha * (m_Toc_box0(i, j) - m_Toc(i, j))); // in m^3/s

      // average the temperature, salinity and overturning over the entire box1
      // this is used as input for box 2
      // (here we sum up)
      lcounter_edge_of_box1_vector[shelf_id]++;
      lmean_salinity_box1_vector[shelf_id] += m_Soc(i, j);
      lmean_temperature_box1_vector[shelf_id] += m_Toc(i, j); // in Kelvin
      lmean_meltrate_box1_vector[shelf_id] += m_basalmeltrate_shelf(i, j);
      lmean_overturning_box1_vector[shelf_id] += m_overturning(i, j);

      // in situ pressure melting point
      m_T_pressure_melting(i, j) = cc.as * m_Soc(i, j) + cc.bs - cc.cs * pressure; //  in Kelvin

    } else { // i.e., not GL_box
      m_basalmeltrate_shelf(i, j) = 0.0;
    }
  }


  // average the temperature, salinity and overturning over box1
  // (here we divide)
  for (int shelf_id = 0; shelf_id < m_numberOfShelves; shelf_id++) {
    double counter_edge_of_box1_vector = 0.0;

    counter_edge_of_box1_vector                  = GlobalSum(m_grid->com, lcounter_edge_of_box1_vector[shelf_id]);
    m_mean_salinity_boundary_vector[shelf_id]    = GlobalSum(m_grid->com, lmean_salinity_box1_vector[shelf_id]);
    m_mean_temperature_boundary_vector[shelf_id] = GlobalSum(m_grid->com, lmean_temperature_box1_vector[shelf_id]);
    m_mean_overturning_box1_vector[shelf_id]     = GlobalSum(m_grid->com, lmean_overturning_box1_vector[shelf_id]);

    if (counter_edge_of_box1_vector > 0.0) {
      m_mean_salinity_boundary_vector[shelf_id] =
          m_mean_salinity_boundary_vector[shelf_id] / counter_edge_of_box1_vector;
      m_mean_temperature_boundary_vector[shelf_id] =
          m_mean_temperature_boundary_vector[shelf_id] / counter_edge_of_box1_vector;
      m_mean_overturning_box1_vector[shelf_id] = m_mean_overturning_box1_vector[shelf_id] / counter_edge_of_box1_vector;
    } else {
      // This means that there is no cells in box 1
      m_mean_salinity_boundary_vector[shelf_id]    = 0.0;
      m_mean_temperature_boundary_vector[shelf_id] = 0.0;
      m_mean_overturning_box1_vector[shelf_id]     = 0.0;
    }

    // print input values for box 2
    m_log->message(5, "  %d: cnt=%.0f, sal=%.3f, temp=%.3f, over=%.1e \n", shelf_id, counter_edge_of_box1_vector,
                   m_mean_salinity_boundary_vector[shelf_id], m_mean_temperature_boundary_vector[shelf_id],
                   m_mean_overturning_box1_vector[shelf_id]);
  }

  countHelpterm = GlobalSum(m_grid->com, lcountHelpterm);
  if (countHelpterm > 0) {
    m_log->message(2, "PICO ocean warning: square-root argument for temperature calculation "
                      "has been negative in %.0f cases!\n",
                   countHelpterm);
  }
}

//! Compute the basal melt for each ice shelf cell in boxes other than box1

//! Here are the core physical equations of the PICO model:
//! We here calculate basal melt rate, ambient ocean temperature and salinity.
//! Overturning is only calculated for box 1.
//! We calculate the average values over box i as input for box i+1.

void Pico::calculate_basal_melt_other_boxes(const Constants &cc) {

  m_log->message(5, "starting calculate_basal_melt_other_boxes routine\n");

  int nBoxes = static_cast<int>(round(m_numberOfBoxes + 1));

  const IceModelVec2S &ice_thickness = *m_grid->variables().get_2d_scalar("land_ice_thickness");

  IceModelVec::AccessList list{
    &ice_thickness, &m_shelf_mask, &m_ocean_box_mask, &m_T_star, &m_Toc_box0, &m_Toc, &m_Soc_box0, &m_Soc,
    &m_overturning, &m_basalmeltrate_shelf, &m_T_pressure_melting
  };

  m_ocean_box_mask.update_ghosts();

  // Iterate over all Boxes i for i > 1
  // box number = numberOfBoxes+1 is used as identifier for Beckmann Goose calculation
  // for cells with missing input and excluded in loop here, i.e. boxi <nBoxes.
  for (int boxi = 2; boxi < nBoxes; ++boxi) {

    m_log->message(5, "computing basal melt rate, temperature and salinity for box i = %d \n", boxi);

    double countGl0 = 0, lcountGl0 = 0;

    // averages over the current box, input for the subsequent box
    std::vector<double> lmean_salinity_boxi_vector(m_numberOfShelves);    // in psu
    std::vector<double> lmean_temperature_boxi_vector(m_numberOfShelves); // in Kelvin
    std::vector<double> lcounter_edge_of_boxi_vector(m_numberOfShelves);

    for (int shelf_id = 0; shelf_id < m_numberOfShelves; shelf_id++) {
      lcounter_edge_of_boxi_vector[shelf_id]  = 0.0;
      lmean_salinity_boxi_vector[shelf_id]    = 0.0;
      lmean_temperature_boxi_vector[shelf_id] = 0.0;
    }

    // for box i compute the melt rates.

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      int shelf_id = m_shelf_mask(i, j);

      if (m_ocean_box_mask(i, j) == boxi && shelf_id > 0.0) {

        double area_boxi, mean_salinity_in_boundary, mean_temperature_in_boundary, mean_overturning_in_box1;

        // get the input from previous box (is from box 1 if boxi=2)
        // overturning is only solved in box 1 and same for other boxes
        // temperature and salinity boundary values will be updated at the
        // end of this routine
        mean_salinity_in_boundary    = m_mean_salinity_boundary_vector[shelf_id];
        mean_temperature_in_boundary = m_mean_temperature_boundary_vector[shelf_id]; // Kelvin
        mean_overturning_in_box1     = m_mean_overturning_box1_vector[shelf_id];

        // if there are no boundary values from the box before
        if (mean_salinity_in_boundary == 0 || mean_overturning_in_box1 == 0 || mean_temperature_in_boundary == 0) {

          // set mask to Beckmann Goose identifier, will be handled in calculate_basal_melt_missing_cells
          m_ocean_box_mask(i, j) = m_numberOfBoxes + 1;
          // flag to print warning later
          lcountGl0 += 1;

        } else {

          // solve the SIMPLE physical model equations for boxes with boxi > 1
          // pressure in dbar, 1dbar = 10000 Pa = 1e4 kg m-1 s-2
          const double pressure = cc.rhoi * cc.earth_grav * ice_thickness(i, j) * 1e-4;
          m_T_star(i, j) =
              cc.a * mean_salinity_in_boundary + cc.b - cc.c * pressure - mean_temperature_in_boundary; // in Kelvin

          //FIXME this assumes rectangular cell areas, adjust with real areas from projection
          area_boxi = (counter_boxes[shelf_id][boxi] * m_dx * m_dy);

          // compute melt rates
          double g1 = area_boxi * m_gamma_T;
          double g2 = g1 / (cc.nu * cc.lambda);

          // temperature for Box i > 1
          m_Toc(i, j) =
              mean_temperature_in_boundary +
              g1 * m_T_star(i, j) / (mean_overturning_in_box1 + g1 - g2 * cc.a * mean_salinity_in_boundary); // K

          // salinity for Box i > 1
          m_Soc(i, j) =
              mean_salinity_in_boundary -
              mean_salinity_in_boundary * (mean_temperature_in_boundary - m_Toc(i, j)) / (cc.nu * cc.lambda); // psu

          // potential pressure melting pointneeded to calculate thermal driving
          // using coefficients for potential temperature
          double potential_pressure_melting_point = cc.a * m_Soc(i, j) + cc.b - cc.c * pressure;

          // basal melt rate for Box i > 1
          m_basalmeltrate_shelf(i, j) =
              (-m_gamma_T / (cc.nu * cc.lambda)) * (potential_pressure_melting_point - m_Toc(i, j)); // in m/s

          // in situ pressure melting point in Kelvin
          m_T_pressure_melting(i, j) = cc.as * m_Soc(i, j) + cc.bs - cc.cs * pressure;

          // average the temperature, salinity over the entire box i
          // this is used as input for box i+1
          // (here we sum up)
          lcounter_edge_of_boxi_vector[shelf_id]++;
          lmean_salinity_boxi_vector[shelf_id] += m_Soc(i, j);
          lmean_temperature_boxi_vector[shelf_id] += m_Toc(i, j);
        }
      } // no else-case, since  calculate_basal_melt_box1() and calculate_basal_melt_missing_cells() cover all other cases and we would overwrite those results here.
    }

    // average the temperature and salinity over box i
    // (here we divide)
    for (int shelf_id = 0; shelf_id < m_numberOfShelves; shelf_id++) {
      // overturning should not be changed, fixed from box 1
      double counter_edge_of_boxi_vector        = 0.0;
      counter_edge_of_boxi_vector               = GlobalSum(m_grid->com, lcounter_edge_of_boxi_vector[shelf_id]);
      m_mean_salinity_boundary_vector[shelf_id] = GlobalSum(m_grid->com, lmean_salinity_boxi_vector[shelf_id]);
      m_mean_temperature_boundary_vector[shelf_id] =
          GlobalSum(m_grid->com, lmean_temperature_boxi_vector[shelf_id]); // in Kelvin

      if (counter_edge_of_boxi_vector > 0.0) {
        m_mean_salinity_boundary_vector[shelf_id] =
            m_mean_salinity_boundary_vector[shelf_id] / counter_edge_of_boxi_vector;
        m_mean_temperature_boundary_vector[shelf_id] =
            m_mean_temperature_boundary_vector[shelf_id] / counter_edge_of_boxi_vector; // in Kelvin
      } else {
        // This means that there is no cell in box i
        m_mean_salinity_boundary_vector[shelf_id]    = 0.0;
        m_mean_temperature_boundary_vector[shelf_id] = 0.0;
      }

      m_log->message(5, "  %d: cnt=%.0f, sal=%.3f, temp=%.3f, over=%.1e \n", shelf_id, counter_edge_of_boxi_vector,
                     m_mean_salinity_boundary_vector[shelf_id], m_mean_temperature_boundary_vector[shelf_id],
                     m_mean_overturning_box1_vector[shelf_id]);
    } // shelves

    countGl0 = GlobalSum(m_grid->com, lcountGl0);
    if (countGl0 > 0) {
      m_log->message(2, "PICO ocean WARNING: box %d, no boundary data from previous box in %.0f case(s)!\n"
                        "switching to Beckmann Goose (2003) meltrate calculation\n",
                     boxi, countGl0);
    }

  } // boxi
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

void Pico::calculate_basal_melt_missing_cells(const Constants &cc) {

  m_log->message(5, "starting calculate_basal_melt_missing_cells routine\n");

  const IceModelVec2S &ice_thickness = *m_grid->variables().get_2d_scalar("land_ice_thickness");

  IceModelVec::AccessList list{
    &ice_thickness, &m_shelf_mask, &m_ocean_box_mask, &m_Toc_box0, &m_Toc, &m_Soc_box0, &m_Soc,
    &m_overturning, &m_basalmeltrate_shelf, &m_T_pressure_melting
  };

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    int shelf_id = m_shelf_mask(i, j);

    // mainly at the boundary of computational domain,
    // or through erroneous basin mask
    if (shelf_id == 0) {
      m_basalmeltrate_shelf(i, j) = 0.0;
    }

    // cell with missing data identifier numberOfBoxes+1, as set in routines before
    if ((shelf_id > 0) && (m_ocean_box_mask(i, j) == (m_numberOfBoxes + 1))) {

      m_Toc(i, j) = m_Toc_box0(i, j); // in Kelvin
      m_Soc(i, j) = m_Soc_box0(i, j); // in psu

      // in dbar, 1dbar = 10000 Pa = 1e4 kg m-1 s-2
      const double pressure = cc.rhoi * cc.earth_grav * ice_thickness(i, j) * 1e-4;

      // potential pressure melting point needed to calculate thermal driving
      // using coefficients for potential temperature
      // these are different to the ones used in POConstantPIK
      double potential_pressure_melting_point = cc.a * m_Soc(i, j) + cc.b - cc.c * pressure;

      double heatflux = cc.meltFactor * cc.rhow * cc.c_p_ocean * cc.default_gamma_T *
                        (m_Toc(i, j) - potential_pressure_melting_point); // in W/m^2

      m_basalmeltrate_shelf(i, j) = heatflux / (cc.latentHeat * cc.rhoi); // in m s-1

      // in situ pressure melting point in Kelvin
      // this will be the temperature boundary condition at the ice at the shelf base
      m_T_pressure_melting(i, j) = cc.as * m_Soc(i, j) + cc.bs - cc.cs * pressure;
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