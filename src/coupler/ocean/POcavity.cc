 // Copyright (C) 2012-2016 Ricarda Winkelmann, Ronja Reese, Torsten Albrecht
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

#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>

#include "POcavity.hh"
#include "base/util/IceGrid.hh"
#include "base/util/PISMVars.hh"
#include "base/util/iceModelVec.hh"
#include "base/util/Mask.hh"
#include "base/util/PISMConfigInterface.hh"

namespace pism {
namespace ocean {

Cavity::Constants::Constants(const Config &config) {

  default_numberOfBasins = 20; // standard value for Antarctic basin mask
  default_numberOfBoxes = 5; // maximum number of boxes (applies for big ice shelves)

  continental_shelf_depth = -800; // threshold between deep ocean and continental shelf

  T_dummy = -1.5 + 273.15; // value for ocean temperature around Antarctica if no other data available FIXME Check
  S_dummy = 34.5; // value for ocean salinity around Antarctica if no other data available FIXME Check

  earth_grav = config.get_double("constants.standard_gravity");
  rhoi       = config.get_double("constants.ice.density");
  rhow       = config.get_double("constants.sea_water.density");
  rho_star   = 1033;                  // kg/m^3
  nu         = rhoi / rho_star;       // no unit

  latentHeat = config.get_double("constants.fresh_water.latent_heat_of_fusion"); //Joule / kg
  c_p_ocean  = 3974.0;       // J/(K*kg), specific heat capacity of ocean mixed layer
  lambda     = latentHeat / c_p_ocean;   // °C, NOTE K vs °C

  // Valus for linearized potential freezing point (from Xylar Asay-Davis, should be in Asay-Davis et al 2016, but not correct in there )
  a          = -0.0572;       // K/psu
  b          = 0.0788 + 273.15;       // K
  c          = 7.77e-4;      // K/dbar

  // in-situ pressure melting point from Jenkins et al. 2010 paper
  as          = -0.0573;       // K/psu
  bs          = 0.0832 + 273.15;       // K
  cs          = 7.53e-4;      // K/dbar

  // in-situ pressure melting point from Olbers & Hellmer 2010 paper
  // as          = -0.057;       // K/psu
  // bs          = 0.0832 + 273.15;       // K
  // cs          = 7.64e-4;      // K/dbar

  alpha      = 7.5e-5;       // 1/K
  beta       = 7.7e-4;       // 1/psu

  default_gamma_T    = 2e-5;        // m/s FIXME check!
  default_overturning_coeff    = 1e6;         // kg−1 s−1 FIXME check!

  // for shelf cells where normal box model is not calculated,
  // used in calculate_basal_melt_missing_cells(), compare POConstantPIK
  // m/s, thermal exchange velocity for Beckmann-Goose parameterization
  // this is a different meltFactor as in POConstantPIK
  meltFactor   = 0.01;

}


const int Cavity::maskfloating = MASK_FLOATING;
const int Cavity::maskocean    = MASK_ICE_FREE_OCEAN;
const int Cavity::maskgrounded = MASK_GROUNDED;

// used in IdentifyMask
const int Cavity::imask_inner        = 2;
const int Cavity::imask_outer        = 0;
const int Cavity::imask_exclude      = 1;
const int Cavity::imask_unidentified = -1;



Cavity::Cavity(IceGrid::ConstPtr g)
  : PGivenClimate<OceanModifier,OceanModel>(g, NULL) {

  m_option_prefix   = "-ocean_cavity";

  // will be de-allocated by the parent's destructor
  m_theta_ocean    = new IceModelVec2T;
  m_salinity_ocean = new IceModelVec2T;

  m_fields["theta_ocean"]     = m_theta_ocean;
  m_fields["salinity_ocean"]  = m_salinity_ocean;

  process_options();

  exicerises_set = options::Bool("-exclude_icerises", "exclude ice rises in ocean cavity model"); // FIXME set always?

  std::map<std::string, std::string> standard_names;
  set_vec_parameters(standard_names);

  Mx = m_grid->Mx(),
  My = m_grid->My(),
  dx = m_grid->dx(),
  dy = m_grid->dy();

  m_theta_ocean->create(m_grid, "theta_ocean");
  m_theta_ocean->set_attrs("climate_forcing",
                           "absolute potential temperature of the adjacent ocean",
                           "Kelvin", "");

  m_salinity_ocean->create(m_grid, "salinity_ocean");
  m_salinity_ocean->set_attrs("climate_forcing",
                                   "salinity of the adjacent ocean",
                                   "g/kg", "");

  m_shelfbtemp.create(m_grid, "shelfbtemp", WITHOUT_GHOSTS);
  m_shelfbtemp.set_attrs("climate_forcing",
                              "absolute temperature at ice shelf base",
                              "Kelvin", "");
  m_variables.push_back(&m_shelfbtemp);

  m_shelfbmassflux.create(m_grid, "shelfbmassflux", WITHOUT_GHOSTS);
  m_shelfbmassflux.set_attrs("climate_forcing",
                                  "ice mass flux from ice shelf base (positive flux is loss from ice shelf)",
                                  "kg m-2 s-1", "");
  m_shelfbmassflux.metadata().set_string("glaciological_units", "kg m-2 year-1");
  m_variables.push_back(&m_shelfbmassflux);


  cbasins.create(m_grid, "basins", WITH_GHOSTS);
  cbasins.set_attrs("climate_forcing","mask determines basins for ocean cavity model",
                    "", "");
  m_variables.push_back(&cbasins);

  // mask to identify the ocean boxes
  ocean_box_mask.create(m_grid, "ocean_box_mask", WITH_GHOSTS);
  ocean_box_mask.set_attrs("model_state", "mask displaying ocean box model grid","", "");
  m_variables.push_back(&ocean_box_mask);

  // mask to identify the ice rises
  icerise_mask.create(m_grid, "icerise_mask", WITH_GHOSTS);
  icerise_mask.set_attrs("model_state", "mask displaying ice rises","", "");
  m_variables.push_back(&icerise_mask);

  // mask displaying continental shelf - region where mean salinity and ocean temperature is calculated
  ocean_contshelf_mask.create(m_grid, "ocean_contshelf_mask", WITH_GHOSTS);
  ocean_contshelf_mask.set_attrs("model_state", "mask displaying ocean region for parameter input","", "");
  m_variables.push_back(&ocean_contshelf_mask);

  // mask displaying open ocean - ice-free regions below sea-level except 'holes' in ice shelves
  ocean_mask.create(m_grid, "ocean_mask", WITH_GHOSTS);
  ocean_mask.set_attrs("model_state", "mask displaying open ocean","", "");
  m_variables.push_back(&ocean_mask);

  // mask with distance (in boxes) to grounding line
  DistGL.create(m_grid, "DistGL", WITH_GHOSTS);
  DistGL.set_attrs("model_state", "mask displaying distance to grounding line","", "");
  m_variables.push_back(&DistGL);

  // mask with distance (in boxes) to ice front
  DistIF.create(m_grid, "DistIF", WITH_GHOSTS);
  DistIF.set_attrs("model_state", "mask displaying distance to ice shelf calving front","", "");
  m_variables.push_back(&DistIF);

  // computed salinity in ocean boxes
  Soc.create(m_grid, "Soc", WITHOUT_GHOSTS);
  Soc.set_attrs("model_state", "ocean salinity field","", "ocean salinity field");  //NOTE unit=psu
  m_variables.push_back(&Soc);

  // salinity input for box 1
  Soc_box0.create(m_grid, "Soc_box0", WITHOUT_GHOSTS);
  Soc_box0.set_attrs("model_state", "ocean base salinity field","", "ocean base salinity field");  //NOTE unit=psu
  m_variables.push_back(&Soc_box0);

  // computed temperature in ocean boxes
  Toc.create(m_grid, "Toc", WITHOUT_GHOSTS);
  Toc.set_attrs("model_state", "ocean temperature field","K", "ocean temperature field");
  m_variables.push_back(&Toc);

  // temperature input for box 1
  Toc_box0.create(m_grid, "Toc_box0", WITHOUT_GHOSTS);
  Toc_box0.set_attrs("model_state", "ocean base temperature","K", "ocean base temperature");
  m_variables.push_back(&Toc_box0);

  // in ocean box i: T_star = aS_{i-1} + b -c p_i - T_{i-1} with T_{-1} = Toc_box0 and S_{-1}=Soc_box0
  // FIXME convert to internal field
  T_star.create(m_grid, "T_star", WITHOUT_GHOSTS);
  T_star.set_attrs("model_state", "T_star field","degree C", "T_star field");
  m_variables.push_back(&T_star);

  overturning.create(m_grid, "overturning", WITHOUT_GHOSTS);
  overturning.set_attrs("model_state", "cavity overturning","m^3 s-1", "cavity overturning"); // no CF standard_name?
  m_variables.push_back(&overturning);

  basalmeltrate_shelf.create(m_grid, "basalmeltrate_shelf", WITHOUT_GHOSTS);
  basalmeltrate_shelf.set_attrs("model_state", "SIMPEL sub-shelf melt rate", "m/s",
                                "SIMPEL sub-shelf melt rate");
  //FIXME unit in field is kg m-2 a-1, but the written unit is m per a
  basalmeltrate_shelf.metadata().set_string("glaciological_units", "m year-1");
  m_variables.push_back(&basalmeltrate_shelf);

  // TODO: this may be initialized to NA, it should only have valid values below ice shelves.
  T_pressure_melting.create(m_grid, "T_pressure_melting", WITHOUT_GHOSTS);
  T_pressure_melting.set_attrs("model_state", "pressure melting temperature at ice shelf base",
                        "Kelvin", "pressure melting temperature at ice shelf base"); // no CF standard_name? // This is the in-situ pressure melting point
  m_variables.push_back(&T_pressure_melting);


  // Initialize this early so that we can check the validity of the "basins" mask read from a file
  // in Cavity::init_impl(). This number is hard-wired, so I don't think it matters that it did not
  // come from Cavity::Constants.
  numberOfBasins = 20;
}

Cavity::~Cavity() {
  // empty
}

void Cavity::init_impl() {

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  m_log->message(2, "* Initializing the Potsdam Cavity Model for the ocean ...\n");

  m_theta_ocean->init(m_filename, m_bc_period, m_bc_reference_time);
  m_salinity_ocean->init(m_filename, m_bc_period, m_bc_reference_time);

  cbasins.regrid(m_filename, CRITICAL);

  m_log->message(4, "SIMPEL basin min=%f,max=%f\n",cbasins.min(),cbasins.max());

  Constants cc(*m_config);
  initBasinsOptions(cc);

  round_basins();

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
  if (m_theta_ocean->get_n_records() == 1 &&
      m_salinity_ocean->get_n_records() == 1) {
        update(m_grid->ctx()->time()->current(), 0); // dt is irrelevant
  }

}

void Cavity::shelf_base_temperature_impl(IceModelVec2S &result) const {
  result.copy_from(m_shelfbtemp);
}

void Cavity::shelf_base_mass_flux_impl(IceModelVec2S &result) const {
  result.copy_from(m_shelfbmassflux);
}

void Cavity::sea_level_elevation_impl(double &result) const {
  result = m_sea_level;
}

void Cavity::melange_back_pressure_fraction_impl(IceModelVec2S &result) const {
  result.set(0.0);
}

//! initialize SIMPEL model variables, can be user-defined.

//! numberOfBasins: number of drainage basins for SIMPEL model
//!                 FIXME: we should infer that from the read-in basin mask
//! numberOfBoxes: maximum number of ocean boxes for SIMPEL model
//!                for smaller shelves, the model may use less.
//! gamma_T: turbulent heat exchange coefficient for ice-ocean boundary layer
//! overturning_coeff: coefficient that scales strength of overturning circulation
//! continental_shelf_depth: threshold for definition of continental shelf area
//!                          area shallower than threshold is used for ocean input

void Cavity::initBasinsOptions(const Constants &cc) {

  m_log->message(5, "starting initBasinOptions\n");

  numberOfBasins = cc.default_numberOfBasins;
  numberOfBasins = options::Integer("-number_of_basins",
                                    "number of drainage basins for SIMPEL model",
                                    numberOfBasins);
  numberOfBoxes = cc.default_numberOfBoxes;
  numberOfBoxes = options::Integer("-number_of_boxes",
                                    "number of ocean boxes for SIMPEL model",
                                    numberOfBoxes);

  Toc_box0_vec.resize(numberOfBasins);
  Soc_box0_vec.resize(numberOfBasins);

  counter_boxes.resize(numberOfBasins, std::vector<double>(2,0));
  mean_salinity_boundary_vector.resize(numberOfBasins);
  mean_temperature_boundary_vector.resize(numberOfBasins);
  mean_overturning_box1_vector.resize(numberOfBasins);

  gamma_T = cc.default_gamma_T;
  gamma_T = options::Real("-gamma_T","gamma_T for ocean cavity model",gamma_T); // meter per second

  overturning_coeff = cc.default_overturning_coeff;
  overturning_coeff = options::Real("-overturning_coeff",
                                    "overturning_coeff for ocean cavity model",overturning_coeff);

  m_log->message(2, "     Using %d drainage basins and values: \n"
                                "     gamma_T= %.2e, overturning_coeff = %.2e... \n"
                                 , numberOfBasins, gamma_T, overturning_coeff);

  continental_shelf_depth = cc.continental_shelf_depth;
  options::Real cont_shelf_depth("-continental_shelf_depth",
                                 "continental shelf depth for ocean cavity model",
                                 (double)continental_shelf_depth);

  if (cont_shelf_depth.is_set()) {
    m_log->message(2,
    "  Depth of continental shelf for computation of temperature and salinity input\n"
    "  is set for whole domain to continental_shelf_depth=%.0f meter\n",
    continental_shelf_depth);
    continental_shelf_depth=cont_shelf_depth;
  }

}

void Cavity::update_impl(double my_t, double my_dt) {

  // Make sure that sea water salinity and sea water potential
  // temperature fields are up to date:
  update_internal(my_t, my_dt);

  m_theta_ocean->average(m_t, m_dt);
  m_salinity_ocean->average(m_t, m_dt);

  Constants cc(*m_config);

  // prepare ocean input temperature and salinity per basin
  identifyMASK(ocean_contshelf_mask,"ocean_continental_shelf");
  compute_ocean_input_per_basin(cc);

  // define the ocean boxes below the ice shelves
  if (exicerises_set) {
    identifyMASK(icerise_mask,"icerises");}
  identifyMASK(ocean_mask,"ocean");
  round_basins();
  compute_distances();
  identify_ocean_box_mask(cc);

  set_ocean_input_fields(cc);

  //basal melt rates underneath ice shelves
  calculate_basal_melt_box1(cc);
  calculate_basal_melt_other_boxes(cc);
  calculate_basal_melt_missing_cells(cc);  //Assumes that mass flux is proportional to the shelf-base heat flux.

  m_shelfbtemp.copy_from(T_pressure_melting); // in-situ freezing point at the ice shelf base
  //
  basalmeltrate_shelf.scale(cc.rhoi);
  m_shelfbmassflux.copy_from(basalmeltrate_shelf); //TODO Check if scaling with ice density
}

// To be used solely in round_basins()
double Cavity::most_frequent_element(const std::vector<double> &v)
  {   // Precondition: v is not empty
      std::map<double, double> frequencyMap;
      int maxFrequency = 0;
      double mostFrequentElement = 0;
      for (double x : v)
      {
          double f = ++frequencyMap[x];
          if (f > maxFrequency)
          {
              maxFrequency = f;
              mostFrequentElement = x;
          }
      }

      return mostFrequentElement;
  }

//! Round non-integer basin mask values to integers.

//! Basin mask can have non-integer values from PISM regridding for points that lie at
//! basin boundaries.
//! Find such point here and set them to the integer value that is most frequent next to it.
void Cavity::round_basins() {

  // FIXME: THIS routine should be applied once in init, and roundbasins should
  // be stored as field (assumed the basins do not change with time).

  double id_fractional;
  std::vector<double> neighbours = {0,0,0,0};

  IceModelVec::AccessList list;
  list.add(cbasins);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // do not consider domain boundaries (they should be far from the shelves.)
    if ((i==0) | (j==0) | (i>(Mx-2)) | (j>(My-2))){
      id_fractional = 0.;
    } else {
      id_fractional = (cbasins)(i,j);
      neighbours[0] = (cbasins)(i+1,j+1);
      neighbours[1] = (cbasins)(i-1,j+1);
      neighbours[2] = (cbasins)(i-1,j-1);
      neighbours[3] = (cbasins)(i+1,j-1);

      // check if this is an interpolated number:
      // first condition: not an integer
      // second condition: has no neighbour with same value
      if ((id_fractional != round(id_fractional)) ||
          ((id_fractional != neighbours[0]) &&
          (id_fractional != neighbours[1]) &&
          (id_fractional != neighbours[2]) &&
          (id_fractional != neighbours[3]))){

        double most_frequent_neighbour = most_frequent_element(neighbours);
        (cbasins)(i,j) = most_frequent_neighbour;
        // m_log->message(2, "most frequent: %f at %d,%d\n",most_frequent_neighbour,i,j);
      }
    }

  }
}

//! Create masks that indicate ocean on continental shelf, ice rises as well as open ocean.

//! ocean_continental_shelf: ocean on the continental shelf without detached submarine islands
//! icerises: grounded ice not connected to the main ice body
//! ocean: ocean without holes in ice shelves, extends beyond continental shelf
//! We here use a search algorithm, starting at the center or the boundary of the domain.
//! We iteratively look for regions which satisfy one of the three types named above.
void Cavity::identifyMASK(IceModelVec2S &inputmask, std::string masktype) {

  m_log->message(5, "starting identifyMASK routine\n");

  // Assume that the center of the domain belongs to main ice body.
  int seed_x = (Mx - 1)/2,
      seed_y = (My - 1)/2;

  double linner_identified = 0.0,
         all_inner_identified = 1.0,
         previous_step_identified = 0.0;

  const IceModelVec2CellType &m_mask = *m_grid->variables().get_2d_cell_type("mask");
  const IceModelVec2S *topg = m_grid->variables().get_2d_scalar("bedrock_altitude");

  IceModelVec::AccessList list;
  list.add(inputmask);
  list.add(m_mask);
  list.add(*topg);

  inputmask.set(imask_unidentified);

  // Find starting points for iteration.
  if ((masktype=="ocean_continental_shelf" || masktype=="icerises") && (seed_x >= m_grid->xs()) && (seed_x < m_grid->xs()+m_grid->xm()) && (seed_y >= m_grid->ys())&& (seed_y < m_grid->ys()+m_grid->ym())){
    inputmask(seed_x,seed_y)=imask_inner;
  }
  else if (masktype=="ocean"){
    //assume that any point on the domain boundary belongs to the open ocean
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();
      if ((i==0) | (j==0) | (i>(Mx-2)) | (j>(My-2))){
        inputmask(i,j)=imask_inner;
      }
    }
  }

  // Iteratively find region which satisfies condition for coninental shelf ocean,
  // ice rise or open ocean.
  int iteration_round = 0;
  while(all_inner_identified > previous_step_identified){

    iteration_round+=1;
    previous_step_identified = all_inner_identified;

    for (Points p(*m_grid); p; p.next()) {

      const int i = p.i(), j = p.j();
      bool masktype_condition = false;

      if (masktype=="ocean_continental_shelf"){
        masktype_condition = (m_mask(i,j)!=maskocean || (*topg)(i,j) >= continental_shelf_depth);}
      else if (masktype=="icerises"){
        masktype_condition = (m_mask(i,j)==maskgrounded);
      }
      else if (masktype=="ocean"){
        masktype_condition = (m_mask(i,j)==maskocean);
      }

      if (masktype_condition && inputmask(i,j)==imask_unidentified &&
        (inputmask(i,j+1)==imask_inner || inputmask(i,j-1)==imask_inner ||
         inputmask(i+1,j)==imask_inner || inputmask(i-1,j)==imask_inner)){
         inputmask(i,j)=imask_inner;
         linner_identified+=1;
      }
      else if (masktype_condition == false){
        inputmask(i,j)=imask_outer;
      }
    }

    inputmask.update_ghosts();

    all_inner_identified = GlobalSum(m_grid->com, linner_identified);

  }

  // Set all unidentified grid cells to value for excluded areas (ice rises
  // or submarine islands)
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (inputmask(i,j)==imask_unidentified){
      inputmask(i,j)=imask_exclude;
    }

    if (masktype=="ocean_continental_shelf"){ //exclude ice covered parts
      if (m_mask(i,j)!=maskocean && inputmask(i,j) == imask_inner){
        inputmask(i,j) = imask_outer;
      }
    }
  }

}


//! Compute temperature and salinity input from ocean data by averaging.

//! We average over ocean_contshelf_mask for each basin.
//! We use dummy ocean data if no such average can be calculated.
//!

void Cavity::compute_ocean_input_per_basin(const Constants &cc) {

  m_log->message(5, "starting compute_ocean_input_per_basin routine \n");

  std::vector<double> lm_count(numberOfBasins); //count cells to take mean over for each basin
  std::vector<double> m_count(numberOfBasins);
  std::vector<double> lm_Sval(numberOfBasins); //add salinity for each basin
  std::vector<double> lm_Tval(numberOfBasins); //add temperature for each basin
  std::vector<double> m_Tval(numberOfBasins);
  std::vector<double> m_Sval(numberOfBasins);

 // initalize to zero per basin
  for(int shelf_id=0;shelf_id<numberOfBasins;shelf_id++){
    m_count[shelf_id]=0.0;
    lm_count[shelf_id]=0.0;
    lm_Sval[shelf_id]=0.0;
    lm_Tval[shelf_id]=0.0;
    m_Tval[shelf_id]=0.0;
    m_Sval[shelf_id]=0.0;
  }

  IceModelVec::AccessList list;
  list.add(*m_theta_ocean);
  list.add(*m_salinity_ocean);
  list.add(cbasins);
  list.add(ocean_contshelf_mask);

  // compute the sum for each basin
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (ocean_contshelf_mask(i,j) == imask_inner ){
      int shelf_id =(cbasins)(i,j);
      lm_count[shelf_id]+=1;
      lm_Sval[shelf_id]+=(*m_salinity_ocean)(i,j);
      lm_Tval[shelf_id]+=(*m_theta_ocean)(i,j);
    }

  }

  // Divide by number of grid cells if more than zero cells belong to the basin.
  // if no ocean_contshelf_mask values intersect with the basin, m_count is zero.
  // in such case, use dummy temperature and salinity. This could happen, for
  // example, if the ice shelf front advances beyond the continental shelf break.
  for(int shelf_id=0;shelf_id<numberOfBasins;shelf_id++) {

    m_count[shelf_id] = GlobalSum(m_grid->com, lm_count[shelf_id]);
    m_Sval[shelf_id] = GlobalSum(m_grid->com, lm_Sval[shelf_id]);
    m_Tval[shelf_id] = GlobalSum(m_grid->com, lm_Tval[shelf_id]);

    // if basin is not dummy basin 0 or there are no ocean cells in this basin to take the mean over.
    // FIXME: the following warning occurs once at initialization before input is available.
    // Please ignore this very first warning for now.
    if(shelf_id>0 && m_count[shelf_id]==0){
      m_log->message(2, "SIMPEL ocean WARNING: basin %d contains no cells with ocean data on continental shelf\n"
                        "(no values with ocean_contshelf_mask=2).\n"
                        "No mean salinity or temperature values are computed, instead using\n"
                        "the standard values T_dummy =%.3f, S_dummy=%.3f.\n"
                        "This might bias your basal melt rates, check your input data carefully.\n",
                        shelf_id, cc.T_dummy, cc.S_dummy);
      Toc_box0_vec[shelf_id] = cc.T_dummy;
      Soc_box0_vec[shelf_id] = cc.S_dummy;
    } else {
      m_Sval[shelf_id] = m_Sval[shelf_id] / m_count[shelf_id];
      m_Tval[shelf_id] = m_Tval[shelf_id] / m_count[shelf_id];

      Toc_box0_vec[shelf_id]=m_Tval[shelf_id];
      Soc_box0_vec[shelf_id]=m_Sval[shelf_id];
      m_log->message(5, "  %d: temp =%.3f, salinity=%.3f\n", shelf_id, Toc_box0_vec[shelf_id], Soc_box0_vec[shelf_id]);
    }
  }

}


//! Compute for each ice shelf cell distance to grounding line and ice front

//! DistGL: distance to grounding line
//! DistIF: distance to calving front
//! Ice holes within the shelf are treated like ice shelf cells,
//! if exicerises_set, also ice rises are treated like ice shelf cells.

void Cavity::compute_distances() {

  m_log->message(5, "starting compute_distances routine\n");

  double currentLabelGL = 1; // to find DistGL, 1 if floating and directly adjacent to a grounded cell
  double currentLabelIF = 1; // to find DistIF, 1 if floating and directly adjacent to an ocean cell

  double global_continue_loop = 1;
  double local_continue_loop  = 0;

  const IceModelVec2CellType &m_mask = *m_grid->variables().get_2d_cell_type("mask");

  IceModelVec::AccessList list;
  list.add(m_mask);
  list.add(DistIF);
  list.add(cbasins);
  list.add(DistGL);
  list.add(ocean_mask);

  if (exicerises_set) { list.add(icerise_mask); }

  DistGL.set(0);
  DistIF.set(0);

  // Find the grounding line and the ice front and
  // set DistGL to 1 if ice shelf cell is next to the grounding line,
  // set DistIF to 1 if ice shelf cell is next to the calving front.
  // Ice holes within the shelf are treated like ice shelf cells,
  // if exicerises_set, also ice rises are treated like ice shelf cells.
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    bool condition;
    if (exicerises_set) {
      condition = (m_mask(i,j)==maskfloating || icerise_mask(i,j)==imask_exclude || ocean_mask(i,j)==imask_exclude);
    }
    else {
      condition = (m_mask(i,j)==maskfloating || ocean_mask(i,j)==imask_exclude);
    }

    if (condition) { //if this is an ice shelf cell (or an ice rise) or a hole in an ice shelf

      // label the shelf cells adjacent to the grounding line with DistGL = 1
      bool neighbor_to_land;
      if (exicerises_set) {
        neighbor_to_land = (  icerise_mask(i,j+1)==imask_inner || icerise_mask(i,j-1)==imask_inner ||
          icerise_mask(i+1,j)==imask_inner || icerise_mask(i-1,j)==imask_inner ||
           icerise_mask(i+1,j+1)==imask_inner || icerise_mask(i+1,j-1)==imask_inner ||
           icerise_mask(i-1,j+1)==imask_inner || icerise_mask(i-1,j-1)==imask_inner );
      } else {
        neighbor_to_land = (  m_mask(i,j+1)<maskfloating || m_mask(i,j-1)<maskfloating ||
           m_mask(i+1,j)<maskfloating || m_mask(i-1,j)<maskfloating ||
          m_mask(i+1,j+1)<maskfloating || m_mask(i+1,j-1)<maskfloating ||
          m_mask(i-1,j+1)<maskfloating || m_mask(i-1,j-1)<maskfloating );
      }

      if (neighbor_to_land ){
        // i.e. there is a grounded neighboring cell (which is not ice rise!)
        DistGL(i,j) = currentLabelGL;
      } // no else

      // label the shelf cells adjacent to the calving front with DistIF = 1,
      // we do not need to exclude ice rises in this case.
      bool neighbor_to_ocean;
      neighbor_to_ocean = (ocean_mask(i,j+1)==imask_inner || ocean_mask(i,j-1)==imask_inner || ocean_mask(i+1,j)==imask_inner || ocean_mask(i-1,j)==imask_inner);

      if (neighbor_to_ocean) {
        DistIF(i,j) = currentLabelIF;
      }

    }
  }

  DistGL.update_ghosts();
  DistIF.update_ghosts();

  // DistGL calculation: Derive the distance from the grounding line for
  // all ice shelf cells iteratively.
  // Ice holes within the shelf are treated like ice shelf cells,
  // if exicerises_set, also ice rises are treated like ice shelf cells.
  global_continue_loop = 1;
  while( global_continue_loop !=0 ) {

    local_continue_loop = 0;

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      bool condition; // this cell is floating or an hole in the ice shelf (or an ice rise)
      if (exicerises_set) {
        condition = (m_mask(i,j)==maskfloating || icerise_mask(i,j)==imask_exclude || ocean_mask(i,j)==imask_exclude);
      }
      else {
        condition = (m_mask(i,j)==maskfloating || ocean_mask(i,j)==imask_exclude);
      }

      if ( condition && DistGL(i,j)==0 &&
        (DistGL(i,j+1)==currentLabelGL || DistGL(i,j-1)==currentLabelGL ||
        DistGL(i+1,j)==currentLabelGL || DistGL(i-1,j)==currentLabelGL) ) {
        // i.e. this is an shelf cell with no distance assigned yet and with a neighbor that has a distance assigned
          DistGL(i,j) = currentLabelGL+1;
          local_continue_loop = 1;
      } //if

    } // for

    currentLabelGL++;
    DistGL.update_ghosts();

    global_continue_loop = GlobalMax(m_grid->com, local_continue_loop);

  } // while: find DistGL

  // DistIF calculation: Derive the distance from the calving front for
  // all ice shelf cells iteratively.
  // Ice holes within the shelf are treated like ice shelf cells,
  // if exicerises_set, also ice rises are treated like ice shelf cells.
  global_continue_loop = 1; // start loop
  while( global_continue_loop !=0  ) {

    local_continue_loop = 0;

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      bool condition; // this cell is floating or an hole in the ice shelf (or an ice rise)
      if (exicerises_set) {
        condition = (m_mask(i,j)==maskfloating || icerise_mask(i,j)==imask_exclude || ocean_mask(i,j)==imask_exclude);
      }
      else {
        condition = (m_mask(i,j)==maskfloating || ocean_mask(i,j)==imask_exclude);
      }

      if ( condition && DistIF(i,j)==0 &&
        (DistIF(i,j+1)==currentLabelIF || DistIF(i,j-1)==currentLabelIF ||
        DistIF(i+1,j)==currentLabelIF || DistIF(i-1,j)==currentLabelIF) ) {
        // i.e. this is an shelf cell with no distance assigned yet and with a neighbor that has a distance assigned
          DistIF(i,j)=currentLabelIF+1;
          local_continue_loop = 1;
      } //if

    } // for

    currentLabelIF++;
    DistIF.update_ghosts();
    global_continue_loop = GlobalMax(m_grid->com, local_continue_loop);

  } // while: find DistIF

}


//! Compute the ocean_box_mask

//! Determine number of boxes for each basin based on max(DistGL).
//! Use a relative distance to the grounding line determine the ocean_box_mask
//! Finally, compute the extent of each ocean box in each basin.

void Cavity::identify_ocean_box_mask(const Constants &cc) {

  m_log->message(5, "starting identify_ocean_box_mask routine\n");

  // Find the maximal DistGL and DistIF for each basin
  std::vector<double> max_distGL(numberOfBasins);
  std::vector<double> max_distIF(numberOfBasins);
  std::vector<double> lmax_distGL(numberOfBasins);
  std::vector<double> lmax_distIF(numberOfBasins);

  const IceModelVec2CellType &m_mask = *m_grid->variables().get_2d_cell_type("mask");

  for(int shelf_id=0;shelf_id<numberOfBasins;shelf_id++){ max_distGL[shelf_id]=0.0; max_distIF[shelf_id]=0.0;lmax_distGL[shelf_id]=0.0; lmax_distIF[shelf_id]=0.0;}

  IceModelVec::AccessList list;
  list.add(cbasins);
  list.add(DistGL);
  list.add(DistIF);
  list.add(ocean_box_mask);
  list.add(m_mask);

  for (Points p(*m_grid); p; p.next()) {
  const int i = p.i(), j = p.j();
    int shelf_id = (cbasins)(i,j);

    if ( DistGL(i,j)> lmax_distGL[shelf_id] ) {
      lmax_distGL[shelf_id] = DistGL(i,j);
    }
    if ( DistIF(i,j)> lmax_distIF[shelf_id] ) {
      lmax_distIF[shelf_id] = DistIF(i,j);
    }
  }


  for (int l=0;l<numberOfBasins;l++){
    max_distGL[l] = GlobalMax(m_grid->com, lmax_distGL[l]);
    max_distIF[l] = GlobalMax(m_grid->com, lmax_distIF[l]);
  }

  // Compute the number of boxes for each basin
  // based on maximum distance between calving front and grounding line (in DistGL)
  // this is done by interpolating between nmin=1 and nmax=numberOfBoxes
  // this will be equal to numberOfBoxes for a 'large' ice shelf

  std::vector<int> lnumberOfBoxes_perBasin(numberOfBasins);

  int n_min = 1; //
  double max_distGL_ref = 500000; // meter //FIXME make this an input parameter
  double zeta = 0.5; // hard coded for now

  for (int l=0;l<numberOfBasins;l++){
    lnumberOfBoxes_perBasin[l] = 0;
    // FIXME: this is only correct for same dx and dy spacing.
    lnumberOfBoxes_perBasin[l] = n_min + static_cast<int>(
        round(pow((max_distGL[l]*dx/max_distGL_ref), zeta) *(numberOfBoxes-n_min)));
    lnumberOfBoxes_perBasin[l] = PetscMin(lnumberOfBoxes_perBasin[l],cc.default_numberOfBoxes);
    m_log->message(5, "lnumberOfBoxes[%d]=%d \n", l, lnumberOfBoxes_perBasin[l]);
  }


  // Define the ocean boxes in ocean_box_mask
  // this is based on the relative distance to the grounding line (computed from DistGL and DistIF)
  // and the number of boxes for the basin
  ocean_box_mask.set(0);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m_mask(i,j)==maskfloating && DistGL(i,j)>0 && DistIF(i,j)>0 && ocean_box_mask(i,j)==0){
      int shelf_id = (cbasins)(i,j);
      int n = lnumberOfBoxes_perBasin[shelf_id];
      // relative distance between grounding line and ice front
      double r = DistGL(i,j)*1.0/(DistGL(i,j)*1.0+DistIF(i,j)*1.0);

      for(int k=0;k<n;++k){

        // define the ocean_box_mask using rule (n-k)/n< (1-r)**2 <(n-k+1)/n
        // FIXME: is there a more elegant way to ensure float?
        if ( ((n*1.0-k*1.0-1.0)/(n*1.0) <= pow((1.0-r),2)) && (pow((1.0-r), 2) <= (n*1.0-k*1.0)/n*1.0) ){


          // ensure that boxnumber of a cell cannot be bigger than the distance to the grounding line
          if (DistGL(i,j) < k+1) {
            ocean_box_mask(i,j) = DistGL(i,j);
          // if smaller or equal, use default case: set to current box number
          } else{
            ocean_box_mask(i,j) = k+1;
          }
        }//if

      } //for
    }
  } // for

  // set all floating cells which have no ocean_box_mask value to numberOfBoxes+1.
  // For these, beckmann-goose melting will be applied, see calculate_basal_melt_missing_cells
  // those are the cells which are not reachable from grounding line or calving front.
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    if (m_mask(i,j)==maskfloating && ocean_box_mask(i,j)==0){ // floating
      ocean_box_mask(i,j) = numberOfBoxes + 1;
    }

  }

  // Compute the number of cells per box and basin and save to counter_boxes.
  const int nBoxes = numberOfBoxes+2;
  std::vector<std::vector<int> > lcounter_boxes(numberOfBasins, std::vector<int>(nBoxes));

  for (int shelf_id=0;shelf_id<numberOfBasins;shelf_id++){
    for (int l=0;l<nBoxes;l++){
      lcounter_boxes[shelf_id][l]=0;
    }
  }

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    int box_id = static_cast<int>(round(ocean_box_mask(i,j)));
    if (box_id > 0){ // floating
      int shelf_id = (cbasins)(i,j);
      lcounter_boxes[shelf_id][box_id]++;
    }
  }

  for (int shelf_id=0;shelf_id<numberOfBasins;shelf_id++){
    counter_boxes[shelf_id].resize(nBoxes);
    for (int l=0;l<nBoxes;l++){
      counter_boxes[shelf_id][l] = GlobalSum(m_grid->com, lcounter_boxes[shelf_id][l]);
    }
  }

}


//! Set ocean ocean input from box 0 as boundary condition for box 1.

//! Set ocean temperature and salinity (Toc_box0, Soc_box0)
//! from box 0 (in front of the ice shelf) as boundary condition for
//! box 1, which is the ocean box adjacent to the grounding line.
//! Toc_box0 and Soc_box0 were computed in function compute_ocean_input_per_basin.
//! We enforce that Toc_box0 is always at least the local pressure melting point.

void Cavity::set_ocean_input_fields(const Constants &cc) {

  m_log->message(5, "starting set_ocean_input_fields routine\n");

  const IceModelVec2S *ice_thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");
  const IceModelVec2CellType &m_mask = *m_grid->variables().get_2d_cell_type("mask");

  IceModelVec::AccessList list;
  list.add(*ice_thickness);
  list.add(cbasins);
  list.add(Soc_box0);
  list.add(Toc_box0);
  list.add(Toc);
  list.add(m_mask);

  double counterTpmp=0.0,
         lcounterTpmp = 0.0;

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // make sure all temperatures are zero at the beginning of each timestep
    Toc(i,j) = 273.15; // in K //FIXME delete?
    Toc_box0(i,j) = 273.15;  // in K
    Soc_box0(i,j) = 0.0; // in psu


    if (m_mask(i,j)==maskfloating){
      int shelf_id = (cbasins)(i,j);
      Toc_box0(i,j) = Toc_box0_vec[shelf_id];
      Soc_box0(i,j) =  Soc_box0_vec[shelf_id];

      // pressure in dbar, 1dbar = 10000 Pa = 1e4 kg m-1 s-2
      const double pressure = cc.rhoi * cc.earth_grav * (*ice_thickness)(i,j) * 1e-4;
      const double T_pmt = cc.a*Soc_box0(i,j) + cc.b - cc.c*pressure; // in Kelvin, here potential freezing point

      // temperature input for grounding line box should not be below pressure melting point
      if (  Toc_box0(i,j) < T_pmt ) {
        // Setting Toc_box0 a little higher than T_pmt ensures that later equations are well solvable.
        Toc_box0(i,j) = T_pmt + 0.001 ;
        lcounterTpmp+=1;
      }

    } // end if herefloating
  }

    counterTpmp = GlobalSum(m_grid->com, lcounterTpmp);
    if (counterTpmp > 0) {
      m_log->message(2, "SIMPEL ocean warning: temperature has been below pressure melting temperature in %.0f cases,\n"
                        "setting it to pressure melting temperature\n", counterTpmp);
    }

}


//! Compute the basal melt for each ice shelf cell in box 1

//! Here are the core physical equations of the SIMPEL model (for box1):
//! We here calculate basal melt rate, ambient ocean temperature and salinity
//! and overturning within box1. We calculate the average values at the boundary
//! between box 1 and box 2 as input for box 2.

void Cavity::calculate_basal_melt_box1(const Constants &cc) {

  m_log->message(5, "starting basal calculate_basal_melt_box1 routine\n");

  std::vector<double> lcounter_edge_of_box1_vector(numberOfBasins);
  std::vector<double> lmean_salinity_box1_vector(numberOfBasins);
  std::vector<double> lmean_temperature_box1_vector(numberOfBasins);
  std::vector<double> lmean_meltrate_box1_vector(numberOfBasins);
  std::vector<double> lmean_overturning_box1_vector(numberOfBasins);

  for (int shelf_id=0;shelf_id<numberOfBasins;shelf_id++){
    lcounter_edge_of_box1_vector[shelf_id]=0.0;
    lmean_salinity_box1_vector[shelf_id]=0.0;
    lmean_temperature_box1_vector[shelf_id]=0.0;
    lmean_meltrate_box1_vector[shelf_id]=0.0;
    lmean_overturning_box1_vector[shelf_id]=0.0;
  }

  const IceModelVec2S *ice_thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  IceModelVec::AccessList list;
  list.add(*ice_thickness);
  list.add(cbasins);
  list.add(ocean_box_mask);
  list.add(T_star);
  list.add(Toc_box0);
  list.add(Toc);
  list.add(Soc_box0);
  list.add(Soc);
  list.add(overturning);
  list.add(basalmeltrate_shelf);
  list.add(T_pressure_melting);

  double countHelpterm=0,
         lcountHelpterm=0;

  ocean_box_mask.update_ghosts();


  // basal melt rate, ambient temperature and salinity and overturning calculation
  // for each box1 grid cell.
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    int shelf_id = (cbasins)(i,j);

    // Make sure everything is at default values at the beginning of each timestep
    T_star(i,j) = 0.0; // in Kelvin
    Toc(i,j) = 273.15; // in Kelvin
    Soc(i,j) = 0.0; // in psu

    basalmeltrate_shelf(i,j) = 0.0;
    overturning(i,j) = 0.0;

    if ((ocean_box_mask(i,j) == 1) && (shelf_id > 0.0)){

      // pressure in dbar, 1dbar = 10000 Pa = 1e4 kg m-1 s-2
      const double pressure = cc.rhoi * cc.earth_grav * (*ice_thickness)(i,j) * 1e-4;
      T_star(i,j) = cc.a*Soc_box0(i,j) + cc.b - cc.c*pressure - Toc_box0(i,j); // in Kelvin

      //FIXME this assumes rectangular cell areas, adjust with real areas from projection
      double area_box1 = (counter_boxes[shelf_id][1] * dx * dy);

      double g1 = area_box1 * gamma_T ;
      double s1 = Soc_box0(i,j) / (cc.nu*cc.lambda);

      // These are the coefficients for solving the quadratic temperature equation
      // trough the p-q formula.
      double p_coeff = g1/( overturning_coeff*cc.rho_star * ( cc.beta * s1 - cc.alpha) ); // in 1 / (1/K) = K
      double q_coeff = (g1*T_star(i,j)) /
                        ( overturning_coeff*cc.rho_star * (  cc.beta * s1 - cc.alpha) ); // in K / (1/K) = K^2

      // This can only happen if T_star > 0.25*p_coeff, in particular T_star > 0
      // which can only happen for values of Toc_box0 close to the local pressure melting point
      if ((0.25*PetscSqr(p_coeff) - q_coeff) < 0.0) {

        m_log->message(5,
          "SIMPEL ocean WARNING: negative square root argument at %d, %d\n"
          "probably because of positive T_star=%f \n"
          "Not aborting, but setting square root to 0... \n",
          i, j, T_star(i,j));

        q_coeff=0.25*PetscSqr(p_coeff);
        lcountHelpterm+=1;
      }

      // temperature for box 1, p-q formula
      Toc(i,j) = Toc_box0(i,j) - ( -0.5*p_coeff + sqrt(0.25*PetscSqr(p_coeff) -q_coeff) ); // in Kelvin
      // salinity for box 1
      Soc(i,j) = Soc_box0(i,j) - (Soc_box0(i,j) / (cc.nu*cc.lambda)) * (Toc_box0(i,j) - Toc(i,j));  // in psu

      // potential pressure melting pointneeded to calculate thermal driving
      // using coefficients for potential temperature
      double potential_pressure_melting_point = cc.a*Soc(i,j) + cc.b - cc.c*pressure;

      // basal melt rate for box 1
      basalmeltrate_shelf(i,j) = (-gamma_T/(cc.nu*cc.lambda)) * (
      potential_pressure_melting_point - Toc(i,j));  // in m/s

      overturning(i,j) = overturning_coeff*cc.rho_star* (cc.beta*(Soc_box0(i,j)-Soc(i,j)) -
                            cc.alpha*(Toc_box0(i,j)-Toc(i,j))); // in m^3/s

      // average the temperature, salinity and overturning at the boundary between box1 and box2
      // this is used as input for box 2
      // using the values at the boundary helps smoothing discontinuities between boxes
      // (here we sum up)
      if (ocean_box_mask(i-1,j)==2 || ocean_box_mask(i+1,j)==2 ||
          ocean_box_mask(i,j-1)==2 || ocean_box_mask(i,j+1)==2){
        lcounter_edge_of_box1_vector[shelf_id]++;
        lmean_salinity_box1_vector[shelf_id] += Soc(i,j);
        lmean_temperature_box1_vector[shelf_id] += Toc(i,j); // in Kelvin
        lmean_meltrate_box1_vector[shelf_id] += basalmeltrate_shelf(i,j);
        lmean_overturning_box1_vector[shelf_id] += overturning(i,j);

      } // no else-case necessary since all variables are set to zero at the beginning of this routine

      // in situ pressure melting point
      T_pressure_melting(i,j) = cc.as*Soc(i,j) + cc.bs - cc.cs*pressure; //  in Kelvin

    }else { // i.e., not GL_box
        basalmeltrate_shelf(i,j) = 0.0;
    }
  }


  // average the temperature, salinity and overturning at the boundary between box1 and box2
  // (here we divide)
  for(int shelf_id=0;shelf_id<numberOfBasins;shelf_id++) {
    double counter_edge_of_box1_vector=0.0;

    counter_edge_of_box1_vector = GlobalSum(m_grid->com, lcounter_edge_of_box1_vector[shelf_id]);
    mean_salinity_boundary_vector[shelf_id] = GlobalSum(m_grid->com, lmean_salinity_box1_vector[shelf_id]);
    mean_temperature_boundary_vector[shelf_id] = GlobalSum(m_grid->com, lmean_temperature_box1_vector[shelf_id]);
    mean_overturning_box1_vector[shelf_id] = GlobalSum(m_grid->com, lmean_overturning_box1_vector[shelf_id]);

    if (counter_edge_of_box1_vector>0.0){
      mean_salinity_boundary_vector[shelf_id] = mean_salinity_boundary_vector[shelf_id]/counter_edge_of_box1_vector;
      mean_temperature_boundary_vector[shelf_id] = mean_temperature_boundary_vector[shelf_id]/counter_edge_of_box1_vector;
      mean_overturning_box1_vector[shelf_id] = mean_overturning_box1_vector[shelf_id]/counter_edge_of_box1_vector;
    } else {
      // This means that there is no [cell from the box 1 neighboring a cell from the box 2,
      // not necessarily that there is no box 1.
      mean_salinity_boundary_vector[shelf_id]=0.0;
      mean_temperature_boundary_vector[shelf_id]=0.0;
      mean_overturning_box1_vector[shelf_id]=0.0;
    }

    // print values at boundary between box 1 and box 2.
    m_log->message(5, "  %d: cnt=%.0f, sal=%.3f, temp=%.3f, over=%.1e \n", shelf_id,counter_edge_of_box1_vector,mean_salinity_boundary_vector[shelf_id],mean_temperature_boundary_vector[shelf_id],mean_overturning_box1_vector[shelf_id]) ;
  }

    countHelpterm = GlobalSum(m_grid->com, lcountHelpterm);
    if (countHelpterm > 0) {
      m_log->message(2, "SIMPEL ocean warning: square-root argument for temperature calculation "
                        "has been negative in %.0f cases!\n", countHelpterm);
    }

}

//! Compute the basal melt for each ice shelf cell in boxes other than box1

//! Here are the core physical equations of the SIMPEL model:
//! We here calculate basal melt rate, ambient ocean temperature and salinity.
//! Overturning is only calculated for box 1.
//! We calculate the average values at the boundary between box i and box i+1 as input for box i+1.

void Cavity::calculate_basal_melt_other_boxes(const Constants &cc) {

  m_log->message(5, "starting calculate_basal_melt_other_boxes routine\n");

  int nBoxes = static_cast<int>(round(numberOfBoxes+1));

  const IceModelVec2S *ice_thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  IceModelVec::AccessList list;
  list.add(*ice_thickness);
  list.add(cbasins);
  list.add(ocean_box_mask);
  list.add(T_star);
  list.add(Toc_box0);
  list.add(Toc);
  list.add(Soc_box0);
  list.add(Soc);
  list.add(overturning);
  list.add(basalmeltrate_shelf);
  list.add(T_pressure_melting);
  ocean_box_mask.update_ghosts();

  // Iterate over all Boxes i for i > 1
  // box number = numberOfBoxes+1 is used as identifier for Beckmann Goose calculation
  // for cells with missing input and excluded in loop here, i.e. boxi <nBoxes.
  for (int boxi=2; boxi <nBoxes; ++boxi) {

    m_log->message(5, "computing basal melt rate, temperature and salinity for box i = %d \n", boxi);

    double countGl0=0,
           lcountGl0=0;

    // averages at boundary between the current box and the subsequent box
    std::vector<double> lmean_salinity_boxi_vector(numberOfBasins); // in psu
    std::vector<double> lmean_temperature_boxi_vector(numberOfBasins); // in Kelvin
    std::vector<double> lcounter_edge_of_boxi_vector(numberOfBasins);

    for (int shelf_id=0;shelf_id<numberOfBasins;shelf_id++){
      lcounter_edge_of_boxi_vector[shelf_id] =0.0;
      lmean_salinity_boxi_vector[shelf_id]   =0.0;
      lmean_temperature_boxi_vector[shelf_id]=0.0;
    }

    // for box i compute the melt rates.

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      int shelf_id = (cbasins)(i,j);

      if (ocean_box_mask(i,j)==boxi && shelf_id > 0.0){

        double area_boxi, mean_salinity_in_boundary,
               mean_temperature_in_boundary,
               mean_overturning_in_box1;

        // get the input from previous box (is from box 1 if boxi=2)
        // overturning is only solved in box 1 and same for other boxes
        // temperature and salinity boundary values will be updated at the
        // end of this routine
        mean_salinity_in_boundary     = mean_salinity_boundary_vector[shelf_id];
        mean_temperature_in_boundary  = mean_temperature_boundary_vector[shelf_id]; // Kelvin
        mean_overturning_in_box1     = mean_overturning_box1_vector[shelf_id];

        // if there are no boundary values from the box before
        if (mean_salinity_in_boundary==0 || mean_overturning_in_box1==0 ||
            mean_temperature_in_boundary==0) {

          // set mask to Beckmann Goose identifier, will be handled in calculate_basal_melt_missing_cells
          ocean_box_mask(i,j) = numberOfBoxes+1;
          // flag to print warning later
          lcountGl0+=1;

        } else {

          // solve the SIMPLE physical model equations for boxes with boxi > 1
          // pressure in dbar, 1dbar = 10000 Pa = 1e4 kg m-1 s-2
          const double pressure = cc.rhoi * cc.earth_grav * (*ice_thickness)(i,j) * 1e-4;
          T_star(i,j) = cc.a*mean_salinity_in_boundary + cc.b - cc.c*pressure - mean_temperature_in_boundary;  // in Kelvin

          //FIXME this assumes rectangular cell areas, adjust with real areas from projection
          area_boxi = (counter_boxes[shelf_id][boxi] * dx * dy);

          // compute melt rates
          double g1 = area_boxi*gamma_T;
          double g2 = g1 / (cc.nu*cc.lambda);

          // temperature for Box i > 1
          Toc(i,j) = mean_temperature_in_boundary + g1 * T_star(i,j)/(mean_overturning_in_box1 + g1 - g2*cc.a*mean_salinity_in_boundary); // K

          // salinity for Box i > 1
          Soc(i,j) = mean_salinity_in_boundary - mean_salinity_in_boundary * (mean_temperature_in_boundary - Toc(i,j))/(cc.nu*cc.lambda); // psu

          // potential pressure melting pointneeded to calculate thermal driving
          // using coefficients for potential temperature
          double potential_pressure_melting_point = cc.a*Soc(i,j) + cc.b - cc.c*pressure;

          // basal melt rate for Box i > 1
          basalmeltrate_shelf(i,j) = (-gamma_T/(cc.nu*cc.lambda)) * (potential_pressure_melting_point - Toc(i,j)); // in m/s

          // in situ pressure melting point in Kelvin
          T_pressure_melting(i,j) = cc.as*Soc(i,j) + cc.bs - cc.cs*pressure;

          // average the temperature, salinity and overturning at the boundary between box i and box i+1
          // this is used as input for box i+1
          // using the values at the boundary helps smoothing discontinuities between boxes
          // (here we sum up)
          if (ocean_box_mask(i-1,j)==(boxi+1) || ocean_box_mask(i+1,j)==(boxi+1) || ocean_box_mask(i,j-1)==(boxi+1) || ocean_box_mask(i,j+1)==(boxi+1)){
            lcounter_edge_of_boxi_vector[shelf_id]++;
            lmean_salinity_boxi_vector[shelf_id] += Soc(i,j);
            lmean_temperature_boxi_vector[shelf_id] += Toc(i,j);
          } // no else-case necessary since all variables are set to zero at the beginning of this routine
        }
      } // no else-case, since  calculate_basal_melt_box1() and calculate_basal_melt_missing_cells() cover all other cases and we would overwrite those results here.
    }

    // average the temperature, salinity and overturning at the boundary between boxi and box i+1
    // (here we divide)
    for (int shelf_id=0;shelf_id<numberOfBasins;shelf_id++) {
      // overturning should not be changed, fixed from box 1
      double counter_edge_of_boxi_vector=0.0;
      counter_edge_of_boxi_vector = GlobalSum(m_grid->com, lcounter_edge_of_boxi_vector[shelf_id]);
      mean_salinity_boundary_vector[shelf_id] = GlobalSum(m_grid->com, lmean_salinity_boxi_vector[shelf_id]);
      mean_temperature_boundary_vector[shelf_id] = GlobalSum(m_grid->com, lmean_temperature_boxi_vector[shelf_id]); // in Kelvin

      if (counter_edge_of_boxi_vector>0.0){
        mean_salinity_boundary_vector[shelf_id] = mean_salinity_boundary_vector[shelf_id]/counter_edge_of_boxi_vector;
        mean_temperature_boundary_vector[shelf_id] = mean_temperature_boundary_vector[shelf_id]/counter_edge_of_boxi_vector; // in Kelvin
      } else {
        // This means that there is no cell in box i neighboring a cell from the box i+1,
        // not necessarily that there is no box i.
        mean_salinity_boundary_vector[shelf_id]=0.0; mean_temperature_boundary_vector[shelf_id]=0.0;
      }

      m_log->message(5, "  %d: cnt=%.0f, sal=%.3f, temp=%.3f, over=%.1e \n", shelf_id,counter_edge_of_boxi_vector,
                        mean_salinity_boundary_vector[shelf_id],mean_temperature_boundary_vector[shelf_id],
                        mean_overturning_box1_vector[shelf_id]) ;
    } // basins

    countGl0 = GlobalSum(m_grid->com, lcountGl0);
    if (countGl0 > 0) {
      m_log->message(2, "SIMPEL ocean WARNING: box %d, no boundary data from previous box in %.0f case(s)!\n"
                        "switching to Beckmann Goose (2003) meltrate calculation\n",
                        boxi,countGl0);
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

void Cavity::calculate_basal_melt_missing_cells(const Constants &cc) {

  m_log->message(5, "starting calculate_basal_melt_missing_cells routine\n");

  const IceModelVec2S *ice_thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  IceModelVec::AccessList list;
  list.add(*ice_thickness);
  list.add(cbasins);
  list.add(ocean_box_mask);
  list.add(Toc_box0);
  list.add(Toc);
  list.add(Soc_box0);
  list.add(Soc);
  list.add(overturning);
  list.add(basalmeltrate_shelf); // in m/s
  list.add(T_pressure_melting);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    int shelf_id = (cbasins)(i,j);

    // mainly at the boundary of computational domain,
    // or through erroneous basin mask
    if (shelf_id == 0) {

      basalmeltrate_shelf(i,j) = 0.0;

    }

    // cell with missing data identifier numberOfBoxes+1, as set in routines before
    if ( (shelf_id > 0) && (ocean_box_mask(i,j)==(numberOfBoxes+1)) ) {

      Toc(i,j) = Toc_box0(i,j); // in Kelvin
      Soc(i,j) = Soc_box0(i,j); // in psu

      // in dbar, 1dbar = 10000 Pa = 1e4 kg m-1 s-2
      const double pressure = cc.rhoi * cc.earth_grav * (*ice_thickness)(i,j) * 1e-4;

      // potential pressure melting point needed to calculate thermal driving
      // using coefficients for potential temperature
      // these are different to the ones used in POConstantPIK
      double potential_pressure_melting_point = cc.a*Soc(i,j) + cc.b - cc.c*pressure;

      double heatflux = cc.meltFactor * cc.rhow * cc.c_p_ocean * cc.default_gamma_T *
                         (Toc(i,j) - potential_pressure_melting_point);  // in W/m^2

      basalmeltrate_shelf(i,j) = heatflux / (cc.latentHeat * cc.rhoi); // in m s-1

      // in situ pressure melting point in Kelvin
      // this will be the temperature boundary condition at the ice at the shelf base
      T_pressure_melting(i,j) =  cc.as*Soc(i,j) + cc.bs - cc.cs*pressure;

    }

  }

}

} // end of namespace ocean
} // end of namespace pism
