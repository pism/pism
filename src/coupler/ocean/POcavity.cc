 // Copyright (C) 2012-2017 Ricarda Winkelmann, Ronja Reese, Torsten Albrecht
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
// Please cite this model as 
//
// Antarctic sub-shelf melt rates via PICO
// R. Reese, T. Albrecht, M. Mengel, X. Asay-Davis and R. Winkelmann 
// The Cryosphere (2018) 
//
// and 
//
// A box model of circulation and melting in ice shelf caverns
// D. Olbers & H. Hellmer
// Ocean Dynamics (2010), Volume 60, Issue 1, pp 141–153
// DOI: 10.1007/s10236-009-0252-z


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

  T_dummy = -1.5 + 273.15; // used as ocean temperature around Antarctica if no other data available (cold conditions)
  S_dummy = 34.7; // used as ocean salinity around Antarctica if no other data available (cold conditions)

  earth_grav = config.get_double("constants.standard_gravity");
  rhoi       = config.get_double("constants.ice.density");
  rhow       = config.get_double("constants.sea_water.density");
  rho_star   = 1033; // kg/m^3
  nu          = rhoi / rhow; // no unit

  latentHeat = config.get_double("constants.fresh_water.latent_heat_of_fusion"); //Joule / kg
  c_p_ocean  = 3974.0;       // J/(K*kg), specific heat capacity of ocean mixed layer
  lambda     = latentHeat / c_p_ocean;   // °C

  // Valus for linearized potential freezing point 
  a          = -0.0572;       // K/psu
  b          = 0.0788 + 273.15;       // K
  c          = 7.77e-4;      // K/dbar

  // in-situ pressure melting point from Jenkins et al. 2010 paper
  as          = -0.0573;       // K/psu
  bs          = 0.0832 + 273.15;       // K
  cs          = 7.53e-4;      // K/dbar

  alpha      = 7.5e-5;       // 1/K
  beta       = 7.7e-4;       // 1/psu

  default_gamma_T    = 2e-5;        // m/s, best fit value from PICO description paper 
  default_overturning_coeff    = 1e6;         // kg−1 s−1, best fit value from PICO description paper  

  // for shelf cells where PICO is not calculated, use Beckmann-Gooosse 2003
  // see also calculate_basal_melt_missing_cells() and POConstantPIK
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

  exicerises_set = options::Bool("-exclude_icerises", "exclude ice rises in ocean cavity model"); 

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

  // mask to identify ice shelves
  shelf_mask.create(m_grid, "shelf_mask", WITH_GHOSTS);
  shelf_mask.set_attrs("model_state", "mask displaying ice shelves","", "");
  m_variables.push_back(&shelf_mask);

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

  // mask displaying subglacial lakes - floating regions with no connection to the ocean
  lake_mask.create(m_grid, "lake_mask", WITH_GHOSTS);
  lake_mask.set_attrs("model_state", "mask displaying subglacial lakes","", "");
  m_variables.push_back(&lake_mask);

  // mask with distance (in boxes) to grounding line
  DistGL.create(m_grid, "DistGL", WITH_GHOSTS);
  DistGL.set_attrs("model_state", "mask displaying distance to grounding line","", "");
  m_variables.push_back(&DistGL);

  // mask with distance (in boxes) to ice front
  DistIF.create(m_grid, "DistIF", WITH_GHOSTS);
  DistIF.set_attrs("model_state", "mask displaying distance to ice shelf calving front","", "");
  m_variables.push_back(&DistIF);

  // salinity in ocean boxes
  Soc.create(m_grid, "Soc", WITHOUT_GHOSTS);
  Soc.set_attrs("model_state", "ocean salinity field","", "ocean salinity field");  // psu
  m_variables.push_back(&Soc);

  // salinity input for box 1
  Soc_box0.create(m_grid, "Soc_box0", WITHOUT_GHOSTS);
  Soc_box0.set_attrs("model_state", "ocean base salinity field","", "ocean base salinity field");  // psu
  m_variables.push_back(&Soc_box0);

  // temperature in ocean boxes
  Toc.create(m_grid, "Toc", WITHOUT_GHOSTS);
  Toc.set_attrs("model_state", "ocean temperature field","K", "ocean temperature field");
  m_variables.push_back(&Toc);

  // temperature input for box 1
  Toc_box0.create(m_grid, "Toc_box0", WITHOUT_GHOSTS);
  Toc_box0.set_attrs("model_state", "ocean base temperature","K", "ocean base temperature");
  m_variables.push_back(&Toc_box0);

  T_star.create(m_grid, "T_star", WITHOUT_GHOSTS);
  T_star.set_attrs("model_state", "T_star field","degree C", "T_star field");
  m_variables.push_back(&T_star);

  overturning.create(m_grid, "overturning", WITHOUT_GHOSTS);
  overturning.set_attrs("model_state", "cavity overturning","m^3 s-1", "cavity overturning"); 
  m_variables.push_back(&overturning);

  basalmeltrate_shelf.create(m_grid, "basalmeltrate_shelf", WITHOUT_GHOSTS);
  basalmeltrate_shelf.set_attrs("model_state", "PICO sub-shelf melt rate", "m/s",
                                "PICO sub-shelf melt rate");
  basalmeltrate_shelf.metadata().set_string("glaciological_units", "m year-1");
  basalmeltrate_shelf.write_in_glaciological_units = true;
  m_variables.push_back(&basalmeltrate_shelf);

  // in-situ pressure melting point
  T_pressure_melting.create(m_grid, "T_pressure_melting", WITHOUT_GHOSTS);
  T_pressure_melting.set_attrs("model_state", "pressure melting temperature at ice shelf base",
                        "Kelvin", "pressure melting temperature at ice shelf base"); 
  m_variables.push_back(&T_pressure_melting);

  numberOfBasins = 20;
  numberOfShelves = numberOfBasins; 
}

Cavity::~Cavity() {
  // empty
}

void Cavity::init_impl() {

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  m_log->message(2, "* Initializing the Potsdam Ice-shelf Cavity mOdel for the ocean ...\n");

  m_theta_ocean->init(m_filename, m_bc_period, m_bc_reference_time);
  m_salinity_ocean->init(m_filename, m_bc_period, m_bc_reference_time);

  cbasins.regrid(m_filename, CRITICAL);

  m_log->message(5, "PICO basin min=%f,max=%f\n",cbasins.min(),cbasins.max());

  Constants cc(*m_config);
  initBasinsOptions(cc);

  round_basins();

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


void Cavity::define_model_state_impl(const PIO &output) const {
  
  cbasins.define(output);
  shelf_mask.define(output);
  ocean_box_mask.define(output);
  icerise_mask.define(output);
  ocean_contshelf_mask.define(output);
  ocean_mask.define(output);
  lake_mask.define(output);
  DistGL.define(output);
  DistIF.define(output);
  Soc.define(output);
  Soc_box0.define(output);
  Toc.define(output);
  Toc_box0.define(output);
  T_star.define(output);
  overturning.define(output);
  basalmeltrate_shelf.define(output);
  T_pressure_melting.define(output);

  OceanModel::define_model_state_impl(output);
}

void Cavity::write_model_state_impl(const PIO &output) const {
  
  cbasins.write(output);
  shelf_mask.write(output);
  ocean_box_mask.write(output);
  icerise_mask.write(output);
  ocean_contshelf_mask.write(output);
  ocean_mask.write(output);
  lake_mask.write(output);
  DistGL.write(output);
  DistIF.write(output);
  Soc.write(output);
  Soc_box0.write(output);
  Toc.write(output);
  Toc_box0.write(output);
  T_star.write(output);
  overturning.write(output);
  basalmeltrate_shelf.write(output);
  T_pressure_melting.write(output);

  OceanModel::define_model_state_impl(output);
}


//! initialize PICO model variables, can be user-defined.

//! numberOfBasins: number of drainage basins
//! numberOfBoxes: maximum number of ocean boxes applied to large shelves, 
//!                for smaller shelves, the model may use less.
//! gamma_T: turbulent heat exchange coefficient for ice-ocean boundary layer
//! overturning_coeff: coefficient that scales strength of overturning circulation
//! continental_shelf_depth: threshold for definition of continental shelf area
//!                          area shallower than threshold is used for ocean input

void Cavity::initBasinsOptions(const Constants &cc) {

  m_log->message(5, "starting initBasinOptions\n");

  numberOfBasins = cc.default_numberOfBasins;
  numberOfBasins = options::Integer("-number_of_basins",
                                    "number of drainage basins for PICO model",
                                    numberOfBasins);
  numberOfBoxes = cc.default_numberOfBoxes;
  numberOfBoxes = options::Integer("-number_of_boxes",
                                    "number of ocean boxes for PICO model",
                                    numberOfBoxes);

  Toc_box0_vec.resize(numberOfBasins);
  Soc_box0_vec.resize(numberOfBasins);
  counter_boxes.resize(numberOfShelves, std::vector<double>(2,0));
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

  // Update sea water salinity and sea water potential
  // temperature fields
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
  identifyMASK(lake_mask,"lakes");
  identify_shelf_mask();
  round_basins();
  compute_distances();
  identify_ocean_box_mask(cc);

  set_ocean_input_fields(cc);

  //basal melt rates underneath ice shelves
  calculate_basal_melt_box1(cc);
  calculate_basal_melt_other_boxes(cc);
  calculate_basal_melt_missing_cells(cc);  

  m_shelfbtemp.copy_from(T_pressure_melting); // in-situ freezing point at the ice shelf base
  m_shelfbmassflux.copy_from(basalmeltrate_shelf); 
  m_shelfbmassflux.scale(cc.rhoi);
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

      if ((id_fractional != round(id_fractional)) ||
          ((id_fractional != neighbours[0]) &&
          (id_fractional != neighbours[1]) &&
          (id_fractional != neighbours[2]) &&
          (id_fractional != neighbours[3]))){

        double most_frequent_neighbour = most_frequent_element(neighbours);
        (cbasins)(i,j) = most_frequent_neighbour;
      }
    }

  }
}

//! Create masks that indicate ocean on continental shelf, ice rises as well as open ocean.

//! ocean_continental_shelf: ocean on the continental shelf without detached submarine islands
//! icerises: grounded ice not connected to the main ice body
//! ocean: ocean without holes in ice shelves, extends beyond continental shelf
//! lakes: subglacial lakes without access to the ocean
//! We here start at the center or the boundary of the domain and 
//! we iteratively look for regions which satisfy one of the three types named above.
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
  else if (masktype=="ocean" || masktype=="lakes"){
    //assume that some point on the domain boundary belongs to the open ocean
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
      else if (masktype=="lakes"){
        masktype_condition = (m_mask(i,j)==maskocean || m_mask(i,j)==maskfloating);
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



//! Create mask that indicates indicidual ice shelves

// Please use the latest PISM version with improved code.
void Cavity::identify_shelf_mask() {

  m_log->message(5, "starting identify_shelf_mask routine \n");

  const IceModelVec2CellType &m_mask = *m_grid->variables().get_2d_cell_type("mask");

  IceModelVec::AccessList list;
  list.add(shelf_mask);
  list.add(m_mask);
  list.add(lake_mask);
  if (exicerises_set) { list.add(icerise_mask); list.add(ocean_mask); }

  shelf_mask.set(0);

  std::vector<double> labels_counter (Mx*My,0);
  std::vector<double> labels_counter_global (Mx*My,0);

  double global_continue_loop = 1;
  double local_continue_loop  = 0;

  while( global_continue_loop !=0 ) { 
    local_continue_loop = 0;

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      bool condition;
      if (exicerises_set) { 
        condition = ( (m_mask(i,j)==maskfloating  && lake_mask(i,j)!=1) || icerise_mask(i,j)==imask_exclude );
      }
      else {
        condition = (m_mask(i,j)==maskfloating && lake_mask(i,j)!=1);
      }

      if (condition){
        
        if (shelf_mask(i,j)==0){ 
          
          if (shelf_mask(i-1,j)>0 && (shelf_mask(i+1,j)>0 && shelf_mask(i+1,j)>=shelf_mask(i-1,j) || shelf_mask(i+1,j)==0) && (shelf_mask(i,j-1)>0 && shelf_mask(i,j-1)>=shelf_mask(i-1,j) || shelf_mask(i,j-1)==0) && (shelf_mask(i,j+1)>0 && shelf_mask(i,j+1)>=shelf_mask(i-1,j) || shelf_mask(i,j+1)==0)){
            shelf_mask(i,j) = shelf_mask(i-1,j);
            labels_counter[static_cast<int>(shelf_mask(i-1,j))]++;
            local_continue_loop = 1;
          }
          if (shelf_mask(i+1,j)>0 && (shelf_mask(i-1,j)>0 && shelf_mask(i-1,j)>shelf_mask(i+1,j) || shelf_mask(i-1,j)==0) && (shelf_mask(i,j-1)>0 && shelf_mask(i,j-1)>=shelf_mask(i+1,j) || shelf_mask(i,j-1)==0) && (shelf_mask(i,j+1)>0 && shelf_mask(i,j+1)>=shelf_mask(i+1,j) || shelf_mask(i,j+1)==0) ){
            shelf_mask(i,j) = shelf_mask(i+1,j);
            labels_counter[static_cast<int>(shelf_mask(i+1,j))]++;
            local_continue_loop = 1;
          }
          if (shelf_mask(i,j-1)>0 && (shelf_mask(i-1,j)>0 && shelf_mask(i-1,j)>shelf_mask(i,j-1) || shelf_mask(i-1,j)==0) && (shelf_mask(i+1,j)>0 && shelf_mask(i+1,j)>shelf_mask(i,j-1) || shelf_mask(i+1,j)==0) && (shelf_mask(i,j+1)>0 && shelf_mask(i,j+1)>=shelf_mask(i,j-1) || shelf_mask(i,j+1)==0)){
            shelf_mask(i,j) = shelf_mask(i,j-1);
            labels_counter[static_cast<int>(shelf_mask(i,j-1))]++;
            local_continue_loop = 1;
          }
          if (shelf_mask(i,j+1)>0 && (shelf_mask(i-1,j)>0 && shelf_mask(i-1,j)>shelf_mask(i,j+1) || shelf_mask(i-1,j)==0) && (shelf_mask(i+1,j)>0 && shelf_mask(i+1,j)>shelf_mask(i,j+1) || shelf_mask(i+1,j)==0) && (shelf_mask(i,j-1)>0 && shelf_mask(i,j-1)>shelf_mask(i,j+1) || shelf_mask(i,j-1)==0)){
            shelf_mask(i,j) = shelf_mask(i,j+1);
            labels_counter[static_cast<int>(shelf_mask(i,j+1))]++;
            local_continue_loop = 1;
          }
          if (shelf_mask(i-1,j)==0 && shelf_mask(i+1,j)==0 && shelf_mask(i,j-1)==0 && shelf_mask(i,j+1)==0 ) {
            shelf_mask(i,j) = Mx*i + j; 
            labels_counter[static_cast<int>(shelf_mask(i,j))]++;
            local_continue_loop = 1;
          }
        } else { 
          
          if ( shelf_mask(i-1,j)>0 && shelf_mask(i-1,j)<shelf_mask(i,j) ){
            labels_counter[static_cast<int>(shelf_mask(i,j))]--;
            shelf_mask(i,j) = shelf_mask(i-1,j);
            labels_counter[static_cast<int>(shelf_mask(i-1,j))]++;
            local_continue_loop = 1;
          }          
          if ( shelf_mask(i+1,j)>0 && shelf_mask(i+1,j)<shelf_mask(i,j) ){
            labels_counter[static_cast<int>(shelf_mask(i,j))]--;
            shelf_mask(i,j) = shelf_mask(i+1,j);
            labels_counter[static_cast<int>(shelf_mask(i+1,j))]++;
            local_continue_loop = 1;
          } 
          if ( shelf_mask(i,j-1)>0 && shelf_mask(i,j-1)<shelf_mask(i,j) ){
            labels_counter[static_cast<int>(shelf_mask(i,j))]--;
            shelf_mask(i,j) = shelf_mask(i,j-1);
            labels_counter[static_cast<int>(shelf_mask(i,j-1))]++;
            local_continue_loop = 1;
          } 
          if ( shelf_mask(i,j+1)>0 && shelf_mask(i,j+1)<shelf_mask(i,j) ){
            labels_counter[static_cast<int>(shelf_mask(i,j))]--;
            shelf_mask(i,j) = shelf_mask(i,j+1);
            labels_counter[static_cast<int>(shelf_mask(i,j+1))]++;
            local_continue_loop = 1;
          } 
          if ( shelf_mask(i-1,j-1)>0 && shelf_mask(i-1,j-1)<shelf_mask(i,j) ){
            labels_counter[static_cast<int>(shelf_mask(i,j))]--;
            shelf_mask(i,j) = shelf_mask(i-1,j-1);
            labels_counter[static_cast<int>(shelf_mask(i-1,j-1))]++;
            local_continue_loop = 1;
          } 
          if ( shelf_mask(i-1,j+1)>0 && shelf_mask(i-1,j+1)<shelf_mask(i,j) ){
            labels_counter[static_cast<int>(shelf_mask(i,j))]--;
            shelf_mask(i,j) = shelf_mask(i-1,j+1);
            labels_counter[static_cast<int>(shelf_mask(i-1,j+1))]++;
            local_continue_loop = 1;
          } 
          if ( shelf_mask(i+1,j-1)>0 && shelf_mask(i+1,j-1)<shelf_mask(i,j) ){
            labels_counter[static_cast<int>(shelf_mask(i,j))]--;
            shelf_mask(i,j) = shelf_mask(i+1,j-1);
            labels_counter[static_cast<int>(shelf_mask(i+1,j-1))]++;
            local_continue_loop = 1;
          } 
          if ( shelf_mask(i+1,j+1)>0 && shelf_mask(i+1,j+1)<shelf_mask(i,j) ){
            labels_counter[static_cast<int>(shelf_mask(i,j))]--;
            shelf_mask(i,j) = shelf_mask(i+1,j+1);
            labels_counter[static_cast<int>(shelf_mask(i+1,j+1))]++;
            local_continue_loop = 1;
          } 
        } // if else
      } // if 
    } // 
    
    shelf_mask.update_ghosts();
    global_continue_loop = GlobalMax(m_grid->com, local_continue_loop);

  } // while 
  
  
  for (int k=0;k<Mx*My;k++){ labels_counter_global[k] = GlobalSum(m_grid->com, labels_counter[k]);}
  
  double new_label_current = 1;
  std::vector<double> new_labels (Mx*My,0); 
  
  for (int k=0;k<Mx*My;k++){
    if (labels_counter_global[k] != 0){
      new_labels[k] = new_label_current;
      new_label_current++;
    } // no else case, skip 
  }

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    int label = static_cast<int>(shelf_mask(i,j));   
    shelf_mask(i,j) = new_labels[label];
  }

  numberOfShelves = new_label_current; 
  m_log->message(5, "Number of shelves = %d\n", numberOfShelves-1); // internally calculated with +1

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
  for(int basin_id=0;basin_id<numberOfBasins;basin_id++){
    m_count[basin_id]=0.0;
    lm_count[basin_id]=0.0;
    lm_Sval[basin_id]=0.0;
    lm_Tval[basin_id]=0.0;
    m_Tval[basin_id]=0.0;
    m_Sval[basin_id]=0.0;
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
      int basin_id =(cbasins)(i,j);
      lm_count[basin_id]+=1;
      lm_Sval[basin_id]+=(*m_salinity_ocean)(i,j);
      lm_Tval[basin_id]+=(*m_theta_ocean)(i,j);
    }

  }

  // Divide by number of grid cells if more than zero cells belong to the basin.
  // if no ocean_contshelf_mask values intersect with the basin, m_count is zero.
  // in such case, use dummy temperature and salinity. This could happen, for
  // example, if the ice shelf front advances beyond the continental shelf break.
  for(int basin_id=0;basin_id<numberOfBasins;basin_id++) {

    m_count[basin_id] = GlobalSum(m_grid->com, lm_count[basin_id]);
    m_Sval[basin_id] = GlobalSum(m_grid->com, lm_Sval[basin_id]);
    m_Tval[basin_id] = GlobalSum(m_grid->com, lm_Tval[basin_id]);

    // if basin is not dummy basin 0 or there are no ocean cells in this basin to take the mean over.
    // Please ignore this very first warning (during the initialization step)
    if(basin_id>0 && m_count[basin_id]==0){
      m_log->message(2, "PICO ocean WARNING: basin %d contains no cells with ocean data on continental shelf\n"
                        "(no values with ocean_contshelf_mask=2).\n"
                        "No mean salinity or temperature values are computed, instead using\n"
                        "the standard values T_dummy =%.3f, S_dummy=%.3f.\n"
                        "This might bias your basal melt rates, check your input data carefully.\n",
                        basin_id, cc.T_dummy, cc.S_dummy);
      Toc_box0_vec[basin_id] = cc.T_dummy;
      Soc_box0_vec[basin_id] = cc.S_dummy;
    } else {
      m_Sval[basin_id] = m_Sval[basin_id] / m_count[basin_id];
      m_Tval[basin_id] = m_Tval[basin_id] / m_count[basin_id];

      Toc_box0_vec[basin_id]=m_Tval[basin_id];
      Soc_box0_vec[basin_id]=m_Sval[basin_id];
      m_log->message(5, "  %d: temp =%.3f, salinity=%.3f\n", basin_id, Toc_box0_vec[basin_id], Soc_box0_vec[basin_id]);
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
        // i.e. there is a grounded neighboring cell (which is not ice rise)
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
  global_continue_loop = 1; 
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
  std::vector<double> max_distGL(numberOfShelves); 
  std::vector<double> max_distIF(numberOfShelves); 
  std::vector<double> lmax_distGL(numberOfShelves); 
  std::vector<double> lmax_distIF(numberOfShelves); 

  double lmax_distGL_ref = 0.0; 
  double max_distGL_ref = 0.0; 

  const IceModelVec2CellType &m_mask = *m_grid->variables().get_2d_cell_type("mask");

  for(int shelf_id=0;shelf_id<numberOfShelves;shelf_id++){ max_distGL[shelf_id]=0.0; max_distIF[shelf_id]=0.0;lmax_distGL[shelf_id]=0.0; lmax_distIF[shelf_id]=0.0;}

  IceModelVec::AccessList list;
  list.add(DistGL);
  list.add(DistIF);
  list.add(ocean_box_mask);
  list.add(lake_mask);
  list.add(m_mask);
  list.add(shelf_mask);

  for (Points p(*m_grid); p; p.next()) {
  const int i = p.i(), j = p.j();
    int shelf_id = shelf_mask(i,j);
    if ( DistGL(i,j)> lmax_distGL[shelf_id] ) {
      lmax_distGL[shelf_id] = DistGL(i,j);
    }
    if ( DistIF(i,j)> lmax_distIF[shelf_id] ) {
      lmax_distIF[shelf_id] = DistIF(i,j);
    }
    if (DistGL(i,j)>lmax_distGL_ref){
      lmax_distGL_ref = DistGL(i,j);
    }
  }

  for (int l=0;l<numberOfShelves;l++){
    max_distGL[l] = GlobalMax(m_grid->com, lmax_distGL[l]);
    max_distIF[l] = GlobalMax(m_grid->com, lmax_distIF[l]);
  }
  max_distGL_ref = GlobalMax(m_grid->com, lmax_distGL_ref);  

  // Compute the number of boxes for each basin
  // based on maximum distance between calving front and grounding line (in DistGL)
  // this is done by interpolating between nmin=1 and nmax=numberOfBoxes
  // this will be equal to numberOfBoxes for a 'large' ice shelf

  std::vector<int> lnumberOfBoxes_perShelf(numberOfShelves); 

  int n_min = 1; 
  double zeta = 0.5; 

  for (int l=0;l<numberOfShelves;l++){
    lnumberOfBoxes_perShelf[l] = 0;
    lnumberOfBoxes_perShelf[l] = n_min + static_cast<int>( 
    		round(pow((max_distGL[l]/max_distGL_ref), zeta) *(numberOfBoxes-n_min))); 
    lnumberOfBoxes_perShelf[l] = PetscMin(lnumberOfBoxes_perShelf[l],cc.default_numberOfBoxes);
    m_log->message(5, "lnumberOfBoxes[%d]=%d \n", l, lnumberOfBoxes_perShelf[l]);
  }


  // Define the ocean boxes in ocean_box_mask
  // this is based on the relative distance to the grounding line (computed from DistGL and DistIF)
  // and the number of boxes for the basin
  ocean_box_mask.set(0);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m_mask(i,j)==maskfloating && DistGL(i,j)>0 && DistIF(i,j)>0 && ocean_box_mask(i,j)==0){
      int shelf_id = shelf_mask(i,j);
      int n = lnumberOfBoxes_perShelf[shelf_id];
      // relative distance between grounding line and ice front
      double r = DistGL(i,j)*1.0/(DistGL(i,j)*1.0+DistIF(i,j)*1.0);

      for(int k=0;k<n;++k){

        // define the ocean_box_mask using rule (n-k)/n< (1-r)**2 <(n-k+1)/n
        if ( ((n*1.0-k*1.0-1.0)/(n*1.0) <= pow((1.0-r),2)) && (pow((1.0-r), 2) <= (n*1.0-k*1.0)/n*1.0) ){

          // ensure that boxnumber of a cell cannot be bigger than the distance to the grounding line
          if (DistGL(i,j) < k+1) {
            ocean_box_mask(i,j) = DistGL(i,j);
          // if smaller or equal, set to current box number
          } else{
            ocean_box_mask(i,j) = k+1;
          }
        }//if

      } //for
    }
  } // for

  // set all floating cells which have no ocean_box_mask value to numberOfBoxes+1.
  // do not apply this to cells that are subglacial lakes, 
  // since they are not accessible for ocean waters and hence should not be treated in this model
  // For these, Beckmann-Goosse (2003) melting will be applied, see calculate_basal_melt_missing_cells
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    if (m_mask(i,j)==maskfloating && ocean_box_mask(i,j)==0 && lake_mask(i,j)!=1){ // floating, no sub-glacial lake
      ocean_box_mask(i,j) = numberOfBoxes + 1;
    }

  }

  // Compute the number of cells per box and basin and save to counter_boxes
  const int nBoxes = numberOfBoxes+2;

  counter_boxes.resize(numberOfShelves, std::vector<double>(2,0));
  std::vector<std::vector<double> > lcounter_boxes(numberOfShelves, std::vector<double>(nBoxes));

  for (int shelf_id=0;shelf_id<numberOfShelves;shelf_id++){
    for (int l=0;l<nBoxes;l++){ 
      lcounter_boxes[shelf_id][l]=0;
    }
  }

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    int box_id = static_cast<int>(round(ocean_box_mask(i,j)));
    if (box_id > 0){ // floating
      int shelf_id = shelf_mask(i,j);
      lcounter_boxes[shelf_id][box_id]++;
    }
  }

  for (int shelf_id=0;shelf_id<numberOfShelves;shelf_id++){
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
  list.add(shelf_mask);


  // compute for each shelf the number of cells  within each basin
  std::vector<std::vector<double> > lcounter_shelf_cells_in_basin(numberOfShelves, std::vector<double>(numberOfBasins));
  std::vector<std::vector<double> > counter_shelf_cells_in_basin(numberOfShelves, std::vector<double>(numberOfBasins)); 

  // compute the number of all shelf cells
  std::vector<double> lcounter_shelf_cells(numberOfShelves);  
  std::vector<double> counter_shelf_cells(numberOfShelves);  

  for (int shelf_id=0;shelf_id<numberOfShelves;shelf_id++){
    lcounter_shelf_cells[shelf_id] =0;
    for (int basin_id=0;basin_id<numberOfBasins;basin_id++){ 
      lcounter_shelf_cells_in_basin[shelf_id][basin_id]=0;
    }
  }
  
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    int shelf_id = shelf_mask(i,j);
    int basin_id = (cbasins)(i,j);
    lcounter_shelf_cells_in_basin[shelf_id][basin_id]++;  
    lcounter_shelf_cells[shelf_id]++;
  }
  
  for (int shelf_id=0;shelf_id<numberOfShelves;shelf_id++){
    counter_shelf_cells[shelf_id] = GlobalSum(m_grid->com, lcounter_shelf_cells[shelf_id]);
    for (int basin_id=0;basin_id<numberOfBasins;basin_id++){
      counter_shelf_cells_in_basin[shelf_id][basin_id] = GlobalSum(m_grid->com, lcounter_shelf_cells_in_basin[shelf_id][basin_id]);
    }
  }


  // now set temp and salinity in box0:
  double counterTpmp=0.0,
         lcounterTpmp = 0.0;

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    Toc(i,j) = 273.15; // in K 
    Toc_box0(i,j) = 0.0;  // in K
    Soc_box0(i,j) = 0.0; // in psu


    if (m_mask(i,j)==maskfloating && shelf_mask(i,j)>0){ // shelf_mask = 0 for lakes 
      int shelf_id = shelf_mask(i,j);
      
      // weighted input depending on the number of shelf cells in each basin
      for (int basin_id=1;basin_id<numberOfBasins;basin_id++){ 
          Toc_box0(i,j) += Toc_box0_vec[basin_id]*counter_shelf_cells_in_basin[shelf_id][basin_id]/counter_shelf_cells[shelf_id];
          Soc_box0(i,j) += Soc_box0_vec[basin_id]*counter_shelf_cells_in_basin[shelf_id][basin_id]/counter_shelf_cells[shelf_id];
      }
      
      // pressure in dbar, 1dbar = 10000 Pa = 1e4 kg m-1 s-2
      const double pressure = cc.rhoi * cc.earth_grav * (*ice_thickness)(i,j) * 1e-4;
      const double T_pmt = cc.a*Soc_box0(i,j) + cc.b - cc.c*pressure; // in Kelvin, potential freezing point

      // temperature input for grounding line box should not be below pressure melting point
      if (  Toc_box0(i,j) < T_pmt ) {
        // Setting Toc_box0 slightly higher than T_pmt ensures that later equations are numerically solvable
        Toc_box0(i,j) = T_pmt + 0.001 ;
        lcounterTpmp+=1;
      }

    } // end if herefloating
  }

    counterTpmp = GlobalSum(m_grid->com, lcounterTpmp);
    if (counterTpmp > 0) {
      m_log->message(2, "PICO ocean warning: temperature has been below pressure melting temperature in %.0f cases,\n"
                        "setting it to pressure melting temperature\n", counterTpmp);
    }

}


//! Compute the basal melt for each ice shelf cell in box 1

//! Here are the core physical equations of the PICO model (for box1):
//! We here calculate basal melt rate, ambient ocean temperature and salinity
//! and overturning within box1. We calculate the average values in box 1 as input for box 2.

void Cavity::calculate_basal_melt_box1(const Constants &cc) {

  m_log->message(5, "starting basal calculate_basal_melt_box1 routine\n");

  std::vector<double> lcounter_edge_of_box1_vector(numberOfShelves);
  std::vector<double> lmean_salinity_box1_vector(numberOfShelves);
  std::vector<double> lmean_temperature_box1_vector(numberOfShelves);
  std::vector<double> lmean_meltrate_box1_vector(numberOfShelves);
  std::vector<double> lmean_overturning_box1_vector(numberOfShelves);

  for (int shelf_id=0;shelf_id<numberOfShelves;shelf_id++){
    lcounter_edge_of_box1_vector[shelf_id]=0.0;
    lmean_salinity_box1_vector[shelf_id]=0.0;
    lmean_temperature_box1_vector[shelf_id]=0.0;
    lmean_meltrate_box1_vector[shelf_id]=0.0;
    lmean_overturning_box1_vector[shelf_id]=0.0;
  }

  mean_salinity_boundary_vector.resize(numberOfShelves);
  mean_temperature_boundary_vector.resize(numberOfShelves);
  mean_overturning_box1_vector.resize(numberOfShelves);

  const IceModelVec2S *ice_thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  IceModelVec::AccessList list;
  list.add(*ice_thickness);
  list.add(ocean_box_mask);
  list.add(T_star);
  list.add(Toc_box0);
  list.add(Toc);
  list.add(Soc_box0);
  list.add(Soc);
  list.add(overturning);
  list.add(basalmeltrate_shelf);
  list.add(T_pressure_melting);
  list.add(shelf_mask);

  double countHelpterm=0,
         lcountHelpterm=0;

  ocean_box_mask.update_ghosts();


  // basal melt rate, ambient temperature and salinity and overturning calculation
  // for each box1 grid cell.
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    int shelf_id = shelf_mask(i,j);

    T_star(i,j) = 0.0; // in Kelvin
    Toc(i,j) = 273.15; // in Kelvin
    Soc(i,j) = 0.0; // in psu

    basalmeltrate_shelf(i,j) = 0.0;
    overturning(i,j) = 0.0;
    T_pressure_melting(i,j)= 0.0;

    if ((ocean_box_mask(i,j) == 1) && (shelf_id > 0.0)){

      // pressure in dbar, 1dbar = 10000 Pa = 1e4 kg m-1 s-2
      const double pressure = cc.rhoi * cc.earth_grav * (*ice_thickness)(i,j) * 1e-4;
      T_star(i,j) = cc.a*Soc_box0(i,j) + cc.b - cc.c*pressure - Toc_box0(i,j); // in Kelvin

      double area_box1 = (counter_boxes[shelf_id][1] * dx * dy);

      double g1 = area_box1 * gamma_T ;
      double s1 = Soc_box0(i,j) / (cc.nu*cc.lambda);

      // These are the coefficients for solving the quadratic temperature equation
      double p_coeff = g1/( overturning_coeff*cc.rho_star * ( cc.beta * s1 - cc.alpha) ); // in 1 / (1/K) = K
      double q_coeff = (g1*T_star(i,j)) /
                        ( overturning_coeff*cc.rho_star * (  cc.beta * s1 - cc.alpha) ); // in K / (1/K) = K^2

      // This can only happen if T_star > 0.25*p_coeff, in particular T_star > 0
      // which can only happen for values of Toc_box0 close to the local pressure melting point
      if ((0.25*PetscSqr(p_coeff) - q_coeff) < 0.0) {

        m_log->message(5,
          "PICO ocean WARNING: negative square root argument at %d, %d\n"
          "probably because of positive T_star=%f \n"
          "Not aborting, but setting square root to 0... \n",
          i, j, T_star(i,j));

        q_coeff=0.25*PetscSqr(p_coeff);
        lcountHelpterm+=1;
      }

      // temperature for box 1
      Toc(i,j) = Toc_box0(i,j) - ( -0.5*p_coeff + sqrt(0.25*PetscSqr(p_coeff) -q_coeff) ); // in Kelvin
      // salinity for box 1
      Soc(i,j) = Soc_box0(i,j) - (Soc_box0(i,j) / (cc.nu*cc.lambda)) * (Toc_box0(i,j) - Toc(i,j));  // in psu

      // potential pressure melting point 
      double potential_pressure_melting_point = cc.a*Soc(i,j) + cc.b - cc.c*pressure;

      // basal melt rate for box 1
      basalmeltrate_shelf(i,j) = (-gamma_T/(cc.nu*cc.lambda)) * (
      potential_pressure_melting_point - Toc(i,j));  // in m/s

      overturning(i,j) = overturning_coeff*cc.rho_star* (cc.beta*(Soc_box0(i,j)-Soc(i,j)) -
                            cc.alpha*(Toc_box0(i,j)-Toc(i,j))); // in m^3/s

      lcounter_edge_of_box1_vector[shelf_id]++;
      lmean_salinity_box1_vector[shelf_id] += Soc(i,j);
      lmean_temperature_box1_vector[shelf_id] += Toc(i,j); // in Kelvin
      lmean_meltrate_box1_vector[shelf_id] += basalmeltrate_shelf(i,j);
      lmean_overturning_box1_vector[shelf_id] += overturning(i,j);

      // in situ pressure melting point
      T_pressure_melting(i,j) = cc.as*Soc(i,j) + cc.bs - cc.cs*pressure; //  in Kelvin

    }else { // not GL_box
        basalmeltrate_shelf(i,j) = 0.0;
    }
  }


  // average the temperature, salinity and overturning over box1
  // (here we divide)
  for(int shelf_id=0;shelf_id<numberOfShelves;shelf_id++) {
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
      mean_salinity_boundary_vector[shelf_id]=0.0;
      mean_temperature_boundary_vector[shelf_id]=0.0;
      mean_overturning_box1_vector[shelf_id]=0.0;
    }

    m_log->message(5, "  %d: cnt=%.0f, sal=%.3f, temp=%.3f, over=%.1e \n", shelf_id,counter_edge_of_box1_vector,mean_salinity_boundary_vector[shelf_id],mean_temperature_boundary_vector[shelf_id],mean_overturning_box1_vector[shelf_id]) ;
  }

    countHelpterm = GlobalSum(m_grid->com, lcountHelpterm);
    if (countHelpterm > 0) {
      m_log->message(2, "PICO ocean warning: square-root argument for temperature calculation "
                        "has been negative in %.0f cases!\n", countHelpterm);
    }

}

//! Compute the basal melt for each ice shelf cell in boxes other than box1

//! Here are the core physical equations of the PICO model:
//! We here calculate basal melt rate, ambient ocean temperature and salinity.
//! Overturning is only calculated for box 1.
//! We calculate the average valuesin box i as input for box i+1.

void Cavity::calculate_basal_melt_other_boxes(const Constants &cc) {

  m_log->message(5, "starting calculate_basal_melt_other_boxes routine\n");

  int nBoxes = static_cast<int>(round(numberOfBoxes+1));

  const IceModelVec2S *ice_thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  IceModelVec::AccessList list;
  list.add(*ice_thickness);
  list.add(ocean_box_mask);
  list.add(T_star);
  list.add(Toc_box0);
  list.add(Toc);
  list.add(Soc_box0);
  list.add(Soc);
  list.add(overturning);
  list.add(basalmeltrate_shelf);
  list.add(T_pressure_melting);
  list.add(shelf_mask);
  ocean_box_mask.update_ghosts();

  // Iterate over all Boxes i for i > 1
  // box number = numberOfBoxes+1 is used as identifier for Beckmann Goose calculation
  // for cells with missing input and excluded in loop here
  for (int boxi=2; boxi <nBoxes; ++boxi) {

    m_log->message(5, "computing basal melt rate, temperature and salinity for box i = %d \n", boxi);

    double countGl0=0,
           lcountGl0=0;

    std::vector<double> lmean_salinity_boxi_vector(numberOfShelves); // in psu
    std::vector<double> lmean_temperature_boxi_vector(numberOfShelves); // in Kelvin
    std::vector<double> lcounter_edge_of_boxi_vector(numberOfShelves);

    for (int shelf_id=0;shelf_id<numberOfShelves;shelf_id++){
      lcounter_edge_of_boxi_vector[shelf_id] =0.0;
      lmean_salinity_boxi_vector[shelf_id]   =0.0;
      lmean_temperature_boxi_vector[shelf_id]=0.0;
    }

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      int shelf_id = shelf_mask(i,j);

      if (ocean_box_mask(i,j)==boxi && shelf_id > 0.0){

        double area_boxi, mean_salinity_in_boundary,
               mean_temperature_in_boundary,
               mean_overturning_in_box1;

        // get the input from previous box (is from box 1 if boxi=2)
        // overturning is only solved in box 1 and same for other boxes
        // temperature and salinity boundary values will be updated at the
        // end of this routine
        mean_salinity_in_boundary     = mean_salinity_boundary_vector[shelf_id]; //psu
        mean_temperature_in_boundary  = mean_temperature_boundary_vector[shelf_id]; // Kelvin
        mean_overturning_in_box1     = mean_overturning_box1_vector[shelf_id];

        // if there are no boundary values from the box before
        if (mean_salinity_in_boundary==0 || mean_overturning_in_box1==0 ||
            mean_temperature_in_boundary==0) {

          // set mask to Beckmann Goose identifier, will be handled in calculate_basal_melt_missing_cells
          ocean_box_mask(i,j) = numberOfBoxes+1;
          lcountGl0+=1;

        } else {

          // solve the PICO physical model equations for boxes with boxi > 1
          // pressure in dbar, 1dbar = 10000 Pa = 1e4 kg m-1 s-2
          const double pressure = cc.rhoi * cc.earth_grav * (*ice_thickness)(i,j) * 1e-4;
          T_star(i,j) = cc.a*mean_salinity_in_boundary + cc.b - cc.c*pressure - mean_temperature_in_boundary;  // in Kelvin

          area_boxi = (counter_boxes[shelf_id][boxi] * dx * dy);

          double g1 = area_boxi*gamma_T;
          double g2 = g1 / (cc.nu*cc.lambda);

          // temperature for Box i > 1
          Toc(i,j) = mean_temperature_in_boundary + g1 * T_star(i,j)/(mean_overturning_in_box1 + g1 - g2*cc.a*mean_salinity_in_boundary); // K

          // salinity for Box i > 1
          Soc(i,j) = mean_salinity_in_boundary - mean_salinity_in_boundary * (mean_temperature_in_boundary - Toc(i,j))/(cc.nu*cc.lambda); // psu

          // potential pressure melting point 
          double potential_pressure_melting_point = cc.a*Soc(i,j) + cc.b - cc.c*pressure;

          // basal melt rate for Box i > 1
          basalmeltrate_shelf(i,j) = (-gamma_T/(cc.nu*cc.lambda)) * (potential_pressure_melting_point - Toc(i,j)); // in m/s

          // in situ pressure melting point in Kelvin
          T_pressure_melting(i,j) = cc.as*Soc(i,j) + cc.bs - cc.cs*pressure;

          lcounter_edge_of_boxi_vector[shelf_id]++;
          lmean_salinity_boxi_vector[shelf_id] += Soc(i,j);
          lmean_temperature_boxi_vector[shelf_id] += Toc(i,j);
          
        }
      } // no else-case, since calculate_basal_melt_box1() and calculate_basal_melt_missing_cells() cover all other cases 
    }

    // average the temperature, salinity in boxi 
    // (here we divide)
    for (int shelf_id=0;shelf_id<numberOfShelves;shelf_id++) {
      double counter_edge_of_boxi_vector=0.0;
      counter_edge_of_boxi_vector = GlobalSum(m_grid->com, lcounter_edge_of_boxi_vector[shelf_id]);
      mean_salinity_boundary_vector[shelf_id] = GlobalSum(m_grid->com, lmean_salinity_boxi_vector[shelf_id]);
      mean_temperature_boundary_vector[shelf_id] = GlobalSum(m_grid->com, lmean_temperature_boxi_vector[shelf_id]); // in Kelvin

      if (counter_edge_of_boxi_vector>0.0){
        mean_salinity_boundary_vector[shelf_id] = mean_salinity_boundary_vector[shelf_id]/counter_edge_of_boxi_vector;
        mean_temperature_boundary_vector[shelf_id] = mean_temperature_boundary_vector[shelf_id]/counter_edge_of_boxi_vector; // in Kelvin
      } else {
        mean_salinity_boundary_vector[shelf_id]=0.0; mean_temperature_boundary_vector[shelf_id]=0.0;
      }

      m_log->message(5, "  %d: cnt=%.0f, sal=%.3f, temp=%.3f, over=%.1e \n", shelf_id,counter_edge_of_boxi_vector,
                        mean_salinity_boundary_vector[shelf_id],mean_temperature_boundary_vector[shelf_id],
                        mean_overturning_box1_vector[shelf_id]) ;
    } // shelves

    countGl0 = GlobalSum(m_grid->com, lcountGl0);
    if (countGl0 > 0) {
      m_log->message(2, "PICO ocean WARNING: box %d, no boundary data from previous box in %.0f case(s)!\n"
                        "switching to Beckmann Goosse (2003) meltrate calculation\n",
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

void Cavity::calculate_basal_melt_missing_cells(const Constants &cc) {

  m_log->message(5, "starting calculate_basal_melt_missing_cells routine\n");

  const IceModelVec2S *ice_thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  IceModelVec::AccessList list;
  list.add(*ice_thickness);
  list.add(ocean_box_mask);
  list.add(Toc_box0);
  list.add(Toc);
  list.add(Soc_box0);
  list.add(Soc);
  list.add(overturning);
  list.add(basalmeltrate_shelf); // in m/s
  list.add(T_pressure_melting);
  list.add(shelf_mask);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    int shelf_id = shelf_mask(i,j);

    if (shelf_id == 0) {
      basalmeltrate_shelf(i,j) = 0.0;
    }

    // cell with missing data identifier numberOfBoxes+1, as set in routines before
    if ( (shelf_id > 0) && (ocean_box_mask(i,j)==(numberOfBoxes+1)) ) {

      Toc(i,j) = Toc_box0(i,j); // in Kelvin
      Soc(i,j) = Soc_box0(i,j); // in psu

      // in dbar, 1dbar = 10000 Pa = 1e4 kg m-1 s-2
      const double pressure = cc.rhoi * cc.earth_grav * (*ice_thickness)(i,j) * 1e-4;

      // potential pressure melting point 
      double potential_pressure_melting_point = cc.a*Soc(i,j) + cc.b - cc.c*pressure;

      double heatflux = cc.meltFactor * cc.rhow * cc.c_p_ocean * cc.default_gamma_T *
                         (Toc(i,j) - potential_pressure_melting_point);  // in W/m^2

      basalmeltrate_shelf(i,j) = heatflux / (cc.latentHeat * cc.rhoi); // in m s-1

      // in situ pressure melting point in Kelvin
      T_pressure_melting(i,j) =  cc.as*Soc(i,j) + cc.bs - cc.cs*pressure;

    }

  }

}

} // end of namespace ocean
} // end of namespace pism
