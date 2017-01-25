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

  a          = -0.057;       // K/psu
  b          = 0.0832 + 273.15;       // K
  c          = 7.64e-4;      // K/dbar

  alpha      = 7.5e-5;       // 1/K
  beta       = 7.7e-4;       // 1/psu

  default_gamma_T    = 1e-6;        // m/s FIXME check!
  default_overturning_coeff    = 5e6;         // kg−1 s−1 FIXME check!

  // for shelf cells where normal box model is not calculated,
  // used in basalMeltRateMissingCells(), compare POConstantPIK
  gamma_T_o    = 1.0e-4; //config.get("gamma_T"); //1e-4;
  // m/s, thermal exchange velocity for Beckmann-Goose parameterization
  meltFactor   = 0.002;     // FIXME add to pism_config, check value
  meltSalinity = 35.0;
  b2           = 0.0939;

}


const int Cavity::box1           = 1;       // ocean box covering the grounding line region
const int Cavity::box2           = 2;       // ocean box neighboring the box 1
// other boxes are covered by boxi

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
  BOXMODELmask.create(m_grid, "BOXMODELmask", WITH_GHOSTS);
  BOXMODELmask.set_attrs("model_state", "mask displaying ocean box model grid","", "");
  m_variables.push_back(&BOXMODELmask);

  // mask to identify the ice rises
  ICERISESmask.create(m_grid, "ICERISESmask", WITH_GHOSTS);
  ICERISESmask.set_attrs("model_state", "mask displaying ice rises","", "");
  m_variables.push_back(&ICERISESmask);

  // mask displaying continental shelf - region where mean salinity and ocean temperature is calculated
  OCEANMEANmask.create(m_grid, "OCEANMEANmask", WITH_GHOSTS);
  OCEANMEANmask.set_attrs("model_state", "mask displaying ocean region for parameter input","", "");
  m_variables.push_back(&OCEANMEANmask);

  // mask displaying open ocean - ice-free regions below sea-level except 'holes' in ice shelves
  OCEANmask.create(m_grid, "OCEANmask", WITH_GHOSTS);
  OCEANmask.set_attrs("model_state", "mask displaying open ocean","", "");
  m_variables.push_back(&OCEANmask);

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
  Soc_base.create(m_grid, "Soc_base", WITHOUT_GHOSTS);
  Soc_base.set_attrs("model_state", "ocean base salinity field","", "ocean base salinity field");  //NOTE unit=psu
  m_variables.push_back(&Soc_base);

  // computed temperature in ocean boxes
  Toc.create(m_grid, "Toc", WITHOUT_GHOSTS);
  Toc.set_attrs("model_state", "ocean temperature field","K", "ocean temperature field");
  m_variables.push_back(&Toc);

  // temperature input for box 1
  Toc_base.create(m_grid, "Toc_base", WITHOUT_GHOSTS);
  Toc_base.set_attrs("model_state", "ocean base temperature","K", "ocean base temperature");
  m_variables.push_back(&Toc_base);

  // in ocean box i: T_star = aS_{i-1} + b -c p_i - T_{i-1} with T_{-1} = Toc_base and S_{-1}=Soc_base
  // FIXME convert to internal field
  T_star.create(m_grid, "T_star", WITHOUT_GHOSTS);
  T_star.set_attrs("model_state", "T_star field","degree C", "T_star field");
  m_variables.push_back(&T_star);

  overturning.create(m_grid, "overturning", WITHOUT_GHOSTS);
  overturning.set_attrs("model_state", "cavity overturning","m^3 s-1", "cavity overturning"); // no CF standard_name?
  m_variables.push_back(&overturning);

  basalmeltrate_shelf.create(m_grid, "basal melt rate from ocean box model", WITHOUT_GHOSTS);
  basalmeltrate_shelf.set_attrs("climate_state", "basal melt rate from ocean box model", "m/s", "");
  //FIXME unit in field is kg m-2 a-1, but the written unit is m per a
  basalmeltrate_shelf.metadata().set_string("glaciological_units", "m year-1");
  m_variables.push_back(&basalmeltrate_shelf);

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

  Range basins_range = cbasins.range();

  if (basins_range.min < 0 or basins_range.max > numberOfBasins - 1) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "Some basin numbers in %s read from %s are invalid:"
                                  "allowed range is [0, %d], found [%d, %d]",
                                  cbasins.get_name().c_str(), m_filename.c_str(),
                                  numberOfBasins - 1,
                                  (int)basins_range.min, (int)basins_range.max);
  }

  // read time-independent data right away:
  if (m_theta_ocean->get_n_records() == 1 &&
      m_salinity_ocean->get_n_records() == 1) {
        update(m_grid->ctx()->time()->current(), 0); // dt is irrelevant
  }

}

void Cavity::add_vars_to_output_impl(const std::string &keyword, std::set<std::string> &result) {
  PGivenClimate<OceanModifier,OceanModel>::add_vars_to_output_impl(keyword, result);

  // add variable here and in define_variables_impl
  // if you want it to appear in snapshots
  if (keyword != "none" && keyword != "small") {
    result.insert(m_shelfbtemp.get_name());
    result.insert(m_shelfbmassflux.get_name());
  }
}

void Cavity::define_variables_impl(const std::set<std::string> &vars,
                                           const PIO &nc, IO_Type nctype) {

  PGivenClimate<OceanModifier,OceanModel>::define_variables_impl(vars, nc, nctype);

  for (unsigned int k = 0; k < m_variables.size(); ++k) {
    IceModelVec *v = m_variables[k];
    std::string name = v->metadata().get_string("short_name");
    if (set_contains(vars, name)) {
      v->define(nc, nctype);
    }
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

void Cavity::write_variables_impl(const std::set<std::string> &vars, const PIO& nc) {

  PGivenClimate<OceanModifier,OceanModel>::write_variables_impl(vars, nc);

  for (unsigned int k = 0; k < m_variables.size(); ++k) {
    IceModelVec *v = m_variables[k];
    std::string name = v->metadata().get_string("short_name");
    if (set_contains(vars, name)) {
      v->write(nc);
    }
  }

}

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

  Toc_base_vec.resize(numberOfBasins);
  Soc_base_vec.resize(numberOfBasins);

  counter_boxes.resize(numberOfBasins, std::vector<double>(2,0));
  mean_salinity_boundary_vector.resize(numberOfBasins);
  mean_temperature_boundary_vector.resize(numberOfBasins);
  mean_overturning_GLbox_vector.resize(numberOfBasins);

  gamma_T = cc.default_gamma_T;
  gamma_T = options::Real("-gamma_T","gamma_T for ocean cavity model",gamma_T);

  overturning_coeff = cc.default_overturning_coeff;
  overturning_coeff = options::Real("-overturning_coeff",
                                    "overturning_coeff for ocean cavity model",overturning_coeff);

  m_log->message(2, "     Using %d drainage basins and values: \n"
                                "     gamma_T= %.2e, overturning_coeff = %.2e... \n"
                                 , numberOfBasins, gamma_T, overturning_coeff);

  continental_shelf_depth = cc.continental_shelf_depth;
  options::Real cont_shelf_depth("-continental_shelf_depth",
                                 "continental shelf depth for ocean cavity model",
                                 continental_shelf_depth);

  if (cont_shelf_depth.is_set()) {
    m_log->message(2,
    "  Depth of continental shelf for computation of temperature and salinity input\n"
    "  is set for whole domain to continental_shelf_depth=%.0f meter\n",
    continental_shelf_depth);
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
  identifyMASK(OCEANMEANmask,"ocean");
  computeOCEANMEANS(cc);

  //geometry of ice shelves and input values of temperature and salinity
  if (exicerises_set) {
    identifyMASK(ICERISESmask,"icerises");}

  identifyMASK(OCEANmask,"openocean");
  extentOfIceShelves();
  identifyBOXMODELmask();
  oceanTemperature(cc);
  m_shelfbtemp.copy_from(Toc);

  //basal melt rates underneath ice shelves
  basalMeltRateGroundingLineBox(cc);
  basalMeltRateOtherBoxes(cc); // TODO Diese Routinen woanders aufrufen (um Dopplung zu vermeiden)
  basalMeltRateMissingCells(cc);  //Assumes that mass flux is proportional to the shelf-base heat flux.
  //const double secpera=31556926.0;
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

//! Round basin mask non integer values to an integral value of the next neighbor
void Cavity::round_basins() {

  //FIXME: THIS routine should be applied once in init, and roundbasins should be stored as field (assumed the basins do not change with time).

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

//! Identify
//!   ocean:    identify ocean up to continental shelf without detached submarine islands regions
//!   icerises: identify grounded regions without detached ice rises
//!   openocean: identify ocean without holes in ice shelves

void Cavity::identifyMASK(IceModelVec2S &inputmask, std::string masktype) {

  m_log->message(5, "starting identifyMASK routine\n");

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

  if ((masktype=="ocean" || masktype=="icerises") && (seed_x >= m_grid->xs()) && (seed_x < m_grid->xs()+m_grid->xm()) && (seed_y >= m_grid->ys())&& (seed_y < m_grid->ys()+m_grid->ym())){
    inputmask(seed_x,seed_y)=imask_inner;
  }
  else if (masktype=="openocean"){
    //only consider domain bounradies:
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();
      if ((i==0) | (j==0) | (i>(Mx-2)) | (j>(My-2))){
        inputmask(i,j)=imask_inner;
      }
    }
  }

  int iteration_round = 0;
  // find inner region first
  while(all_inner_identified > previous_step_identified){

    iteration_round+=1;
    previous_step_identified = all_inner_identified;

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

      bool masktype_condition = false;
      if (masktype=="ocean"){
        masktype_condition = (m_mask(i,j)!=maskocean || (*topg)(i,j) >= continental_shelf_depth);}
      else if (masktype=="icerises"){
        masktype_condition = (m_mask(i,j)==maskgrounded);
      }
      else if (masktype=="openocean"){
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

  // TODO: Not sure if we have to reinitialize m_mask and inputmask here.
  // set value for excluded areas (ice rises or submarine islands)

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (inputmask(i,j)==imask_unidentified){
      inputmask(i,j)=imask_exclude;
    }

    if (masktype=="ocean"){ //exclude ice covered parts
      if (m_mask(i,j)!=maskocean && inputmask(i,j) == imask_inner){
        inputmask(i,j) = imask_outer;
      }
    }

  }

}


//! When ocean_given is set compute mean salinity and temperature in each basin.
void Cavity::computeOCEANMEANS(const Constants &cc) {

  m_log->message(5, "starting computeOCEANMEANS routine \n");

  std::vector<double> lm_count(numberOfBasins); //count cells to take mean over for each basin
  std::vector<double> m_count(numberOfBasins);
  std::vector<double> lm_Sval(numberOfBasins); //add salinity for each basin
  std::vector<double> lm_Tval(numberOfBasins); //add temperature for each basin
  std::vector<double> m_Tval(numberOfBasins);
  std::vector<double> m_Sval(numberOfBasins);

  for(int k=0;k<numberOfBasins;k++){
    m_count[k]=0.0;
    lm_count[k]=0.0;
    lm_Sval[k]=0.0;
    lm_Tval[k]=0.0;
    m_Tval[k]=0.0;
    m_Sval[k]=0.0;
  }

  IceModelVec::AccessList list;
  list.add(*m_theta_ocean);
  list.add(*m_salinity_ocean);
  list.add(cbasins);
  list.add(OCEANMEANmask);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (OCEANMEANmask(i,j) == imask_inner ){
      int shelf_id =(cbasins)(i,j);
      lm_count[shelf_id]+=1;
      lm_Sval[shelf_id]+=(*m_salinity_ocean)(i,j);
      lm_Tval[shelf_id]+=(*m_theta_ocean)(i,j);
    } //if

  }


  for(int k=0;k<numberOfBasins;k++) {

    m_count[k] = GlobalSum(m_grid->com, lm_count[k]);
    m_Sval[k] = GlobalSum(m_grid->com, lm_Sval[k]);
    m_Tval[k] = GlobalSum(m_grid->com, lm_Tval[k]);

    //if basin is not dummy basin 0 or there are no ocean cells in this basin to take the mean over.
    //FIXME: the following error occurs once at initialization before input is available.
    // Please ignore this very first warning for now.
    if(k>0 && m_count[k]==0){
      m_log->message(2, "SIMPEL ocean WARNING: basin %d contains no cells with OCEANMEANmask=2.\n"
                        "no mean salinity or temperature values are computed, using\n"
                        "the standard values T_dummy =%.3f, S_dummy=%.3f.\n", k, cc.T_dummy, cc.S_dummy);
      Toc_base_vec[k] = cc.T_dummy;
      Soc_base_vec[k] = cc.S_dummy;
    } else {
      m_Sval[k] = m_Sval[k] / m_count[k];
      m_Tval[k] = m_Tval[k] / m_count[k];

      Toc_base_vec[k]=m_Tval[k];
      Soc_base_vec[k]=m_Sval[k];
      m_log->message(5, "  %d: temp =%.3f, salinity=%.3f\n", k, Toc_base_vec[k], Soc_base_vec[k]);
    }
  }

}


//! Compute the extent of the ice shelves of each basin/region (i.e. counter) and
//  compute for each ice shelf cell the distance to the grounding line (i.e. DistGL) and the calving front (i.e. DistIF)


void Cavity::extentOfIceShelves() {

  m_log->message(5, "starting extentOfIceShelves routine\n");

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
  list.add(OCEANmask);

	if (exicerises_set) { list.add(ICERISESmask); }

	DistGL.set(0);
	DistIF.set(0);

	// find the grounding line and the ice front
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    bool condition;
    if (exicerises_set) {
      condition = (m_mask(i,j)==maskfloating || ICERISESmask(i,j)==imask_exclude || OCEANmask(i,j)==imask_exclude);
    }
		else {
      condition = (m_mask(i,j)==maskfloating || OCEANmask(i,j)==imask_exclude);
    }

    if (condition) { //if this is a ice shelf cell (or an ice rise) or a hole in an ice shelf

			// label the shelf cells adjacent to the grounding line with DistGL = 1
			bool neighbor_to_land;
			if (exicerises_set) {
				neighbor_to_land = (  ICERISESmask(i,j+1)==imask_inner || ICERISESmask(i,j-1)==imask_inner ||
					ICERISESmask(i+1,j)==imask_inner || ICERISESmask(i-1,j)==imask_inner ||
 					ICERISESmask(i+1,j+1)==imask_inner || ICERISESmask(i+1,j-1)==imask_inner ||
 					ICERISESmask(i-1,j+1)==imask_inner || ICERISESmask(i-1,j-1)==imask_inner );
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
      neighbor_to_ocean = (OCEANmask(i,j+1)==imask_inner || OCEANmask(i,j-1)==imask_inner || OCEANmask(i+1,j)==imask_inner || OCEANmask(i-1,j)==imask_inner);
      //(m_mask(i,j+1)==maskocean || m_mask(i,j-1)== maskocean || m_mask(i+1,j)==maskocean || m_mask(i-1,j)==maskocean)
			if (neighbor_to_ocean) {
				DistIF(i,j) = currentLabelIF;
			}// no else

		}
	}

	DistGL.update_ghosts();
	DistIF.update_ghosts();

  // Find DistGL for all shelf cells
  // FIXME: Do we want to take compute DistGL using four direct neigbors or
  //        also diagonal-neighbor (some points might not be reached otherwise)?

  global_continue_loop = 1;
  while( global_continue_loop !=0 ) {

    local_continue_loop = 0;

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      bool condition; // this cell is floating or an hole in the ice shelf (or an ice rise)
      if (exicerises_set) {
        condition = (m_mask(i,j)==maskfloating || ICERISESmask(i,j)==imask_exclude || OCEANmask(i,j)==imask_exclude);
      }
      else {
        condition = (m_mask(i,j)==maskfloating || OCEANmask(i,j)==imask_exclude);
      }

      if ( condition && DistGL(i,j)==0 &&
        (DistGL(i,j+1)==currentLabelGL || DistGL(i,j-1)==currentLabelGL ||
        DistGL(i+1,j)==currentLabelGL || DistGL(i-1,j)==currentLabelGL // ||
        //DistGL(i+1,j+1)==currentLabelGL || DistGL(i+1,j-1)==currentLabelGL ||
        //DistGL(i-1,j+1)==currentLabelGL || DistGL(i-1,j-1)==currentLabelGL
        ) ) { // i.e. this is an shelf cell with no distance assigned yet and with a neighbor that has a distance assigned
          DistGL(i,j) = currentLabelGL+1;
          local_continue_loop = 1;
      } //if

    } // for

    currentLabelGL++;
    DistGL.update_ghosts();

    global_continue_loop = GlobalMax(m_grid->com, local_continue_loop);

  } // while: find DistGL

  // Find DistIF for all shelf cells
  // FIXME: Do we want to take compute DistIF using four direct neigbors or
  //        also diagonal-neighbor (some points might not be reached otherwise)?


  global_continue_loop = 1; // start loop
  while( global_continue_loop !=0  ) {

    local_continue_loop = 0;

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      bool condition; // this cell is floating or an hole in the ice shelf (or an ice rise)
      if (exicerises_set) {
        condition = (m_mask(i,j)==maskfloating || ICERISESmask(i,j)==imask_exclude || OCEANmask(i,j)==imask_exclude);
      }
      else {
        condition = (m_mask(i,j)==maskfloating || OCEANmask(i,j)==imask_exclude);
      }

      if ( condition && DistIF(i,j)==0 &&
        (DistIF(i,j+1)==currentLabelIF || DistIF(i,j-1)==currentLabelIF ||
        DistIF(i+1,j)==currentLabelIF || DistIF(i-1,j)==currentLabelIF // ||
        //DistIF(i+1,j+1)==currentLabelIF || DistIF(i+1,j-1)==currentLabelIF ||
        //DistIF(i-1,j+1)==currentLabelIF || DistIF(i-1,j-1)==currentLabelIF
        ) ) { // i.e. this is an shelf cell with no distance assigned yet and with a neighbor that has a distance assigned
          DistIF(i,j)=currentLabelIF+1;
          local_continue_loop = 1;
      } //if

    } // for


    currentLabelIF++;
    DistIF.update_ghosts();

    global_continue_loop = GlobalMax(m_grid->com, local_continue_loop);

  } // while: find DistIF

}


//! Compute the BOXMODELmask based on DistGL and DistIF, calculate the extent of each box in each region

void Cavity::identifyBOXMODELmask() {

  m_log->message(5, "starting identifyBOXMODELmask routine\n");

  // Find the maximal DistGL and DistIF
  // FIXME! this could already be done in routine where DistGL and DistIF are computed
  std::vector<double> max_distGL(numberOfBasins);
  std::vector<double> max_distIF(numberOfBasins);
  std::vector<double> lmax_distGL(numberOfBasins);
  std::vector<double> lmax_distIF(numberOfBasins);

  const IceModelVec2CellType &m_mask = *m_grid->variables().get_2d_cell_type("mask");

  for(int k=0;k<numberOfBasins;k++){ max_distGL[k]=0.0; max_distIF[k]=0.0;lmax_distGL[k]=0.0; lmax_distIF[k]=0.0;}

  IceModelVec::AccessList list;
  list.add(cbasins);
  list.add(DistGL);
  list.add(DistIF);
  list.add(BOXMODELmask);
  list.add(m_mask);

  for (Points p(*m_grid); p; p.next()) {
  const int i = p.i(), j = p.j();
    int shelf_id = (cbasins)(i,j);

    if ( DistGL(i,j)> lmax_distGL[shelf_id] ) {
      lmax_distGL[shelf_id] = DistGL(i,j);
    } //if
    if ( DistIF(i,j)> lmax_distIF[shelf_id] ) {
      lmax_distIF[shelf_id] = DistIF(i,j);
    } //if
  } // for


  for (int l=0;l<numberOfBasins;l++){
    max_distGL[l] = GlobalMax(m_grid->com, lmax_distGL[l]);
    max_distIF[l] = GlobalMax(m_grid->com, lmax_distIF[l]);
  }



  // Define the number of boxes for each basin
  std::vector<int> lnumberOfBoxes_perBasin(numberOfBasins);

  int n_min = 1; //
  double max_distGL_ref = 500000; // meter
  double zeta = 0.5;

  for (int l=0;l<numberOfBasins;l++){
    lnumberOfBoxes_perBasin[l] = 0;
    //ATTENTION, this is only correct for same dx and dy spacing.
    // Otherwise, we need to change the calculation of DistGL and DistIF
    lnumberOfBoxes_perBasin[l] = n_min + static_cast<int>(
        round(pow((max_distGL[l]*dx/max_distGL_ref), zeta) *(numberOfBoxes-n_min)));
    m_log->message(5, "lnumberOfBoxes[%d]=%d \n", l, lnumberOfBoxes_perBasin[l]);
  }

  // Define the BOXMODELmask

  // IceModelVec::AccessList list;
  // list.add(cbasins);
  // list.add(DistGL);
  // list.add(DistIF);


  BOXMODELmask.set(0);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m_mask(i,j)==maskfloating && DistGL(i,j)>0 && DistIF(i,j)>0 && BOXMODELmask(i,j)==0){
      int shelf_id = (cbasins)(i,j);
      int n = lnumberOfBoxes_perBasin[shelf_id];
      double r = DistGL(i,j)*1.0/(DistGL(i,j)*1.0+DistIF(i,j)*1.0); // relative distance between grounding line and ice front

      for(int k=0;k<n;++k){

        // First variant to define the BOXMODELmask using a rule like k/n< (1-r)**2 <k+1/n
        // this rule is motivated by splitting a half-circle into halfcircles of same area and using 1-r like some kind of radius

        if (  ((n*1.0-k*1.0-1.0)/(n*1.0) <= pow((1.0-r),2)) && (pow((1.0-r), 2) <= (n*1.0-k*1.0)/n*1.0) ){ // FIXME do we need to multiply by 1.0 here?
          if (DistGL(i,j) < k+1) {
            BOXMODELmask(i,j) = DistGL(i,j); // the boxnumber of a cell cannot be bigger then the distance to the grounding line //FIXME Discuss!!!
          } else{
          BOXMODELmask(i,j) = k+1;
          }
        }//if
        /*
        // Second variant to define the BOXMODELmask using a rule like k/n < r**0.5 < k+1/n
        if (  ((k*1.0)/(n*1.0) <= pow(r,0.5)) && (pow(r, 0.5) <= (k*1.0+1.0)/n*1.0) ){ // FIXME do we need to multiply by 1.0 here?
          if (DistGL(i,j) < k+1) {
            BOXMODELmask(i,j) = DistGL(i,j); // the boxnumber of a cell cannot be bigger then the distance to the grounding line //FIXME Discuss!!!
          } else{
          BOXMODELmask(i,j) = k+1;
          }
        }//if */

      } //for
    }
  } // for


  // set all floating cells which have no BOXMODELmask value as numberOfBoxes+1 -> beckmann-goose for melting
  // those are the cells which are not reachable from GL or IF //FIXME does that make sense?


  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    if (m_mask(i,j)==maskfloating && BOXMODELmask(i,j)==0){ // floating
      BOXMODELmask(i,j) = numberOfBoxes + 1;
    }

  }

  // Compute the number of cells per box and basin. Later: Include this in the loop above to save time...
  const int nBoxes = numberOfBoxes+2;
  std::vector<std::vector<int> > lcounter_boxes(
    numberOfBasins, std::vector<int>(nBoxes));
  for (int k=0;k<numberOfBasins;k++){
    for (int l=0;l<nBoxes;l++){
      lcounter_boxes[k][l]=0;
    }
  }

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    int box_id = static_cast<int>(round(BOXMODELmask(i,j)));
    if (box_id > 0){ // floating
      int shelf_id = (cbasins)(i,j);
      lcounter_boxes[shelf_id][box_id]++;
    }
  }

  for (int k=0;k<numberOfBasins;k++){
    counter_boxes[k].resize(nBoxes);
    for (int l=0;l<nBoxes;l++){
      counter_boxes[k][l] = GlobalSum(m_grid->com, lcounter_boxes[k][l]);
    }
  }

}



/*!
Set ocean temperature in Box 0, used as boundary condition for Box 1
*/


void Cavity::oceanTemperature(const Constants &cc) {

  m_log->message(5, "starting oceanTemperature routine\n");

  const IceModelVec2S *ice_thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");
  const IceModelVec2CellType &m_mask = *m_grid->variables().get_2d_cell_type("mask");

  IceModelVec::AccessList list;
  list.add(*ice_thickness);
  list.add(cbasins);
  list.add(Soc_base);
  list.add(Toc_base);
  list.add(Toc);
  list.add(m_mask);

  double counterTpmp=0.0,
         lcounterTpmp = 0.0;

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // make sure all temperatures are zero at the beginning of each timestep
    Toc(i,j) = 273.15; // in K
    Toc_base(i,j) = 273.15;  // in K
    Soc_base(i,j) = 0.0; // in psu


    if (m_mask(i,j)==maskfloating){
      int shelf_id = (cbasins)(i,j);
      Toc_base(i,j) = Toc_base_vec[shelf_id];
      Soc_base(i,j) =  Soc_base_vec[shelf_id];

      //! salinity and temperature for grounding line box
      if ( Soc_base(i,j) == 0.0 || Toc_base_vec[shelf_id] == 273.15 ) { //FIXME is there a reason that Toc and Soc are different?
        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "PISM_ERROR: Missing Soc_base and Toc_base for"
                                      "%d, %d, basin %d \n   Aborting... \n", i, j, shelf_id);
      }


      //! temperature input for grounding line box should not be below pressure melting point
      const double pressure = cc.rhoi * cc.earth_grav * (*ice_thickness)(i,j) * 1e-4; // MUST be in dbar  // NOTE 1dbar = 10000 Pa = 1e4 kg m-1 s-2,
      const double T_pmt = cc.a*Soc_base(i,j) + cc.b - cc.c*pressure; // in Kelvin


      if (  Toc_base(i,j) < T_pmt ) {
        //FIXME: only works in serial runs
        m_log->message(2, "SIMPEL ocean WARNING: Toc_base is below the local pressure melting temperature\n"
                      "for %d, %d, basin %d, setting it to pressure melting point \n",i,j,shelf_id);
        Toc_base(i,j) = T_pmt ;
        lcounterTpmp+=1;
      }

    } // end if herefloating
  } // end i
    counterTpmp = GlobalSum(m_grid->com, lcounterTpmp);
    if (counterTpmp > 0) {
      m_log->message(2, "SIMPEL ocean warning: temperature has been below pressure melting temperature in %.0f cases,\n"
                        "setting it to pressure melting temperature\n", counterTpmp);
    }

}


//! Compute the basal melt / refreezing rates for each shelf cell bordering the grounding line box
void Cavity::basalMeltRateGroundingLineBox(const Constants &cc) {

  m_log->message(5, "starting basal basalMeltRateGroundingLineBox routine\n");

  std::vector<double> lcounter_edge_of_GLbox_vector(numberOfBasins);
  std::vector<double> lmean_salinity_GLbox_vector(numberOfBasins);
  std::vector<double> lmean_temperature_GLbox_vector(numberOfBasins);
  std::vector<double> lmean_meltrate_GLbox_vector(numberOfBasins);
  std::vector<double> lmean_overturning_GLbox_vector(numberOfBasins);

  for (int k=0;k<numberOfBasins;k++){
    lcounter_edge_of_GLbox_vector[k]=0.0;
    lmean_salinity_GLbox_vector[k]=0.0;
    lmean_temperature_GLbox_vector[k]=0.0;
    lmean_meltrate_GLbox_vector[k]=0.0;
    lmean_overturning_GLbox_vector[k]=0.0;
  }

  const IceModelVec2S *ice_thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  IceModelVec::AccessList list;
  list.add(*ice_thickness);
  list.add(cbasins);
  list.add(BOXMODELmask);
  list.add(T_star);
  list.add(Toc_base);
  list.add(Toc);
  list.add(Soc_base);
  list.add(Soc);
  list.add(overturning);
  list.add(basalmeltrate_shelf);


  double countHelpterm=0,
         lcountHelpterm=0;


  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    int shelf_id = (cbasins)(i,j);

    // Make sure everything is at default values at the beginning of each timestep
    T_star(i,j) = 0.0; // in Kelvin
    Toc(i,j) = 273.15; // in Kelvin
    Soc(i,j) = 0.0; // in psu

    basalmeltrate_shelf(i,j) = 0.0;
    overturning(i,j) = 0.0;


    if ((BOXMODELmask(i,j) == box1) && (shelf_id > 0.0)){

      const double pressure = cc.rhoi * cc.earth_grav * (*ice_thickness)(i,j) * 1e-4; // MUST be in dbar  // NOTE 1dbar = 10000 Pa = 1e4 kg m-1 s-2
      // FIXME need to include atmospheric pressure?
      T_star(i,j) = cc.a*Soc_base(i,j) + cc.b - cc.c*pressure - Toc_base(i,j); // in Kelvin

      double g1 = (counter_boxes[shelf_id][box1] * dx * dy) * gamma_T / (overturning_coeff*cc.rho_star);

      double helpterm1 = g1/(cc.beta*(Soc_base(i,j) / (cc.nu*cc.lambda)) - cc.alpha); // in 1 / (1/K) = K
      double helpterm2 = (g1*T_star(i,j)) / (cc.beta*(Soc_base(i,j) / (cc.nu*cc.lambda)) - cc.alpha); // in K / (1/K) = K^2



      // This can only happen if T_star > 0.25 helpterm1, in particular T_star > 0
      // which can only happen for very small values of Toc_base
      if ((0.25*PetscSqr(helpterm1) - helpterm2) < 0.0) {

        m_log->message(5,
          "SIMPEL ocean WARNING: negative square root argument at %d, %d\n"
          "probably because of positive T_star=%f \n"
          "Not aborting, but setting square root to 0... \n",
          i, j, T_star(i,j));

        helpterm2=0.25*PetscSqr(helpterm1);

        lcountHelpterm+=1;
      }

      //! temperature for grounding line box
      Toc(i,j) = Toc_base(i,j) - ( -0.5*helpterm1 + sqrt(0.25*PetscSqr(helpterm1) -helpterm2) ); // in Kelvin
      //! salinity for grounding line box
      Soc(i,j) = Soc_base(i,j) - (Soc_base(i,j) / (cc.nu*cc.lambda)) * (Toc_base(i,j) - Toc(i,j));  // in psu

      //! basal melt rate for grounding line box
      basalmeltrate_shelf(i,j) = (-gamma_T/(cc.nu*cc.lambda)) * (cc.a*Soc(i,j) + cc.b - cc.c*pressure - Toc(i,j));  // in m/s

      //! overturning
      // NOTE Actually, there is of course no overturning-FIELD, it is only a scalar for each shelf.
      // Here, we compute overturning as   MEAN[C1*cc.rho_star* (cc.beta*(Soc_base(i,j)-Soc(i,j)) - cc.alpha*((Toc_base(i,j)-273.15+Toc_anomaly(i,j))-Toc_inCelsius(i,j)))]
      // while in fact it should be   C1*cc.rho_star* (cc.beta*(Soc_base-MEAN[Soc(i,j)]) - cc.alpha*((Toc_base-273.15+Toc_anomaly)-MEAN[Toc_inCelsius(i,j)]))
      // which is the SAME since Soc_base, Toc_base and Toc_anomaly are the same FOR ALL i,j CONSIDERED, so this is just nomenclature!
      overturning(i,j) = overturning_coeff*cc.rho_star* (cc.beta*(Soc_base(i,j)-Soc(i,j)) - cc.alpha*(Toc_base(i,j)-Toc(i,j))); // in m^3/s

      if (BOXMODELmask(i-1,j)==box2 || BOXMODELmask(i+1,j)==box2 || BOXMODELmask(i,j-1)==box2 || BOXMODELmask(i,j+1)==box2){
      // i.e., if this cell is from the GL box and one of the neighbours is from the CF box - It is important to only take the border of the grounding line box
      // to the calving front box into account, because the following mean value will be used to compute the value for the calving front box. I.e., this helps avoiding discontinuities!
        lcounter_edge_of_GLbox_vector[shelf_id]++;
        lmean_salinity_GLbox_vector[shelf_id] += Soc(i,j);
        lmean_temperature_GLbox_vector[shelf_id] += Toc(i,j); // in Kelvin
        lmean_meltrate_GLbox_vector[shelf_id] += basalmeltrate_shelf(i,j);
        lmean_overturning_GLbox_vector[shelf_id] += overturning(i,j);

      } // no else-case necessary since all variables are set to zero at the beginning of this routine

    }else { // i.e., not GL_box
        basalmeltrate_shelf(i,j) = 0.0;
    }
  }

  for(int k=0;k<numberOfBasins;k++) {
    double counter_edge_of_GLbox_vector=0.0;

    counter_edge_of_GLbox_vector = GlobalSum(m_grid->com, lcounter_edge_of_GLbox_vector[k]);
    mean_salinity_boundary_vector[k] = GlobalSum(m_grid->com, lmean_salinity_GLbox_vector[k]);
    mean_temperature_boundary_vector[k] = GlobalSum(m_grid->com, lmean_temperature_GLbox_vector[k]);
    mean_overturning_GLbox_vector[k] = GlobalSum(m_grid->com, lmean_overturning_GLbox_vector[k]);


    if (counter_edge_of_GLbox_vector>0.0){
      mean_salinity_boundary_vector[k] = mean_salinity_boundary_vector[k]/counter_edge_of_GLbox_vector;
      mean_temperature_boundary_vector[k] = mean_temperature_boundary_vector[k]/counter_edge_of_GLbox_vector;
      mean_overturning_GLbox_vector[k] = mean_overturning_GLbox_vector[k]/counter_edge_of_GLbox_vector;
    } else { // This means that there is no [cell from the GLbox neighboring a cell from the CFbox], NOT necessarily that there is no GLbox!
      mean_salinity_boundary_vector[k]=0.0; mean_temperature_boundary_vector[k]=0.0; mean_overturning_GLbox_vector[k]=0.0;
    }

    m_log->message(5, "  %d: cnt=%.0f, sal=%.3f, temp=%.3f, over=%.1e \n", k,counter_edge_of_GLbox_vector,mean_salinity_boundary_vector[k],mean_temperature_boundary_vector[k],mean_overturning_GLbox_vector[k]) ;
  }

    countHelpterm = GlobalSum(m_grid->com, lcountHelpterm);
    if (countHelpterm > 0) {
      m_log->message(2, "SIMPEL ocean warning: square-root argument for temperature calculation "
                        "has been negative in %.0f cases!\n", countHelpterm);
    }

}



//! Iteratively over i=2,..,n: Compute the basal melt / refreezing rates for each shelf cell bordering the Box i

void Cavity::basalMeltRateOtherBoxes(const Constants &cc) { //FIXME rename routine!!

  m_log->message(5, "starting basalMeltRateOtherBoxes routine\n");

  int nBoxes = static_cast<int>(round(numberOfBoxes+1)); // do not include the Beckmann-Goose (=last) Box

  const IceModelVec2S *ice_thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  //! Iterate over all Boxes i for i > 1
  for (int boxi=2; boxi <nBoxes; ++boxi) {
    m_log->message(5, "computing basal melt rate, temperature and salinity for box i = %d \n", boxi);

    double countGl0=0,
           lcountGl0=0;

    std::vector<double> lcounter_edge_of_boxi_vector(numberOfBasins);     // to compute means at boundary for the current box
    std::vector<double> lmean_salinity_boxi_vector(numberOfBasins);
    std::vector<double> lmean_temperature_boxi_vector(numberOfBasins); // in Kelvin

    for (int k=0;k<numberOfBasins;k++){
      lcounter_edge_of_boxi_vector[k] =0.0;
      lmean_salinity_boxi_vector[k]   =0.0;
      lmean_temperature_boxi_vector[k]=0.0;
    }

    // TODO: does this need to be within the loop over boxes?
    // TODO: do we really need all these variables as full fields?
    IceModelVec::AccessList list;
    list.add(*ice_thickness);
    list.add(cbasins);
    list.add(BOXMODELmask);
    list.add(T_star);
    list.add(Toc_base);
    list.add(Toc);
    list.add(Soc_base);
    list.add(Soc);
    list.add(overturning);
    list.add(basalmeltrate_shelf);

    // for box i compute the melt rates.

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      int shelf_id = (cbasins)(i,j);

      if (BOXMODELmask(i,j)==boxi && shelf_id > 0.0){

        double  area_boxi,mean_salinity_in_boundary,mean_temperature_in_boundary,mean_overturning_in_GLbox;

        // FIXME RENAME THESE in GENERAL
        mean_salinity_in_boundary     = mean_salinity_boundary_vector[shelf_id];
        mean_temperature_in_boundary  = mean_temperature_boundary_vector[shelf_id]; // note: in Kelvin, mean over Toc
        mean_overturning_in_GLbox     = mean_overturning_GLbox_vector[shelf_id]; // !!!leave this one with the grounding line box

        const double pressure = cc.rhoi * cc.earth_grav * (*ice_thickness)(i,j) * 1e-4; // MUST be in dbar  // NOTE 1dbar = 10000 Pa = 1e4 kg m-1 s-2
        T_star(i,j) = cc.a*mean_salinity_in_boundary + cc.b - cc.c*pressure - mean_temperature_in_boundary;  // in °C or Kelvin since anomaly


        area_boxi = (counter_boxes[shelf_id][boxi] * dx * dy); //FIXME this assumes rectangular cell areas, adjust with real areas from projection


        if (mean_salinity_in_boundary==0 || mean_overturning_in_GLbox==0 || mean_temperature_in_boundary==273.15) {

          // if there are no boundary values from the box before
          m_log->message(2, "SIMPEL ocean WARNING: No input from box i-1 for box %d at %d,%d, basin=%d \n"
                            "switching to Beckmann Goose (2003) meltrate calculation\n",
                         boxi, i, j, shelf_id);
          BOXMODELmask(i,j) = numberOfBoxes+1;
          lcountGl0+=1;

        } else {
          // compute melt rates

          double g1, g2;
          g1 = (area_boxi*gamma_T);
          g2 = g1 / (cc.nu*cc.lambda);

          //! temperature for Box i > 1
          Toc(i,j) = mean_temperature_in_boundary + g1 * T_star(i,j)/(mean_overturning_in_GLbox + g1 - g2*cc.a*mean_salinity_in_boundary); // K

          //! salinity for Box i > 1
          Soc(i,j) = mean_salinity_in_boundary - mean_salinity_in_boundary * (mean_temperature_in_boundary - Toc(i,j))/(cc.nu*cc.lambda); // psu FIXME: add term in denominator? Then also above in GL box + mean_temperature_in_boundary - Toc_inCelsius(i,j))

          //! basal melt rate for Box i > 1
          basalmeltrate_shelf(i,j) = (-gamma_T/(cc.nu*cc.lambda)) * (cc.a*Soc(i,j) + cc.b - cc.c*pressure - Toc(i,j)); // in m/s


          // compute means at boundary to next box
          if (BOXMODELmask(i-1,j)==(boxi+1) || BOXMODELmask(i+1,j)==(boxi+1) || BOXMODELmask(i,j-1)==(boxi+1) || BOXMODELmask(i,j+1)==(boxi+1)){
            // i.e., if this cell is from the current Box and one of the neighbours is from the next higher box - It is important to only take the border of the current box
            // to the calving front box into account, because the following mean value will be used to compute the value for the calving front box. I.e., this helps avoiding discontinuities!
            lcounter_edge_of_boxi_vector[shelf_id]++;
            lmean_salinity_boxi_vector[shelf_id] += Soc(i,j);
            lmean_temperature_boxi_vector[shelf_id] += Toc(i,j);
          } // no else-case necessary since all variables are set to zero at the beginning of this routine
        }
      } // NOTE NO else-case, since  basalmeltRateGroundingLineBox() and basalMeltRateMissingCells() cover all other cases and we would overwrite those results here.
    }

    for(int k=0;k<numberOfBasins;k++) {
      // NOTE: overturning should not be changed!!!
      double counter_edge_of_boxi_vector=0.0;
      counter_edge_of_boxi_vector = GlobalSum(m_grid->com, lcounter_edge_of_boxi_vector[k]);
      mean_salinity_boundary_vector[k] = GlobalSum(m_grid->com, lmean_salinity_boxi_vector[k]);
      mean_temperature_boundary_vector[k] = GlobalSum(m_grid->com, lmean_temperature_boxi_vector[k]); // in Kelvin

      if (counter_edge_of_boxi_vector>0.0){
        mean_salinity_boundary_vector[k] = mean_salinity_boundary_vector[k]/counter_edge_of_boxi_vector;
        mean_temperature_boundary_vector[k] = mean_temperature_boundary_vector[k]/counter_edge_of_boxi_vector; // in Kelvin
      } else { // This means that there is no [cell from the GLbox neighboring a cell from the CFbox], NOT necessarily that there is no GLbox!
        mean_salinity_boundary_vector[k]=0.0; mean_temperature_boundary_vector[k]=0.0;
      }

      m_log->message(5, "  %d: cnt=%.0f, sal=%.3f, temp=%.3f, over=%.1e \n", k,counter_edge_of_boxi_vector,
                        mean_salinity_boundary_vector[k],mean_temperature_boundary_vector[k],
                        mean_overturning_GLbox_vector[k]) ;
    } // basins

    countGl0 = GlobalSum(m_grid->com, lcountGl0);
    if (countGl0 > 0) {
      m_log->message(2, "SIMPEL ocean WARNING: box %d, no boundary data from previous box in %.0f case(s)!\n"
                        "switching to Beckmann Goose (2003) meltrate calculation\n",
                        boxi,countGl0);
    }

  } // boxi

}


//! Compute the melt rate for all other ice shelves.
void Cavity::basalMeltRateMissingCells(const Constants &cc) {

  m_log->message(5, "starting basalMeltRateMissingCells routine\n");

  const IceModelVec2S *ice_thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  // TODO: do we really need all these variables as full fields?
  IceModelVec::AccessList list;
  list.add(*ice_thickness);
  list.add(cbasins);
  list.add(BOXMODELmask);
  list.add(Toc_base);
  list.add(Toc);
  list.add(overturning);
  list.add(basalmeltrate_shelf);  // NOTE meltrate has units:   J m-2 s-1 / (J kg-1 * kg m-3) = m s-1

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    int shelf_id = (cbasins)(i,j);

    if (shelf_id == 0) { // boundary of computational domain

      basalmeltrate_shelf(i,j) = 0.0;

    } else if (BOXMODELmask(i,j)==(numberOfBoxes+1) ) {

      Toc(i,j) = Toc_base(i,j); // in K, NOTE: Toc_base is already in K, so no (+273.15)
      // default: compute the melt rate from the temperature field according to beckmann_goosse03 (see below)


      const double shelfbaseelev = - (cc.rhoi / cc.rhow) * (*ice_thickness)(i,j);


      //FIXME: for consistency reasons there should be constants a,b,c, gamma_T used
      double T_f = 273.15 + (cc.a*cc.meltSalinity + cc.b2 + cc.c*shelfbaseelev); // add 273.15 to get it in Kelvin... 35 is the salinity

      double heatflux = cc.meltFactor * cc.rhow * cc.c_p_ocean * cc.gamma_T_o * (Toc(i,j) - T_f);  // in W/m^2
      basalmeltrate_shelf(i,j) = heatflux / (cc.latentHeat * cc.rhoi); // in m s-1

    } else if (shelf_id > 0.0) {
      continue; // standard case

    } else { // This must not happen

      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "PISM_ERROR: [rank %d] at %d, %d  -- basins(i,j)=%d causes problems.\n"
                                    "Aborting... \n",m_grid->rank(), i, j, shelf_id);
    }
  }

}

} // end of namespace ocean
} // end of namespace pism
