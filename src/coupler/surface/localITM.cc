// Copyright (C) 2009, 2010, 2011, 2013, 2014, 2015, 2016, 2017 Ed Bueler and Constantine Khroulev and Andy Aschwanden
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

#include <cassert>
#include <ctime>  // for time(), used to initialize random number gen
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>       // M_PI
#include <cmath>                // for erfc() in CalovGreveIntegrand()
#include <algorithm>
#include <iostream>

#include "pism/util/pism_utilities.hh"
#include "pism/util/ConfigInterface.hh"
#include "localITM.hh"
#include "pism/util/IceGrid.hh"

namespace pism {
namespace surface {

LocalMassBalanceITM::Changes::Changes() {
  snow_depth = 0.0;
  melt       = 0.0;
  runoff     = 0.0;
  smb        = 0.0;
}

LocalMassBalanceITM::Melt::Melt() {
  T_melt    = 0.0;
  I_melt    = 0.0;
  c_melt    = 0.0;
  ITM_melt  = 0.0;
}

LocalMassBalanceITM::LocalMassBalanceITM(Config::ConstPtr myconfig, units::System::Ptr system)
  : m_config(myconfig), m_unit_system(system),
    m_seconds_per_day(86400) {
  // empty
}

LocalMassBalanceITM::~LocalMassBalanceITM() {
  // empty
}

std::string LocalMassBalanceITM::method() const {
  return m_method;
}

ITMMassBalance::ITMMassBalance(Config::ConstPtr config, units::System::Ptr system)
  : LocalMassBalanceITM(config, system) {
  precip_as_snow     = m_config->get_flag("surface.pdd.interpret_precip_as_snow");
  Tmin               = m_config->get_number("surface.pdd.air_temp_all_precip_as_snow");
  Tmax               = m_config->get_number("surface.pdd.air_temp_all_precip_as_rain");
  refreeze_ice_melt  = m_config->get_flag("surface.pdd.refreeze_ice_melt");
  pdd_threshold_temp = m_config->get_number("surface.pdd.positive_threshold_temp");
  

  //FIXME use itm instead of pdd! 

  m_method = "insolation temperature melt";
}


/*! \brief Compute the number of points for temperature and
    precipitation time-series.
 */
unsigned int ITMMassBalance::get_timeseries_length(double dt) {
  const unsigned int    NperYear = static_cast<unsigned int>(m_config->get_number("surface.itm.max_evals_per_year"));
  const double dt_years = units::convert(m_unit_system, dt, "seconds", "years");

  return std::max(1U, static_cast<unsigned int>(ceil(NperYear * dt_years)));
}


double ITMMassBalance::CalovGreveIntegrand(double sigma, double TacC) {

  if (sigma == 0) {
    return std::max(TacC, 0.0);
  } else {
    const double Z = TacC / (sqrt(2.0) * sigma);
    return (sigma / sqrt(2.0 * M_PI)) * exp(-Z*Z) + (TacC / 2.0) * erfc(-Z);
  }
}


/* Determine albedo over ice cells by parametrization with pdd temperature*/

double ITMMassBalance::get_albedo_pdd(double T, double S,  int mask_value, bool print){
  double albedo = 0. ;
  double max = 0.78; //m_config->get_number("surface.itm.max_albedo");
  double min = 0.47; //m_config->get_number("surface.itm.min_albedo");
  double slope = -0.037; // m_config->get_number("surface.itm.albedo_slope");
  double Teff = 0;
  if (mask_value == 2 || mask_value == 3){ // mask value for grounded or floating ice
      Teff = CalovGreveIntegrand(S, T - pdd_threshold_temp);
      albedo = max + slope * Teff;
      if (albedo < min ){
        albedo = min; 
      }
  }

  return albedo;
}





/*  Determine albedo over ice cells (no ground or ocean cells respected here)
    if melt >= critical_melt, change albedo to melting snow
    TODO: put constants into config. 
*/
double ITMMassBalance::get_albedo(double melt, double snow_depth, double firn_depth, int mask_value, bool print){
  double albedo = 0. ;
  const double critical_melt = 100.;
  const double critical_melt_ice = 200.;
  const double albedo_dry_snow = 0.8;
  const double albedo_melting_snow = 0.6;
  const double albedo_firn = 0.5;
  const double albedo_ice = 0.4; 
  const double albedo_melting_ice = 0.3;
  const double albedo_land = 0.2;
  const double albedo_ocean = 0.1;

  if (snow_depth > 0.) {
    if (melt < critical_melt) {
      albedo = albedo_dry_snow;
    }
    else {
      albedo = albedo_melting_snow;
    }
  }
  else if (firn_depth > 0.) {
    albedo = albedo_firn;
  }
  else {


    if (mask_value == 4){ // mask value for ice free ocean
      albedo = albedo_ocean;
    }
    if (mask_value == 0){ // mask value for bedrock
      albedo =  albedo_land;
    }
    if (mask_value == 2 || mask_value == 3){ // mask value for grounded or floating ice
      if (melt >= critical_melt_ice){
        albedo = albedo_melting_ice;
      }
      else{
      albedo = albedo_ice;
      }
    }
  }

  return albedo;
}


double ITMMassBalance::get_albedo_melt(double melt, int mask_value, double dtseries, bool print){
  const double ice_density = m_config->get_number("constants.ice.density");
  double albedo =  0.82;
  const double albedo_land = 0.2;
  const double albedo_ocean = 0.1;
  // melt has a unit of meters ice equivalent
  // dtseries has a unit of seconds
  double intersection = 0.82; //FIXME: make those part of the config file
  double slope = -790; 

  if (mask_value == 4){ // mask value for ice free ocean
      albedo = albedo_ocean;
  }
  if (mask_value == 0){ // mask value for bedrock
      albedo =  albedo_land;
  }
  else {
      albedo = intersection + slope * melt * ice_density / (dtseries); //check if this is fine. 
      if (albedo < 0.47){
        albedo = 0.47;
      }
  }

  if(print){
    std::cout << "25. function get albedo melt\n" 
              << "26. melt input  = " << melt << '\n'
              << "27. time step = " << dtseries << '\n'
              << "28. with conversion = " << melt * ice_density   / (dtseries) << '\n'
              << "29. slope = " << slope << '\n'
              << "30. interserction = " << intersection << '\n'
              << "31. ----------------\n"
              << "32. resulting albedo = " << albedo  
              << "33. -------\n";
  }
  return albedo;
}



//! 
/* compute melt by equation (16) from Robinson2010
 * @param dt_series length of the step for the time-series
 * @param T air temperature at time [k]
 * @param insolation at time [k]
 * @param surface_elevation
 * @param albedo which was should be figured by get_albedo (?)
 * @param[out] melt pointer to a pre-allocated array with N-1 elements
 * output in mm water equivalent
 */
ITMMassBalance::Melt ITMMassBalance::calculate_ITM_melt(
  double dt_series,
  const double &S, 
  const double &T, 
  const double &surface_elevation, 
  const double &delta, 
  const double &lat, 
  double &albedo,
  bool print) {



  Melt ITM_melt;

  // double ITM_melt = 0.;
  // double T_melt = 0.;
  // double I_melt = 0.; 

  const double rho_w = m_config->get_number("constants.fresh_water.density");    // mass density of water //default 1000 kg m-3
  const double L_m = m_config->get_number("constants.fresh_water.latent_heat_of_fusion");      // latent heat of ice melting //default 3.34e5 Joule / kg
  double z = surface_elevation;               // surface elevation in meters
  double tau_a = 0.65 +  0.000032 * z;  // transmissivity of the atmosphere, linear fit, plug in values //FIXME parameters in config
  const double itm_c = m_config->get_number("surface.itm.itm_c");
  const double itm_lambda = m_config->get_number("surface.itm.itm_lambda");
  const double bm_temp    = m_config->get_number("surface.itm.background_melting_temp");
  // if background melting is true, use effective pdd temperatures and do not allow melting below background meltin temp. 
  const bool background_melting = m_config->get_flag("surface.itm.background_melting");



  assert(dt_series > 0.0);

  double input_h_0 = (- sin(lat) * sin(delta)) / (cos(lat) * cos(delta));
  double input_h_0_clipped = std::max(-1., std::min(input_h_0, 1.));
  double h_0 = acos(input_h_0_clipped);
  double q_insol;



  if (h_0 == 0. ){
    q_insol = 0;
  }
  else{
    q_insol = 1367. * (h_0 * sin(lat) * sin(delta) + cos(lat) * cos(delta) * sin(h_0)) / M_PI;
  }



  double Teff = 0.;
  if (background_melting){
    Teff = CalovGreveIntegrand(S, T - pdd_threshold_temp);
    if (Teff < 1.e-4){ Teff = 0;}
    ITM_melt.T_melt = dt_series / (rho_w * L_m) * itm_lambda * (Teff);
    if (T < bm_temp){
      ITM_melt.ITM_melt = 0.;
    }
    else{
      ITM_melt.ITM_melt = dt_series / (rho_w * L_m) * (tau_a*(1. - albedo) * q_insol + itm_c + itm_lambda * (Teff ));
    }
  }
  else{
    ITM_melt.ITM_melt = dt_series / (rho_w * L_m) * (tau_a*(1. - albedo) * q_insol + itm_c + itm_lambda * (T - 273.15 ));
    ITM_melt.T_melt = dt_series / (rho_w * L_m) * itm_lambda * (T - 273.15);
  }


  if (print){
    std::cout << "timestep dt = " << dt_series << '\n';
    // std::cout << "input temperature as seen by itm melt \n T = " << T << "\n";
    std::cout << "\n######\n\n";
    std::cout << "surface elevation z = " << z << "\n";
    std::cout << "delta = " << delta << "\n";
    std::cout << "input_h_0 = " << input_h_0 << "\n";
    std::cout << "input_h_0_clipped = " << input_h_0_clipped << "\n";
    std::cout << "h_0 = " << h_0 << "\n";
    std::cout << "lat = " << lat << "\n";
    std::cout << "insolation i = " << q_insol <<"\n";
    std::cout << "swd = " << q_insol * tau_a << "\n";
    std::cout << "albedo a = " <<  albedo << "\n";
    std::cout << "insolation melt = " << tau_a*(1. - albedo) * q_insol << "\n";
    std::cout << "temperature T = " << T << "\n";
    std::cout << "effective temperature Teff = " << Teff << "\n";
    std::cout << "temperature melt = " << itm_lambda * (Teff ) << "\n";
    std::cout << "itm_c = " << itm_c << '\n';
  }


  ITM_melt.I_melt = dt_series / (rho_w * L_m) * (tau_a*(1. - albedo) * q_insol);
  ITM_melt.c_melt = dt_series / (rho_w * L_m) * itm_c;

  if (print){
    std::cout << "in calculate ITM melt: temp melt rate = " << ITM_melt.T_melt  / dt_series << "\n";
    std::cout << "in calculate ITM melt: insol melt rate = " << ITM_melt.I_melt / dt_series << "\n";
    std::cout << "in calculate ITM melt: itm melt rate = " << ITM_melt.ITM_melt / dt_series << "\n";
    std::cout << "in calculate ITM melt: c melt rate = " << ITM_melt.c_melt / dt_series << "\n";
    std::cout << "\n prefactor " << dt_series / (rho_w * L_m) << "\n###### \n";
  }

  return ITM_melt;
}

/*
compute diurnal melt scheme  by equation (6) by Uta Krebs-Kanzow, The Cryosphere, 2018
*/

ITMMassBalance::Melt ITMMassBalance::calculate_ETIM_melt(double dt_series,
                                         const double &S,
                                         const double &T,
                                         const double &surface_elevation,
                                         const double &delta,
                                         const double &lat,
                                         double &albedo, bool print ) {

  Melt ETIM_melt;

  const double rho_w = 1e3;    // mass density of water //FIXME in config
  const double L_m = 3.34e5;      // latent heat of ice melting //FIXME in config
  const double z = surface_elevation;               // surface elevation 
  const double tau_a = 0.65 +  0.000032 * z;  // transmissivity of the atmosphere, linear fit, plug in values //FIXME parameters in config
  const double itm_c = m_config->get_number("surface.itm.itm_c");
  const double itm_lambda = m_config->get_number("surface.itm.itm_lambda");
  const double bm_temp    = m_config->get_number("surface.itm.background_melting_temp");
  // if background melting is true, use effective pdd temperatures and do not allow melting below background meltin temp. 
  const bool background_melting = m_config->get_flag("surface.itm.background_melting");



  const double phi = 17.5 * M_PI / 180.;

  double input_h_phi = ( sin(phi) - sin(lat) * sin(delta)) / (cos(lat) * cos(delta));
  double input_h_phi_clipped = std::max(-1., std::min(input_h_phi, 1.));
  double h_phi = acos(input_h_phi_clipped);
  double quotient_delta_t =  h_phi /M_PI ;
  double q_insol;


  ETIM_melt.transmissivity = tau_a;



  if (h_phi == 0. ){
    q_insol = 0;
  }
  else{
    q_insol = 1367. * (h_phi * sin(lat) * sin(delta) + cos(lat) * cos(delta) * sin(h_phi))  / h_phi;
  }
  
  ETIM_melt.TOA_insol = q_insol; 

  assert(dt_series > 0.0);
  double Teff = 0;
  if (background_melting){
    Teff = CalovGreveIntegrand(S, T - pdd_threshold_temp);
    if (Teff < 1.e-4){ Teff = 0;}

    ETIM_melt.T_melt = quotient_delta_t * dt_series / (rho_w * L_m) * itm_lambda * (Teff);
    
    if (T < bm_temp){
      ETIM_melt.ITM_melt = 0.;
    }
    else{
      ETIM_melt.ITM_melt = quotient_delta_t * dt_series / (rho_w * L_m) * (tau_a*(1. - albedo) * q_insol + itm_c + itm_lambda * (Teff ));
    }
  }
  else{
    ETIM_melt.ITM_melt = quotient_delta_t * dt_series / (rho_w * L_m) * (tau_a*(1. - albedo) * q_insol + itm_c + itm_lambda * (T - 273.15 ));
    ETIM_melt.T_melt = dt_series  * quotient_delta_t / (rho_w * L_m) * itm_lambda * (T - 273.15);

  }
  
  ETIM_melt.I_melt = dt_series / (rho_w * L_m) * (tau_a * (1. - albedo) * q_insol) * quotient_delta_t;
  ETIM_melt.c_melt = dt_series / (rho_w * L_m) * itm_c * quotient_delta_t;

  // TODO: check units!!!!



  if (print){
    std::cout << " 3. ######\n";
    std::cout << " 4. surface elevation z = " << z << "\n";
    std::cout << " 5. delta = " << delta << "\n";
    std::cout << " 6. lat = " << lat << "\n";
    std::cout << " 7. input_h_phi = " << input_h_phi << "\n";
    std::cout << " 8. input_h_phi_clipped = " << input_h_phi_clipped << "\n";
    std::cout << " 9. h_phi = " << h_phi << "\n";
    std::cout << "10. qdeltat = " << quotient_delta_t << "\n";
    std::cout << "11. q_insol = " << q_insol <<"\n";
    std::cout << "12. swd = " << q_insol * tau_a << "\n";
    std::cout << "13. albedo a = " <<  albedo << "\n";
    std::cout << "14. insolation melt = " << tau_a*(1. - albedo) * q_insol << "\n";
    std::cout << "15. temperature T = " << T << "\n";
    std::cout << "16. effective temperature Teff = " << Teff << "\n";
    std::cout << "17. temperature melt = " << itm_lambda * (Teff ) << "\n";
    std::cout << "18. itm_c = " << itm_c << '\n';
    std::cout << "19. prefactor = " << dt_series / (rho_w * L_m) * quotient_delta_t << '\n';
  }

  //T - 273.15));
  //FIXME Teff )); // hier muss irgendwo die Formel (16) from Robinson2010;


  if (print){
    std::cout << "20. in calculate ETIM melt: temp melt rate (m) = " << ETIM_melt.T_melt << "\n";
    std::cout << "21. in calculate ETIM melt: insol melt rate (m) = " << ETIM_melt.I_melt  << "\n";
    std::cout << "22. in calculate ETIM melt: c melt rate (m) = " << ETIM_melt.c_melt  << "\n";
    std::cout << "23. ---> in calculate ETIM melt: total melt rate (m) = " << "bla bla " << "\n";
    std::cout << "24. ###### \n";


  }

  return ETIM_melt;
}


//! \brief Extract snow accumulation from mixed (snow and rain)
//! precipitation using the temperature time-series.
/** Uses the temperature time-series to determine whether the
 * precipitation is snow or rain. Rain is removed entirely from the
 * surface mass balance, and will not be included in the computed
 * runoff, which is meltwater runoff. There is an allowed linear
 * transition for Tmin below which all precipitation is interpreted as
 * snow, and Tmax above which all precipitation is rain (see, e.g.
 * [\ref Hock2005b]).
 *
 * Sets P[i] to the *solid* (snow) accumulation *rate*.
 *
 * @param[in,out] P precipitation rate (array of length N)
 * @param[in] T air temperature (array of length N)
 * @param[in] N array length
 */
void ITMMassBalance::get_snow_accumulationITM(const std::vector<double> &T,
                                           std::vector<double> &P) {


  assert(T.size() == P.size());
  const size_t N = T.size();

  // Following \ref Hock2005b we employ a linear transition from Tmin to Tmax
  for (unsigned int i = 0; i < N; i++) {
    // do not allow negative precipitation
    if (P[i] < 0.0) {
      P[i] = 0.0;
      continue;
    }

    if (precip_as_snow || T[i] <= Tmin) { // T <= Tmin, all precip is snow
      // no change
    } else if (T[i] < Tmax) { // linear transition from Tmin to Tmax
      P[i] *= (Tmax - T[i]) / (Tmax - Tmin);
    } else { // T >= Tmax, all precip is rain -- ignore it
      P[i] = 0.0;
    }
  }

}


double ITMMassBalance::get_refreeze_fraction(const double &T) {
  double refreeze;
  double Tmin_refreeze  = m_config->get_number("surface.itm.air_temp_all_refreeze");
  double Tmax_refreeze  = m_config->get_number("surface.itm.air_temp_no_refreeze");
  if (T <= Tmin_refreeze){refreeze = 1. ;}
  else if ((Tmin_refreeze<  T) and (T <= Tmax_refreeze)){
    refreeze = 1./(Tmin_refreeze - Tmax_refreeze) * T + Tmax_refreeze / (Tmax_refreeze - Tmin_refreeze) ; 
  }
  else {refreeze = 0.;}
  return refreeze;
}


//! \brief Compute the surface mass balance at a location from the amount of 
//! melted snow and the accumulation amount in a time interval.
/*!
 * This is a ITM scheme. The input parameter `melt_conversion_factor` is a conversion between
 * snow melting and ice melting per time period.
 *
 * `accumulation` has units "meter / second".
 *
 * - a fraction of the melted snow and ice refreezes, conceptualized
 *   as superimposed ice, and this is controlled by parameter \c
 *   refreeze_fraction
 *
 * - the excess of 'ITM_melt' is used to melt both the ice that came
 *   from refreeze and then any ice which is already present.
 *
 * The scheme here came from EISMINT-Greenland [\ref RitzEISMINT], but
 * is influenced by R. Hock (personal communication).
 */
ITMMassBalance::Changes ITMMassBalance::step(const double &melt_conversion_factor,
                                             const double &refreeze_fraction,
                                             double thickness,
                                             double ITM_melt,
                                             double old_firn_depth,
                                             double old_snow_depth,
                                             double accumulation, bool print) {


  Changes result;

  if (print){
    // std::cout<<"** this comes from the step function \n read in variables are \n ITM melt = " << ITM_melt << "\n old_firn_depth = " << old_firn_depth << "\n old snow depth = " << old_snow_depth << "\n accumulation = " << accumulation << "\n **** \n";
  }
  double
    firn_depth      = old_firn_depth,
    snow_depth      = old_snow_depth,
    max_snow_melted = ITM_melt,
    firn_melted     = 0.0,
    snow_melted     = 0.0,
    excess_melt     = 0.0;

    // FIXME check if itm melt is read in correctly!!!!!

  // if (old_snow_depth > 0){
  //   std::cout << "initial snow depth " << old_snow_depth << '\n';
  //   std::cout << "snow depth before accumulation " << snow_depth << '\n';
  //   std::cout << "accumulation " << accumulation << '\n';
  //   std::cout << "melt " << ITM_melt << '\n';
  // }
  assert(thickness >= 0);

  // snow depth cannot exceed total thickness
  snow_depth = std::min(snow_depth, thickness);

  assert(snow_depth >= 0);


  // firn depth cannot exceed thickness - snow_depth
  firn_depth = std::min(firn_depth, thickness - snow_depth);

  // if (old_snow_depth > 0){
  //   std::cout << "thickness = " << thickness << " snow depth = " << snow_depth << '\n';
  //   std::cout << "firn depth " << firn_depth << '\n';
  // }

  assert(firn_depth >= 0);

  double ice_thickness = thickness - snow_depth - firn_depth;

  assert(ice_thickness >= 0);

  snow_depth += accumulation;
  // if (old_snow_depth > 0){
  //   std::cout << "snow depth after accumulation " << snow_depth << '\n';
  // }
  
  if (ITM_melt <= 0.0) {            // The "no melt" case.
    snow_melted = 0.0;
    firn_melted = 0.0,
    excess_melt = 0.0;
    if (print) {std::cout << "no melt case \n";}
  } else if (max_snow_melted <= snow_depth) {
    // Some of the snow melted and some is left; in any case, all of
    // the energy available for melt was used up in melting snow.
    snow_melted = max_snow_melted;
    firn_melted = 0.0;
    excess_melt = 0.0;
  } else if (max_snow_melted <= firn_depth + snow_depth) {
    // All of the snow is melted but some firn is left; in any case, all of
    // the energy available for melt was used up in melting snow.
    snow_melted = snow_depth;
    firn_melted = max_snow_melted - snow_melted;
    excess_melt = 0.0;
  } else {
    // All (firn and snow_depth meters) of snow melted. Excess_pddsum is the
    // positive degree days available to melt ice.
    firn_melted = firn_depth;
    snow_melted = snow_depth;
    excess_melt = ITM_melt - (firn_melted + snow_melted) ; // 
  }

  double
    ice_melted              = excess_melt * melt_conversion_factor,
    melt                    = snow_melted + firn_melted + ice_melted,
    ice_created_by_refreeze = 0.0;

  if (print){
    std::cout<< " snow melted  = " << snow_melted << "\n firn melted = " << firn_melted << "\n ice_melted = " << ice_melted <<  
    "\n ITM melt " << ITM_melt << "\n melt " << melt << "\n **** \n";
  }

  // if (old_snow_depth > 0){
  // std::cout << "snow melted " << snow_melted << '\n';
  // }

  if (refreeze_ice_melt) {
    ice_created_by_refreeze = melt * refreeze_fraction;
  } else {
    // Should this only be snow melted?
    ice_created_by_refreeze = (firn_melted + snow_melted) * refreeze_fraction;
  }


  snow_depth = std::max(snow_depth - snow_melted, 0.0);
  // if (old_snow_depth > 0){
  // std::cout << "snow depth after checking " << snow_depth << '\n';
  // }
  firn_depth = std::max(firn_depth - firn_melted, 0.0);
  // FIXME: need to add snow that hasn't melted, is this correct?
  // firn_depth += (snow_depth - snow_melted);
  // Turn firn into ice at X times accumulation
  // firn_depth -= accumulation *  m_config->get_number("surface.pdd.firn_compaction_to_accumulation_ratio");

  

  const double runoff = melt - ice_created_by_refreeze;
  const double smb = accumulation - runoff;

  result.firn_depth = firn_depth - old_firn_depth;
  result.snow_depth = snow_depth - old_snow_depth;
  result.melt       = melt;
  result.runoff     = runoff;
  // result.smb        = accumulation - runoff;
  result.smb        = thickness + smb >= 0 ? smb : -thickness;

  // if (print){
    // std::cout << "ITM melt " << ITM_melt << '\n';
    // std::cout << "accumulation " << accumulation << '\n';
    // if (old_snow_depth > 0){
    // // std::cout << "initial snow depth " << old_snow_depth << '\n';
    // std::cout << "new snow depth " << snow_depth << '\n';
    // std::cout << "snow melted " << snow_melted << '\n';
    // std::cout << "changes snow depth " << result.snow_depth << "\n _______________________ \n";

    // }
    // std::cout << "now at the end of step function \n the melt calculated is ==== " << melt << " \n but melt returned is ==== " << result.melt << "\n ** \n ";
  // }
  assert(thickness + result.smb >= 0);

  return result;


}




} // end of namespace surface
} // end of namespace pism
