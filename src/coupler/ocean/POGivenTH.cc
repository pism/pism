// Copyright (C) 2011, 2012, 2013, 2014 PISM Authors
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

#include "POGivenTH.hh"
#include "IceGrid.hh"
#include "PISMVars.hh"

/* TO DO:
 *
 * Add detailed references ("this implements equations (3)-(7) in Smith and
 * Jones, 2001" and similar).
 *
 */

POGivenTH::POGivenTH(IceGrid &g, const PISMConfig &conf)
  : PGivenClimate<POModifier,PISMOceanModel>(g, conf, NULL)
{
  PetscErrorCode ierr = allocate_POGivenTH(); CHKERRCONTINUE(ierr);
  if (ierr != 0)
    PISMEnd();
}

POGivenTH::~POGivenTH() {
  // empty
}

PetscErrorCode POGivenTH::allocate_POGivenTH() {
  PetscErrorCode ierr;
  option_prefix   = "-ocean_th";

  // will be de-allocated by the parent's destructor
  theta_ocean    = new IceModelVec2T;
  salinity_ocean = new IceModelVec2T;

  m_fields["theta_ocean"]     = theta_ocean;
  m_fields["salinity_ocean"]  = salinity_ocean;

  ierr = process_options(); CHKERRQ(ierr);

  std::map<std::string, std::string> standard_names;
  ierr = set_vec_parameters(standard_names); CHKERRQ(ierr);

  ierr = theta_ocean->create(grid, "theta_ocean", false); CHKERRQ(ierr);
  ierr = theta_ocean->set_attrs("climate_forcing",
                                "absolute potential temperature of the adjacent ocean",
                                "Kelvin", ""); CHKERRQ(ierr);

  ierr = salinity_ocean->create(grid, "salinity_ocean", false); CHKERRQ(ierr);
  ierr = salinity_ocean->set_attrs("climate_forcing",
                                   "salinity of the adjacent ocean",
                                   "g/kg", ""); CHKERRQ(ierr);

  ierr = shelfbtemp.create(grid, "shelfbtemp", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = shelfbtemp.set_attrs("climate_forcing",
                              "absolute temperature at ice shelf base",
                              "Kelvin", ""); CHKERRQ(ierr);

  ierr = shelfbmassflux.create(grid, "shelfbmassflux", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = shelfbmassflux.set_attrs("climate_forcing",
                                  "ice mass flux from ice shelf base (positive flux is loss from ice shelf)",
                                  "kg m-2 s-1", ""); CHKERRQ(ierr);
  ierr = shelfbmassflux.set_glaciological_units("kg m-2 year-1"); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode POGivenTH::init(PISMVars &vars) {
  PetscErrorCode ierr;

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  ierr = verbPrintf(2, grid.com,
                    "* Initializing the 3eqn melting parameterization ocean model\n"
                    "  reading ocean temperature and salinity from a file...\n"); CHKERRQ(ierr);

  ice_thickness = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
  if (!ice_thickness) {SETERRQ(grid.com, 1, "ERROR: ice thickness is not available");}

  ierr = theta_ocean->init(filename, bc_period, bc_reference_time); CHKERRQ(ierr);
  ierr = salinity_ocean->init(filename, bc_period, bc_reference_time); CHKERRQ(ierr);

  // read time-independent data right away:
  if (theta_ocean->get_n_records() == 1 && salinity_ocean->get_n_records() == 1) {
    ierr = update(grid.time->current(), 0); CHKERRQ(ierr); // dt is irrelevant
  }

  return 0;
}

PetscErrorCode POGivenTH::update(double my_t, double my_dt) {

  PetscErrorCode ierr = update_internal(my_t, my_dt); CHKERRQ(ierr);

  ierr = theta_ocean->average(m_t, m_dt); CHKERRQ(ierr);
  ierr = salinity_ocean->average(m_t, m_dt); CHKERRQ(ierr);

  ierr = calc_shelfbtemp_shelfbmassflux(); CHKERRQ(ierr);

  // convert from [m s-1] to [kg m-2 s-1]:
  ierr = shelfbmassflux.scale(config.get("ice_density")); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode POGivenTH::shelf_base_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr = shelfbtemp.copy_to(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode POGivenTH::shelf_base_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr = shelfbmassflux.copy_to(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode POGivenTH::calc_shelfbtemp_shelfbmassflux() {

  PetscErrorCode ierr;

  const double rhoi = config.get("ice_density");
  const double rhow = config.get("sea_water_density");
  const double reference_pressure = 1.01325; // pressure of atmosphere in bar

  double pressure_at_shelf_base, bmeltrate, thetao, temp_base;

  ierr = ice_thickness->begin_access();   CHKERRQ(ierr);
  ierr = theta_ocean->begin_access(); CHKERRQ(ierr);
  ierr = salinity_ocean->begin_access(); CHKERRQ(ierr);
  ierr = shelfbmassflux.begin_access(); CHKERRQ(ierr);
  ierr = shelfbtemp.begin_access(); CHKERRQ(ierr);

  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {

      thetao = (*theta_ocean)(i,j) - 273.15; // convert from Kelvin to Celsius

      btemp_bmelt_3eqn(rhow, rhoi,(*salinity_ocean)(i,j), thetao, (*ice_thickness)(i,j), temp_base, bmeltrate);

      // the ice/ocean boundary layer temperature is seen by PISM as shelfbtemp.
      shelfbtemp(i,j) = temp_base + 273.15; // convert from Celsius to Kelvin

      // FIXME: do we need this "-1"? (note the definition of the sub-shelf mass flux above).
      shelfbmassflux(i,j) = bmeltrate;
    }
  }

  ierr = ice_thickness->end_access(); CHKERRQ(ierr);
  ierr = theta_ocean->end_access(); CHKERRQ(ierr);
  ierr = salinity_ocean->end_access(); CHKERRQ(ierr);
  ierr = shelfbmassflux.end_access(); CHKERRQ(ierr);
  ierr = shelfbtemp.end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode POGivenTH::btemp_bmelt_3eqn(double rhow, double rhoi,
                                           double sal_ocean, double theta_ocean, double zice,
                                           double &temp_base, double &meltrate) {

  // This function solves the three equation model of ice-shelf ocean interaction (Hellmer and Olbers, 1989).
  // Equations are
  // (1) freezing point dependence on salinity and pressure
  // (2) heat conservation in brine layer
  // (3) salinity conservation in brine layer

  // Code derived from BRIOS subroutine iceshelf (which goes back to H.Hellmer's 2D ice shelf model code)
  // and adjusted for use in FESOM by Ralph Timmermann, 16.02.2011
  // adapted for PISM by matthias.mengel@pik-potsdam.de

  // The model is described in section 3 of
  // Hellmer, Hartmut, S. S. Jacobs, and A. Jenkins.
  // "Oceanic erosion of a floating Antarctic glacier in the Amundsen Sea."
  // Ocean, Ice, and Atmosphere: Interactions at the Antarctic continental margin (S Jacobs, R Weiss, eds)
  // Antarctic Research Series, AGU, Washington DC, USA 75 (1998): 319-339.

  // FIXME: uses fixed temperature at ice surface tob=-20 degC to calculate
  //        heat flux from ice/ocean brine boundary layer into ice.
  //        We could do this better by usings PISMs ice temperature or heat flux
  //        calculus at the bottom layers of the ice shelf.

  double rhor, sal_base;
  double ep1,ep2,ep3,ep4,ep4b,ep5;
  double ex1,ex2,ex3,ex4,ex5;
  double sr1,sr2,sf1,sf2,tf1,tf2,tf,sf;

  // coefficients for linearized freezing point equation
  // for in situ temperature
  double ai   = -0.0575;                // [°C/psu] Foldvik&Kvinge (1974)
  double bi   =  0.0901;                // [°C]
  double ci   =  7.61e-4;               // [°C/m]

  // coefficients for linearized freezing point equation
  // for potential temperature

  double ap  =  -0.0575;   // [°C/psu]
  double bp  =   0.0921;   // [°C]
  double cp  =   7.85e-4; // [°C/m]

  double tob=  -20.;                   //temperature at the ice surface
  double cpw =  4180.0;                //Barnier et al. (1995)
  double lhf =  3.33e+5;               //latent heat of fusion
  double atk =  273.15;                //0 deg C in Kelvin
  //FIXME: can use PISMs surface temp for tob?
  double cpi =  152.5+7.122*(atk+tob); //Paterson:"The Physics of Glaciers"

  // Prescribe the turbulent heat and salt transfer coeff. GAT and GAS
  double gat  = 1.00e-4;   //[m/s] RG3417 Default value from Hellmer and Olbers 89
  double gas  = 5.05e-7;   //[m/s] RG3417 Default value from Hellmer and Olbers 89

  // calculate salinity and in situ temperature of ice/ocean boundary layer,
  // by solving a quadratic equation in salinity (sf).
  ep1  = cpw*gat;
  ep2  = cpi*gas;
  ep3  = lhf*gas;
  ep4  = bi-ci*zice;
  ep4b = bp-cp*zice;
  // negative heat flux term in the ice (due to -kappa/D)
  //ex1 = ai*(ep1-ep2);
  ex1 = ai*ep1-ap*ep2;
  ex2 = ep1*(ep4b-theta_ocean)+ep2*(tob+ai*sal_ocean-ep4)-ep3;
  ex3 = sal_ocean*(ep2*(ep4-tob)+ep3);
  ex4 = ex2/ex1;
  ex5 = ex3/ex1;

  sr1 = 0.25*ex4*ex4-ex5;
  sr2 = -0.5*ex4;
  sf1 = sr2+sqrt(sr1);
  tf1 = ai*sf1+ep4;
  sf2 = sr2-sqrt(sr1);
  tf2 = ai*sf2+ep4;

  // sf is solution of quadratic equation in salinity.
  // salinities < 0 psu are not defined, therefore pick the positive of the two solutions.
  if(sf1 > 0.) {
    tf = tf1;
    sf = sf1;
  }else{
    tf = tf2;
    sf = sf2;
  }

  temp_base = tf;
  sal_base  = sf;

  // Calculate
  // density in the boundary layer: rhow
  // and interface pressure pg [dbar]
  // to determine the melting/freezing rate.

  // FIXME: need to calculate water density instead of const value.
  //  call fcn_density(thetao,sal,zice,rho)
  // matthias.mengel: meltrate scales linear with rhow, so
  //                  the error should be not more than 1e-2.

  // Calculate the melting/freezing rate [m/s]
  meltrate = -1*gas*rhow/rhoi*(1.0-sal_ocean/sal_base);

  return 0;
}


PetscErrorCode POGivenTH::sea_level_elevation(double &result) {
  result = sea_level;
  return 0;
}

PetscErrorCode POGivenTH::melange_back_pressure_fraction(IceModelVec2S &result) {
  PetscErrorCode ierr = result.set(0.0); CHKERRQ(ierr);
  return 0;
}
