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

#include <cassert>

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

  double bmeltrate, thetao, temp_base;

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
                                           double sal_ocean, double potential_temperature, double zice,
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

  double sal_base;
  double ep1,ep2,ep3,ep4,ep4b;
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
  ex2 = ep1*(ep4b-potential_temperature)+ep2*(tob+ai*sal_ocean-ep4)-ep3;
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

/** Compute temperature and melt rate at the base of the shelf.
 *
 * Use equations for the heat and salt flux balance at the base of the
 * shelf to compute the temperature at the base of the shelf and the
 * sub-shelf melt rate.
 *
 * @note This model is not applicable in the case of basal freeze-on.
 *
 * Following [@ref HellmerOlbers1989], let @f$ Q_T @f$ be the total heat
 * flux crossing the interface between the shelf base and the ocean,
 * @f$ Q_T^B @f$ be the amount of heat lost by the ocean due to
 * melting of glacial ice, and @f$ Q_T^I @f$ be the conductive flux
 * into the ice column.
 *
 * @f$ Q_{T} @f$ is parameterized by (see [@ref HellmerOlbers1989], equation
 * 10):
 *
 * @f[ Q_{T} = \rho_W\, c_{pW}\, \gamma_{T}\, (T^B - T^W),@f]
 *
 * where @f$ \rho_{W} @f$ is the sea water density, @f$ c_{pW} @f$ is
 * the heat capacity of sea water, and @f$ \gamma_{T} @f$ is a
 * turbulent heat exchange coefficient.
 *
 * We assume that the difference between the basal temperature and
 * adjacent ocean temperature @f$ T^B - T^W @f$ is well approximated
 * by @f$ \Theta_B - \Theta_W, @f$ where @f$ \Theta_{\cdot} @f$ is the
 * corresponding potential temperature.
 *
 * @f$ Q_T^B @f$ is parameterized by (see [@ref HellmerOlbers1989], equation 11):
 *
 * @f[ Q_T^B = \rho_I\, L\, \dot h, @f]
 *
 * where @f$ \rho_I @f$ is the ice density, @f$ L @f$ is the latent
 * heat of fusion, and @f$ \dot h @f$ is the ice thickening rate
 * (equal to minus the melt rate).
 *
 * The conductive flux into the ice column is parameterized by ([@ref
 * Hellmeretal1998], equation 7):
 *
 * @f[ Q_T^I = \rho_I\, c_{pI}\, \kappa\, T_{\text{grad}}, @f]
 *
 * where @f$ \rho_I @f$ is the ice density, @f$ c_{pI} @f$ is the heat
 * capacity of ice, @f$ \kappa @f$ is the ice thermal diffusivity, and
 * @f$ T_{\text{grad}} @f$ is the vertical temperature gradient at the
 * base of a column of ice.
 *
 * We use the parameterization of the temperature gradient from [@ref
 * Hellmeretal1998], equation 13:
 *
 * @f[ T_{\text{grad}} = -\Delta T\, \frac{\dot h}{\kappa}, @f]
 *
 * where @f$ \Delta T @f$ is the difference between the ice
 * temperature at the top of the ice column and its bottom:
 * @f$ \Delta T = T^S - T^B. @f$ With this parameterization, we have
 *
 * @f[ Q_T^I = \rho_I\, c_{pI}\, {\dot h}\, (T^S - T^B). @f]
 *
 * Now, the heat flux balance implies
 *
 * @f[ Q_T = Q_T^B + Q_T^I. @f]
 *
 * For the salt flux balance, we have
 *
 * @f[ Q_S = Q_S^B + Q_S^I, @f]
 *
 * where @f$ Q_S @f$ is the total salt flux across the interface, @f$
 * Q_S^B @f$ is the salt input as a result of freezing, @f$ Q_S^I = 0
 * @f$ is the salt flux due to molecular diffusion of salt through
 * ice.
 *
 * @f$ Q_S @f$ is parameterized by ([@ref Hellmeretal1998], equation 13)
 *
 * @f[ Q_S = \rho_W\, \gamma_S\, (S^B - S^W), @f]
 *
 * where @f$ \gamma_S @f$ is a turbulent salt exchange coefficient,
 * @f$ S_B @f$ is salinity at the shelf base, and @f$ S^W @f$ is the
 * salinity of adjacent ocean.
 *
 * The basal salt flux @f$ Q_S^B @f$ is parameterized by ([@ref
 * Hellmeretal1998], equation 14)
 *
 * @f[ Q_S^B = \rho_I\, S^B\, {\dot h}. @f]
 * 
 * To avoid converting shelf base temperature to shelf base potential
 * temperature and back, we use following parameterizations:
 *
 * @f[ T^{B}(S,h) = a_0\cdot S + a_1 + a_2\cdot h, @f]
 *
 * @f[ \Theta^{B}(S,h) = b_0\cdot S + b_1 + b_2\cdot h, @f]
 *
 * where @f$ S @f$ is salinity and @f$ h @f$ is ice shelf thickness.
 *
 * The parameterization for the basal temperature @f$ T^B(S,h) @f$ is
 * based on [@ref Hellmeretal1998] and [@ref FoldvikKvinge1974].
 *
 * Given this parameterization and a function @f$ \Theta_T^B(T)
 * @f$ one can define @f$ \Theta^B_{*}(S,h) = \Theta_T^B(T^B(S,h) @f$.
 *
 * The parameterization @f$ \Theta^B(S,h) @f$ used here was produced
 * by linearizing @f$ \Theta^B_{*}(S,h) @f$ near the melting point.
 * (The definition of @f$ \Theta_T^B(T), @f$, converting in situ
 * temperature into potential temperature, was adopted from FESOM
 * [@ref Wangetal2013]).
 *
 * Treating ice thickness, sea water salinity, and sea water potential
 * temperature as "known" we can write down a system of equations
 *
 * @f{align*}{
 * Q_T &= Q_T^B + Q_T^I,\\
 * Q_S &= Q_S^B + Q_S^I,\\
 * T^{B}(S,h) &= a_0\cdot S + a_1 + a_2\cdot h,\\
 * \Theta^{B}(S,h) &= b_0\cdot S + b_1 + b_2\cdot h\\
 * @f}
 *
 * and simplify it to produce a quadratic equation for the salinity at the shelf base, @f$ S^B @f$:
 *
 * @f[ A\cdot (S^B)^2 + B\cdot S^B + C = 0 @f]
 *
 * with
 * @f{align*}{
 * A &= a_{0}\,\gamma_S\,c_{pI}-b_{0}\,\gamma_T\,c_{pW}\\
 * B &= \gamma_S\,\left(L-c_{pI}\,\left(T^S+a_{0}\,S^W-a_{2}\,h-a_{1}\right)\right)+
 *      \gamma_T\,c_{pW}\,\left(\Theta^W-b_{2}\,h-b_{1}\right)\\
 * C &= -\gamma_S\,S^W\,\left(L-c_{pI}\,\left(T^S-a_{2}\,h-a_{1}\right)\right)
 * @f}
 *
 * Once @f$ S_B @f$ is found by solving this quadratic equation, we can
 * compute the basal temperature using the parameterization for @f$
 * T^{B}(S,h) @f$.
 *
 * To find the basal melt rate, we solve the salt flux balance
 * equation for @f$ {\dot h}, @f$ obtaining
 *
 * @f[ w_b={{\gamma_S\,\rho_W\,\left(S^W-S^B\right)}\over{\rho_I\,S^B}}. @f]
 * 
 * 
 * @param[in] c physical constants, stored here to avoid looking them up in a double for loop
 * @param[in] sea_water_salinity salinity of the ocean immediately adjacent to the shelf, [g/kg]
 * @param[in] sea_water_potential_temperature potential temperature of the sea water, [degrees Celsius]
 * @param[in] thickness thickness of the ice shelf, [meters]
 * @param[out] shelf_base_temperature_out computed basal temperature, [degrees Celsius]
 * @param[out] shelf_base_melt_rate_out computed basal melt rate, [m/second]
 *
 * @return 0 on success
 */
PetscErrorCode POGivenTH::pointwise_calculation(const POGivenTHConstants &c,
                                                double sea_water_salinity,
                                                double sea_water_potential_temperature,
                                                double thickness,
                                                double *shelf_base_temperature_out,
                                                double *shelf_base_melt_rate_out) {

  assert(sea_water_salinity >= 0.0);
  assert(sea_water_salinity <= 400.0); // that would be saltier than the Dead Sea
  assert(thickness >= 0.0);

  // Coefficients for linearized freezing point equation for in situ
  // temperature:
  //
  // Tb(salinity, thickness) = a[0] * salinity + a[1] + a[2] * thickness
  const double a[3] = {-0.0575, 0.0901, -7.61e-4};

  // Coefficients for linearized freezing point equation for potential
  // temperature
  //
  // Theta_b(salinity, thickness) = b[0] * salinity + b[1] + b[2] * thickness
  const double b[3] = {-0.0575, 0.0921, -7.85e-4};

  // Turbulent heat and salt transfer coefficients:
  const double gamma_t = 1.00e-4;   // [m/s] RG3417 Default value from Hellmer and Olbers 89
  const double gamma_s = 5.05e-7;   // [m/s] RG3417 Default value from Hellmer and Olbers 89

  const double
    c_pI    = c.ice_specific_heat_capacity,
    c_pW    = c.sea_water_specific_heat_capacity,
    L       = c.water_latent_heat_fusion,
    T_S     = c.shelf_top_surface_temperature,
    S_W     = sea_water_salinity,
    Theta_W = sea_water_potential_temperature;

  // We solve a quadratic equation for Sb, the salinity at the shelf
  // base.
  //
  // A*Sb^2 + B*Sb + C = 0
  const double A = a[0] * gamma_s * c_pI - b[0] * gamma_t * c_pW;
  const double B = (gamma_s * (L - c_pI * (T_S + a[0] * S_W - a[2] * thickness - a[1])) +
                    gamma_t * c_pW * (Theta_W - b[2] * thickness - b[1]));
  const double C = -gamma_s * S_W * (L - c_pI * (T_S - a[2] * thickness - a[1]));
  // Find two roots of the equation:
  const double S1 = (-B + sqrt(B*B - 4.0 * A * C)) / (2.0 * A);
  const double S2 = (-B - sqrt(B*B - 4.0 * A * C)) / (2.0 * A);

  // Both roots cannot be negative at the same time
  assert((S1 < 0.0 && S2 < 0.0) == false);

  // pick the positive root
  double basal_salinity = 0.0;
  if (S1 > 0.0) {
    basal_salinity = S1;
  } else {
    basal_salinity = S2;
  }

  assert(basal_salinity >= 0.0);
  assert(basal_salinity <= 400.0); // this would be saltier than the Dead Sea
  assert(basal_salinity <= sea_water_salinity); // ice melt should lower salinity

  if (shelf_base_temperature_out != NULL) {
    *shelf_base_temperature_out = a[0] * basal_salinity + a[1] + a[2] * thickness;
  }

  if (shelf_base_melt_rate_out != NULL) {
    *shelf_base_melt_rate_out = gamma_s * c.sea_water_density * (sea_water_salinity - basal_salinity) / (c.ice_density * basal_salinity);

    // we use an approximation of the temperature gradient at the base
    // of the shelf that is invalid for negative melt rates.
    assert(*shelf_base_melt_rate_out >= 0.0);
  }

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
