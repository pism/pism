// Copyright (C) 2011, 2012, 2014 PISM Authors
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

#ifndef _POGIVENTH_H_
#define _POGIVENTH_H_

#include "PGivenClimate.hh"
#include "POModifier.hh"

namespace pism {

class POGivenTH : public PGivenClimate<POModifier,PISMOceanModel>
{
public:
  POGivenTH(IceGrid &g, const PISMConfig &conf);
  virtual ~POGivenTH();

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(double my_t, double my_dt);

  virtual PetscErrorCode sea_level_elevation(double &result);

  virtual PetscErrorCode shelf_base_temperature(IceModelVec2S &result);
  virtual PetscErrorCode shelf_base_mass_flux(IceModelVec2S &result);

  virtual PetscErrorCode melange_back_pressure_fraction(IceModelVec2S &result);

  class POGivenTHConstants {
  public:
    POGivenTHConstants(const PISMConfig &config);
    //! Coefficients for linearized freezing point equation for in situ
    //! temperature:
    //!
    //! Tb(salinity, thickness) = a[0] * salinity + a[1] + a[2] * thickness
    double a[3];
    //! Coefficients for linearized freezing point equation for potential
    //! temperature
    //!
    //! Theta_b(salinity, thickness) = b[0] * salinity + b[1] + b[2] * thickness
    double b[3];

    //! Turbulent heat transfer coefficient:
    double gamma_T;
    //! Turbulent salt transfer coefficient:
    double gamma_S;

    double shelf_top_surface_temperature;
    double water_latent_heat_fusion;
    double sea_water_density;
    double sea_water_specific_heat_capacity;
    double ice_density;
    double ice_specific_heat_capacity;
    double ice_thermal_diffusivity;
    bool limit_salinity_range;
  };
private:
  IceModelVec2S shelfbtemp, shelfbmassflux;
  IceModelVec2S *ice_thickness;
  IceModelVec2T *theta_ocean, *salinity_ocean;

  PetscErrorCode pointwise_update(const POGivenTHConstants &constants,
                                  double sea_water_salinity,
                                  double sea_water_potential_temperature,
                                  double ice_thickness,
                                  double *shelf_base_temperature_out,
                                  double *shelf_base_melt_rate_out);

  PetscErrorCode subshelf_salinity(const POGivenTHConstants &constants,
                                   double sea_water_salinity,
                                   double sea_water_potential_temperature,
                                   double ice_thickness,
                                   double *shelf_base_salinity);

  PetscErrorCode subshelf_salinity_melt(const POGivenTHConstants &constants,
                                        double sea_water_salinity,
                                        double sea_water_potential_temperature,
                                        double ice_thickness,
                                        double *shelf_base_salinity);

  PetscErrorCode subshelf_salinity_freeze_on(const POGivenTHConstants &constants,
                                             double sea_water_salinity,
                                             double sea_water_potential_temperature,
                                             double ice_thickness,
                                             double *shelf_base_salinity);

  PetscErrorCode subshelf_salinity_diffusion_only(const POGivenTHConstants &constants,
                                                  double sea_water_salinity,
                                                  double sea_water_potential_temperature,
                                                  double ice_thickness,
                                                  double *shelf_base_salinity);

  PetscErrorCode allocate_POGivenTH();

};

} // end of namespace pism

#endif /* _POGIVENTH_H_ */
