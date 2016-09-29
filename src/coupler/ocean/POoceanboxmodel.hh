// Copyright (C) 2012--2015 Ricarda Winkelmann, Ronja Reese and Torsten Albrecht
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

#ifndef _POOCEANBOXMODEL_H_
#define _POOCEANBOXMODEL_H_

#include "PGivenClimate.hh"
#include "POModifier.hh"
#include "Timeseries.hh"


//! A class defining the interface of a PISM ocean model modifier.

//! \brief A class implementing an ocean model.
//! Computes the subshelf melt/refreezing rate based on a simple ocean box model by Olbers & Hellmer (2010).

class POoceanboxmodel : public PGivenClimate<POModifier,PISMOceanModel>
{
public:
  POoceanboxmodel(IceGrid &g, const PISMConfig &conf);
  virtual ~POoceanboxmodel();

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(double my_t, double my_dt);

  virtual void add_vars_to_output(std::string keyword, std::set<std::string> &result); 

  virtual PetscErrorCode define_variables(std::set<std::string> vars, const PIO &nc,
                                          PISM_IO_Type nctype); 

  virtual PetscErrorCode write_variables(std::set<std::string> vars, const PIO& nc);

  virtual PetscErrorCode sea_level_elevation(double &result);

  virtual PetscErrorCode shelf_base_temperature(IceModelVec2S &result);

  virtual PetscErrorCode shelf_base_mass_flux(IceModelVec2S &result);

  virtual PetscErrorCode melange_back_pressure_fraction(IceModelVec2S &result);

  //////////////////////////////////////////////////////////////////////////////////////////


  class POBMConstants {
  public:
    POBMConstants(const PISMConfig &config);

      double        earth_grav,
                    rhoi, rhow, rho_star, nu,
                    latentHeat, c_p_ocean, lambda,
                    a, b, c,
                    alpha, beta;

      double        gamma_T, value_C,
                    T_dummy, S_dummy;

      double        gamma_T_o, meltFactor, meltSalinity, b2;
      double        continental_shelf_depth;

      int           numberOfBasins; 

  };

private:

  IceModelVec2S   shelfbtemp, 
                  shelfbmassflux;

  IceModelVec2T   *theta_ocean, 
                  *salinity_ocean;

  IceModelVec2S   *ice_thickness, 
                  *topg, 
                  *basins;  // not owned by this class 

  IceModelVec2Int *mask;  // not owned by this class


  virtual PetscErrorCode initBasinsOptions(const POBMConstants &constants);
  virtual PetscErrorCode roundBasins();
  virtual PetscErrorCode identifyMASK(IceModelVec2S &inputmask, std::string masktype);
  virtual PetscErrorCode computeOCEANMEANS();
  virtual PetscErrorCode extentOfIceShelves();
  virtual PetscErrorCode identifyBOXMODELmask();
  //virtual PetscErrorCode extendGLBox();
  //virtual PetscErrorCode extendIFBox();
  virtual PetscErrorCode oceanTemperature(const POBMConstants &constants);
  virtual PetscErrorCode basalMeltRateForGroundingLineBox(const POBMConstants &constants);
  virtual PetscErrorCode basalMeltRateForIceFrontBox(const POBMConstants &constants);
  virtual PetscErrorCode basalMeltRateForOtherShelves(const POBMConstants &constants);



  static const int  box_unidentified, 
                    box_noshelf, 
                    box_GL, 
                    box_neighboring, 
                    box_IF, 
                    box_other, 

                    maskfloating, 
                    maskocean, 
                    maskgrounded,

                    imask_inner,
                    imask_outer,
                    imask_exclude,
                    imask_unidentified;

  PetscScalar     counter_box_unidentified, 
                  counter_floating,
                  numberOfBoxes; // number of OBM-Boxes, there is one more box where Beckmann-Goose is computed..

  std::vector<double> Toc_base_vec,
                      Soc_base_vec,
                      gamma_T_star_vec,
                      C_vec,

                      mean_salinity_boundary_vector, //FIXME rename these, used at all boundaries
                      mean_temperature_boundary_vector,
                      mean_meltrate_boundary_vector,
                      mean_overturning_GLbox_vector; // execpt for the overturning...

  std::vector< std::vector<double> >  counter_boxes; 

  IceModelVec2S ICERISESmask, 
                BOXMODELmask,
                OCEANMEANmask, //FIXME delete when development finished
                DistGL,
                DistIF,
                Soc, 
                Soc_base, 
                Toc, 
                Toc_base, 
                Toc_inCelsius, 
                T_star, 
                Toc_anomaly, 
                overturning, 
                heatflux, 
                basalmeltrate_shelf;

  double        gamma_T, value_C,
                T_dummy, S_dummy,
                continental_shelf_depth;

  int      numberOfBasins; 


protected:

  Timeseries *delta_T;
  double delta_T_factor;
  double temp_anomaly;


  bool  ocean_oceanboxmodel_deltaT_set, 
        exicerises_set, 
        roundbasins_set, 
        continental_shelf_depth_set;


PetscErrorCode allocate_POoceanboxmodel();
};

#endif /* _POOCEANBOXMODEL_H_ */