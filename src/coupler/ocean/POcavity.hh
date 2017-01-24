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

#ifndef _POCAVITY_H_
#define _POCAVITY_H_

#include "coupler/util/PGivenClimate.hh"
#include "POModifier.hh"
#include "base/util/IceModelVec2CellType.hh"

namespace pism {
namespace ocean {

class Cavity : public PGivenClimate<OceanModifier,OceanModel> {
public:
  Cavity(IceGrid::ConstPtr g);
  virtual ~Cavity();

  class Constants {
  public:
    Constants(const Config &config);

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

protected:
  virtual void update_impl(double my_t, double my_dt);
  virtual void write_variables_impl(const std::set<std::string> &vars, const PIO& nc);
  virtual void add_vars_to_output_impl(const std::string &keyword, std::set<std::string> &result);
  virtual void define_variables_impl(const std::set<std::string> &vars,
                                     const PIO &nc, IO_Type nctype);
  virtual void init_impl();
  virtual void melange_back_pressure_fraction_impl(IceModelVec2S &result) const;
  virtual void sea_level_elevation_impl(double &result) const;
  virtual void shelf_base_temperature_impl(IceModelVec2S &result) const;
  virtual void shelf_base_mass_flux_impl(IceModelVec2S &result) const;

  std::vector<IceModelVec*> m_variables;

  bool   exicerises_set; // FIXME shouldn't this be always used?

private:
  IceModelVec2S   m_shelfbtemp,
                  m_shelfbmassflux,
                  cbasins, // a basin defines the domain where one box model instance is solved
                  ICERISESmask, 
                  BOXMODELmask,
                  OCEANMEANmask, 
                  OCEANmask,
                  DistGL,
                  DistIF,
                  Soc,
                  Soc_base,
                  Toc,
                  Toc_base,
                  T_star,
                  overturning,
                  basalmeltrate_shelf;

  IceModelVec2T   *m_theta_ocean,
                  *m_salinity_ocean;

  void initBasinsOptions(const Constants &constants);
  void round_basins();
  void identifyMASK(IceModelVec2S &inputmask, std::string masktype);
  void computeOCEANMEANS();
  void extentOfIceShelves();
  void identifyBOXMODELmask();
  void oceanTemperature(const Constants &constants);
  void basalMeltRateGroundingLineBox(const Constants &constants);
  void basalMeltRateOtherBoxes(const Constants &constants);
  void basalMeltRateMissingCells(const Constants &constants);
  double most_frequent_element(const std::vector<double>&);

  static const int  numberOfBoxes, // max number of ocean boxes (reached for big ice shelves)

                    box1, // ocean box covering the grounding line region
                    box2, // ocean box neighboring the box 1, other boxes are covered by boxi

                    maskfloating,
                    maskocean,
                    maskgrounded,

                    imask_inner, // used in IdentifyMask 
                    imask_outer, // used in IdentifyMask 
                    imask_exclude, // used in IdentifyMask 
                    imask_unidentified; // used in IdentifyMask 
                    
  std::vector<double> Toc_base_vec, // temperature input for box 1 per basin
                      Soc_base_vec, // salinity input for box 1 per basin
                      gamma_T_star_vec, // FIXME delete 
                      C_vec, // FIXME delete 
                      mean_salinity_boundary_vector, // salinity input for box i>1 per basin
                      mean_temperature_boundary_vector, // temperature input for box i>1 per basin
                      mean_meltrate_boundary_vector, // mean melt rate in box i-1 as input for box i>1 per basin
                      mean_overturning_GLbox_vector; // mean overturning, computed in box 1, as input for box i>1 per basin

  std::vector< std::vector<double> >  counter_boxes; // matrix containing the number of shelf cells per basin and box
                                                     // used for area calculation

  // standard values are defined in Constants
  // here needed to store custom values from user options.            
  double        gamma_T, value_C, 
                continental_shelf_depth;

  int      numberOfBasins,
           Mx, My, dx, dy;
};

} // end of namespace ocean
} // end of namespace pism

#endif /* _POCAVITY_H_ */
