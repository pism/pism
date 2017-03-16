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
//! \brief Implements the SIMPEL ocean model as submitted to The Cryosphere (March 2017).
//!
//! Generalizes the two dimensional ocean box model of [@ref OlbersHellmer2010] for
//! use in PISM, i.e. three dimensions.
//!
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
                    a, b, c, as, bs, cs,
                    alpha, beta;

      double        default_gamma_T, default_overturning_coeff,
                    T_dummy, S_dummy;

      double        gamma_T_o, meltFactor, meltSalinity, b2;
      double        continental_shelf_depth;

      int           default_numberOfBasins, default_numberOfBoxes;

  };

protected:
  virtual void update_impl(double my_t, double my_dt);
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
                  icerise_mask,
                  ocean_box_mask,
                  ocean_contshelf_mask,
                  ocean_mask,
                  DistGL,
                  DistIF,
                  Soc,
                  Soc_box0,
                  Toc,
                  Toc_box0,
                  T_star,
                  overturning,
                  basalmeltrate_shelf,
                  T_pressure_melting;

  IceModelVec2T   *m_theta_ocean,
                  *m_salinity_ocean;

  void initBasinsOptions(const Constants &constants);
  void round_basins();
  void identifyMASK(IceModelVec2S &inputmask, std::string masktype);
  void compute_ocean_input_per_basin(const Constants &constants);
  void compute_distances();
  void identify_ocean_box_mask(const Constants &constants);
  void set_ocean_input_fields(const Constants &constants);
  void calculate_basal_melt_box1(const Constants &constants);
  void calculate_basal_melt_other_boxes(const Constants &constants);
  void calculate_basal_melt_missing_cells(const Constants &constants);
  double most_frequent_element(const std::vector<double>&);

  static const int  maskfloating,
                    maskocean,
                    maskgrounded,

                    imask_inner, // used in IdentifyMask
                    imask_outer, // used in IdentifyMask
                    imask_exclude, // used in IdentifyMask
                    imask_unidentified; // used in IdentifyMask

  std::vector<double> Toc_box0_vec, // temperature input for box 1 per basin
                      Soc_box0_vec, // salinity input for box 1 per basin
                      mean_salinity_boundary_vector, // salinity input for box i>1 per basin
                      mean_temperature_boundary_vector, // temperature input for box i>1 per basin
                      mean_overturning_box1_vector; // mean overturning, computed in box 1, as input for box i>1 per basin

  std::vector< std::vector<double> >  counter_boxes; // matrix containing the number of shelf cells per basin and box
                                                     // used for area calculation

  // standard values are defined in Constants
  // here needed to store custom values from user options.
  double        gamma_T, overturning_coeff,
                continental_shelf_depth;

  int      numberOfBasins, numberOfBoxes,
           Mx, My, dx, dy;
};

} // end of namespace ocean
} // end of namespace pism

#endif /* _POCAVITY_H_ */
