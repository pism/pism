// Copyright (C) 2012-2016, 2018 Ricarda Winkelmann, Ronja Reese, Torsten Albrecht
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

#ifndef _POPICO_H_
#define _POPICO_H_

#include "CompleteOceanModel.hh"

#include "pism/coupler/util/PGivenClimate.hh"
#include "pism/util/IceModelVec2CellType.hh"

namespace pism {
namespace ocean {
//! \brief Implements the PICO ocean model as submitted to The Cryosphere (March 2017).
//!
//! Generalizes the two dimensional ocean box model of [@ref OlbersHellmer2010] for
//! use in PISM, i.e. three dimensions.
//!
class Pico : public PGivenClimate<CompleteOceanModel, CompleteOceanModel> {
public:
  Pico(IceGrid::ConstPtr g);
  virtual ~Pico();

  class Constants {
  public:
    Constants(const Config &config);

    double earth_grav, rhoi, rhow, rho_star, nu, latentHeat, c_p_ocean, lambda, a, b, c, as, bs, cs, alpha, beta;

    double default_gamma_T, default_overturning_coeff, T_dummy, S_dummy;

    double gamma_T_o, meltFactor, meltSalinity, b2;
    double continental_shelf_depth;

    int default_numberOfBasins, default_numberOfBoxes;
  };

protected:
  void update_impl(double t, double dt);
  void init_impl();

  void define_model_state_impl(const PIO &output) const;
  void write_model_state_impl(const PIO &output) const;
  void test();

  std::map<std::string, Diagnostic::Ptr> diagnostics_impl() const;

  bool m_exicerises_set; // FIXME shouldn't this be always used?

private:
  IceModelVec2S m_cbasins, // a basin defines the domain where one box model instance is solved
      m_icerise_mask, m_ocean_box_mask, m_shelf_mask, m_ocean_contshelf_mask, m_ocean_mask, m_lake_mask, m_DistGL,
      m_DistIF, m_Soc, m_Soc_box0, m_Toc, m_Toc_box0, m_T_star, m_overturning, m_basalmeltrate_shelf,
      m_T_pressure_melting;

  IceModelVec2T *m_theta_ocean, *m_salinity_ocean;

  void initBasinsOptions(const Constants &constants);
  void round_basins();
  void identifyMASK(IceModelVec2S &inputmask, std::string masktype);
  void identify_shelf_mask();
  void compute_ocean_input_per_basin(const Constants &constants);
  void compute_distances();
  void identify_ocean_box_mask(const Constants &constants);
  void set_ocean_input_fields(const Constants &constants);
  void calculate_basal_melt_box1(const Constants &constants);
  void calculate_basal_melt_other_boxes(const Constants &constants);
  void calculate_basal_melt_missing_cells(const Constants &constants);
  double most_frequent_element(const std::vector<double> &);

  enum IdentifyMaskFlags {INNER = 2, OUTER = 0, EXCLUDE = 1, UNIDENTIFIED = -1};

  std::vector<double> m_Toc_box0_vec,     // temperature input for box 1 per basin
      m_Soc_box0_vec,                     // salinity input for box 1 per basin
      m_mean_salinity_boundary_vector,    // salinity input for box i>1 per basin
      m_mean_temperature_boundary_vector, // temperature input for box i>1 per basin
      m_mean_overturning_box1_vector;     // mean overturning, computed in box 1, as input for box i>1 per basin

  std::vector<std::vector<double> > counter_boxes; // matrix containing the number of shelf cells per basin and box
                                                   // used for area calculation

  // standard values are defined in Constants
  // here needed to store custom values from user options.
  double m_gamma_T, m_overturning_coeff, m_continental_shelf_depth;

  int m_numberOfBasins, m_numberOfBoxes, m_numberOfShelves, m_Mx, m_My, m_dx, m_dy;
};


} // end of namespace ocean
} // end of namespace pism

#endif /* _POPICO_H_ */
