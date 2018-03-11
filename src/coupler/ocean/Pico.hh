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

double f_area(double counter_boxes, double m_dx, double m_dy);

struct TocBox1 {
  bool failed;
  double value;
};

class BoxModel {
public:
  BoxModel(const Config &config);

  double pressure(double ice_thickness) const;
  double T_star(double salinity, double temperature, double pressure) const;

  TocBox1 Toc_box1(double area, double T_star, double Soc_box0, double Toc_box0) const;
  double Toc_other_boxes(double area,
                         double temp_in_boundary, double T_star,
                         double overturning, double salinity_in_boundary) const;

  double Soc_box1(double Toc_box0, double Soc_box0, double Toc) const;
  double Soc_other_boxes(double salinity_in_boundary, double temperature_in_boundary, double Toc) const;

  double pot_pressure_melting(double salinity, double pressure) const;
  double pressure_melting(double salinity, double pressure) const;
  double bmelt_rate(double pm_point, double Toc) const;
  double bmelt_rate_beckm_goose(double Toc, double pot_pm_point) const;
  double overturning(double Soc_box0, double Soc, double Toc_box0,
                       double Toc) const;
  double p_coeff(double g1, double s1) const;
  double q_coeff(double g1, double s1, double T_star) const;

  double gamma_T() const;
  double overturning_coeff() const;
  double T_dummy() const;
  double S_dummy() const;
  double ice_density() const;
  double continental_shelf_depth() const;
private:
private:
  double m_gamma_T, m_overturning_coeff, m_T_dummy, m_S_dummy;
  double m_ice_density, m_continental_shelf_depth;


  double m_earth_grav, m_sea_water_density, m_rho_star, m_nu, m_latentHeat, m_c_p_ocean, m_alpha, m_beta;

  double m_lambda;

  // coefficients of the parameterization of the potential temperature
  double m_a_pot, m_b_pot, m_c_pot;

  // coefficients of the parameterization of the in situ temperature
  double m_b_in_situ, m_c_in_situ, m_a_in_situ;

  double m_meltFactor;
};

//! \brief Implements the PICO ocean model as submitted to The Cryosphere (March 2017).
//!
//! Generalizes the two dimensional ocean box model of [@ref OlbersHellmer2010] for
//! use in PISM, i.e. three dimensions.
//!
class Pico : public PGivenClimate<CompleteOceanModel, CompleteOceanModel> {
public:
  Pico(IceGrid::ConstPtr g);
  virtual ~Pico();

protected:
  void update_impl(double t, double dt);
  void init_impl();

  void define_model_state_impl(const PIO &output) const;
  void write_model_state_impl(const PIO &output) const;
  void test();

  std::map<std::string, Diagnostic::Ptr> diagnostics_impl() const;

  bool m_exicerises_set; // FIXME shouldn't this be always used?

private:
  IceModelVec2S m_DistGL,
      m_DistIF, m_Soc, m_Soc_box0, m_Toc, m_Toc_box0, m_T_star, m_overturning, m_basalmeltrate_shelf,
      m_T_pressure_melting;

  // a basin defines the domain where one box model instance is solved
  IceModelVec2Int m_icerise_mask, m_cbasins, m_shelf_mask, m_lake_mask,
    m_ocean_box_mask, m_ocean_mask, m_ocean_contshelf_mask;

  IceModelVec2T *m_theta_ocean, *m_salinity_ocean;

  void initBasinsOptions(const BoxModel &constants);
  void identifyMASK(IceModelVec2S &inputmask, std::string masktype);
  void identify_shelf_mask();
  void compute_ocean_input_per_basin(const BoxModel &cc,
                                     const IceModelVec2Int &basin_mask,
                                     const IceModelVec2Int &continental_shelf_mask,
                                     const IceModelVec2S &salinity_ocean,
                                     const IceModelVec2S &theta_ocean);
  void compute_distances();
  void identify_ocean_box_mask(const BoxModel &constants);
  void set_ocean_input_fields(const IceModelVec2S &ice_thickness,
                              const IceModelVec2CellType &mask,
                              const IceModelVec2Int &m_cbasins,
                              const IceModelVec2Int &m_shelf_mask,
                              const BoxModel &cc,
                              IceModelVec2S &Toc_box0,
                              IceModelVec2S &Soc_box0
                              );
  void calculate_basal_melt_box1(const IceModelVec2S &ice_thickness,
                                 const IceModelVec2Int &shelf_mask,
                                 const IceModelVec2Int &box_mask,
                                 const IceModelVec2S &Toc_box0,
                                 const IceModelVec2S &Soc_box0,
                                 const BoxModel &cc,
                                 IceModelVec2S &T_star,
                                 IceModelVec2S &Toc,
                                 IceModelVec2S &Soc,
                                 IceModelVec2S &basal_melt_rate,
                                 IceModelVec2S &overturning,
                                 IceModelVec2S &T_pressure_melting);
  void calculate_basal_melt_other_boxes(const IceModelVec2S &ice_thickness,
                                        const IceModelVec2Int &shelf_mask,
                                        const BoxModel &cc,
                                        IceModelVec2Int &box_mask,
                                        IceModelVec2S &T_star,
                                        IceModelVec2S &Toc,
                                        IceModelVec2S &Soc,
                                        IceModelVec2S &basal_melt_rate,
                                        IceModelVec2S &T_pressure_melting);
  void calculate_basal_melt_missing_cells(const BoxModel &cc,
                                          const IceModelVec2Int &shelf_mask,
                                          const IceModelVec2Int &box_mask,
                                          const IceModelVec2S &ice_thickness,
                                          const IceModelVec2S &Toc_box0,
                                          const IceModelVec2S &Soc_box0,
                                          IceModelVec2S &Toc,
                                          IceModelVec2S &Soc,
                                          IceModelVec2S &basal_melt_rate,
                                          IceModelVec2S &T_pressure_melting);

  enum IdentifyMaskFlags {INNER = 2, OUTER = 0, EXCLUDE = 1, UNIDENTIFIED = -1};

  std::vector<double> m_Toc_box0_vec,     // temperature input for box 1 per basin
      m_Soc_box0_vec,                     // salinity input for box 1 per basin
      m_mean_salinity_boundary_vector,    // salinity input for box i>1 per basin
      m_mean_temperature_boundary_vector, // temperature input for box i>1 per basin
      m_mean_overturning_box1_vector;     // mean overturning, computed in box 1, as input for box i>1 per basin

  std::vector<std::vector<double> > counter_boxes; // matrix containing the number of shelf cells per basin and box
                                                   // used for area calculation

  int m_numberOfBasins, m_numberOfBoxes, m_numberOfShelves, m_Mx, m_My, m_dx, m_dy;
};

void round_basins(IceModelVec2S &basin_mask);

} // end of namespace ocean
} // end of namespace pism

#endif /* _POPICO_H_ */
