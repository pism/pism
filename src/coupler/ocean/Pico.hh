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

class PicoConstants {
public:
  PicoConstants(const Config &config);

  double earth_grav, rhoi, rhow, rho_star, nu, latentHeat, c_p_ocean, alpha, beta;

  double lambda_;

  // coefficients of the parameterization of the potential temperature
  double a_pot, b_pot, c_pot;

  // coefficients of the parameterization of the in situ temperature
  double b_in_situ, c_in_situ, a_in_situ;

  double gamma_T, overturning_coeff, T_dummy, S_dummy;

  double gamma_T_o, meltFactor, meltSalinity, b2;
  double continental_shelf_depth;

  int default_numberOfBasins, default_numberOfBoxes;
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

  void initBasinsOptions(const PicoConstants &constants);
  void identifyMASK(IceModelVec2S &inputmask, std::string masktype);
  void identify_shelf_mask();
  void compute_ocean_input_per_basin(const PicoConstants &cc,
                                     const IceModelVec2Int &basin_mask,
                                     const IceModelVec2Int &continental_shelf_mask,
                                     const IceModelVec2S &salinity_ocean,
                                     const IceModelVec2S &theta_ocean);
  void compute_distances();
  void identify_ocean_box_mask(const PicoConstants &constants);
  void set_ocean_input_fields(const IceModelVec2S &ice_thickness,
                              const IceModelVec2CellType &mask,
                              const IceModelVec2Int &m_cbasins,
                              const IceModelVec2Int &m_shelf_mask,
                              const PicoConstants &cc,
                              IceModelVec2S &Toc_box0,
                              IceModelVec2S &Soc_box0
                              );
  void calculate_basal_melt_box1(const IceModelVec2S &ice_thickness,
                                 const IceModelVec2Int &shelf_mask,
                                 const IceModelVec2Int &box_mask,
                                 const IceModelVec2S &Toc_box0,
                                 const IceModelVec2S &Soc_box0,
                                 const PicoConstants &cc,
                                 IceModelVec2S &T_star,
                                 IceModelVec2S &Toc,
                                 IceModelVec2S &Soc,
                                 IceModelVec2S &basal_melt_rate,
                                 IceModelVec2S &overturning,
                                 IceModelVec2S &T_pressure_melting);
  void calculate_basal_melt_other_boxes(const IceModelVec2S &ice_thickness,
                                        const IceModelVec2Int &shelf_mask,
                                        const PicoConstants &cc,
                                        IceModelVec2Int &box_mask,
                                        IceModelVec2S &T_star,
                                        IceModelVec2S &Toc,
                                        IceModelVec2S &Soc,
                                        IceModelVec2S &basal_melt_rate,
                                        IceModelVec2S &T_pressure_melting);
  void calculate_basal_melt_missing_cells(const PicoConstants &cc,
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

double f_pressure(const PicoConstants &cc, double ice_thickness);
double f_T_star(const PicoConstants &cc, double salinity, double temperature, double pressure);
double f_area(double counter_boxes, double m_dx, double m_dy);

double f_Toc_box1(const PicoConstants &cc, double area, double T_star, double Soc_box0, double Toc_box0, bool *success);
double f_Toc_other_boxes(const PicoConstants &cc, double area, double gamma_T,
  double temp_in_boundary, double T_star,
    double overturning, double salinity_in_boundary);

double f_Soc_box1(const PicoConstants &cc, double Toc_box0, double Soc_box0, double Toc);
double f_Soc_other_boxes(const PicoConstants &cc, double salinity_in_boundary, double temperature_in_boundary, double Toc);

double f_pot_pressure_melting(const PicoConstants &cc, double salinity, double pressure);
double f_pressure_melting(const PicoConstants &cc, double salinity, double pressure);
double f_bmelt_rate(const PicoConstants &cc, double m_gamma_T, double pm_point, double Toc);
double f_bmelt_rate_beckm_goose(const PicoConstants &cc, double Toc, double pot_pm_point);
double f_overturning(const PicoConstants &cc, double overturning_coeff, double Soc_box0,
    double Soc, double Toc_box0, double Toc);
double f_p_coeff(const PicoConstants &cc, double g1, double overturning_coeff,
    double s1);
double f_q_coeff(const PicoConstants &cc, double g1, double overturning_coeff,
    double s1, double T_star);

} // end of namespace ocean
} // end of namespace pism

#endif /* _POPICO_H_ */
