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
// #include "base/util/iceModelVec.hh"

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

  bool ocean_oceanboxmodel_deltaT_set,
       omeans_set,
       roundbasins_set,
       continental_shelf_depth_set,
       exicerises_set;

private:
  IceModelVec2S   m_shelfbtemp,
                  m_shelfbmassflux,
                  cbasins,
                  ICERISESmask,
                  BOXMODELmask,
                  OCEANMEANmask, 
                  DistGL,
                  DistIF,
                  Soc,
                  Soc_base,
                  Toc,
                  Toc_base,
                  Toc_inCelsius,
                  T_star,
                  overturning,
                  heatflux,
                  basalmeltrate_shelf;

  // const IceModelVec2CellType m_mask;

  IceModelVec2T   *m_theta_ocean,
                  *m_salinity_ocean,
                  *basins;

  // The following are declared within POcavity.cc
  // IceModelVec2S   *ice_thickness;
                  // *topg;  // not owned by this class

  IceModelVec2Int *mask;  // not owned by this class


  void initBasinsOptions(const Constants &constants);
  void round_basins();
  void identifyMASK(IceModelVec2S &inputmask, std::string masktype);
  void computeOCEANMEANS();
  void extentOfIceShelves();
  void identifyBOXMODELmask();
  void oceanTemperature(const Constants &constants);
  void basalMeltRateForGroundingLineBox(const Constants &constants);
  void basalMeltRateForIceFrontBox(const Constants &constants);
  void basalMeltRateForOtherShelves(const Constants &constants);
  double most_frequent_element(const std::vector<double>&);
  // FIXME: move these declarations where they are also be initiated?
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
                    imask_unidentified,
                    numberOfBoxes;

  double counter_box_unidentified,
         counter_floating;
          // number of OBM-Boxes, there is one more box where Beckmann-Goose is computed..
         // numberOfBoxes; // FIXME: Why should these be doubles and not ints?

  std::vector<double> Toc_base_vec,
                      Soc_base_vec,
                      gamma_T_star_vec,
                      C_vec,
                      mean_salinity_boundary_vector, //FIXME rename these, used at all boundaries
                      mean_temperature_boundary_vector,
                      mean_meltrate_boundary_vector,
                      mean_overturning_GLbox_vector; // execpt for the overturning...

  std::vector< std::vector<double> >  counter_boxes;

  double        gamma_T, value_C,
                T_dummy, S_dummy,
                continental_shelf_depth;

  int      numberOfBasins,
           Mx, My, xs, xm, ys, ym, dx, dy;
};

} // end of namespace ocean
} // end of namespace pism

#endif /* _POCAVITY_H_ */
