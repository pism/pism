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

#include "pism/util/IceModelVec2CellType.hh"
#include "pism/util/iceModelVec2T.hh"

namespace pism {
namespace ocean {

class PicoGeometry;
class PicoPhysics;

//! Implements the PICO ocean model as submitted to The Cryosphere (March 2017).
//!
//! Generalizes the two dimensional ocean box model of [@ref OlbersHellmer2010] for
//! use in PISM, i.e. three dimensions.
//!
class Pico : public CompleteOceanModel {
public:
  Pico(IceGrid::ConstPtr g);
  virtual ~Pico();

protected:
  void init_impl(const Geometry &geometry);
  void update_impl(const Geometry &geometry, double t, double dt);
  MaxTimestep max_timestep_impl(double t) const;

  void define_model_state_impl(const PIO &output) const;
  void write_model_state_impl(const PIO &output) const;

  std::map<std::string, Diagnostic::Ptr> diagnostics_impl() const;

private:
  IceModelVec2S m_Soc, m_Soc_box0;
  IceModelVec2S m_Toc, m_Toc_box0, m_T_star;
  IceModelVec2S m_overturning;
  IceModelVec2S m_basal_melt_rate;

  IceModelVec2Int m_basin_mask;

  std::unique_ptr<PicoGeometry> m_geometry;

  IceModelVec2T::Ptr m_theta_ocean, m_salinity_ocean;

  void compute_ocean_input_per_basin(const PicoPhysics &physics,
                                     const IceModelVec2Int &basin_mask,
                                     const IceModelVec2Int &continental_shelf_mask,
                                     const IceModelVec2S &salinity_ocean,
                                     const IceModelVec2S &theta_ocean,
                                     std::vector<double> &temperature,
                                     std::vector<double> &salinity);

  void set_ocean_input_fields(const PicoPhysics &physics,
                              const IceModelVec2S &ice_thickness,
                              const IceModelVec2CellType &mask,
                              const IceModelVec2Int &basin_mask,
                              const IceModelVec2Int &shelf_mask,
                              const std::vector<double> basin_temperature,
                              const std::vector<double> basin_salinity,
                              IceModelVec2S &Toc_box0,
                              IceModelVec2S &Soc_box0);

  void process_box1(const PicoPhysics &physics,
                    const IceModelVec2S &ice_thickness,
                    const IceModelVec2Int &shelf_mask,
                    const IceModelVec2Int &box_mask,
                    const IceModelVec2S &Toc_box0,
                    const IceModelVec2S &Soc_box0,
                    IceModelVec2S &basal_melt_rate,
                    IceModelVec2S &basal_temperature,
                    IceModelVec2S &T_star,
                    IceModelVec2S &Toc,
                    IceModelVec2S &Soc,
                    IceModelVec2S &overturning);

  void process_other_boxes(const PicoPhysics &cc,
                           const IceModelVec2S &ice_thickness,
                           const IceModelVec2Int &shelf_mask,
                           const IceModelVec2Int &box_mask,
                           IceModelVec2S &basal_melt_rate,
                           IceModelVec2S &basal_temperature,
                           IceModelVec2S &T_star,
                           IceModelVec2S &Toc,
                           IceModelVec2S &Soc);

  void beckmann_goosse(const PicoPhysics &physics,
                       const IceModelVec2S &ice_thickness,
                       const IceModelVec2Int &shelf_mask,
                       const IceModelVec2CellType &cell_type,
                       const IceModelVec2S &Toc_box0,
                       const IceModelVec2S &Soc_box0,
                       IceModelVec2S &basal_melt_rate,
                       IceModelVec2S &T_pressure_melting,
                       IceModelVec2S &Toc,
                       IceModelVec2S &Soc);

  void compute_box_average(int box_id,
                           const IceModelVec2S &field,
                           const IceModelVec2Int &shelf_mask,
                           const IceModelVec2Int &box_mask,
                           std::vector<double> &result);

  void compute_box_area(int box_id,
                        const IceModelVec2Int &shelf_mask,
                        const IceModelVec2Int &box_mask,
                        std::vector<double> &result);


  int m_n_basins, m_n_boxes, m_n_shelves;
};

} // end of namespace ocean
} // end of namespace pism

#endif /* _POPICO_H_ */
