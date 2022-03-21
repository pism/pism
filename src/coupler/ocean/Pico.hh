// Copyright (C) 2012-2016, 2018, 2020, 2021, 2022 Ricarda Winkelmann, Ronja Reese, Torsten Albrecht
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

#include "pism/util/array/Forcing.hh"
#include "PicoGeometry.hh"

namespace pism {

namespace ocean {

class PicoPhysics;

//! Implements the PICO ocean model as submitted to The Cryosphere (March 2017).
//!
//! Generalizes the two dimensional ocean box model of [@ref OlbersHellmer2010] for
//! use in PISM, i.e. three dimensions.
//!
class Pico : public CompleteOceanModel {
public:
  Pico(IceGrid::ConstPtr g);
  virtual ~Pico() = default;

protected:
  void init_impl(const Geometry &geometry);
  void update_impl(const Geometry &geometry, double t, double dt);
  MaxTimestep max_timestep_impl(double t) const;

  void define_model_state_impl(const File &output) const;
  void write_model_state_impl(const File &output) const;

  std::map<std::string, Diagnostic::Ptr> diagnostics_impl() const;

private:
  array::Scalar m_Soc, m_Soc_box0;
  array::Scalar m_Toc, m_Toc_box0, m_T_star;
  array::Scalar m_overturning;
  array::Scalar1 m_basal_melt_rate;

  PicoGeometry m_geometry;

  std::shared_ptr<array::Forcing> m_theta_ocean, m_salinity_ocean;

  void compute_ocean_input_per_basin(const PicoPhysics &physics,
                                     const array::Scalar &basin_mask,
                                     const array::Scalar &continental_shelf_mask,
                                     const array::Scalar &salinity_ocean,
                                     const array::Scalar &theta_ocean,
                                     std::vector<double> &temperature,
                                     std::vector<double> &salinity) const;

  void set_ocean_input_fields(const PicoPhysics &physics,
                              const array::Scalar &ice_thickness,
                              const array::CellType1 &mask,
                              const array::Scalar &basin_mask,
                              const array::Scalar &shelf_mask,
                              const std::vector<double> &basin_temperature,
                              const std::vector<double> &basin_salinity,
                              array::Scalar &Toc_box0,
                              array::Scalar &Soc_box0) const;

  void process_box1(const PicoPhysics &physics,
                    const array::Scalar &ice_thickness,
                    const array::Scalar &shelf_mask,
                    const array::Scalar &box_mask,
                    const array::Scalar &Toc_box0,
                    const array::Scalar &Soc_box0,
                    array::Scalar &basal_melt_rate,
                    array::Scalar &basal_temperature,
                    array::Scalar &T_star,
                    array::Scalar &Toc,
                    array::Scalar &Soc,
                    array::Scalar &overturning);

  void process_other_boxes(const PicoPhysics &physics,
                           const array::Scalar &ice_thickness,
                           const array::Scalar &shelf_mask,
                           const array::Scalar &box_mask,
                           array::Scalar &basal_melt_rate,
                           array::Scalar &basal_temperature,
                           array::Scalar &T_star,
                           array::Scalar &Toc,
                           array::Scalar &Soc) const;

  void beckmann_goosse(const PicoPhysics &physics,
                       const array::Scalar &ice_thickness,
                       const array::Scalar &shelf_mask,
                       const array::CellType0 &cell_type,
                       const array::Scalar &Toc_box0,
                       const array::Scalar &Soc_box0,
                       array::Scalar &basal_melt_rate,
                       array::Scalar &basal_temperature,
                       array::Scalar &Toc,
                       array::Scalar &Soc);

  void compute_box_average(int box_id,
                           const array::Scalar &field,
                           const array::Scalar &shelf_mask,
                           const array::Scalar &box_mask,
                           std::vector<double> &result) const;

  void compute_box_area(int box_id,
                        const array::Scalar &shelf_mask,
                        const array::Scalar &box_mask,
                        std::vector<double> &result) const;

  int m_n_basins, m_n_boxes, m_n_shelves;
};

} // end of namespace ocean
} // end of namespace pism

#endif /* _POPICO_H_ */
