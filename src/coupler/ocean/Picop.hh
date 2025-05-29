// Copyright (C) 2012-2016, 2018, 2020, 2021, 2022, 2023, 2025 Ricarda Winkelmann, Ronja Reese, Torsten Albrecht
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

#ifndef _POPICOP_H_
#define _POPICOP_H_

#include "pism/geometry/UNO.hh"
#include "pism/coupler/ocean/CompleteOceanModel.hh"
#include "pism/coupler/ocean/Pico.hh"
#include "pism/coupler/ocean/PicoGeometry.hh"

namespace pism {

namespace ocean {

class PicopPhysics;

//! Implements the PICO ocean model published in The Cryosphere (2018)
//! and adds the Plume model Pelle et al (2019).
//!
class Picop : public CompleteOceanModel {
public:
  Picop(std::shared_ptr<const Grid> g);
  virtual ~Picop() = default;

protected:
  void init_impl(const Geometry &geometry);
  void update_impl(const Inputs &inputs, double t, double dt);
  MaxTimestep max_timestep_impl(double t) const;

  void define_model_state_impl(const File &output) const;
  void write_model_state_impl(const File &output) const;

  std::map<std::string, Diagnostic::Ptr> diagnostics_impl() const;

private:
  
  std::shared_ptr<Pico> m_pico;
  std::unique_ptr<pism::UNO> m_uno;

  array::Scalar1 m_basal_melt_rate;
  array::Scalar1 m_grounding_line_elevation;
  array::Scalar1 m_shelf_base_elevation;
  array::Scalar1 m_slope;
  
  const array::Scalar &m_theta_ocean;
  const array::Scalar &m_salinity_ocean;
  
  // array::Scalar m_work;
  array::Vector m_flow_direction;
  
  void compute_melt_rate(const PicopPhysics &physics,
                         const array::Scalar &T_a,
                         const array::Scalar &S_a,
                         array::Scalar1 &grounding_line_elevation,
                         array::Scalar1 &shelf_base_elevation,
                         array::Scalar1 &slope,
                         array::Scalar1 &basal_melt_rate) const;
  
  void compute_grounding_line_elevation(const Inputs &inputs,
                                        array::Scalar1 &result);
  
  void compute_shelf_base_elevation(const Geometry &geometry,
                                    array::Scalar1 &result) const;

  void compute_slope(const Geometry &geometry,
                     array::Scalar1 &shelf_base_elevation,
                     array::Scalar &slope) const;

};
} // end of namespace ocean
} // end of namespace pism

#endif /* _POPICOP_H_ */
