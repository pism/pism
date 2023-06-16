/* Copyright (C) 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2023 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#ifndef _PBLINGLECLARK_H_
#define _PBLINGLECLARK_H_

#include <memory>               // std::unique_ptr

#include "BedDef.hh"

namespace pism {
namespace bed {

class LingleClarkSerial;

//! A wrapper class around LingleClarkSerial.
class LingleClark : public BedDef {
public:
  LingleClark(IceGrid::ConstPtr g);
  virtual ~LingleClark();

  const array::Scalar& total_displacement() const;

  const array::Scalar& viscous_displacement() const;

  const array::Scalar& elastic_displacement() const;

  const array::Scalar& relief() const;

  void step(const array::Scalar &ice_thickness,
            const array::Scalar &sea_level_elevation,
            double dt);

  array::Scalar::Ptr elastic_load_response_matrix() const;
protected:
  virtual void define_model_state_impl(const File &output) const;
  virtual void write_model_state_impl(const File &output) const;

  DiagnosticList diagnostics_impl() const;

  MaxTimestep max_timestep_impl(double t) const;
  void init_impl(const InputOptions &opts, const array::Scalar &ice_thickness,
                 const array::Scalar &sea_level_elevation);
  void bootstrap_impl(const array::Scalar &bed_elevation,
                      const array::Scalar &bed_uplift,
                      const array::Scalar &ice_thickness,
                      const array::Scalar &sea_level_elevation);
  void update_impl(const array::Scalar &ice_thickness,
                   const array::Scalar &sea_level_elevation,
                   double t, double dt);

  //! Total (viscous and elastic) bed displacement.
  array::Scalar m_total_displacement;

  //! Storage on rank zero. Used to pass the load to the serial deformation model and get
  //! bed displacement back.
  std::shared_ptr<petsc::Vec> m_work0;

  //! Bed relief relative to the bed displacement.
  array::Scalar m_relief;

  //! Ice-equivalent load thickness.
  array::Scalar m_load_thickness;

  //! Serial viscoelastic bed deformation model.
  std::unique_ptr<LingleClarkSerial> m_serial_model;

  //! extended grid for the viscous plate displacement
  std::shared_ptr<IceGrid> m_extended_grid;

  //! Viscous displacement on the extended grid (part of the model state).
  std::shared_ptr<array::Scalar> m_viscous_displacement;
  //! rank 0 storage using the extended grid
  std::shared_ptr<petsc::Vec> m_viscous_displacement0;

  //! Elastic bed displacement (part of the model state)
  array::Scalar m_elastic_displacement;
  //! rank 0 storage for the elastic displacement
  std::shared_ptr<petsc::Vec> m_elastic_displacement0;

  //! time of the last bed deformation update
  double m_t_last;
  //! Update interval in seconds
  double m_update_interval;
  //! Temporal resolution to use when checking whether it's time to update
  double m_t_eps;
  //! Name of the variable used to store the last update time.
  std::string m_time_name;
};

} // end of namespace bed
} // end of namespace pism

#endif /* _PBLINGLECLARK_H_ */
