// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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

#ifndef __BedDef_hh
#define __BedDef_hh

#include "pism/util/Component.hh"

namespace pism {

//! @brief Bed-related models: bed deformation (provide bed elevation
//! and uplift) and (soon) bed erosion.
namespace bed {

double compute_load(double bed, double ice_thickness, double sea_level,
                    double ice_density, double ocean_density);

void accumulate_load(const array::Scalar &bed_elevation, const array::Scalar &ice_thickness,
                     const array::Scalar &sea_level_elevation, double C, array::Scalar &result);

//! PISM bed deformation model (base class).
class BedDef : public Component {
public:
  BedDef(std::shared_ptr<const Grid> g, const std::string &model_name);
  virtual ~BedDef() = default;

  void init(const InputOptions &opts, const array::Scalar &ice_thickness,
            const array::Scalar &sea_level_elevation);
  void bootstrap(const array::Scalar &bed_elevation,
                 const array::Scalar &bed_uplift,
                 const array::Scalar &ice_thickness,
                 const array::Scalar &sea_level_elevation);

  void update(const array::Scalar &ice_thickness,
              const array::Scalar &sea_level_elevation,
              double t, double dt);

  const array::Scalar& bed_elevation() const;
  const array::Scalar& uplift() const;

protected:
  virtual MaxTimestep max_timestep_impl(double t) const;

  virtual void define_model_state_impl(const File &output) const;
  virtual void write_model_state_impl(const File &output) const;

  virtual DiagnosticList diagnostics_impl() const;

  virtual void update_impl(const array::Scalar &load,
                           double t, double dt) = 0;

  virtual void init_impl(const InputOptions &opts, const array::Scalar &ice_thickness,
                         const array::Scalar &sea_level_elevation) = 0;

  virtual void bootstrap_impl(const array::Scalar &bed_elevation,
                              const array::Scalar &bed_uplift,
                              const array::Scalar &ice_thickness,
                              const array::Scalar &sea_level_elevation) = 0;

  static void apply_topg_offset(const std::string &filename, array::Scalar &bed_topography);

  //! current bed elevation
  array::Scalar2 m_topg;

  //! bed elevation at the time of the last update
  array::Scalar m_topg_last;

  array::Scalar m_load;
  array::Scalar m_load_accumulator;

  //! bed uplift rate
  array::Scalar m_uplift;

  //! time of the last bed deformation update
  double m_t_last;
  //! Update interval in seconds
  double m_update_interval;
  //! Temporal resolution to use when checking whether it's time to update
  double m_t_eps;
  //! Name of the variable used to store the last update time.
  std::string m_time_name;

  std::string m_model_name;
};

/*!
 * The do-nothing bed deformation model.
 */
class Null : public BedDef {
public:
  Null(std::shared_ptr<const Grid> g);
protected:
  void update_impl(const array::Scalar &load, double t, double dt);

  void init_impl(const InputOptions &opts, const array::Scalar &ice_thickness,
                 const array::Scalar &sea_level_elevation);

  void bootstrap_impl(const array::Scalar &bed_elevation,
                      const array::Scalar &bed_uplift,
                      const array::Scalar &ice_thickness,
                      const array::Scalar &sea_level_elevation);

  MaxTimestep max_timestep_impl(double t) const;
};

//! Point-wise isostasy bed deformation model.
class PointwiseIsostasy : public BedDef {
public:
  PointwiseIsostasy(std::shared_ptr<const Grid> g);
  virtual ~PointwiseIsostasy() = default;
protected:
  void init_impl(const InputOptions &opts, const array::Scalar &ice_thickness,
                 const array::Scalar &sea_level_elevation);

  void bootstrap_impl(const array::Scalar &bed_elevation,
                      const array::Scalar &bed_uplift,
                      const array::Scalar &ice_thickness,
                      const array::Scalar &sea_level_elevation);

  void update_impl(const array::Scalar &load, double t, double dt);

  //! last ice load (ice-equivalent thickness)
  array::Scalar m_load_last;
};

} // end of namespace bed
} // end of namespace pism

#endif  // __BedDef_hh
