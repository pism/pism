// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022 PISM Authors
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

void compute_load(const array::Scalar &bed_elevation,
                  const array::Scalar &ice_thickness,
                  const array::Scalar &sea_level_elevation,
                  array::Scalar &result);

//! PISM bed deformation model (base class).
class BedDef : public Component {
public:
  BedDef(IceGrid::ConstPtr g);
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
  virtual void define_model_state_impl(const File &output) const;
  virtual void write_model_state_impl(const File &output) const;

  virtual DiagnosticList diagnostics_impl() const;

  virtual void update_impl(const array::Scalar &ice_thickness,
                           const array::Scalar &sea_level_elevation,
                           double t, double dt) = 0;
  virtual void init_impl(const InputOptions &opts, const array::Scalar &ice_thickness,
                         const array::Scalar &sea_level_elevation);
  virtual void bootstrap_impl(const array::Scalar &bed_elevation,
                              const array::Scalar &bed_uplift,
                              const array::Scalar &ice_thickness,
                              const array::Scalar &sea_level_elevation);
  virtual void apply_topg_offset(const std::string &filename);

  void compute_uplift(const array::Scalar &bed, const array::Scalar &bed_last,
                            double dt, array::Scalar &result);
protected:
  const int m_wide_stencil;
  //! current bed elevation
  array::Scalar2 m_topg;

  //! bed elevation at the time of the last update
  array::Scalar2 m_topg_last;

  //! bed uplift rate
  array::Scalar m_uplift;
};

/*!
 * The do-nothing bed deformation model.
 */
class Null : public BedDef {
public:
  Null(IceGrid::ConstPtr g);
protected:
  void update_impl(const array::Scalar &ice_thickness,
                   const array::Scalar &sea_level_elevation,
                   double t, double dt);
  MaxTimestep max_timestep_impl(double t) const;
  void init_impl(const InputOptions &opts, const array::Scalar &ice_thickness,
                 const array::Scalar &sea_level_elevation);
};

//! Point-wise isostasy bed deformation model.
class PointwiseIsostasy : public BedDef {
public:
  PointwiseIsostasy(IceGrid::ConstPtr g);
  virtual ~PointwiseIsostasy() = default;
protected:
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
  //! last ice load (ice-equivalent thickness)
  array::Scalar m_load_last;
};

} // end of namespace bed
} // end of namespace pism

#endif  // __BedDef_hh
