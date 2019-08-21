// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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
#include "pism/util/iceModelVec.hh"

namespace pism {

//! @brief Bed-related models: bed deformation (provide bed elevation
//! and uplift) and (soon) bed erosion.
namespace bed {

double compute_load(double bed, double ice_thickness, double sea_level,
                    double ice_density, double ocean_density);

void compute_load(const IceModelVec2S &bed_elevation,
                  const IceModelVec2S &ice_thickness,
                  const IceModelVec2S &sea_level_elevation,
                  IceModelVec2S &result);

//! PISM bed deformation model (base class).
class BedDef : public Component {
public:
  BedDef(IceGrid::ConstPtr g);
  virtual ~BedDef();

  void init(const InputOptions &opts, const IceModelVec2S &ice_thickness,
            const IceModelVec2S &sea_level_elevation);
  void bootstrap(const IceModelVec2S &bed_elevation,
                 const IceModelVec2S &bed_uplift,
                 const IceModelVec2S &ice_thickness,
                 const IceModelVec2S &sea_level_elevation);

  void update(const IceModelVec2S &ice_thickness,
              const IceModelVec2S &sea_level_elevation,
              double t, double dt);

  const IceModelVec2S& bed_elevation() const;
  const IceModelVec2S& uplift() const;

protected:
  virtual void define_model_state_impl(const PIO &output) const;
  virtual void write_model_state_impl(const PIO &output) const;

  virtual DiagnosticList diagnostics_impl() const;

  virtual void update_impl(const IceModelVec2S &ice_thickness,
                           const IceModelVec2S &sea_level_elevation,
                           double t, double dt) = 0;
  virtual void init_impl(const InputOptions &opts, const IceModelVec2S &ice_thickness,
                         const IceModelVec2S &sea_level_elevation);
  virtual void bootstrap_impl(const IceModelVec2S &bed_elevation,
                              const IceModelVec2S &bed_uplift,
                              const IceModelVec2S &ice_thickness,
                              const IceModelVec2S &sea_level_elevation);
  virtual void apply_topg_offset(const std::string &filename);

  void compute_uplift(const IceModelVec2S &bed, const IceModelVec2S &bed_last,
                            double dt, IceModelVec2S &result);
protected:
  //! time of the last bed deformation update
  double m_t_beddef_last;

  //! current bed elevation
  IceModelVec2S m_topg;

  //! bed elevation at the time of the last update
  IceModelVec2S m_topg_last;

  //! bed uplift rate
  IceModelVec2S m_uplift;
};

/*!
 * The do-nothing bed deformation model.
 */
class Null : public BedDef {
public:
  Null(IceGrid::ConstPtr g);
protected:
  void update_impl(const IceModelVec2S &ice_thickness,
                   const IceModelVec2S &sea_level_elevation,
                   double t, double dt);
  MaxTimestep max_timestep_impl(double t) const;
  void init_impl(const InputOptions &opts, const IceModelVec2S &ice_thickness,
                 const IceModelVec2S &sea_level_elevation);
};

//! Point-wise isostasy bed deformation model.
class PointwiseIsostasy : public BedDef {
public:
  PointwiseIsostasy(IceGrid::ConstPtr g);
  virtual ~PointwiseIsostasy();
protected:
  MaxTimestep max_timestep_impl(double t) const;
  void init_impl(const InputOptions &opts, const IceModelVec2S &ice_thickness,
                 const IceModelVec2S &sea_level_elevation);
  void bootstrap_impl(const IceModelVec2S &bed_elevation,
                      const IceModelVec2S &bed_uplift,
                      const IceModelVec2S &ice_thickness,
                      const IceModelVec2S &sea_level_elevation);
  void update_impl(const IceModelVec2S &ice_thickness,
                   const IceModelVec2S &sea_level_elevation,
                   double t, double dt);
  IceModelVec2S m_load_last;       //!< last ice load (ice-equivalent thickness)
};

} // end of namespace bed
} // end of namespace pism

#endif  // __BedDef_hh
