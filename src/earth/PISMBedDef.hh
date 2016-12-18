// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016 PISM Authors
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

#include "base/util/PISMComponent.hh"
#include "base/util/iceModelVec.hh"

namespace pism {

//! @brief Bed-related models: bed deformation (provide bed elevation
//! and uplift) and (soon) bed erosion.
namespace bed {

//! PISM bed deformation model (base class).
class BedDef : public Component_TS {
public:
  BedDef(IceGrid::ConstPtr g);
  virtual ~BedDef();

  void init();
  void init(const IceModelVec2S &bed_elevation,
            const IceModelVec2S &bed_uplift,
            const IceModelVec2S &ice_thickness);

  using Component_TS::update;
  void update(const IceModelVec2S &ice_thickness,
              double my_t, double my_dt);

  const IceModelVec2S& bed_elevation() const;
  const IceModelVec2S& uplift() const;

  void set_elevation(const IceModelVec2S &input);
  void set_uplift(const IceModelVec2S &input);

protected:
  virtual void define_model_state_impl(const PIO &output) const;
  virtual void write_model_state_impl(const PIO &output) const;

  void update_impl(double my_t, double my_dt);
  virtual void update_with_thickness_impl(const IceModelVec2S &ice_thickness,
                                          double my_t, double my_dt) = 0;
  virtual void init_impl();
  virtual void init_with_inputs_impl(const IceModelVec2S &bed_elevation,
                                     const IceModelVec2S &bed_uplift,
                                     const IceModelVec2S &ice_thickness);
  virtual void apply_topg_offset(const std::string &filename);

  void compute_uplift(double dt_beddef);
protected:
  //! time of the last bed deformation update
  double m_t_beddef_last;

  //! current bed elevation
  IceModelVec2S m_topg;

  //! bed elevation at the beginning of a run
  IceModelVec2S m_topg_initial;

  //! bed elevation at the time of the last update
  IceModelVec2S m_topg_last;

  //! bed uplift rate
  IceModelVec2S m_uplift;
};

class PBNull : public BedDef {
public:
  PBNull(IceGrid::ConstPtr g);
protected:
  void update_with_thickness_impl(const IceModelVec2S &ice_thickness,
                                  double my_t, double my_dt);
  MaxTimestep max_timestep_impl(double t) const;
  void init_impl();
};

//! Pointwide isostasy bed deformation model.
class PBPointwiseIsostasy : public BedDef {
public:
  PBPointwiseIsostasy(IceGrid::ConstPtr g);
  virtual ~PBPointwiseIsostasy();
protected:
  virtual MaxTimestep max_timestep_impl(double t) const;
  virtual void init_impl();
  void update_with_thickness_impl(const IceModelVec2S &ice_thickness,
                                  double my_t, double my_dt);
  IceModelVec2S m_thk_last;       //!< last ice thickness
};

} // end of namespace bed
} // end of namespace pism

#endif  // __BedDef_hh
