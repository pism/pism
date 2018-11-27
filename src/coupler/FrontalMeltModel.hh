// Copyright (C) 2018 Constantine Khroulev and Andy Aschwanden
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

#ifndef __PISMFrontalMeltModel_hh
#define __PISMFrontalMeltModel_hh

#include <memory>

#include "pism/util/Component.hh"
#include "pism/hydrology/Hydrology.hh"

namespace pism {

class IceModelVec2S;
class Geometry;
  
class FrontalMeltInputs {
public:
  FrontalMeltInputs();

  const Geometry *geometry;

  // used by the discharge routing
  const IceModelVec2S *subglacial_water_speed;

};

  //! @brief FrontalMelt models and modifiers: provides frontal
//! melt rate.
namespace frontalmelt {

//! A very rudimentary PISM frontal melt model.
class FrontalMeltModel : public Component {
public:
  // "modifier" constructor
  FrontalMeltModel(IceGrid::ConstPtr g, std::shared_ptr<FrontalMeltModel> input);
  // "model" constructor
  FrontalMeltModel(IceGrid::ConstPtr g);

  virtual ~FrontalMeltModel();

  void init(const Geometry &geometry);
  void bootstrap(const Geometry &geometry);

  void update(const FrontalMeltInputs &inputs, double t, double dt);

  const IceModelVec2S& frontal_melt_rate() const;

protected:
  virtual void init_impl(const Geometry &geometry);
  virtual void bootstrap_impl(const Geometry &geometry);
  // provides default (pass-through) implementations for "modifiers"
  virtual void update_impl(const FrontalMeltInputs &inputs, double t, double dt);
  virtual MaxTimestep max_timestep_impl(double t) const;
  virtual void define_model_state_impl(const PIO &output) const;
  virtual void write_model_state_impl(const PIO &output) const;

  virtual DiagnosticList diagnostics_impl() const;
  virtual TSDiagnosticList ts_diagnostics_impl() const;

  virtual const IceModelVec2S& frontal_melt_rate_impl() const;

protected:
  std::shared_ptr<FrontalMeltModel> m_input_model;

  static IceModelVec2S::Ptr allocate_frontal_melt_rate(IceGrid::ConstPtr g);
  
};

} // end of namespace frontalmelt
} // end of namespace pism

#endif  // __PISMFrontalMeltModel_hh
