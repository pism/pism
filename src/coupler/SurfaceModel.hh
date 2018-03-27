// Copyright (C) 2008-2018 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
// Gudfinna Adalgeirsdottir and Andy Aschwanden
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

#ifndef __PISMSurfaceModel_hh
#define __PISMSurfaceModel_hh

/*!
 * This file should contain the class definition and nothing else.
 * Implementations should go in separate files.
 */

#include "pism/util/Component.hh"

namespace pism {

namespace atmosphere {
class AtmosphereModel;
}

class Geometry;
class IceModelVec2S;

//! @brief Surface models and modifiers: provide top-surface
//! temperature, mass flux, liquid water fraction, mass and thickness of the surface
//! layer.
namespace surface {

//! \brief The interface of PISM's surface models.
class SurfaceModel : public Component {
public:
  SurfaceModel(IceGrid::ConstPtr g);
  SurfaceModel(IceGrid::ConstPtr g, std::shared_ptr<SurfaceModel> input);
  SurfaceModel(IceGrid::ConstPtr g, std::shared_ptr<atmosphere::AtmosphereModel> atmosphere);

  virtual ~SurfaceModel();

  void init(const Geometry &geometry);

  // the interface:
  void update(const Geometry &geometry, double t, double dt);

  const IceModelVec2S& layer_mass() const;
  const IceModelVec2S& layer_thickness() const;
  const IceModelVec2S& liquid_water_fraction() const;
  const IceModelVec2S& mass_flux() const;
  const IceModelVec2S& temperature() const;

protected:

  virtual const IceModelVec2S& layer_mass_impl() const;
  virtual const IceModelVec2S& layer_thickness_impl() const;
  virtual const IceModelVec2S& liquid_water_fraction_impl() const;
  virtual const IceModelVec2S& mass_flux_impl() const;
  virtual const IceModelVec2S& temperature_impl() const;

  virtual void init_impl(const Geometry &geometry);
  virtual void update_impl(const Geometry &geometry, double t, double dt);

  virtual void define_model_state_impl(const PIO &output) const;
  virtual void write_model_state_impl(const PIO &output) const;

  virtual MaxTimestep max_timestep_impl(double my_t) const;

  virtual DiagnosticList diagnostics_impl() const;
  virtual TSDiagnosticList ts_diagnostics_impl() const;

  static IceModelVec2S::Ptr allocate_layer_mass(IceGrid::ConstPtr grid);
  static IceModelVec2S::Ptr allocate_layer_thickness(IceGrid::ConstPtr grid);
  static IceModelVec2S::Ptr allocate_liquid_water_fraction(IceGrid::ConstPtr grid);
  static IceModelVec2S::Ptr allocate_mass_flux(IceGrid::ConstPtr grid);
  static IceModelVec2S::Ptr allocate_temperature(IceGrid::ConstPtr grid);

protected:
  IceModelVec2S::Ptr m_liquid_water_fraction;
  IceModelVec2S::Ptr m_layer_mass;
  IceModelVec2S::Ptr m_layer_thickness;

  std::shared_ptr<SurfaceModel> m_input_model;
  std::shared_ptr<atmosphere::AtmosphereModel> m_atmosphere;
};

} // end of namespace surface
} // end of namespace pism

#endif  // __PISMSurfaceModel_hh
