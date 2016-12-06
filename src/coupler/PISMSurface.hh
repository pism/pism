// Copyright (C) 2008-2016 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
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

#include "base/util/PISMComponent.hh"

namespace pism {

namespace atmosphere {
class AtmosphereModel;
}

class IceModelVec2S;

//! @brief Surface models and modifiers: provide top-surface
//! temperature, mass flux, liquid water fraction, mass and thickness of the surface
//! layer.
namespace surface {

//! \brief The interface of PISM's surface models.
class SurfaceModel : public Component_TS {
public:
  SurfaceModel(IceGrid::ConstPtr g);
  virtual ~SurfaceModel();

  void init();

  void attach_atmosphere_model(atmosphere::AtmosphereModel *input);

  // the interface:
  void ice_surface_mass_flux(IceModelVec2S &result) const;

  void ice_surface_temperature(IceModelVec2S &result) const;
  void ice_surface_liquid_water_fraction(IceModelVec2S &result) const;

  void mass_held_in_surface_layer(IceModelVec2S &result) const;
  void surface_layer_thickness(IceModelVec2S &result) const;
protected:
  virtual void init_impl();

  virtual void attach_atmosphere_model_impl(atmosphere::AtmosphereModel *input);

  virtual void define_model_state_impl(const PIO &output) const;
  virtual void write_model_state_impl(const PIO &output) const;

  virtual MaxTimestep max_timestep_impl(double my_t) const;

  virtual void surface_layer_thickness_impl(IceModelVec2S &result) const;
  virtual void mass_held_in_surface_layer_impl(IceModelVec2S &result) const;

  virtual void ice_surface_temperature_impl(IceModelVec2S &result) const = 0;
  virtual void ice_surface_liquid_water_fraction_impl(IceModelVec2S &result) const;

  virtual void ice_surface_mass_flux_impl(IceModelVec2S &result) const = 0;

  virtual std::map<std::string, Diagnostic::Ptr> diagnostics_impl() const;
  virtual std::map<std::string, TSDiagnostic::Ptr> ts_diagnostics_impl() const;
protected:
  atmosphere::AtmosphereModel *m_atmosphere;
};

/*! @brief Climatic mass balance */
class PS_climatic_mass_balance : public Diag<SurfaceModel>
{
public:
  PS_climatic_mass_balance(const SurfaceModel *m);
protected:
  IceModelVec::Ptr compute_impl();
};

/*! @brief Ice surface temperature. */
class PS_ice_surface_temp : public Diag<SurfaceModel>
{
public:
  PS_ice_surface_temp(const SurfaceModel *m);
protected:
  IceModelVec::Ptr compute_impl();
};

/*! @brief Ice liquid water fraction at the ice surface. */
class PS_liquid_water_fraction : public Diag<SurfaceModel>
{
public:
  PS_liquid_water_fraction(const SurfaceModel *m);
protected:
  IceModelVec::Ptr compute_impl();
};

/*! @brief Mass of the surface layer (snow and firn). */
class PS_surface_layer_mass : public Diag<SurfaceModel>
{
public:
  PS_surface_layer_mass(const SurfaceModel *m);
protected:
  IceModelVec::Ptr compute_impl();
};

/*! @brief Surface layer (snow and firn) thickness. */
class PS_surface_layer_thickness : public Diag<SurfaceModel>
{
public:
  PS_surface_layer_thickness(const SurfaceModel *m);
protected:
  IceModelVec::Ptr compute_impl();
};

} // end of namespace surface
} // end of namespace pism

#endif  // __PISMSurfaceModel_hh

