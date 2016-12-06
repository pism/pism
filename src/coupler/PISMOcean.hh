// Copyright (C) 2008-2011, 2013, 2014, 2015, 2016 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
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

#ifndef __PISMOceanModel_hh
#define __PISMOceanModel_hh

#include "base/util/PISMComponent.hh"

namespace pism {
class IceModelVec2S;

//! @brief Ocean models and modifiers: provide sea level elevation,
//! melange back pressure, shelf base mass flux and shelf base
//! temperature.
namespace ocean {
//! A very rudimentary PISM ocean model.
class OceanModel : public Component_TS {
public:
  OceanModel(IceGrid::ConstPtr g);
  virtual ~OceanModel();

  void init();

  double sea_level_elevation() const;
  void shelf_base_temperature(IceModelVec2S &result) const;
  void shelf_base_mass_flux(IceModelVec2S &result) const;
  void melange_back_pressure_fraction(IceModelVec2S &result) const;
protected:
  virtual void init_impl() = 0;

  virtual std::map<std::string, Diagnostic::Ptr> diagnostics_impl() const;

  virtual void melange_back_pressure_fraction_impl(IceModelVec2S &result) const ;
  virtual void shelf_base_mass_flux_impl(IceModelVec2S &result) const = 0;
  virtual void shelf_base_temperature_impl(IceModelVec2S &result) const = 0;
  virtual void sea_level_elevation_impl(double &result) const = 0;

protected:
  double m_sea_level;
};

/*! @brief Sea level elevation. */
class PO_sea_level : public Diag<OceanModel>
{
public:
  PO_sea_level(const OceanModel *m);
protected:
  IceModelVec::Ptr compute_impl();
};

/*! @brief Shelf base temperature. */
class PO_shelf_base_temperature : public Diag<OceanModel>
{
public:
  PO_shelf_base_temperature(const OceanModel *m);
protected:
  IceModelVec::Ptr compute_impl();
};

/*! @brief Shelf base mass flux. */
class PO_shelf_base_mass_flux : public Diag<OceanModel>
{
public:
  PO_shelf_base_mass_flux(const OceanModel *m);
protected:
  IceModelVec::Ptr compute_impl();
};

/*! @brief Melange back pressure fraction. */
class PO_melange_back_pressure_fraction : public Diag<OceanModel>
{
public:
  PO_melange_back_pressure_fraction(const OceanModel *m);
protected:
  IceModelVec::Ptr compute_impl();
};

} // end of namespace ocean
} // end of namespace pism

#endif  // __PISMOceanModel_hh
