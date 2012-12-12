// Copyright (C) 2008-2012 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
// Gudfinna Adalgeirsdottir and Andy Aschwanden
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

#include "PISMSurface.hh"
#include "PISMAtmosphere.hh"
#include "PIO.hh"
#include "PISMVars.hh"
#include "PISMTime.hh"
#include "IceGrid.hh"
#include "pism_options.hh"

///// PISMSurfaceModel base class:

PISMSurfaceModel::PISMSurfaceModel(IceGrid &g, const NCConfigVariable &conf)
  : PISMComponent_TS(g, conf) {
  atmosphere = NULL;
}

PISMSurfaceModel::~PISMSurfaceModel() {
  delete atmosphere;
}

void PISMSurfaceModel::get_diagnostics(map<string, PISMDiagnostic*> &dict) {
  if (atmosphere)
    atmosphere->get_diagnostics(dict);
}

void PISMSurfaceModel::attach_atmosphere_model(PISMAtmosphereModel *input) {
  if (atmosphere != NULL) {
    delete atmosphere;
  }
  atmosphere = input;
}

PetscErrorCode PISMSurfaceModel::init(PISMVars &vars) {
  PetscErrorCode ierr;

  if (atmosphere == NULL)
    SETERRQ(grid.com, 1, "PISMSurfaceModel::init(PISMVars &vars): atmosphere == NULL");

  ierr = atmosphere->init(vars); CHKERRQ(ierr);

  return 0;
}

//! \brief Returns mass held in the surface layer.
/*!
 * Basic surface models currently implemented in PISM do not model the mass of
 * the surface layer.
 */
PetscErrorCode PISMSurfaceModel::mass_held_in_surface_layer(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = result.set(0.0); CHKERRQ(ierr);

  return 0;
}

//! \brief Returns thickness of the surface layer. Used to compute surface
//! elevation as a sum of elevation of the top surface of the ice and surface
//! layer (firn, etc) thickness.
/*!
 * Basic surface models currently implemented in PISM do not model surface
 * layer thickness.
 */
PetscErrorCode PISMSurfaceModel::surface_layer_thickness(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = result.set(0.0); CHKERRQ(ierr);

  return 0;
}

//! \brief Returns the liquid water fraction of the ice at the top ice surface.
/*!
 * Most PISM surface models return 0.
 */
PetscErrorCode PISMSurfaceModel::ice_surface_liquid_water_fraction(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = result.set(0.0); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PISMSurfaceModel::define_variables(set<string> vars, const PIO &nc, PISM_IO_Type nctype) {
  PetscErrorCode ierr;

  if (atmosphere != NULL) {
    ierr = atmosphere->define_variables(vars, nc, nctype); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PISMSurfaceModel::write_variables(set<string> vars, const PIO &nc) {
  PetscErrorCode ierr;

  if (atmosphere != NULL) {
    ierr = atmosphere->write_variables(vars, nc); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PISMSurfaceModel::max_timestep(PetscReal my_t, PetscReal &my_dt, bool &restrict) {
  PetscErrorCode ierr;

  if (atmosphere != NULL) {
    ierr = atmosphere->max_timestep(my_t, my_dt, restrict); CHKERRQ(ierr);
  } else {
    my_dt = -1;
    restrict = false;
  }

  return 0;
}

void PISMSurfaceModel::add_vars_to_output(string keyword, map<string,NCSpatialVariable> &result) {
  if (atmosphere != NULL) {
    atmosphere->add_vars_to_output(keyword, result);
  }
}
