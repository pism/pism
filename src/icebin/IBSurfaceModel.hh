// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 PISM Authors
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

#ifndef _IBSURFACEMODEL_H_
#define _IBSURFACEMODEL_H_

#include "base/util/iceModelVec.hh"
#include "coupler/PISMAtmosphere.hh"
#include "coupler/PISMSurface.hh"

namespace pism {
namespace icebin {

//! \brief A class implementing a constant-in-time surface model for the surface mass balance.
//!
//! Reads data from a PISM input file.
//!
//! Ice surface temperature is parameterized as in PISM-IBSurfaceModel, using a latitude
//! and surface elevation-dependent formula.

class IBSurfaceModel : public pism::surface::SurfaceModel {
public:
  IBSurfaceModel(IceGrid::ConstPtr g);

protected:
  virtual void init_impl();
  virtual void attach_atmosphere_model_impl(atmosphere::AtmosphereModel *input);
  virtual void ice_surface_mass_flux_impl(IceModelVec2S &result);
  virtual void ice_surface_temperature_impl(IceModelVec2S &result);
  virtual void ice_surface_liquid_water_fraction_impl(IceModelVec2S &result);
  virtual MaxTimestep max_timestep_impl(double t);
  virtual void update_impl(double my_t, double my_dt);
  virtual void get_diagnostics_impl(std::map<std::string, Diagnostic::Ptr> &dict,
                                    std::map<std::string, TSDiagnostic::Ptr> &ts_dict);
  virtual void write_variables_impl(const std::set<std::string> &vars, const PIO &nc);
  virtual void add_vars_to_output_impl(const std::string &keyword, std::set<std::string> &result);
  virtual void define_variables_impl(const std::set<std::string> &vars, const PIO &nc, IO_Type nctype);

protected:
  bool _initialized;
  std::string m_input_file;

  // Organized list of the variables...
  std::vector<std::pair<std::string, pism::IceModelVec2S *>> vecs;

  // Used internally
  void create(pism::IceModelVec2S &vec, std::string const &name);


public:
  // ------ See icebin/contracts/modele_pism.cpp
  // Inputs from IceBin


  // Mass of ice being transferred GCM --> Ice Model
  pism::IceModelVec2S massxfer; // [kg m-2 s-1]
  // Enthalpy of ice being transferred Stieglitz --> Icebin
  pism::IceModelVec2S enthxfer; // [J m-2 s-1]

  // GCM's idea of energy transfer into ice sheet.
  // Used to compute mass/energy budget
  pism::IceModelVec2S deltah;

  // Temperature of the Dirichlet B.C.
  pism::IceModelVec2S ice_top_bc_temp;
  // Water content of the Dirichlet B.C.
  pism::IceModelVec2S ice_top_bc_wc;
};

} // end of namespace surface
} // end of namespace pism

#endif /* _IBSURFACEMODEL_H_ */
