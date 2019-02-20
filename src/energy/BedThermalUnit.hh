// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019 Ed Bueler and Constantine Khroulev
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

#ifndef _PISMBEDTHERMALUNIT_H_
#define _PISMBEDTHERMALUNIT_H_

#include "pism/util/Component.hh"
#include "pism/util/iceModelVec3Custom.hh"

#include "pism/util/Diagnostic.hh"

namespace pism {

class Vars;

//! @brief Energy balance models and utilities.
namespace energy {

// Vertical grid information for BTU_Full.
struct BTUGrid {
  BTUGrid(Context::ConstPtr ctx);
  static BTUGrid FromOptions(Context::ConstPtr ctx);

  unsigned int Mbz;             // number of vertical levels
  double Lbz;                   // depth of the bed thermal layer
};

//! Given the temperature of the top of the bedrock, for the duration of one time-step, provides upward geothermal flux at that interface at the end of the time-step.
/*!
  The geothermal flux actually applied to the base of an ice sheet is dependent, over time,
  on the temperature of the basal ice itself.  The purpose of a bedrock thermal layer
  in an ice sheet model is to implement this dependency by using a physical model
  for the temperature within that layer, the upper lithosphere.  Because the
  upper part of the lithosphere stores or releases energy into the ice,
  the typical lithosphere geothermal flux rate is not the same thing as the
  geothermal flux applied to the base of the ice.  This issue has long been
  recognized by ice sheet modelers [%e.g. \ref RitzFabreLetreguilly].

  For instance, suppose the ice sheet is in a balanced state in which the geothermal
  flux deep in the crust is equal to the heat flux into the ice base.  If the
  near-surface ice cools from this state then, because the ice temperature gradient
  is now greater in magnitude, between the warm bedrock and the cooler ice, the ice
  will for some period receive more than the deep geothermal flux rate. Similarly,
  if the ice warms from the balanced state then the temperature difference with
  the bedrock has become smaller and the magnitude of the ice basal heat flux will
  be less than the deep geothermal rate.

  We regard the lithosphere geothermal flux rate, which is applied in this model
  to the base of the bedrock thermal layer, as a time-independent quantity.  This
  concept is the same as in all published ice sheet models, to our knowledge.

  Because the relevant layer of bedrock below an ice sheet is typically shallow,
  modeling the bedrock temperature is quite simple.
  Let \f$T_b(t,x,y,z)\f$ be the temperature of the bedrock layer, for elevations
  \f$-L_b \le z \le 0\f$.  In this routine, \f$z=0\f$ refers to the top of the
  bedrock, the ice/bedrock interface.  (Note \f$z=0\f$ is the base of the ice in
  IceModel, and thus a different location if ice is floating.)
  Let \f$G\f$ be the lithosphere geothermal flux rate, namely the PISM input
  variable `bheatflx`; see Related Page \ref std_names .  Let \f$k_b\f$
  = `bedrock_thermal_conductivity` in pism_config.cdl) be the constant thermal
  conductivity of the upper lithosphere.  In these terms the actual
  upward heat flux into the ice/bedrock interface is the quantity,
  \f[G_0 = -k_b \frac{\partial T_b}{\partial z}.\f]
  This is the \e output of the method top_heat_flux() in this class.

  The evolution equation solved in this class, for which a timestep is done by the
  update() method, is the standard 1D heat equation
  \f[\rho_b c_b \frac{\partial T_b}{\partial t} = k_b \frac{\partial^2 T_b}{\partial z^2}\f]
  where \f$\rho_b\f$ = `bedrock_thermal_density` and \f$c_b\f$ =
  `bedrock_thermal_specific_heat_capacity` in pism_config.cdl.

  If `n_levels` >= 3 then everything is the general case.  The lithospheric temperature
  in `temp` is saved in files as `litho_temp`.  The top_heat_flux()
  method uses second-order differencing to compute the values of \f$G_0\f$.

  If `n_levels` <= 1 then this object becomes very simplified: there is no internal
  state in IceModelVec3 temp.  The update() and allocate() methods are null,
  and the top_heat_flux() method does nothing other than to copy the
  field \f$G\f$ = `bheatflx` into `result`.

  If `n_levels` == 2 then everything is the general case except that 
  top_heat_flux() method uses first-order differencing to compute the
  values of \f$G_0\f$.
*/
class BedThermalUnit : public Component {
public:

  static BedThermalUnit* FromOptions(IceGrid::ConstPtr g,
                                     Context::ConstPtr ctx);

  BedThermalUnit(IceGrid::ConstPtr g);

  virtual ~BedThermalUnit();

  typedef std::shared_ptr<BedThermalUnit> Ptr;
  typedef std::shared_ptr<const BedThermalUnit> ConstPtr;

  void init(const InputOptions &opts);

  //! Return the upward heat flux through the top surface of the bedrock thermal layer.
  const IceModelVec2S& flux_through_top_surface() const;

  //! Return the upward heat flux through the bottom surface of the bedrock thermal layer.
  const IceModelVec2S& flux_through_bottom_surface() const;

  void update(const IceModelVec2S &bedrock_top_temperature,
              double t, double dt);

  double vertical_spacing() const;
  double depth() const;

  unsigned int Mz() const;

protected:
  virtual void initialize_bottom_surface_flux();

  virtual void init_impl(const InputOptions &opts);

  virtual void update_impl(const IceModelVec2S &bedrock_top_temperature,
                           double t, double dt) = 0;

  virtual double vertical_spacing_impl() const = 0;
  virtual double depth_impl() const = 0;
  virtual unsigned int Mz_impl() const = 0;

  virtual void define_model_state_impl(const PIO &output) const;
  virtual void write_model_state_impl(const PIO &output) const;

  virtual DiagnosticList diagnostics_impl() const;
protected:
  //! upward heat flux through the bottom surface of the bed thermal layer
  IceModelVec2S m_bottom_surface_flux;

  //! upward heat flux through the top surface of the bed thermal layer
  IceModelVec2S m_top_surface_flux;
};

class BTU_geothermal_flux_at_ground_level : public Diag<BedThermalUnit> {
public:
  BTU_geothermal_flux_at_ground_level(const BedThermalUnit *m);
protected:
  virtual IceModelVec::Ptr compute_impl() const;
};

} // end of namespace energy
} // end of namespace pism

#endif /* _PISMBEDTHERMALUNIT_H_ */

