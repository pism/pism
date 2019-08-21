/* Copyright (C) 2016, 2017, 2019 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef BTU_FULL_H
#define BTU_FULL_H

#include "BedThermalUnit.hh"
#include "pism/util/Context.hh"

namespace pism {
namespace energy {

class BedrockColumn;

//! @brief Given the temperature of the top of the bedrock, for the duration of one time-step,
//! provides upward geothermal flux at that interface at the end of the time-step.
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
  This is the \e output of the method flux_through_top_surface() in this class.

  The evolution equation solved in this class, for which a timestep is done by the
  update() method, is the standard 1D heat equation
  \f[\rho_b c_b \frac{\partial T_b}{\partial t} = k_b \frac{\partial^2 T_b}{\partial z^2}\f]
  where \f$\rho_b\f$ = `bedrock_thermal_density` and \f$c_b\f$ =
  `bedrock_thermal_specific_heat_capacity` in pism_config.cdl.

  If 3 or more levels are used then everything is the general case. The lithospheric temperature in
  `temp` is saved in files as `litho_temp`. The flux_through_top_surface() method uses second-order
  differencing to compute the values of \f$G_0\f$.

  If 2 levels are used then everything is the general case except that flux_through_top_surface()
  method uses first-order differencing to compute the values of \f$G_0\f$.
*/
class BTU_Full : public BedThermalUnit {
public:
  BTU_Full(IceGrid::ConstPtr g, const BTUGrid &vertical_grid);
  virtual ~BTU_Full();

  //! Bedrock thermal layer temperature field.
  const IceModelVec3Custom& temperature() const;

protected:
  virtual void bootstrap(const IceModelVec2S &bedrock_top_temperature);

  virtual void init_impl(const InputOptions &opts);

  virtual double vertical_spacing_impl() const;
  virtual double depth_impl() const;
  virtual unsigned int Mz_impl() const;

  virtual MaxTimestep max_timestep_impl(double my_t) const;

  using BedThermalUnit::update_impl;
  virtual void update_impl(const IceModelVec2S &bedrock_top_temperature,
                           double t, double dt);

  virtual void define_model_state_impl(const PIO &output) const;
  virtual void write_model_state_impl(const PIO &output) const;
protected:
  //! bedrock thermal layer temperature, in degrees Kelvin; part of state; uses equally-spaced
  //! layers.
  IceModelVec3Custom m_temp;

  //! bedrock thermal conductivity
  double m_k;
  //! diffusivity of the heat flow within the bedrock layer
  double m_D;

  //! number of vertical levels within the bedrock
  unsigned int m_Mbz;
  //! thickness of the bedrock layer, in meters
  double m_Lbz;

  //! true if the model needs to "bootstrap" the temperature field during the first time step
  bool m_bootstrapping_needed;

  void update_flux_through_top_surface();

  std::shared_ptr<BedrockColumn> m_column;
};

} // end of namespace energy
} // end of namespace pism


#endif /* BTU_FULL_H */
