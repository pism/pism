// Copyright (C) 2012-2018 PISM Authors
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

#ifndef _PISMHYDROLOGY_H_
#define _PISMHYDROLOGY_H_

#include "pism/util/iceModelVec.hh"
#include "pism/util/Component.hh"
#include "pism/util/iceModelVec2T.hh"

namespace pism {

class IceModelVec2CellType;

//! @brief Sub-glacial hydrology models and related diagnostics.
namespace hydrology {

class Inputs {
public:
  Inputs();

  // modeling domain (set to NULL in whole-ice-sheet configurations)
  const IceModelVec2Int      *no_model_mask;
  // geometry
  const IceModelVec2CellType *cell_type;
  const IceModelVec2S        *ice_thickness;
  const IceModelVec2S        *bed_elevation;
  // hydrological inputs
  const IceModelVec2S        *surface_input_rate;
  const IceModelVec2S        *basal_melt_rate;
  const IceModelVec2S        *ice_sliding_speed;
};

//! \brief The PISM subglacial hydrology model interface.
/*!
  This is a virtual base class.

  The purpose of this class and its derived classes is to provide
  \code
  subglacial_water_thickness()
  subglacial_water_pressure()
  till_water_thickness()
  \endcode
  These correspond to state variables \f$W\f$, \f$P\f$, and \f$W_{\text{till}}\f$
  in [\ref BuelervanPeltDRAFT], though not all derived classes of Hydrology
  have all of them as state variables.

  Additional modeled fields, for diagnostic purposes, are
  \code
  overburden_pressure(IceModelVec2S &result)
  wall_melt(IceModelVec2S &result)
  \endcode

  This interface is appropriate to subglacial hydrology models which track a
  two-dimensional water layer with a well-defined thickness and pressure at each
  map-plane location.  The methods subglacial_water_thickness() and
  subglacial_water_pressure() return amount and pressure of the *transportable*
  water, that is, the subglacial water which moves along a modeled hydraulic
  head gradient, in contrast to the water stored in the till.

  The transportable water moves through a subglacial morphology which is not
  determined in this base class.

  The Hydrology models have separate, but potentially-coupled, water
  which is held in local till storage.  Thus the
  transportable water (bwat) and till water (tillwat) thicknesses are different.
  Published models with till storage include [\ref BBssasliding, \ref SchoofTill,
  \ref TrufferEchelmeyerHarrison2001, \ref Tulaczyketal2000b, \ref vanderWeletal2013].

  The till water thickness is can be used, via the theory of
  [\ref Tulaczyketal2000], to compute an effective pressure for the water in the
  pore spaces of the till, which can then be used by the Mohr-Coulomb criterion
  to provide a yield stress.  Class MohrCoulombYieldStress does this
  calculation.  Here in Hydrology only the till water thickness tillwat is
  computed.

  Hydrology is a timestepping component.  Because of the
  short physical timescales associated to liquid water moving under a glacier,
  Hydrology (and derived) classes generally take many substeps in PISM's major
  ice dynamics time steps.  Thus when an update() method in a Hydrology
  class is called it will advance its internal time to the new goal t+dt
  using its own internal time steps.

  Generally Hydrology classes use the ice geometry, the basal melt
  rate, and the basal sliding velocity in determining the evolution of the
  hydrology state variables.  Note that the basal melt rate is an
  energy-conservation-derived field and the basal-sliding velocity is derived
  from the solution of a stress balance.  The basal melt rate and
  sliding velocity fields therefore generally come from IceModel and
  StressBalance, respectively.

  Additional, time-dependent and spatially-variable water input to the basal
  layer, taken directly from a file, is possible too.

  Ice geometry and energy fields are normally treated as constant in time
  during the update() call for the interval [t,t+dt].  Thus the coupling is
  one-way during the update() call.
*/
class Hydrology : public Component {
public:
  Hydrology(IceGrid::ConstPtr g);
  virtual ~Hydrology();

  void restart(const PIO &input_file, int record);

  void bootstrap(const PIO &input_file,
                 const IceModelVec2S &ice_thickness);

  void initialize(const IceModelVec2S &W_till,
                  const IceModelVec2S &W,
                  const IceModelVec2S &P);

  void update(double t, double dt, const Inputs& inputs);

  const IceModelVec2S& till_water_thickness() const;
  const IceModelVec2S& subglacial_water_thickness() const;
  const IceModelVec2S& overburden_pressure() const;
  const IceModelVec2S& total_input_rate() const;

  const IceModelVec2S& mass_change() const;
  const IceModelVec2S& mass_change_at_grounded_margin() const;
  const IceModelVec2S& mass_change_at_grounding_line() const;
  const IceModelVec2S& mass_change_at_domain_boundary() const;
  const IceModelVec2S& mass_change_due_to_conservation_error() const;
  const IceModelVec2S& mass_change_due_to_input() const;
  const IceModelVec2S& mass_change_due_to_lateral_flow() const;

protected:
  virtual void restart_impl(const PIO &input_file, int record);

  virtual void bootstrap_impl(const PIO &input_file,
                              const IceModelVec2S &ice_thickness);

  virtual void initialize_impl(const IceModelVec2S &W_till,
                               const IceModelVec2S &W,
                               const IceModelVec2S &P);

  virtual void update_impl(double t, double dt, const Inputs& inputs) = 0;
  virtual std::map<std::string, Diagnostic::Ptr> diagnostics_impl() const;

  virtual void define_model_state_impl(const PIO &output) const;
  virtual void write_model_state_impl(const PIO &output) const;

  void compute_overburden_pressure(const IceModelVec2S &ice_thickness,
                                   IceModelVec2S &result) const;

  void compute_input_rate(const IceModelVec2CellType &mask,
                          const IceModelVec2S &basal_melt_rate,
                          const IceModelVec2S *surface_input_rate,
                          IceModelVec2S &result);
protected:
  //! effective thickness of basal water stored in till
  IceModelVec2S m_Wtill;

  //! effective thickness of transportable basal water
  IceModelVec2S m_W;

  //! overburden pressure
  IceModelVec2S m_Pover;

  // total input rate, combining basal melt rate and the input from the surface
  IceModelVec2S m_input_rate;

  // change due to flow for the current hydrology time step
  IceModelVec2S m_flow_change_incremental;

  // changes in water thickness
  //
  // these quantities are re-set to zero at the beginning of the PISM time step
  IceModelVec2S m_conservation_error_change;
  IceModelVec2S m_grounded_margin_change;
  IceModelVec2S m_grounding_line_change;
  IceModelVec2S m_input_change;
  IceModelVec2S m_no_model_mask_change;
  IceModelVec2S m_total_change;
  IceModelVec2S m_flow_change;

  // when we update the water amounts, careful mass accounting at the boundary
  // is needed
  void enforce_bounds(const IceModelVec2CellType &cell_type,
                      const IceModelVec2Int *no_model_mask,
                      double max_thickness,
                      IceModelVec2S &water_thickness,
                      IceModelVec2S &grounded_margin_change,
                      IceModelVec2S &grounding_line_change,
                      IceModelVec2S &conservation_error_change,
                      IceModelVec2S &no_model_mask_change);
private:
  virtual void initialization_message() const = 0;
};

void check_bounds(const IceModelVec2S& W, double W_max);

} // end of namespace hydrology
} // end of namespace pism

#endif /* _PISMHYDROLOGY_H_ */
