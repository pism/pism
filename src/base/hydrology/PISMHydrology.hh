// Copyright (C) 2012-2016 PISM Authors
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

#include "base/util/iceModelVec.hh"
#include "base/util/PISMComponent.hh"

namespace pism {

class IceModelVec2T;

namespace stressbalance {
class StressBalance;
}

//! @brief Sub-glacial hydrology models and related diagnostics.
namespace hydrology {

//! \brief The PISM subglacial hydrology model interface.
/*!
  This is a virtual base class.

  The purpose of this class and its derived classes is to provide
  \code
  subglacial_water_thickness(IceModelVec2S &result)
  subglacial_water_pressure(IceModelVec2S &result)
  till_water_thickness(IceModelVec2S &result)
  \endcode
  These correspond to state variables \f$W\f$, \f$P\f$, and \f$W_{\text{til}}\f$
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

  Hydrology is a timestepping component (Component_TS).  Because of the
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
class Hydrology : public Component_TS {
public:
  Hydrology(IceGrid::ConstPtr g);
  virtual ~Hydrology();

  virtual void init();

  friend class Hydrology_hydrobmelt;
  friend class Hydrology_hydroinput;

  // all Hydrology models have a Wtil state variable, which this returns
  virtual void till_water_thickness(IceModelVec2S &result) const;

  // this diagnostic method returns the standard shallow approximation
  virtual void overburden_pressure(IceModelVec2S &result) const;

  // this diagnostic method returns zero in the base class
  virtual void wall_melt(IceModelVec2S &result) const;

  // these methods MUST be implemented in the derived class
  virtual void subglacial_water_thickness(IceModelVec2S &result) const = 0;
  virtual void subglacial_water_pressure(IceModelVec2S &result) const = 0;

protected:
  virtual void update_impl(double icet, double icedt) = 0;
  virtual std::map<std::string, Diagnostic::Ptr> diagnostics_impl() const;

  virtual void define_model_state_impl(const PIO &output) const;
  virtual void write_model_state_impl(const PIO &output) const;

  virtual void get_input_rate(double hydro_t, double hydro_dt, IceModelVec2S &result);
  virtual void check_Wtil_bounds();
protected:
  // this model's state
  IceModelVec2S m_Wtil;      // effective thickness of till
  // this model's workspace
  IceModelVec2S m_total_input, m_bmelt_local;

  bool m_hold_bmelt;

  IceModelVec2T *m_inputtobed;// time dependent input of water to bed, in addition to bmelt
  unsigned int m_inputtobed_period;      // in years
  double m_inputtobed_reference_time; // in seconds

};


//! The PISM minimal model has till, but water that exceeds the capacity of the till is not conserved. There is no model for lateral transport.
/*!
  This is the minimum functional derived class.  It updates till water thickness.
  It implements a version of the "undrained plastic bed" model of [\ref Tulaczyketal2000b],
  but with non-conserved drainage.

  It has no transportable water and subglacial_water_thickness() returns zero.

  This model can give no meaningful report on conservation errors, and thus it
  does not use the TSDiag objects used by mass-conserving derived classes.

  This talk illustrates a "till-can" metaphor applicable to this model:
  http://www2.gi.alaska.edu/snowice/glaciers/iceflow/bueler-igs-fairbanks-june2012.pdf
*/
class NullTransport : public Hydrology {
public:
  NullTransport(IceGrid::ConstPtr g);
  virtual ~NullTransport();

  virtual void init();

  //! Sets result to 0.
  virtual void subglacial_water_thickness(IceModelVec2S &result) const;

  //! Returns the overburden pressure in hope it is harmless.
  virtual void subglacial_water_pressure(IceModelVec2S &result) const;

protected:
  virtual MaxTimestep max_timestep_impl(double t) const;
  //! Solves an implicit step of a highly-simplified ODE.
  virtual void update_impl(double icet, double icedt);

  void diffuse_till_water(double dt);

private:
  double m_diffuse_tillwat;
  double m_diffusion_time;
  double m_diffusion_distance;
  double m_tillwat_max;
  double m_tillwat_decay_rate;

  IceModelVec2S m_Wtil_old;
};


//! \brief A subglacial hydrology model which assumes water pressure
//! equals overburden pressure.
/*!
  This PISM hydrology model has lateral motion of subglacial water and which
  conserves the water mass.  Further documentation is in [\ref BuelervanPeltDRAFT].

  The water velocity is along the steepest descent route for the hydraulic
  potential.  This potential is (mostly) a function of ice sheet geometry,
  because the water pressure is set to the overburden pressure, a simplified but
  well-established model [\ref Shreve1972].  However, the water layer thickness
  is also a part of the hydraulic potential because it is actually the potential
  of the *top* of the water layer.

  This (essential) model has been used for finding locations of subglacial lakes
  [\ref Siegertetal2009, \ref Livingstoneetal2013].  Subglacial lakes occur
  at local minima of the hydraulic potential.  If water builds up significantly
  (e.g. thickness of 10s of meters or more) then in the model here the resulting
  lakes diffuse instead of becoming infinitely deep.  Thus we avoid delta
  functions of water thickness at the minima of the hydraulic potential in this
  well-posed model.

  This model should generally be tested using static ice geometry first, i.e.
  using option -no_mass.

  The state space includes both the till water effective thickness \f$W_{til}\f$,
  which is in Hydrology, and the transportable water layer thickness \f$W\f$.

  For more complete modeling where the water pressure is determined by a
  physical model for the opening and closing of cavities, and where the state
  space includes a nontrivial pressure variable, see hydrology::Distributed.

  There is an option `-hydrology_null_strip` `X` which produces a strip of
  `X` km around the edge of the computational domain.  In that strip the water flow
  velocity is set to zero.  The water amount is also reset to zero at the end
  of each time step in this strip (in an accounted way).

  As noted this is the minimal model which has a lateral water flux.  This flux is
  \f[ \mathbf{q} = - K \nabla \psi = \mathbf{V} W - D \nabla W \f]
  where \f$\psi\f$ is the hydraulic potential
  \f[ \psi = P + \rho_w g (b + W). \f]
  The generalized conductivity \f$K\f$ is nontrivial and it generally also
  depends on the water thickness:
  \f[ K = k W^{\alpha-1} |\nabla (P+\rho_w g b)|^{\beta-2}. \f]

  This model contains enough information (enough modeled fields) so that we can
  compute the wall melt generated by dissipating the gravitational
  potential energy in the moving, presumably turbulent, subglacial water.  If we
  suppose that this heat is dissipated immediately as melt on the
  cavity/conduit walls then we get a formula for a wall melt contribution.  (This
  is in addition to the `bmelt` field coming from conserving energy in the flowing
  ice.)  See wall_melt().  At this time the wall melt is diagnostic only and does
  not add to the water amount W; such an addition is generally unstable.
*/
class Routing : public Hydrology {
public:
  Routing(IceGrid::ConstPtr g);
  virtual ~Routing();

  virtual void init();

  friend class MCHydrology_ice_free_land_loss_cumulative;
  friend class MCHydrology_ice_free_land_loss;
  friend class MCHydrology_ocean_loss_cumulative;
  friend class MCHydrology_ocean_loss;
  friend class MCHydrology_negative_thickness_gain_cumulative;
  friend class MCHydrology_negative_thickness_gain;
  friend class MCHydrology_null_strip_loss_cumulative;
  friend class MCHydrology_null_strip_loss;

  virtual void wall_melt(IceModelVec2S &result) const;

  virtual void subglacial_water_thickness(IceModelVec2S &result) const;

  virtual void subglacial_water_pressure(IceModelVec2S &result) const;

protected:
  virtual void update_impl(double icet, double icedt);

  virtual MaxTimestep max_timestep_impl(double t) const;
  virtual std::map<std::string, Diagnostic::Ptr> diagnostics_impl() const;
  virtual std::map<std::string, TSDiagnostic::Ptr> ts_diagnostics_impl() const;

  virtual void define_model_state_impl(const PIO &output) const;
  virtual void write_model_state_impl(const PIO &output) const;

protected:
  // this model's state
  IceModelVec2S m_W;      // water layer thickness
  // this model's auxiliary variables
  IceModelVec2Stag m_V,   // components are
  //   V(i,j,0) = u(i,j) = east-edge  centered x-component of water velocity
  //   V(i,j,1) = v(i,j) = north-edge centered y-component of water velocity
    m_Wstag,// edge-centered (staggered) W values (averaged from regular)
    m_K,// edge-centered (staggered) values of nonlinear conductivity
    m_Q;// edge-centered (staggered) advection fluxes
  // this model's workspace variables
  IceModelVec2S m_Wnew, m_Wtilnew, m_Pover;
  mutable IceModelVec2S m_R;

  double m_stripwidth; // width in m of strip around margin where V and W are set to zero;
  // if negative then the strip mechanism is inactive inactive

  virtual void init_bwat();

  // when we update the water amounts, careful mass accounting at the boundary
  // is needed; we update the new thickness variable, a temporary during update
  virtual void boundary_mass_changes(IceModelVec2S &newthk,
                                     double &icefreelost, double &oceanlost,
                                     double &negativegain, double &nullstriplost);

  double m_ice_free_land_loss_cumulative,
         m_ocean_loss_cumulative,
         m_negative_thickness_gain_cumulative,
         m_null_strip_loss_cumulative;

  virtual void check_water_thickness_nonnegative(IceModelVec2S &thk);

  virtual void water_thickness_staggered(IceModelVec2Stag &result);
  virtual void subglacial_hydraulic_potential(IceModelVec2S &result);

  virtual void conductivity_staggered(IceModelVec2Stag &result, double &maxKW);
  virtual void velocity_staggered(IceModelVec2Stag &result) const;
  friend class Routing_bwatvel;  // needed because bwatvel diagnostic needs protected velocity_staggered()
  virtual void advective_fluxes(IceModelVec2Stag &result);

  virtual void adaptive_for_W_evolution(double t_current, double t_end, double maxKW,
                                        double &dt_result,
                                        double &maxV_result, double &maxD_result,
                                        double &dtCFL_result, double &dtDIFFW_result);

  void raw_update_W(double hdt);
  void raw_update_Wtil(double hdt);
protected:
  double m_dx, m_dy;
};

//! \brief The PISM subglacial hydrology model for a distributed linked-cavity system.
/*!
  This class implements the model documented in [\ref BuelervanPeltDRAFT].

  Unlike hydrology::Routing, the water pressure \f$P\f$ is a state variable, and there
  are modeled mechanisms for cavity geometry evolution, including creep closure
  and opening through sliding ("cavitation").  Because of cavitation, this model
  needs access to a StressBalance object.   Background references for this kind of
  model includes especially [\ref Kamb1987, \ref Schoofetal2012], but see also
  [\ref Hewitt2011, \ref Hewittetal2012, \ref Hewitt2013].

  In addition to the actions within the null strip taken by hydrology::Routing,
  this model also sets the staggered grid values of the gradient of the hydraulic
  potential to zero if either regular grid neighbor is in the null strip.
*/
class Distributed : public Routing {
public:
  Distributed(IceGrid::ConstPtr g, stressbalance::StressBalance *sb);
  virtual ~Distributed();

  virtual void init();

  friend class Distributed_hydrovelbase_mag;

  virtual void subglacial_water_pressure(IceModelVec2S &result) const;

protected:
  virtual void update_impl(double icet, double icedt);

  virtual std::map<std::string, Diagnostic::Ptr> diagnostics_impl() const;
  virtual std::map<std::string, TSDiagnostic::Ptr> ts_diagnostics_impl() const;

  virtual void define_model_state_impl(const PIO &output) const;
  virtual void write_model_state_impl(const PIO &output) const;

  virtual void init_bwp();

  virtual void check_P_bounds(bool enforce_upper);

  virtual void update_velbase_mag(IceModelVec2S &result);
  virtual void P_from_W_steady(IceModelVec2S &result);

  virtual void adaptive_for_WandP_evolution(double t_current, double t_end, double maxKW,
                                            double &dt_result,
                                            double &maxV_result, double &maxD_result,
                                            double &PtoCFLratio);
protected:
  // this model's state, in addition to what is in hydrology::Routing
  IceModelVec2S m_P;      //!< water pressure
  // this model's auxiliary variables, in addition ...
  IceModelVec2S m_psi,    //!< hydraulic potential
    m_velbase_mag,  //!< sliding speed of overlying ice
    m_Pnew;   //!< pressure during update
  bool m_hold_velbase_mag;

  // need to get basal sliding velocity (thus speed):
  stressbalance::StressBalance* m_stressbalance;
};

} // end of namespace hydrology
} // end of namespace pism

#endif /* _PISMHYDROLOGY_H_ */

