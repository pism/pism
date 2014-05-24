// Copyright (C) 2012-2014 PISM Authors
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

#include <assert.h>

#include "iceModelVec.hh"
#include "iceModelVec2T.hh"
#include "PISMComponent.hh"
#include "PISMStressBalance.hh"

namespace pism {
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
  subglacial_water_pressure() return amount and pressure.  This subglacial water
  is *transportable*, that is, it moves along a modeled hydraulic head gradient.
  Background references for such models include [\ref FlowersClarke2002_theory,
  \ref Hewittetal2012, \ref Schoofetal2012, \ref Hewitt2013].

  These models always have a separate, but potentially-coupled, amount of water
  which is held in local till storage.  It is important to note that the
  transportable water (bwat) and till water (tillwat) thicknesses are different.
  Published models with till storage include [\ref BBssasliding, \ref SchoofTill,
  \ref TrufferEchelmeyerHarrison2001, \ref Tulaczyketal2000b].

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
  Hydrology(IceGrid &g, const Config &conf);
  virtual ~Hydrology();

  virtual PetscErrorCode init(Vars &vars);

  virtual void get_diagnostics(std::map<std::string, Diagnostic*> &dict,
                               std::map<std::string, TSDiagnostic*> &ts_dict);
  friend class Hydrology_hydrobmelt;
  friend class Hydrology_hydroinput;

  // in the base class these only add/define/write tillwat
  virtual void add_vars_to_output(const std::string &keyword, std::set<std::string> &result);
  virtual PetscErrorCode define_variables(const std::set<std::string> &vars, const PIO &nc,
                                          IO_Type nctype);
  virtual PetscErrorCode write_variables(const std::set<std::string> &vars, const PIO &nc);

  // all Hydrology models have a Wtil state variable, which this returns
  virtual PetscErrorCode till_water_thickness(IceModelVec2S &result);

  // this diagnostic method returns the standard shallow approximation
  virtual PetscErrorCode overburden_pressure(IceModelVec2S &result);

  // this diagnostic method returns zero in the base class
  virtual PetscErrorCode wall_melt(IceModelVec2S &result);

  // these methods MUST be implemented in the derived class
  virtual PetscErrorCode subglacial_water_thickness(IceModelVec2S &result) = 0;
  virtual PetscErrorCode subglacial_water_pressure(IceModelVec2S &result) = 0;
  virtual PetscErrorCode update(double icet, double icedt) = 0;

protected:
  // this model's state
  IceModelVec2S Wtil;      // effective thickness of till
  // this model's workspace
  IceModelVec2S total_input, bmelt_local;

  // pointers into IceModel; these describe the ice sheet and the source
  IceModelVec2S *thk,   // ice thickness
    *bed,   // bed elevation (not all models need this)
    *cellarea, // projection-dependent area of each cell, used in mass reporting
    *bmelt; // ice sheet basal melt rate
  IceModelVec2Int *mask;// floating, grounded, etc. mask

  bool hold_bmelt;

  IceModelVec2T *inputtobed;// time dependent input of water to bed, in addition to bmelt
  unsigned int inputtobed_period;      // in years
  double inputtobed_reference_time; // in seconds

  Vars *variables;

  virtual PetscErrorCode get_input_rate(double hydro_t, double hydro_dt, IceModelVec2S &result);

  virtual PetscErrorCode check_Wtil_bounds();
};


//! The PISM minimal model has till in a "can". Water that overflows
//! the can is not conserved. There is no model for lateral transport.
/*!
  This is the minimum functional derived class.  It updates till water thickness.

  It has no transportable water and subglacial_water_thickness() returns zero.

  This model can give no meaningful report on conservation errors.

  Here is a talk which illustrates the "till-can" metaphor:
  http://www2.gi.alaska.edu/snowice/glaciers/iceflow/bueler-igs-fairbanks-june2012.pdf
*/
class NullTransportHydrology : public Hydrology {
public:
  NullTransportHydrology(IceGrid &g, const Config &conf);
  virtual ~NullTransportHydrology();

  virtual PetscErrorCode init(Vars &vars);

  //! Sets result to 0.
  virtual PetscErrorCode subglacial_water_thickness(IceModelVec2S &result);

  //! Returns the overburden pressure in hope it is harmless.
  virtual PetscErrorCode subglacial_water_pressure(IceModelVec2S &result);

  //! Solves an implicit step of a highly-simplified ODE.
  virtual PetscErrorCode update(double icet, double icedt);
};


//! \brief A subglacial hydrology model which assumes water pressure
//! equals overburden pressure.
/*!
  This is the minimal PISM hydrology model that has lateral motion of
  subglacial water and which conserves the water mass.  It was promised
  as a PISM addition in in Bueler's talk at IGS 2012 Fairbanks:
  http://www2.gi.alaska.edu/snowice/glaciers/iceflow/bueler-igs-fairbanks-june2012.pdf

  The water velocity is along the steepest descent route for the hydraulic
  potential.  This potential is (mostly) a function of ice sheet geometry,
  because the water pressure is set to the overburden pressure, a simplified but
  well-established model [\ref Shreve1972].  However, the water layer thickness
  is also a part of the hydraulic potential because it is actually the potential
  of the top of the water layer.

  This (essential) model has been used for finding locations of subglacial lakes
  [\ref Siegertetal2009, \ref Livingstoneetal2013TCD].  Subglacial lakes occur
  at local minima of the hydraulic potential.  If water builds up significantly
  (e.g. thickness of 10s of meters or more) then in the model here the resulting
  lakes diffuse instead of becoming infinitely deep.  Thus we avoid delta
  functions of water thickness at the minima of the hydraulic potential in this
  well-posed model.

  This model should generally be tested using static ice geometry first, i.e.
  using option -no_mass.

  Use option `-report_mass_accounting` to see stdout reports which balance the
  books on this model.

  The state space includes both the till water effective thickness \f$W_{til}\f$,
  which is in Hydrology, and the transportable water layer thickness \f$W\f$.

  For more complete modeling where the water pressure is determined by a
  physical model for the opening and closing of cavities, and where the state
  space includes a nontrivial pressure variable, see DistributedHydrology.

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
class RoutingHydrology : public Hydrology {
public:
  RoutingHydrology(IceGrid &g, const Config &conf);
  virtual ~RoutingHydrology();

  virtual PetscErrorCode init(Vars &vars);

  virtual void add_vars_to_output(const std::string &keyword, std::set<std::string> &result);
  virtual PetscErrorCode define_variables(const std::set<std::string> &vars, const PIO &nc,
                                          IO_Type nctype);
  virtual PetscErrorCode write_variables(const std::set<std::string> &vars, const PIO &nc);

  virtual void get_diagnostics(std::map<std::string, Diagnostic*> &dict,
                               std::map<std::string, TSDiagnostic*> &ts_dict);

  virtual PetscErrorCode wall_melt(IceModelVec2S &result);

  virtual PetscErrorCode subglacial_water_thickness(IceModelVec2S &result);

  virtual PetscErrorCode subglacial_water_pressure(IceModelVec2S &result);

  virtual PetscErrorCode update(double icet, double icedt);

protected:
  // this model's state
  IceModelVec2S W;      // water layer thickness
  // this model's auxiliary variables
  IceModelVec2Stag V,   // components are
  //   V(i,j,0) = u(i,j) = east-edge  centered x-component of water velocity
  //   V(i,j,1) = v(i,j) = north-edge centered y-component of water velocity
    Wstag,// edge-centered (staggered) W values (averaged from regular)
    Kstag,// edge-centered (staggered) values of nonlinear conductivity
    Qstag;// edge-centered (staggered) advection fluxes
  // this model's workspace variables
  IceModelVec2S Wnew, Wtilnew, Pover, R;

  double stripwidth; // width in m of strip around margin where V and W are set to zero;
  // if negative then the strip mechanism is inactive inactive

  PetscErrorCode allocate();
  virtual PetscErrorCode init_bwat(Vars &vars);

  // when we update the water amounts, careful mass accounting at the
  // boundary is needed; we update the new thickness variable, typically a
  // temporary during the update
  bool report_mass_accounting;
  virtual PetscErrorCode boundary_mass_changes(IceModelVec2S &newthk,
                                               double &icefreelost, double &oceanlost,
                                               double &negativegain, double &nullstriplost);

  virtual PetscErrorCode check_water_thickness_nonnegative(IceModelVec2S &thk);

  virtual PetscErrorCode water_thickness_staggered(IceModelVec2Stag &result);
  virtual PetscErrorCode subglacial_hydraulic_potential(IceModelVec2S &result);

  virtual PetscErrorCode conductivity_staggered(IceModelVec2Stag &result, double &maxKW);
  virtual PetscErrorCode velocity_staggered(IceModelVec2Stag &result);
  friend class RoutingHydrology_bwatvel;  // needed because bwatvel diagnostic needs protected velocity_staggered()
  virtual PetscErrorCode advective_fluxes(IceModelVec2Stag &result);

  virtual PetscErrorCode adaptive_for_W_evolution(
                                                  double t_current, double t_end, double maxKW,
                                                  double &dt_result,
                                                  double &maxV_result, double &maxD_result,
                                                  double &dtCFL_result, double &dtDIFFW_result);

  PetscErrorCode raw_update_W(double hdt);
  PetscErrorCode raw_update_Wtil(double hdt);
};


//! \brief The PISM subglacial hydrology model for a distributed linked-cavity system.
/*!
  This implements the new Bueler & van Pelt model documented at the repo (currently
  private):
  https://github.com/bueler/hydrolakes
  Unlike RoutingHydrology, the water pressure P is a state variable, and there
  are modeled mechanisms for cavity geometry evolution, including creep closure
  and opening through sliding ("cavitation").  Because of cavitation, this model
  needs access to a StressBalance object.

  In addition to the actions within the null strip taken by RoutingHydrology,
  this model also sets the staggered grid values of the gradient of the hydraulic
  potential to zero if either regular grid neighbor is in the null strip.
*/
class DistributedHydrology : public RoutingHydrology {
public:
  DistributedHydrology(IceGrid &g, const Config &conf, StressBalance *sb);
  virtual ~DistributedHydrology();

  virtual PetscErrorCode init(Vars &vars);

  virtual void add_vars_to_output(const std::string &keyword, std::set<std::string> &result);
  virtual void get_diagnostics(std::map<std::string, Diagnostic*> &dict,
                               std::map<std::string, TSDiagnostic*> &ts_dict);
  virtual PetscErrorCode define_variables(const std::set<std::string> &vars, const PIO &nc,
                                          IO_Type nctype);
  virtual PetscErrorCode write_variables(const std::set<std::string> &vars, const PIO &nc);

  virtual PetscErrorCode update(double icet, double icedt);

  virtual PetscErrorCode subglacial_water_pressure(IceModelVec2S &result);

protected:
  // this model's state, in addition to what is in RoutingHydrology
  IceModelVec2S P;      //!< water pressure
  // this model's auxiliary variables, in addition ...
  IceModelVec2S psi,    //!< hydraulic potential
    velbase_mag,  //!< sliding speed of overlying ice
    Pnew;   //!< pressure during update

  // need to get basal sliding velocity (thus speed):
  StressBalance* stressbalance;

  PetscErrorCode allocate_pressure();
  virtual PetscErrorCode init_bwp(Vars &vars);

  virtual PetscErrorCode check_P_bounds(bool enforce_upper);

  virtual PetscErrorCode update_velbase_mag(IceModelVec2S &result);
  virtual PetscErrorCode P_from_W_steady(IceModelVec2S &result);

  virtual PetscErrorCode adaptive_for_WandP_evolution(
                                                      double t_current, double t_end, double maxKW,
                                                      double &dt_result,
                                                      double &maxV_result, double &maxD_result,
                                                      double &PtoCFLratio);
};

} // end of namespace pism

#endif /* _PISMHYDROLOGY_H_ */

