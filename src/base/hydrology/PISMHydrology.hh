// Copyright (C) 2012-2013 PISM Authors
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


//! \brief The PISM subglacial hydrology model interface.
/*!
This is a virtual base class.

The purpose of this class and its derived classes is to provide
\code
  subglacial_water_thickness(IceModelVec2S &result)
  subglacial_water_pressure(IceModelVec2S &result)
  englacial_water_thickness(IceModelVec2S &result)
  till_water_thickness(IceModelVec2S &result)
  till_water_pressure(IceModelVec2S &result)
\endcode

Additional modeled fields, for diagnostic purposes, are
\code
  overburden_pressure(IceModelVec2S &result)
  wall_melt(IceModelVec2S &result)
\endcode

This interface is specific to subglacial hydrology models which track a
two-dimensional water layer with a well-defined thickness and pressure at each
map-plane location.  The methods subglacial_water_thickness() and
subglacial_water_pressure() return amount and pressure.  This subglacial water
is *transportable*, that is, it moves along a modeled hydraulic head gradient.
For more information see [\ref vanPeltBuelerDRAFT].  Background references for
such models include [\ref FlowersClarke2002_theory, \ref Hewittetal2012,
\ref Schoofetal2012].

These models alway have a separate, but potentially coupled, amount of water
which is held in local till storage.  Generally the transportable water (bwat)
and till water (tillwat) thicknesses are different.  Also, generally the
tranportable water (bwp) till water (tillwp) pressures are different.
References for models with till storage include [\ref BBssasliding,
\ref SchoofTill, \ref TrufferEchelmeyerHarrison, \ref Tulaczyketal2000b].

The till water pressure, not the tranportable water pressure, is used by the
Mohr-Coulomb criterion to provide a yield stress.

The base class does not implement the evolution of the till water thickness.
It does not report a till water pressure.

These models also either track the amount of englacial water, in a manner
which allows computation of an effective thickness and which is returned by
englacial_water_thickness(), or they lack the mechanism and
englacial_water_thickness() returns zero.  A reference for such a model with
englacial storage is [\ref Bartholomausetal2011].

PISMHydrology is a timestepping component (PISMComponent_TS).  Because of the
short physical timescales associated to liquid water moving under a glacier,
PISMHydrology derived classes generally take many substeps in PISM's major
ice dynamics time steps.  Thus when an update() method in a PISMHydrology
derived class is called it will advance its internal time to the new goal t+dt
using its own internal time steps.

Generally PISMHydrology and derived classes use the ice geometry, the basal melt
rate, and the basal sliding velocity.  Note that the basal melt rate is an
energy-conservation-derived field.  These fields generally
come from IceModel and PISMStressBalance.  Additionally, time-dependent
and spatially-variable water input to the basal layer, taken directly from a
file, is possible too.  Potentially PISMSurfaceModel could supply such a
quantity.

Ice geometry and energy fields are normally treated as constant in time
during the update() call for the interval [t,t+dt].  Thus the coupling is
one-way during the update() call.
 */
class PISMHydrology : public PISMComponent_TS {
public:
  PISMHydrology(IceGrid &g, const NCConfigVariable &conf);
  virtual ~PISMHydrology() {}

  virtual PetscErrorCode init(PISMVars &vars);

  virtual PetscErrorCode regrid(IceModelVec2S &myvar);

  virtual void get_diagnostics(map<string, PISMDiagnostic*> &dict);
  friend class PISMHydrology_hydroinput;

  virtual PetscErrorCode max_timestep(PetscReal my_t, PetscReal &my_dt, bool &restrict_dt);

  // in the base class these only add/define/write tillwat
  virtual void add_vars_to_output(string keyword, set<string> &result);
  virtual PetscErrorCode define_variables(set<string> vars, const PIO &nc,PISM_IO_Type nctype);
  virtual PetscErrorCode write_variables(set<string> vars, const PIO &nc);

  // all PISMHydrology models have a Wtil state variable, which this returns
  virtual PetscErrorCode till_water_thickness(IceModelVec2S &result);

  // this exists in the base class and sets result = 0:
  virtual PetscErrorCode englacial_water_thickness(IceModelVec2S &result);

  // this diagnostic method returns the standard shallow approximation
  virtual PetscErrorCode overburden_pressure(IceModelVec2S &result);

  // this diagnostic method returns zero in the base class
  virtual PetscErrorCode wall_melt(IceModelVec2S &result);

  // these methods MUST be implemented in the derived class
  virtual PetscErrorCode subglacial_water_thickness(IceModelVec2S &result) = 0;
  virtual PetscErrorCode subglacial_water_pressure(IceModelVec2S &result) = 0;
  virtual PetscErrorCode till_water_pressure(IceModelVec2S &result) = 0;
  using PISMComponent_TS::update;
  virtual PetscErrorCode update(PetscReal icet, PetscReal icedt) = 0;

protected:
  // this model's state
  IceModelVec2S Wtil;      // effective thickness of till
  // this model's workspace
  IceModelVec2S total_input;

  // pointers into IceModel; these describe the ice sheet and the source
  IceModelVec2S *thk,   // ice thickness
                *bed,   // bed elevation (not all models need this)
                *cellarea, // projection-dependent area of each cell, used in mass reporting
                *bmelt; // ice sheet basal melt rate
  IceModelVec2Int *mask;// floating, grounded, etc. mask
  PISMVars *variables;

  // for time dependent input of water to bed (in addition to bmelt)
  IceModelVec2T *inputtobed;
  PetscReal     inputtobed_period, inputtobed_reference_time;
  virtual PetscErrorCode get_input_rate(
                            PetscReal hydro_t, PetscReal hydro_dt, IceModelVec2S &result);

  virtual PetscErrorCode check_Wtil_bounds();
};


//! The PISM minimal model has till in a "can".  Water that overflows the can is not conserved.  Thus there is no true hydrology, i.e. no model for transport.
/*!
This is the minimum functional derived class.  It updates till water thickness.
It returns a simple model for the pressure of the water stored in till.

It has no tranportable water so subglacial_water_thickness() returns zero.

The method subglacial_water_pressure() is trivialized: it returns overburden
pressure.

This model can give no meaningful report on conservation errors.

Here is a talk which illustrates the "till-can" metaphor:
  http://www2.gi.alaska.edu/snowice/glaciers/iceflow/bueler-igs-fairbanks-june2012.pdf
 */
class PISMNullTransportHydrology : public PISMHydrology {
public:
  PISMNullTransportHydrology(IceGrid &g, const NCConfigVariable &conf)
    : PISMHydrology(g, conf) {}
  virtual ~PISMNullTransportHydrology() {}

  // sets result = 0
  virtual PetscErrorCode subglacial_water_thickness(IceModelVec2S &result);

  // sets result = overburden pressure
  virtual PetscErrorCode subglacial_water_pressure(IceModelVec2S &result);

  // sets result = Bueler&Brown version of pressure of till water
  virtual PetscErrorCode till_water_pressure(IceModelVec2S &result);

  // solves an implicit step of a highly-simplified ODE
  virtual PetscErrorCode update(PetscReal icet, PetscReal icedt);
};


//! \brief A subglacial hydrology model which assumes water pressure equals overburden pressure.
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
functions at the minima of the hydraulic potential.  This model is a
well-posed PDE which finds subglacial lakes.

This model should generally be tested using static ice geometry first, i.e.
using option -no_mass.

Use option `-report_mass_accounting` to see stdout reports which balance the
books on this model.

The state space includes both the till water effective thickness \f$W_{til}\f$,
which is in PISMHydrology, and the transportable water layer thickness \f$W\f$.

For more complete modeling where the water pressure is determined by a
physical model for the opening and closing of cavities, and where the state
space includes a nontrivial pressure variable, see PISMDistributedHydrology.

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
    \f[ K = k W^{\alpha-1} |\nabla (P+\rho_w g b)|^{\beta - 2}. \f]

This model contains enough information (modeled fields) so that we can
compute the wall melt generated by dissipating the gravitational
potential energy in the moving, presumably turbulent, subglacial water.  If we
suppose that this heat is dissipated immediately as melt on the
cavity/conduit walls then we get a formula for a wall melt contribution.  (This
is in addition to the `bmelt` field coming from conserving energy in the flowing
ice.)  See wall_melt().  At this time the wall melt is diagnostic only and does
not add to the water amount in the till, because that addition is clearly
unstable.
 */
class PISMRoutingHydrology : public PISMHydrology {
public:
  PISMRoutingHydrology(IceGrid &g, const NCConfigVariable &conf);
  virtual ~PISMRoutingHydrology() {}

  virtual PetscErrorCode init(PISMVars &vars);

  virtual void add_vars_to_output(string keyword, set<string> &result);
  virtual PetscErrorCode define_variables(set<string> vars, const PIO &nc,PISM_IO_Type nctype);
  virtual PetscErrorCode write_variables(set<string> vars, const PIO &nc);

  virtual void get_diagnostics(map<string, PISMDiagnostic*> &dict);

  virtual PetscErrorCode update(PetscReal icet, PetscReal icedt);

  virtual PetscErrorCode subglacial_water_thickness(IceModelVec2S &result);
  virtual PetscErrorCode subglacial_water_pressure(IceModelVec2S &result);
  virtual PetscErrorCode subglacial_hydraulic_potential(IceModelVec2S &result);
  virtual PetscErrorCode wall_melt(IceModelVec2S &result);

  virtual PetscErrorCode velocity_staggered(IceModelVec2Stag &result);

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
  IceModelVec2S Wnew, Pover, R;

  PetscReal stripwidth; // width in m of strip around margin where V and W are set to zero;
                        // if negative then the strip mechanism is inactive inactive

  virtual PetscErrorCode allocate();
  virtual PetscErrorCode init_bwat(PISMVars &vars, bool i_set, bool bootstrap_set);

  // when we update the transportable water, careful mass accounting at the
  // boundary is needed; we update Wnew and so model state (Wtil or W) is not touched
  bool report_mass_accounting;
  virtual PetscErrorCode boundary_mass_changes(IceModelVec2S &Wnew,
             PetscReal &icefreelost, PetscReal &oceanlost,
             PetscReal &negativegain, PetscReal &nullstriplost);

  virtual PetscErrorCode check_W_nonnegative();
  virtual PetscErrorCode water_thickness_staggered(IceModelVec2Stag &result);

  virtual PetscErrorCode conductivity_staggered(IceModelVec2Stag &result, PetscReal &maxKW);
  virtual PetscErrorCode advective_fluxes(IceModelVec2Stag &result);

  virtual PetscErrorCode adaptive_for_W_evolution(
            PetscReal t_current, PetscReal t_end, PetscReal maxKW,
            PetscReal &dt_result,
            PetscReal &maxV_result, PetscReal &maxD_result,
            PetscReal &dtCFL_result, PetscReal &dtDIFFW_result);

  PetscErrorCode raw_update_W(PetscReal hdt);

  inline bool in_null_strip(PetscInt i, PetscInt j) {
    if (stripwidth < 0.0) return false;
    return ((grid.x[i] <= grid.x[0] + stripwidth) || (grid.x[i] >= grid.x[grid.Mx-1] - stripwidth)
            || (grid.y[j] <= grid.y[0] + stripwidth) || (grid.y[j] >= grid.y[grid.My-1] - stripwidth));
  }
};


//! \brief The PISM subglacial hydrology model for a distributed linked-cavity system.
/*!
This implements the new van Pelt & Bueler model documented at the repo (currently
private):
  https://github.com/bueler/hydrolakes
Unlike PISMRoutingHydrology, the water pressure P is a state variable, and there
are modeled mechanisms for cavity geometry evolution, including creep closure
and opening through sliding ("cavitation").  Because of cavitation, this model
needs access to a PISMStressBalance object.

In addition to the actions within the null strip taken by PISMRoutingHydrology,
this model also sets the staggered grid values of the gradient of the hydraulic
potential to zero if either regular grid neighbor is in the null strip.
 */
class PISMDistributedHydrology : public PISMRoutingHydrology {
public:
  PISMDistributedHydrology(IceGrid &g, const NCConfigVariable &conf, PISMStressBalance *sb);
  virtual ~PISMDistributedHydrology() {}

  virtual PetscErrorCode init(PISMVars &vars);

  virtual void add_vars_to_output(string keyword, set<string> &result);
  virtual void get_diagnostics(map<string, PISMDiagnostic*> &dict);
  virtual PetscErrorCode define_variables(set<string> vars, const PIO &nc,PISM_IO_Type nctype);
  virtual PetscErrorCode write_variables(set<string> vars, const PIO &nc);

  virtual PetscErrorCode update(PetscReal icet, PetscReal icedt);

  virtual PetscErrorCode subglacial_water_pressure(IceModelVec2S &result);
  virtual PetscErrorCode englacial_water_thickness(IceModelVec2S &result);

protected:
  // this model's state, in addition to what is in PISMRoutingHydrology
  IceModelVec2S Wen,    // englacial water thickness
                P;      // water pressure
  // this model's auxiliary variables, in addition ...
  IceModelVec2S psi,    // hydraulic potential
                cbase,  // sliding speed of overlying ice
                Pnew;   // pressure during update

  // need to get basal sliding velocity (thus speed):
  PISMStressBalance* stressbalance;

  virtual PetscErrorCode allocate_englacial();
  virtual PetscErrorCode allocate_pressure();

  virtual PetscErrorCode check_Wen_nonnegative();
  virtual PetscErrorCode check_P_bounds(bool enforce_upper);

  virtual PetscErrorCode update_cbase(IceModelVec2S &result);
  virtual PetscErrorCode P_from_W_steady(IceModelVec2S &result);

  virtual PetscErrorCode adaptive_for_WandP_evolution(
                           PetscReal t_current, PetscReal t_end, PetscReal maxKW,
                           PetscReal &dt_result,
                           PetscReal &maxV_result, PetscReal &maxD_result,
                           PetscReal &PtoCFLratio);

  virtual PetscErrorCode update_englacial_storage(
                               IceModelVec2S &myPnew, IceModelVec2S &Wnew_tot);
};

#endif /* _PISMHYDROLOGY_H_ */

