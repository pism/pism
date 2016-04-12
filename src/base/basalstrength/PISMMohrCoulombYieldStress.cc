// Copyright (C) 2004--2016 PISM Authors
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

#include <cmath>
#include <cassert>
#include <gsl/gsl_math.h>

#include "PISMMohrCoulombYieldStress.hh"
#include "base/hydrology/PISMHydrology.hh"
#include "base/util/IceGrid.hh"
#include "base/util/Mask.hh"
#include "base/util/PISMVars.hh"
#include "base/util/error_handling.hh"
#include "base/util/io/PIO.hh"
#include "base/util/pism_options.hh"
#include "base/util/MaxTimestep.hh"
#include "base/util/pism_utilities.hh"
#include "base/util/IceModelVec2CellType.hh"

namespace pism {

//! \file PISMMohrCoulombYieldStress.cc  Process model which computes pseudo-plastic yield stress for the subglacial layer.
/*! \file PISMMohrCoulombYieldStress.cc
The output variable of this submodel is `tauc`, the pseudo-plastic yield stress
field that is used in the ShallowStressBalance objects.  This quantity is
computed by the Mohr-Coulomb criterion [\ref SchoofTill], but using an empirical
relation between the amount of water in the till and the effective pressure
of the overlying glacier resting on the till [\ref Tulaczyketal2000].

The "dry" strength of the till is a state variable which is private to
the submodel, namely `tillphi`.  Its initialization is nontrivial: either the
`-topg_to_phi`  heuristic is used or inverse modeling can be used.  (In the
latter case `tillphi` can be read-in at the beginning of the run.  Currently
`tillphi` does not evolve during the run.)

This submodel uses a pointer to a Hydrology instance to get the till (pore)
water amount, the effective water layer thickness.  The effective pressure is
derived from this.  Then the effective pressure is combined with tillphi to
compute an updated `tauc` by the Mohr-Coulomb criterion.

This submodel is inactive in floating areas.
*/


MohrCoulombYieldStress::MohrCoulombYieldStress(IceGrid::ConstPtr g,
                                               hydrology::Hydrology *hydro)
  : YieldStress(g) {

  m_topg_to_phi = false;
  m_tauc_to_phi = false;

  m_hydrology = hydro;

  unsigned int stencil_width = m_config->get_double("grid_max_stencil_width");

  m_till_phi.create(m_grid, "tillphi", WITH_GHOSTS, stencil_width);
  m_till_phi.set_attrs("model_state",
                       "friction angle for till under grounded ice sheet",
                       "degrees", "");
  m_till_phi.set_time_independent(true);
  // in this model; need not be time-independent in general

  // internal working space; stencil width needed because redundant computation
  // on overlaps
  m_tillwat.create(m_grid, "tillwat_for_MohrCoulomb",
                   WITH_GHOSTS, stencil_width);
  m_tillwat.set_attrs("internal",
                      "copy of till water thickness held by MohrCoulombYieldStress",
                      "m", "");
  bool addtransportable = m_config->get_boolean("tauc_add_transportable_water");
  if (addtransportable == true) {
    m_bwat.create(m_grid, "bwat_for_MohrCoulomb", WITHOUT_GHOSTS);
    m_bwat.set_attrs("internal",
                     "copy of transportable water thickness held by MohrCoulombYieldStress",
                     "m", "");
  }
  m_Po.create(m_grid, "overburden_pressure_for_MohrCoulomb",
              WITH_GHOSTS, stencil_width);
  m_Po.set_attrs("internal",
                 "copy of overburden pressure held by MohrCoulombYieldStress",
                 "Pa", "");
}

MohrCoulombYieldStress::~MohrCoulombYieldStress() {
  // empty
}


//! Initialize the pseudo-plastic till mechanical model.
/*!
The pseudo-plastic till basal resistance model is governed by this power law
equation,
    \f[ \tau_b = - \frac{\tau_c}{|\mathbf{U}|^{1-q} U_{\mathtt{th}}^q} \mathbf{U}, \f]
where \f$\tau_b=(\tau_{(b)x},\tau_{(b)y})\f$ is the basal shear stress and
\f$U=(u,v)\f$ is the sliding velocity.

We call the scalar field \f$\tau_c(t,x,y)\f$ the \e yield \e stress even when
the power \f$q\f$ is not zero; when that power is zero the formula describes
a plastic material with an actual yield stress.  The constant
\f$U_{\mathtt{th}}\f$ is the \e threshhold \e speed, and \f$q\f$ is the \e pseudo
\e plasticity \e exponent.  The current class computes this yield stress field.
See also IceBasalResistancePlasticLaw::drag().

The strength of the saturated till material, the yield stress, is modeled by a
Mohr-Coulomb relation [\ref Paterson, \ref SchoofStream],
    \f[   \tau_c = c_0 + (\tan \varphi) N_{til}, \f]
where \f$N_{til}\f$ is the effective pressure of the glacier on the mineral
till.

The determination of the till friction angle \f$\varphi(x,y)\f$  is important.
It is assumed in this default model to be a time-independent factor which
describes the strength of the unsaturated "dry" (mineral) till material.  Thus
it is assumed to change more slowly than the till water pressure, and it follows
that it changes more slowly than the yield stress and the basal shear stress.

Option `-topg_to_phi` causes call to topg_to_phi() at the beginning of the run.
This determines the map of \f$\varphi(x,y)\f$.  If this option is note given,
the current method leaves `tillphi` unchanged, and thus either in its
read-in-from-file state or with a default constant value from the config file.
*/
void MohrCoulombYieldStress::init_impl() {
  {
    std::string hydrology_tillwat_max = "hydrology_tillwat_max";
    bool till_is_present = m_config->get_double(hydrology_tillwat_max) > 0.0;

    if (not till_is_present) {
      throw RuntimeError::formatted("The Mohr-Coulomb yield stress model cannot be used without till.\n"
                                    "Reset %s or choose a different yield stress model.",
                                    hydrology_tillwat_max.c_str());
    }
  }

  {
    const std::string flag_name = "tauc_add_transportable_water";
    hydrology::Routing *hydrology_routing = dynamic_cast<hydrology::Routing*>(m_hydrology);
    if (m_config->get_boolean(flag_name) == true && hydrology_routing == NULL) {
      throw RuntimeError::formatted("Flag %s is set.\n"
                                    "Thus the Mohr-Coulomb yield stress model needs a hydrology::Routing\n"
                                    "(or derived like hydrology::Distributed) object with transportable water.\n"
                                    "The current Hydrology instance is not suitable.  Set flag\n"
                                    "%s to 'no' or choose a different yield stress model.",
                                    flag_name.c_str(), flag_name.c_str());
    }
  }

  m_log->message(2, "* Initializing the default basal yield stress model...\n");

  options::Real
    plastic_phi("-plastic_phi", "constant in space till friction angle",
                m_config->get_double("default_till_phi"));

  options::RealList
    topg_to_phi_option("-topg_to_phi",
                       "Turn on, and specify, the till friction angle parameterization"
                       " based on bedrock elevation (topg)");

  bool bootstrap = false;
  int start = 0;
  std::string filename;
  bool use_input_file = find_pism_input(filename, bootstrap, start);

  if (topg_to_phi_option.is_set() and plastic_phi.is_set()) {
    throw RuntimeError("only one of -plastic_phi and -topg_to_phi is allowed.");
  }

  if (topg_to_phi_option.is_set()) {

    m_log->message(2,
                 "  option -topg_to_phi seen; creating tillphi map from bed elev ...\n");

    if (use_input_file) {

      PIO nc(m_grid->com, "guess_mode");

      nc.open(filename, PISM_READONLY);
      bool tillphi_present = nc.inq_var(m_till_phi.metadata().get_string("short_name"));
      nc.close();

      if (tillphi_present) {
        m_log->message(2,
                     "PISM WARNING: -topg_to_phi computation will override the '%s' field\n"
                     "              present in the input file '%s'!\n",
                     m_till_phi.metadata().get_string("short_name").c_str(), filename.c_str());
      }
    }

    // note option -topg_to_phi will be read again to get comma separated array of parameters
    m_topg_to_phi = true;

    double phi_min  = topg_to_phi_option[0];
    double phi_max  = topg_to_phi_option[1];
    double topg_min = topg_to_phi_option[2];
    double topg_max = topg_to_phi_option[3];

    m_log->message(2,
                 "  till friction angle (phi) is piecewise-linear function of bed elev (topg):\n"
                 "            /  %5.2f                                 for   topg < %.f\n"
                 "      phi = |  %5.2f + (topg - (%.f)) * (%.2f / %.f)   for   %.f < topg < %.f\n"
                 "            \\  %5.2f                                 for   %.f < topg\n",
                 phi_min, topg_min,
                 phi_min, topg_min, phi_max - phi_min, topg_max - topg_min, topg_min, topg_max,
                 phi_max, topg_max);

  } else if (use_input_file) {
    if (bootstrap) {
      m_till_phi.regrid(filename, OPTIONAL,
                        m_config->get_double("bootstrapping_tillphi_value_no_var"));
    } else {
      m_till_phi.read(filename, start);
    }
  } else {
    // Use the default value *or* the value set using the -plastic_phi
    // command-line option.
    m_till_phi.set(plastic_phi);
  }

  // regrid if requested, regardless of how initialized
  regrid("MohrCoulombYieldStress", m_till_phi);

  options::String tauc_to_phi_file("-tauc_to_phi",
                                   "Turn on, and specify, the till friction angle computation"
                                   " which uses basal yield stress (tauc) and the rest of the model state",
                                   "", options::ALLOW_EMPTY);

  if (tauc_to_phi_file.is_set()) {

    if (not tauc_to_phi_file->empty()) {
      // "-tauc_to_phi filename.nc" is given
      m_tauc.regrid(tauc_to_phi_file, CRITICAL);
    } else {
      // "-tauc_to_phi" is given (without a file name); assume that tauc has to
      // be present in an input file
      if (bootstrap) {
        m_tauc.regrid(filename, CRITICAL);
      } else {
        m_tauc.read(filename, start);
      }
    }

    m_log->message(2,
                 "  Will compute till friction angle (tillphi) as a function"
                 " of the yield stress (tauc)...\n");

    m_tauc_to_phi = true;

  } else {
    m_tauc.set(0.0);
  }

  // ensure that update() computes tauc at the beginning of the run:
  m_t = m_dt = GSL_NAN;
}

MaxTimestep MohrCoulombYieldStress::max_timestep_impl(double t) {
  (void) t;
  return MaxTimestep();
}

void MohrCoulombYieldStress::set_till_friction_angle(const IceModelVec2S &input) {
  m_till_phi.copy_from(input);
}

void MohrCoulombYieldStress::add_vars_to_output_impl(const std::string &/*keyword*/,
                                                     std::set<std::string> &result) {
  result.insert("tillphi");
}


void MohrCoulombYieldStress::define_variables_impl(const std::set<std::string> &vars,
                                                   const PIO &nc, IO_Type nctype) {
  if (set_contains(vars, "tillphi")) {
    m_till_phi.define(nc, nctype);
  }
}


void MohrCoulombYieldStress::write_variables_impl(const std::set<std::string> &vars,
                                                  const PIO &nc) {
  if (set_contains(vars, "tillphi")) {
    m_till_phi.write(nc);
  }
}


//! Update the till yield stress for use in the pseudo-plastic till basal stress
//! model.  See also IceBasalResistancePlasticLaw.
/*!
Updates yield stress  @f$ \tau_c @f$  based on modeled till water layer thickness
from a Hydrology object.  We implement the Mohr-Coulomb criterion allowing
a (typically small) till cohesion  @f$ c_0 @f$
and by expressing the coefficient as the tangent of a till friction angle
 @f$ \varphi @f$ :
    @f[   \tau_c = c_0 + (\tan \varphi) N_{til}. @f]
See [@ref Paterson] table 8.1 regarding values.

The effective pressure on the till is empirically-related
to the amount of water in the till.  We use this formula derived from
[@ref Tulaczyketal2000] and documented in [@ref BuelervanPeltDRAFT]:

@f[ N_{til} = \min\left\{P_o, N_0 \left(\frac{\delta P_o}{N_0}\right)^s 10^{(e_0/C_c) (1 - s)}\right\} @f]

where  @f$ s = W_{til} / W_{til}^{max} @f$,  @f$ W_{til}^{max} @f$ =`hydrology_tillwat_max`,
@f$ \delta @f$ =`till_effective_fraction_overburden`,  @f$ P_o @f$  is the
overburden pressure,  @f$ N_0 @f$ =`till_reference_effective_pressure` is a
reference effective pressure,   @f$ e_0 @f$ =`till_reference_void_ratio` is the void ratio
at the reference effective pressure, and  @f$ C_c @f$ =`till_compressibility_coefficient`
is the coefficient of compressibility of the till.  Constants  @f$ N_0, e_0, C_c @f$  are
found by [@ref Tulaczyketal2000] from laboratory experiments on samples of
till.

If `tauc_add_transportable_water` is yes then @f$ s @f$ in the above formula
becomes @f$ s = (W + W_{til}) / W_{til}^{max} @f$,
that is, the water amount is the sum @f$ W+W_{til} @f$.  This only works
if @f$ W @f$ is present, that is, if `hydrology` points to a
hydrology::Routing (or derived class thereof).
 */
void MohrCoulombYieldStress::update_impl(double my_t, double my_dt) {

  if (m_topg_to_phi) {
    this->topg_to_phi();
    m_topg_to_phi = false;
  }

  if (m_tauc_to_phi) {
    this->tauc_to_phi();
    m_tauc_to_phi = false;
  }

  if ((fabs(my_t - m_t) < 1e-12) &&
      (fabs(my_dt - m_dt) < 1e-12)) {
    return;
  }

  m_t = my_t;
  m_dt = my_dt;
  // this model does no internal time-stepping

  bool slipperygl       = m_config->get_boolean("tauc_slippery_grounding_lines"),
       addtransportable = m_config->get_boolean("tauc_add_transportable_water");

  const double high_tauc   = m_config->get_double("high_tauc"),
               tillwat_max = m_config->get_double("hydrology_tillwat_max"),
               c0          = m_config->get_double("till_cohesion"),
               N0          = m_config->get_double("till_reference_effective_pressure"),
               e0overCc    = m_config->get_double("till_reference_void_ratio")
                                / m_config->get_double("till_compressibility_coefficient"),
               delta       = m_config->get_double("till_effective_fraction_overburden"),
               tlftw       = m_config->get_double("till_log_factor_transportable_water");

  hydrology::Routing* hydrowithtransport = dynamic_cast<hydrology::Routing*>(m_hydrology);
  if (m_hydrology) {
    m_hydrology->till_water_thickness(m_tillwat);
    m_hydrology->overburden_pressure(m_Po);
    if (addtransportable == true) {
        assert(hydrowithtransport != NULL);
        hydrowithtransport->subglacial_water_thickness(m_bwat);
    }
  }

  const IceModelVec2CellType &mask           = *m_grid->variables().get_2d_cell_type("mask");
  const IceModelVec2S        &bed_topography = *m_grid->variables().get_2d_scalar("bedrock_altitude");

  IceModelVec::AccessList list;
  if (addtransportable == true) {
    list.add(m_bwat);
  }
  list.add(m_tillwat);
  list.add(m_till_phi);
  list.add(m_tauc);
  list.add(mask);
  list.add(bed_topography);
  list.add(m_Po);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.ocean(i, j)) {
      m_tauc(i, j) = 0.0;
    } else if (mask.ice_free(i, j)) {
      m_tauc(i, j) = high_tauc;  // large yield stress if grounded and ice-free
    } else { // grounded and there is some ice
      // user can ask that marine grounding lines get special treatment
      const double sea_level = 0.0; // FIXME: get sea-level from correct PISM source
      double water = m_tillwat(i,j); // usual case
      if (slipperygl == true &&
          bed_topography(i,j) <= sea_level &&
          (mask.next_to_floating_ice(i,j) || mask.next_to_ice_free_ocean(i,j))) {
        water = tillwat_max;
      } else if (addtransportable == true) {
        water = m_tillwat(i,j) + tlftw * log(1.0 + m_bwat(i,j) / tlftw);
      }
      double s    = water / tillwat_max,
        Ntil = N0 * pow(delta * m_Po(i,j) / N0, s) * pow(10.0, e0overCc * (1.0 - s));
      Ntil = std::min(m_Po(i,j), Ntil);

      m_tauc(i, j) = c0 + Ntil * tan((M_PI/180.0) * m_till_phi(i, j));
    }
  }

  m_tauc.update_ghosts();
}

//! Computes the till friction angle phi as a piecewise linear function of bed elevation, according to user options.
/*!
Computes the till friction angle \f$\phi(x,y)\f$ at a location as the following
increasing, piecewise-linear function of the bed elevation \f$b(x,y)\f$.  Let
        \f[ M = (\phi_{\text{max}} - \phi_{\text{min}}) / (b_{\text{max}} - b_{\text{min}}) \f]
be the slope of the nontrivial part.  Then
        \f[ \phi(x,y) = \begin{cases}
                \phi_{\text{min}}, & b(x,y) \le b_{\text{min}}, \\
                \phi_{\text{min}} + (b(x,y) - b_{\text{min}}) \,M,
                                  &  b_{\text{min}} < b(x,y) < b_{\text{max}}, \\
                \phi_{\text{max}}, & b_{\text{max}} \le b(x,y), \end{cases} \f]
where \f$\phi_{\text{min}}=\f$`phi_min`, \f$\phi_{\text{max}}=\f$`phi_max`,
\f$b_{\text{min}}=\f$`topg_min`, \f$b_{\text{max}}=\f$`topg_max`.

The default values are vaguely suitable for Antarctica.  See src/pism_config.cdl.
*/
void MohrCoulombYieldStress::topg_to_phi() {

  double phi_min  = m_config->get_double("till_topg_to_phi_phi_min"),
         phi_max  = m_config->get_double("till_topg_to_phi_phi_max"),
         topg_min = m_config->get_double("till_topg_to_phi_topg_min"),
         topg_max = m_config->get_double("till_topg_to_phi_topg_max");

  options::RealList option("-topg_to_phi",
                           "Turn on, and specify, the till friction angle parameterization"
                           " based on bedrock elevation (topg)");

  if (option.is_set() and option->size() != 4) {
    throw RuntimeError::formatted("invalid -topg_to_phi arguments: has to be a list"
                                  " of 4 numbers, got %d", (int)option->size());
  }

  if (option.is_set()) {
    phi_min  = option[0];
    phi_max  = option[1];
    topg_min = option[2];
    topg_max = option[3];
  }

  if (phi_min >= phi_max) {
    throw RuntimeError("invalid -topg_to_phi arguments: phi_min < phi_max is required");
  }

  if (topg_min >= topg_max) {
    throw RuntimeError("invalid -topg_to_phi arguments: topg_min < topg_max is required");
  }

  const IceModelVec2S &bed_topography = *m_grid->variables().get_2d_scalar("bedrock_altitude");

  double slope = (phi_max - phi_min) / (topg_max - topg_min);

  IceModelVec::AccessList list;
  list.add(bed_topography);
  list.add(m_till_phi);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    const double bed = bed_topography(i, j);

    if (bed <= topg_min) {
      m_till_phi(i, j) = phi_min;
    } else if (bed >= topg_max) {
      m_till_phi(i, j) = phi_max;
    } else {
      m_till_phi(i, j) = phi_min + (bed - topg_min) * slope;
    }
  }

  // communicate ghosts so that the tauc computation can be performed locally
  // (including ghosts of tauc, that is)
  m_till_phi.update_ghosts();
}


void MohrCoulombYieldStress::tauc_to_phi() {
  const double c0 = m_config->get_double("till_cohesion"),
    N0            = m_config->get_double("till_reference_effective_pressure"),
    e0overCc      = (m_config->get_double("till_reference_void_ratio") /
                     m_config->get_double("till_compressibility_coefficient")),
    delta         = m_config->get_double("till_effective_fraction_overburden"),
    tillwat_max   = m_config->get_double("hydrology_tillwat_max");

  assert(m_hydrology != NULL);

  m_hydrology->till_water_thickness(m_tillwat);
  m_hydrology->overburden_pressure(m_Po);

  const IceModelVec2CellType &mask = *m_grid->variables().get_2d_cell_type("mask");

  IceModelVec::AccessList list;
  list.add(mask);
  list.add(m_tauc);
  list.add(m_tillwat);
  list.add(m_Po);
  list.add(m_till_phi);

  // make sure that we have enough ghosts:
  unsigned int GHOSTS = m_till_phi.get_stencil_width();
  assert(mask.get_stencil_width()      >= GHOSTS);
  assert(m_tauc.get_stencil_width()    >= GHOSTS);
  assert(m_tillwat.get_stencil_width() >= GHOSTS);
  assert(m_Po.get_stencil_width()      >= GHOSTS);

  for (PointsWithGhosts p(*m_grid, GHOSTS); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.ocean(i, j)) {
      // no change
    } else if (mask.ice_free(i, j)) {
      // no change
    } else { // grounded and there is some ice
      double s    = m_tillwat(i,j) / tillwat_max,
        Ntil = N0 * pow(delta * m_Po(i,j) / N0, s) * pow(10.0, e0overCc * (1.0 - s));
      Ntil = std::min(m_Po(i,j), Ntil);
      m_till_phi(i, j) = 180.0/M_PI * atan((m_tauc(i, j) - c0) / Ntil);
    }
  }
}

} // end of namespace pism
