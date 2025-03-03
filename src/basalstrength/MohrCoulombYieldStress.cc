// Copyright (C) 2004--2023, 2025 PISM Authors
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

#include "pism/basalstrength/MohrCoulombYieldStress.hh"
#include "pism/basalstrength/MohrCoulombPointwise.hh"

#include "pism/util/Grid.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/File.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/Time.hh"
#include "pism/util/array/CellType.hh"
#include "pism/util/array/Forcing.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/coupler/util/options.hh" // ForcingOptions

namespace pism {

//! \file MohrCoulombYieldStress.cc  Process model which computes pseudo-plastic yield stress for the subglacial layer.

/*! \file MohrCoulombYieldStress.cc
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

The effective pressure is derived from the till (pore) water amount (the effective water
layer thickness). Then the effective pressure is combined with tillphi to compute an
updated `tauc` by the Mohr-Coulomb criterion.

This submodel is inactive in floating areas.
*/


/*!
The pseudo-plastic till basal resistance model is governed by this power law
equation,
    @f[ \tau_b = - \frac{\tau_c}{|\mathbf{U}|^{1-q} U_{\mathtt{th}}^q} \mathbf{U}, @f]
where @f$\tau_b=(\tau_{(b)x},\tau_{(b)y})@f$ is the basal shear stress and
@f$U=(u,v)@f$ is the sliding velocity.

We call the scalar field @f$\tau_c(t,x,y)@f$ the *yield stress* even when
the power @f$q@f$ is not zero; when that power is zero the formula describes
a plastic material with an actual yield stress.  The constant
@f$U_{\mathtt{th}}@f$ is the *threshold speed*, and @f$q@f$ is the *pseudo*
*plasticity exponent*.  The current class computes this yield stress field.
See also IceBasalResistancePlasticLaw::drag().

The strength of the saturated till material, the yield stress, is modeled by a
Mohr-Coulomb relation [\ref Paterson, \ref SchoofStream],
    @f[   \tau_c = c_0 + (\tan \varphi) N_{till}, @f]
where @f$N_{till}@f$ is the effective pressure of the glacier on the mineral
till.

The determination of the till friction angle @f$\varphi(x,y)@f$  is important.
It is assumed in this default model to be a time-independent factor which
describes the strength of the unsaturated "dry" (mineral) till material.  Thus
it is assumed to change more slowly than the till water pressure, and it follows
that it changes more slowly than the yield stress and the basal shear stress.

Option `-topg_to_phi` causes call to topg_to_phi() at the beginning of the run.
This determines the map of @f$\varphi(x,y)@f$.  If this option is note given,
the current method leaves `tillphi` unchanged, and thus either in its
read-in-from-file state or with a default constant value from the config file.
*/
MohrCoulombYieldStress::MohrCoulombYieldStress(std::shared_ptr<const Grid> grid)
  : YieldStress(grid),
  m_till_phi(m_grid, "tillphi") {

  m_name = "Mohr-Coulomb yield stress model";

  m_till_phi.metadata()
      .long_name("friction angle for till under grounded ice sheet")
      .units("degrees")
      .set_time_independent(true);

  // in this model; need not be time-independent in general

  {
    std::string hydrology_tillwat_max = "hydrology.tillwat_max";
    bool till_is_present = m_config->get_number(hydrology_tillwat_max) > 0.0;

    if (not till_is_present) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "The Mohr-Coulomb yield stress model cannot be used without till.\n"
                                    "Reset %s or choose a different yield stress model.",
                                    hydrology_tillwat_max.c_str());
    }
  }

  auto delta_file = m_config->get_string("basal_yield_stress.mohr_coulomb.delta.file");

  if (not delta_file.empty()) {
    ForcingOptions opt(*m_grid->ctx(), "basal_yield_stress.mohr_coulomb.delta");

    unsigned int buffer_size = m_config->get_number("input.forcing.buffer_size");

    File file(m_grid->com, opt.filename, io::PISM_NETCDF3, io::PISM_READONLY);

    m_delta = std::make_shared<array::Forcing>(m_grid,
                                          file,
                                          "mohr_coulomb_delta",
                                          "", // no standard name
                                          buffer_size,
                                          opt.periodic, LINEAR);
    m_delta->metadata()
        .long_name("minimum effective pressure on till as a fraction of overburden pressure")
        .units("1");
  }
}

void MohrCoulombYieldStress::restart_impl(const File &input_file, int record) {
  m_basal_yield_stress.read(input_file, record);
  m_till_phi.read(input_file, record);
}


//! Initialize the pseudo-plastic till mechanical model.
void MohrCoulombYieldStress::bootstrap_impl(const File &input_file,
                                            const YieldStressInputs &inputs) {

  auto tauc_to_phi_file = m_config->get_string("basal_yield_stress.mohr_coulomb.tauc_to_phi.file");

  if (not tauc_to_phi_file.empty()) {
    m_basal_yield_stress.regrid(tauc_to_phi_file, io::Default::Nil());

    m_log->message(2,
                   "  Will compute till friction angle (tillphi) as a function"
                   " of the yield stress (tauc)...\n");

    till_friction_angle(m_basal_yield_stress,
                        *inputs.till_water_thickness,
                        inputs.geometry->ice_thickness,
                        inputs.geometry->cell_type,
                        m_till_phi);

  } else if (m_config->get_flag("basal_yield_stress.mohr_coulomb.topg_to_phi.enabled")) {

    m_log->message(2, "  creating till friction angle map from bed elevation...\n");

    if (input_file.variable_exists(m_till_phi.metadata().get_name())) {
      // till_phi is present in the input file
      m_log->message(2,
                     "PISM WARNING: -topg_to_phi computation will override the '%s' field\n"
                     "              present in the input file '%s'!\n",
                     m_till_phi.metadata().get_name().c_str(), input_file.name().c_str());
    }

    till_friction_angle(inputs.geometry->bed_elevation, m_till_phi);

  } else {
    double till_phi_default = m_config->get_number("basal_yield_stress.mohr_coulomb.till_phi_default");
    m_till_phi.regrid(input_file, io::Default(till_phi_default));
  }

  finish_initialization(inputs);
}

void MohrCoulombYieldStress::init_impl(const YieldStressInputs &inputs) {
  double till_phi_default = m_config->get_number("basal_yield_stress.mohr_coulomb.till_phi_default");
  m_till_phi.set(till_phi_default);

  finish_initialization(inputs);
}

/*!
 * Finish initialization after bootstrapping or initializing using constants.
 */
void MohrCoulombYieldStress::finish_initialization(const YieldStressInputs &inputs) {
  // regrid if requested, regardless of how initialized
  regrid(name(), m_till_phi);

  if (m_delta) {
    ForcingOptions opt(*m_grid->ctx(), "basal_yield_stress.mohr_coulomb.delta");

    m_delta->init(opt.filename, opt.periodic);
  }

  // We use a short time step length because we can get away with it here, but we can
  // probably do better...
  this->update(inputs, time().current(), 1.0 /* one second time step */);
}

MaxTimestep MohrCoulombYieldStress::max_timestep_impl(double t) const {
  (void) t;

  if (m_delta) {
    auto dt = m_delta->max_timestep(t);

    if (dt.finite()) {
      return MaxTimestep(dt.value(), name());
    }
  }

  return MaxTimestep(name());
}

void MohrCoulombYieldStress::set_till_friction_angle(const array::Scalar &input) {
  m_till_phi.copy_from(input);
}

void MohrCoulombYieldStress::define_model_state_impl(const File &output) const {
  m_basal_yield_stress.define(output, io::PISM_DOUBLE);
  m_till_phi.define(output, io::PISM_DOUBLE);
}

void MohrCoulombYieldStress::write_model_state_impl(const File &output) const {
  m_basal_yield_stress.write(output);
  m_till_phi.write(output);
}

//! Update the till yield stress for use in the pseudo-plastic till basal stress
//! model.  See also IceBasalResistancePlasticLaw.
/*!
Updates yield stress  @f$ \tau_c @f$  based on modeled till water layer thickness
from a Hydrology object.  We implement the Mohr-Coulomb criterion allowing
a (typically small) till cohesion  @f$ c_0 @f$
and by expressing the coefficient as the tangent of a till friction angle
 @f$ \varphi @f$ :
    @f[   \tau_c = c_0 + (\tan \varphi) N_{till}. @f]
See [@ref Paterson] table 8.1 regarding values.

The effective pressure on the till is empirically-related
to the amount of water in the till.  We use this formula derived from
[@ref Tulaczyketal2000] and documented in [@ref BuelervanPeltDRAFT]:

@f[ N_{till} = \min\left\{P_o, N_0 \left(\frac{\delta P_o}{N_0}\right)^s 10^{(e_0/C_c) (1 - s)}\right\} @f]

where  @f$ s = W_{till} / W_{till}^{max} @f$,  @f$ W_{till}^{max} @f$ =`hydrology_tillwat_max`,
@f$ \delta @f$ =`basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden`,  @f$ P_o @f$  is the
overburden pressure,  @f$ N_0 @f$ =`basal_yield_stress.mohr_coulomb.till_reference_effective_pressure` is a
reference effective pressure,   @f$ e_0 @f$ =`basal_yield_stress.mohr_coulomb.till_reference_void_ratio` is the void ratio
at the reference effective pressure, and  @f$ C_c @f$ =`basal_yield_stress.mohr_coulomb.till_compressibility_coefficient`
is the coefficient of compressibility of the till.  Constants  @f$ N_0, e_0, C_c @f$  are
found by [@ref Tulaczyketal2000] from laboratory experiments on samples of
till.

If `basal_yield_stress.add_transportable_water` is yes then @f$ s @f$ in the above formula
becomes @f$ s = (W + W_{till}) / W_{till}^{max} @f$,
that is, the water amount is the sum @f$ W+W_{till} @f$.
 */
void MohrCoulombYieldStress::update_impl(const YieldStressInputs &inputs,
                                         double t, double dt) {
  (void) t;
  (void) dt;

  bool slippery_grounding_lines = m_config->get_flag("basal_yield_stress.slippery_grounding_lines"),
       add_transportable_water  = m_config->get_flag("basal_yield_stress.add_transportable_water");

  const double
    ice_density      = m_config->get_number("constants.ice.density"),
    standard_gravity = m_config->get_number("constants.standard_gravity");

  const double high_tauc  = m_config->get_number("basal_yield_stress.ice_free_bedrock"),
               W_till_max = m_config->get_number("hydrology.tillwat_max"),
               delta      = m_config->get_number("basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden"),
               tlftw      = m_config->get_number("basal_yield_stress.mohr_coulomb.till_log_factor_transportable_water");

  MohrCoulombPointwise mc(m_config);

  const array::Scalar
    &W_till        = *inputs.till_water_thickness,
    &W_subglacial  = *inputs.subglacial_water_thickness,
    &ice_thickness = inputs.geometry->ice_thickness;

  const auto &cell_type      = inputs.geometry->cell_type;
  const auto &bed_topography = inputs.geometry->bed_elevation;
  const auto &sea_level      = inputs.geometry->sea_level_elevation;

  array::AccessScope list{&W_till, &m_till_phi, &m_basal_yield_stress, &cell_type,
                               &bed_topography, &sea_level, &ice_thickness};

  if (add_transportable_water) {
    list.add(W_subglacial);
  }

  if (m_delta) {
    m_delta->update(t, dt);
    m_delta->average(t, dt);
    list.add(*m_delta);
  }

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (cell_type.ice_free(i, j)) {
      m_basal_yield_stress(i, j) = high_tauc;  // large yield stress if ice-free
    } else { // grounded and there is some ice

      // user can ask that marine grounding lines get special treatment
      double water = W_till(i,j); // usual case

      if (slippery_grounding_lines and
          bed_topography(i, j) <= sea_level(i, j) and
          (cell_type.next_to_floating_ice(i, j) or cell_type.next_to_ice_free_ocean(i, j))) {
        water = W_till_max;
      } else if (add_transportable_water) {
        water = W_till(i, j) + tlftw * log(1.0 + W_subglacial(i, j) / tlftw);
      }

      double P_overburden = ice_density * standard_gravity * ice_thickness(i, j);

      m_basal_yield_stress(i, j) = mc.yield_stress(m_delta ? (*m_delta)(i, j) : delta,
                                                   P_overburden, water, m_till_phi(i, j));
    }
  }

  m_basal_yield_stress.update_ghosts();
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
void MohrCoulombYieldStress::till_friction_angle(const array::Scalar &bed_topography,
                                                 array::Scalar &result) {

  const double
    phi_min  = m_config->get_number("basal_yield_stress.mohr_coulomb.topg_to_phi.phi_min"),
    phi_max  = m_config->get_number("basal_yield_stress.mohr_coulomb.topg_to_phi.phi_max"),
    topg_min = m_config->get_number("basal_yield_stress.mohr_coulomb.topg_to_phi.topg_min"),
    topg_max = m_config->get_number("basal_yield_stress.mohr_coulomb.topg_to_phi.topg_max");

  m_log->message(2,
                 "  till friction angle (phi) is piecewise-linear function of bed elev (topg):\n"
                 "            /  %5.2f                                 for   topg < %.f\n"
                 "      phi = |  %5.2f + (topg - (%.f)) * (%.2f / %.f)   for   %.f < topg < %.f\n"
                 "            \\  %5.2f                                 for   %.f < topg\n",
                 phi_min, topg_min,
                 phi_min, topg_min, phi_max - phi_min, topg_max - topg_min, topg_min, topg_max,
                 phi_max, topg_max);

  if (phi_min >= phi_max) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "invalid -topg_to_phi arguments: phi_min < phi_max is required");
  }

  if (topg_min >= topg_max) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "invalid -topg_to_phi arguments: topg_min < topg_max is required");
  }

  const double slope = (phi_max - phi_min) / (topg_max - topg_min);

  array::AccessScope list{&bed_topography, &result};

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();
    const double bed = bed_topography(i, j);

    if (bed <= topg_min) {
      result(i, j) = phi_min;
    } else if (bed >= topg_max) {
      result(i, j) = phi_max;
    } else {
      result(i, j) = phi_min + (bed - topg_min) * slope;
    }
  }

  // communicate ghosts so that the tauc computation can be performed locally
  // (including ghosts of tauc, that is)
  result.update_ghosts();
}

/*!
 * Compute the till friction angle in grounded areas using available basal yield stress,
 * till water thickness, and overburden pressure.
 *
 * This is the inverse of the formula used by `update_impl()`.
 */
void MohrCoulombYieldStress::till_friction_angle(const array::Scalar &basal_yield_stress,
                                                 const array::Scalar &till_water_thickness,
                                                 const array::Scalar &ice_thickness,
                                                 const array::CellType &cell_type,
                                                 array::Scalar &result) {

  MohrCoulombPointwise mc(m_config);

  double
    ice_density      = m_config->get_number("constants.ice.density"),
    standard_gravity = m_config->get_number("constants.standard_gravity");

  double
    delta = m_config->get_number("basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden");

  const array::Scalar
    &W_till = till_water_thickness;

  array::AccessScope list{&cell_type, &basal_yield_stress, &W_till, &ice_thickness, &result};

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (cell_type.ocean(i, j) or cell_type.ice_free(i, j)) {
      // no change
    } else { // grounded and there is some ice
      double P_overburden = ice_density * standard_gravity * ice_thickness(i, j);

      result(i, j) = mc.till_friction_angle(delta, P_overburden, W_till(i, j), basal_yield_stress(i, j));
    }
  }

  result.update_ghosts();
}

DiagnosticList MohrCoulombYieldStress::diagnostics_impl() const {
  return combine({{"tillphi", Diagnostic::wrap(m_till_phi)}},
                 YieldStress::diagnostics_impl());
}

} // end of namespace pism
