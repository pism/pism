// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 Constantine Khroulev and Ed Bueler
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

#include "StressBalance.hh"
#include "ShallowStressBalance.hh"
#include "SSB_Modifier.hh"
#include "pism/coupler/OceanModel.hh"
#include "pism/util/EnthalpyConverter.hh"
#include "pism/rheology/FlowLaw.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/Mask.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Vars.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/Profiling.hh"
#include "pism/util/IceModelVec2CellType.hh"
#include "pism/util/Time.hh"
#include "pism/geometry/Geometry.hh"

namespace pism {
namespace stressbalance {

Inputs::Inputs() {
  geometry = NULL;
  new_bed_elevation = true;

  basal_melt_rate       = NULL;
  melange_back_pressure = NULL;
  fracture_density      = NULL;
  basal_yield_stress    = NULL;

  enthalpy = NULL;
  age      = NULL;

  bc_mask   = NULL;
  bc_values = NULL;

  no_model_mask              = NULL;
  no_model_ice_thickness     = NULL;
  no_model_surface_elevation = NULL;
}

/*!
 * Save stress balance inputs to a file (for debugging).
 */
void Inputs::dump(const char *filename) const {
  if (not geometry) {
    return;
  }

  Context::ConstPtr ctx = geometry->ice_thickness.grid()->ctx();
  Config::ConstPtr config = ctx->config();

  PIO output(ctx->com(), config->get_string("output.format"), filename, PISM_READWRITE_MOVE);

  config->write(output);

  io::define_time(output, *ctx);
  io::append_time(output, config->get_string("time.dimension_name"), ctx->time()->current());

  {
    geometry->latitude.write(output);
    geometry->longitude.write(output);

    geometry->bed_elevation.write(output);
    geometry->sea_level_elevation.write(output);

    geometry->ice_thickness.write(output);
    geometry->ice_area_specific_volume.write(output);

    geometry->cell_type.write(output);
    geometry->cell_grounded_fraction.write(output);
    geometry->ice_surface_elevation.write(output);
  }

  if (basal_melt_rate) {
    basal_melt_rate->write(output);
  }

  if (melange_back_pressure) {
    melange_back_pressure->write(output);
  }

  if (fracture_density) {
    fracture_density->write(output);
  }

  if (basal_yield_stress) {
    basal_yield_stress->write(output);
  }

  if (enthalpy) {
    enthalpy->write(output);
  }

  if (age) {
    age->write(output);
  }

  if (bc_mask) {
    bc_mask->write(output);
  }

  if (bc_values) {
    bc_values->write(output);
  }

  if (no_model_mask) {
    no_model_mask->write(output);
  }

  if (no_model_ice_thickness) {
    no_model_ice_thickness->write(output);
  }

  if (no_model_surface_elevation) {
    no_model_surface_elevation->write(output);
  }
}

StressBalance::StressBalance(IceGrid::ConstPtr g,
                             ShallowStressBalance *sb,
                             SSB_Modifier *ssb_mod)
  : Component(g), m_shallow_stress_balance(sb), m_modifier(ssb_mod) {

  // allocate the vertical velocity field:
  m_w.create(m_grid, "wvel_rel", WITHOUT_GHOSTS);
  m_w.set_attrs("diagnostic",
                "vertical velocity of ice, relative to base of ice directly below",
                "m s-1", "");
  m_w.set_time_independent(false);
  m_w.metadata().set_string("glaciological_units", "m year-1");

  m_strain_heating.create(m_grid, "strain_heating", WITHOUT_GHOSTS);
  m_strain_heating.set_attrs("internal",
                             "rate of strain heating in ice (dissipation heating)",
                             "W m-3", "");
}

StressBalance::~StressBalance() {
  delete m_shallow_stress_balance;
  delete m_modifier;
}

//! \brief Initialize the StressBalance object.
void StressBalance::init() {
  m_shallow_stress_balance->init();
  m_modifier->init();
}

//! \brief Performs the shallow stress balance computation.
void StressBalance::update(const Inputs &inputs, bool full_update) {

  const Profiling &profiling = m_grid->ctx()->profiling();

  try {
    profiling.begin("stress_balance.shallow");
    m_shallow_stress_balance->update(inputs, full_update);
    profiling.end("stress_balance.shallow");

    profiling.begin("stress_balance.modifier");
    m_modifier->update(m_shallow_stress_balance->velocity(),
                       inputs, full_update);
    profiling.end("stress_balance.modifier");

    if (full_update) {
      const IceModelVec3 &u = m_modifier->velocity_u();
      const IceModelVec3 &v = m_modifier->velocity_v();

      profiling.begin("stress_balance.strain_heat");
      this->compute_volumetric_strain_heating(inputs);
      profiling.end("stress_balance.strain_heat");

      profiling.begin("stress_balance.vertical_velocity");
      this->compute_vertical_velocity(inputs.geometry->cell_type,
                                      u, v, inputs.basal_melt_rate, m_w);
      profiling.end("stress_balance.vertical_velocity");

      m_cfl_3d = ::pism::max_timestep_cfl_3d(inputs.geometry->ice_thickness,
                                             inputs.geometry->cell_type,
                                             m_modifier->velocity_u(),
                                             m_modifier->velocity_v(),
                                             m_w);
    }

    m_cfl_2d = ::pism::max_timestep_cfl_2d(inputs.geometry->ice_thickness,
                                           inputs.geometry->cell_type,
                                           m_shallow_stress_balance->velocity());
  }
  catch (RuntimeError &e) {
    e.add_context("updating the stress balance");
    throw;
  }
}

CFLData StressBalance::max_timestep_cfl_2d() const {
  return m_cfl_2d;
}

CFLData StressBalance::max_timestep_cfl_3d() const {
  return m_cfl_3d;
}

const IceModelVec2V& StressBalance::advective_velocity() const {
  return m_shallow_stress_balance->velocity();
}

const IceModelVec2Stag& StressBalance::diffusive_flux() const {
  return m_modifier->diffusive_flux();
}

double StressBalance::max_diffusivity() const {
  return m_modifier->max_diffusivity();
}

const IceModelVec3& StressBalance::velocity_u() const {
  return m_modifier->velocity_u();
}

const IceModelVec3& StressBalance::velocity_v() const {
  return m_modifier->velocity_v();
}

const IceModelVec3& StressBalance::velocity_w() const {
  return m_w;
}

const IceModelVec2S& StressBalance::basal_frictional_heating() const {
  return m_shallow_stress_balance->basal_frictional_heating();
}

const IceModelVec3& StressBalance::volumetric_strain_heating() const {
  return m_strain_heating;
}

void StressBalance::compute_2D_stresses(const IceModelVec2V &velocity,
                                        const IceModelVec2S &hardness,
                                        const IceModelVec2CellType &cell_type,
                                        IceModelVec2 &result) const {
  m_shallow_stress_balance->compute_2D_stresses(velocity, hardness, cell_type, result);
}

//! Compute vertical velocity using incompressibility of the ice.
/*!
The vertical velocity \f$w(x,y,z,t)\f$ is the velocity *relative to the
location of the base of the ice column*.  That is, the vertical velocity
computed here is identified as \f$\tilde w(x,y,s,t)\f$ in the page
[]@ref vertchange.

Thus \f$w<0\f$ here means that that
that part of the ice is getting closer to the base of the ice, and so on.
The slope of the bed (i.e. relative to the geoid) and/or the motion of the
bed (i.e. from bed deformation) do not affect the vertical velocity.

In fact the following statement is exactly true if the basal melt rate is zero:
the vertical velocity at a point in the ice is positive (negative) if and only
if the average horizontal divergence of the horizontal velocity, in the portion
of the ice column below that point, is negative (positive).
In particular, because \f$z=0\f$ is the location of the base of the ice
always, the only way to have \f$w(x,y,0,t) \ne 0\f$ is to have a basal melt
rate.

Incompressibility itself says
   \f[ \nabla\cdot\mathbf{U} + \frac{\partial w}{\partial z} = 0. \f]
This is immediately equivalent to the integral
   \f[ w(x,y,z,t) = - \int_{b(x,y,t)}^{z} \nabla\cdot\mathbf{U}\,d\zeta
                           + w_b(x,y,t). \f]
Here the value \f$w_b(x,y,t)\f$ is either zero or the negative of the basal melt rate
according to the value of the flag `geometry.update.use_basal_melt_rate`.

The vertical integral is computed by the trapezoid rule.
 */
void StressBalance::compute_vertical_velocity(const IceModelVec2CellType &mask,
                                              const IceModelVec3 &u,
                                              const IceModelVec3 &v,
                                              const IceModelVec2S *basal_melt_rate,
                                              IceModelVec3 &result) {

  const bool use_upstream_fd = m_config->get_string("stress_balance.vertical_velocity_approximation") == "upstream";

  IceModelVec::AccessList list{&u, &v, &mask, &result};

  if (basal_melt_rate) {
    list.add(*basal_melt_rate);
  }

  const std::vector<double> &z = m_grid->z();
  const unsigned int Mz = m_grid->Mz();

  const double
    dx = m_grid->dx(),
    dy = m_grid->dy();

  std::vector<double> u_x_plus_v_y(Mz);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double *w_ij = result.get_column(i,j);

    const double
      *u_w  = u.get_column(i-1,j),
      *u_ij = u.get_column(i,j),
      *u_e  = u.get_column(i+1,j);
    const double
      *v_s  = v.get_column(i,j-1),
      *v_ij = v.get_column(i,j),
      *v_n  = v.get_column(i,j+1);

    double
      west  = 1.0,
      east  = 1.0,
      south = 1.0,
      north = 1.0;
    double
      D_x = 0,                  // 1/(dx), 1/(2dx), or 0
      D_y = 0;                  // 1/(dy), 1/(2dy), or 0

    // Switch between second-order centered differences in the interior and
    // first-order one-sided differences at ice margins.

    // x-derivative
    {
      // use basal velocity to determine FD direction ("upwind" when it's clear, centered when it's
      // not)
      if (use_upstream_fd) {
        const double
          uw = 0.5 * (u_w[0] + u_ij[0]),
          ue = 0.5 * (u_ij[0] + u_e[0]);

        if (uw > 0.0 and ue >= 0.0) {
          west = 1.0;
          east = 0.0;
        } else if (uw <= 0.0 and ue < 0.0) {
          west = 0.0;
          east = 1.0;
        } else {
          west = 1.0;
          east = 1.0;
        }
      }

      if ((mask.icy(i,j) and mask.ice_free(i+1,j)) or (mask.ice_free(i,j) and mask.icy(i+1,j))) {
        east = 0;
      }
      if ((mask.icy(i,j) and mask.ice_free(i-1,j)) or (mask.ice_free(i,j) and mask.icy(i-1,j))) {
        west = 0;
      }

      if (east + west > 0) {
        D_x = 1.0 / (dx * (east + west));
      } else {
        D_x = 0.0;
      }
    }

    // y-derivative
    {
      // use basal velocity to determine FD direction ("upwind" when it's clear, centered when it's
      // not)
      if (use_upstream_fd) {
        const double
          vs = 0.5 * (v_s[0] + v_ij[0]),
          vn = 0.5 * (v_ij[0] + v_n[0]);

        if (vs > 0.0 and vn >= 0.0) {
          south = 1.0;
          north = 0.0;
        } else if (vs <= 0.0 and vn < 0.0) {
          south = 0.0;
          north = 1.0;
        } else {
          south = 1.0;
          north = 1.0;
        }
      }

      if ((mask.icy(i,j) and mask.ice_free(i,j+1)) or (mask.ice_free(i,j) and mask.icy(i,j+1))) {
        north = 0;
      }
      if ((mask.icy(i,j) and mask.ice_free(i,j-1)) or (mask.ice_free(i,j) and mask.icy(i,j-1))) {
        south = 0;
      }

      if (north + south > 0) {
        D_y = 1.0 / (dy * (north + south));
      } else {
        D_y = 0.0;
      }
    }

    // compute u_x + v_y using a vectorizable loop
    for (unsigned int k = 0; k < Mz; ++k) {
      double
        u_x = D_x * (west  * (u_ij[k] - u_w[k]) + east  * (u_e[k] - u_ij[k])),
        v_y = D_y * (south * (v_ij[k] - v_s[k]) + north * (v_n[k] - v_ij[k]));
      u_x_plus_v_y[k] = u_x + v_y;
    }

    // at the base: include the basal melt rate
    if (basal_melt_rate != NULL) {
      w_ij[0] = - (*basal_melt_rate)(i,j);
    } else {
      w_ij[0] = 0.0;
    }

    // within the ice and above:
    for (unsigned int k = 1; k < Mz; ++k) {
      const double dz = z[k] - z[k-1];

      w_ij[k] = w_ij[k - 1] - (0.5 * dz) * (u_x_plus_v_y[k] + u_x_plus_v_y[k - 1]);
    }
  }
}

/**
 * This function computes \f$D^2\f$ defined by
 *
 * \f[ 2D^2 = D_{ij} D_{ij}\f]
 * or
 * \f[
 * D^2 = \frac{1}{2}\,\left(\frac{1}{2}\,(v_{z})^2 + (v_{y} + u_{x})^2 +
 *       (v_{y})^2 + \frac{1}{2}\,(v_{x} + u_{y})^2 + \frac{1}{2}\,(u_{z})^2 +
 *       (u_{x})^2\right)
 * \f]
 *
 * (note the use of the summation convention). Here \f$D_{ij}\f$ is the
 * strain rate tensor. See
 * StressBalance::compute_volumetric_strain_heating() for details.
 *
 * @param u_x,u_y,u_z partial derivatives of \f$u\f$, the x-component of the ice velocity
 * @param v_x,v_y,v_z partial derivatives of \f$v\f$, the y-component of the ice velocity
 *
 * @return \f$D^2\f$, where \f$D\f$ is defined above.
 */
static inline double D2(double u_x, double u_y, double u_z, double v_x, double v_y, double v_z) {
  return 0.5 * (PetscSqr(u_x + v_y) + u_x*u_x + v_y*v_y + 0.5 * (PetscSqr(u_y + v_x) + u_z*u_z + v_z*v_z));
}

/**
  \brief Computes the volumetric strain heating using horizontal
  velocity.

  Following the notation used in [\ref BBssasliding], let \f$u\f$ be a
  three-dimensional *vector* velocity field. Then the strain rate
  tensor \f$D_{ij}\f$ is defined by

  \f[ D_{ij} = \frac 12 \left(\diff{u_{i}}{x_{j}} + \diff{u_{j}}{x_{i}} \right), \f]

  Where \f$i\f$ and \f$j\f$ range from \f$1\f$ to \f$3\f$.

  The flow law in the viscosity form states

  \f[ \tau_{ij} = 2 \eta D_{ij}, \f]

  and the nonlinear ice viscosity satisfies

  \f[ 2 \eta = B(T) D^{(1/n) - 1}. \f]

  Here \f$D^{2}\f$ is defined by \f$2D^{2} = D_{ij}D_{ij}\f$ (using the
  summation convention) and \f$B(T) = A(T)^{-1/n}\f$ is the ice hardness.

  Now the volumetric strain heating is

  \f[ \Sigma = \sum_{i,j=1}^{3}D_{ij}\tau_{ij} = 2 B(T) D^{(1/n) + 1}. \f]

  We use an *approximation* of \f$D_{ij}\f$ common in shallow ice models:

  - we assume that horizontal derivatives of the vertical velocity are
    much smaller than \f$z\f$ derivatives horizontal velocity
    components \f$u\f$ and \f$v\f$. (We drop \f$w_x\f$ and \f$w_y\f$
    terms in \f$D_{ij}\f$.)

  - we use the incompressibility of ice to approximate \f$w_z\f$:

  \f[ w_z = - (u_x + v_y). \f]

  Requires ghosts of `u` and `v` velocity components and uses the fact
  that `u` and `v` above the ice are filled using constant
  extrapolation.

  Resulting field does not have ghosts.

  Below is the *Maxima* code that produces the expression evaluated by D2().

       derivabbrev : true;
       U : [u, v, w]; X : [x, y, z]; depends(U, X);
       gradef(w, x, 0); gradef(w, y, 0);
       gradef(w, z, -(diff(u, x) + diff(v, y)));
       d[i,j] := 1/2 * (diff(U[i], X[j]) + diff(U[j], X[i]));
       D : genmatrix(d, 3, 3), ratsimp, factor;
       tex('D = D);
       tex('D^2 = 1/2 * mat_trace(D . D));

  @return 0 on success
 */
void StressBalance::compute_volumetric_strain_heating(const Inputs &inputs) {
  PetscErrorCode ierr;

  const rheology::FlowLaw &flow_law = *m_shallow_stress_balance->flow_law();
  EnthalpyConverter::Ptr EC = m_shallow_stress_balance->enthalpy_converter();

  const IceModelVec3
    &u = m_modifier->velocity_u(),
    &v = m_modifier->velocity_v();

  const IceModelVec2S &thickness = inputs.geometry->ice_thickness;
  const IceModelVec3  *enthalpy  = inputs.enthalpy;

  const IceModelVec2CellType &mask = inputs.geometry->cell_type;

  double
    enhancement_factor = flow_law.enhancement_factor(),
    n = flow_law.exponent(),
    exponent = 0.5 * (1.0 / n + 1.0),
    e_to_a_power = pow(enhancement_factor,-1.0/n);

  IceModelVec::AccessList list{&mask, enthalpy, &m_strain_heating, &thickness, &u, &v};

  const std::vector<double> &z = m_grid->z();
  const unsigned int Mz = m_grid->Mz();
  std::vector<double> depth(Mz), pressure(Mz), hardness(Mz);

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double H = thickness(i, j);
      int ks = m_grid->kBelowHeight(H);
      const double
        *u_ij, *u_w, *u_n, *u_e, *u_s,
        *v_ij, *v_w, *v_n, *v_e, *v_s;
      double *Sigma;
      const double *E_ij;

      double west = 1, east = 1, south = 1, north = 1,
        D_x = 0,                // 1/(dx), 1/(2dx), or 0
        D_y = 0;                // 1/(dy), 1/(2dy), or 0

      // x-derivative
      {
        if ((mask.icy(i,j) and mask.ice_free(i+1,j)) or (mask.ice_free(i,j) and mask.icy(i+1,j))) {
          east = 0;
        }
        if ((mask.icy(i,j) and mask.ice_free(i-1,j)) or (mask.ice_free(i,j) and mask.icy(i-1,j))) {
          west = 0;
        }

        if (east + west > 0) {
          D_x = 1.0 / (m_grid->dx() * (east + west));
        } else {
          D_x = 0.0;
        }
      }

      // y-derivative
      {
        if ((mask.icy(i,j) and mask.ice_free(i,j+1)) or (mask.ice_free(i,j) and mask.icy(i,j+1))) {
          north = 0;
        }
        if ((mask.icy(i,j) and mask.ice_free(i,j-1)) or (mask.ice_free(i,j) and mask.icy(i,j-1))) {
          south = 0;
        }

        if (north + south > 0) {
          D_y = 1.0 / (m_grid->dy() * (north + south));
        } else {
          D_y = 0.0;
        }
      }

      u_ij = u.get_column(i,     j);
      u_w  = u.get_column(i - 1, j);
      u_e  = u.get_column(i + 1, j);
      u_s  = u.get_column(i,     j - 1);
      u_n  = u.get_column(i,     j + 1);

      v_ij = v.get_column(i,     j);
      v_w  = v.get_column(i - 1, j);
      v_e  = v.get_column(i + 1, j);
      v_s  = v.get_column(i,     j - 1);
      v_n  = v.get_column(i,     j + 1);

      E_ij = enthalpy->get_column(i, j);
      Sigma = m_strain_heating.get_column(i, j);

      for (int k = 0; k <= ks; ++k) {
        depth[k] = H - z[k];
      }

      // pressure added by the ice (i.e. pressure difference between the
      // current level and the top of the column)
      EC->pressure(depth, ks, pressure); // FIXME issue #15

      flow_law.hardness_n(E_ij, &pressure[0], ks + 1, &hardness[0]);

      for (int k = 0; k <= ks; ++k) {
        double dz;

        double u_z = 0.0, v_z = 0.0,
          u_x = D_x * (west  * (u_ij[k] - u_w[k]) + east  * (u_e[k] - u_ij[k])),
          u_y = D_y * (south * (u_ij[k] - u_s[k]) + north * (u_n[k] - u_ij[k])),
          v_x = D_x * (west  * (v_ij[k] - v_w[k]) + east  * (v_e[k] - v_ij[k])),
          v_y = D_y * (south * (v_ij[k] - v_s[k]) + north * (v_n[k] - v_ij[k]));

        if (k > 0) {
          dz = z[k+1] - z[k-1];
          u_z = (u_ij[k+1] - u_ij[k-1]) / dz;
          v_z = (v_ij[k+1] - v_ij[k-1]) / dz;
        } else {
          // use one-sided differences for u_z and v_z on the bottom level
          dz = z[1] - z[0];
          u_z = (u_ij[1] - u_ij[0]) / dz;
          v_z = (v_ij[1] - v_ij[0]) / dz;
        }

        Sigma[k] = 2.0 * e_to_a_power * hardness[k] * pow(D2(u_x, u_y, u_z, v_x, v_y, v_z), exponent);
      } // k-loop

      int remaining_levels = Mz - (ks + 1);
      if (remaining_levels > 0) {
        ierr = PetscMemzero(&Sigma[ks+1],
                            remaining_levels*sizeof(double));
        PISM_CHK(ierr, "PetscMemzero");
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

std::string StressBalance::stdout_report() const {
  return m_shallow_stress_balance->stdout_report() + m_modifier->stdout_report();
}

const ShallowStressBalance* StressBalance::shallow() const {
  return m_shallow_stress_balance;
}

const SSB_Modifier* StressBalance::modifier() const {
  return m_modifier;
}


void StressBalance::define_model_state_impl(const PIO &output) const {
  m_shallow_stress_balance->define_model_state(output);
  m_modifier->define_model_state(output);
}

void StressBalance::write_model_state_impl(const PIO &output) const {
  m_shallow_stress_balance->write_model_state(output);
  m_modifier->write_model_state(output);
}

//! \brief Compute eigenvalues of the horizontal, vertically-integrated strain rate tensor.
/*!
Calculates all components \f$D_{xx}, D_{yy}, D_{xy}=D_{yx}\f$ of the
vertically-averaged strain rate tensor \f$D\f$ [\ref SchoofStream].  Then computes
the eigenvalues `result(i,j,0)` = (maximum eigenvalue), `result(i,j,1)` = (minimum
eigenvalue).  Uses the provided thickness to make decisions (PIK) about computing
strain rates near calving front.

Note that `result(i,j,0)` >= `result(i,j,1)`, but there is no necessary relation between
the magnitudes, and either principal strain rate could be negative or positive.

Result can be used in a calving law, for example in eigencalving (PIK).

Note: strain rates will be derived from SSA velocities, using ghosts when
necessary. Both implementations (SSAFD and SSAFEM) call
update_ghosts() to ensure that ghost values are up to date.
 */
void compute_2D_principal_strain_rates(const IceModelVec2V &V,
                                       const IceModelVec2CellType &mask,
                                       IceModelVec2 &result) {

  using mask::ice_free;

  IceGrid::ConstPtr grid = result.grid();
  double    dx = grid->dx(), dy = grid->dy();

  if (result.ndof() != 2) {
    throw RuntimeError(PISM_ERROR_LOCATION, "result.dof() == 2 is required");
  }

  IceModelVec::AccessList list{&V, &mask, &result};

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.ice_free(i,j)) {
      result(i,j,0) = 0.0;
      result(i,j,1) = 0.0;
      continue;
    }

    StarStencil<int> m = mask.int_star(i,j);
    StarStencil<Vector2> U = V.star(i,j);

    // strain in units s-1
    double u_x = 0, u_y = 0, v_x = 0, v_y = 0,
      east = 1, west = 1, south = 1, north = 1;

    // Computes u_x using second-order centered finite differences written as
    // weighted sums of first-order one-sided finite differences.
    //
    // Given the cell layout
    // *----n----*
    // |         |
    // |         |
    // w         e
    // |         |
    // |         |
    // *----s----*
    // east == 0 if the east neighbor of the current cell is ice-free. In
    // this case we use the left- (west-) sided difference.
    //
    // If both neighbors in the east-west (x) direction are ice-free the
    // x-derivative is set to zero (see u_x, v_x initialization above).
    //
    // Similarly in other directions.
    if (ice_free(m.e)) {
      east = 0;
    }
    if (ice_free(m.w)) {
      west = 0;
    }
    if (ice_free(m.n)) {
      north = 0;
    }
    if (ice_free(m.s)) {
      south = 0;
    }

    if (west + east > 0) {
      u_x = 1.0 / (dx * (west + east)) * (west * (U.ij.u - U[West].u) + east * (U[East].u - U.ij.u));
      v_x = 1.0 / (dx * (west + east)) * (west * (U.ij.v - U[West].v) + east * (U[East].v - U.ij.v));
    }

    if (south + north > 0) {
      u_y = 1.0 / (dy * (south + north)) * (south * (U.ij.u - U[South].u) + north * (U[North].u - U.ij.u));
      v_y = 1.0 / (dy * (south + north)) * (south * (U.ij.v - U[South].v) + north * (U[North].v - U.ij.v));
    }

    const double A = 0.5 * (u_x + v_y),  // A = (1/2) trace(D)
      B   = 0.5 * (u_x - v_y),
      Dxy = 0.5 * (v_x + u_y),  // B^2 = A^2 - u_x v_y
      q   = sqrt(PetscSqr(B) + PetscSqr(Dxy));
    result(i,j,0) = A + q;
    result(i,j,1) = A - q; // q >= 0 so e1 >= e2

  }
}

} // end of namespace stressbalance
} // end of namespace pism
