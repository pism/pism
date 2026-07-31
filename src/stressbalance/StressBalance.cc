// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026 Constantine Khroulev and Ed Bueler
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
#include <memory>
#include <vector>
#include <petsc.h>

#include "pism/stressbalance/StressBalance.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/rheology/FlowLaw.hh"
#include "pism/util/Config.hh"
#include "pism/util/Context.hh"
#include "pism/util/EnthalpyConverter.hh"
#include "pism/util/Grid.hh"
#include "pism/util/Mask.hh"
#include "pism/util/Time.hh"
#include "pism/util/array/CellType.hh"
#include "pism/util/array/Vector.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/SynchronousOutputWriter.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/basalstrength/basal_resistance.hh"

namespace pism {
namespace stressbalance {

Inputs::Inputs() {
  geometry          = nullptr;
  new_bed_elevation = true;

  averaged_hardness     = nullptr;
  water_column_pressure = nullptr;
  fracture_density      = nullptr;
  basal_yield_stress    = nullptr;

  enthalpy = nullptr;
  age      = nullptr;

  bc_mask   = nullptr;
  bc_values = nullptr;

  no_model_mask              = nullptr;
  no_model_ice_thickness     = nullptr;
  no_model_surface_elevation = nullptr;
}

/*!
 * Save stress balance inputs to a file (for debugging).
 */
void Inputs::dump(const char *filename) const {
  if (geometry == nullptr) {
    return;
  }

  auto grid   = geometry->ice_thickness.grid();
  auto ctx    = grid->ctx();
  auto config = ctx->config();

  VariableMetadata mapping{ "mapping", ctx->unit_system() };
  auto writer = std::make_shared<SynchronousOutputWriter>(ctx->com(), *config);
  writer->initialize({}, true);

  OutputFile output(writer, filename);

  io::write_config(*config, "pism_config", output);

  auto time = ctx->time();
  output.define_variable(time->metadata());
  output.append_time(time->current());

  std::vector<const array::Array *> geom = { &geometry->bed_elevation,
                                             &geometry->sea_level_elevation,
                                             &geometry->ice_thickness,
                                             &geometry->ice_area_specific_volume,
                                             &geometry->cell_type,
                                             &geometry->cell_grounded_fraction,
                                             &geometry->ice_surface_elevation };
  if (grid->has_longitude_latitude()) {
    geom.push_back(&grid->longitude());
    geom.push_back(&grid->latitude());
  }

  // define
  for (const auto * vec : geom) {
    for (const auto &var : vec->all_metadata()) {
      output.define_variable(var);
    }
  }
  // write
  for (const auto * vec : geom) {
    vec->write(output);
  }

  const array::Array *optional[] = { water_column_pressure,
                                     fracture_density,
                                     basal_yield_stress,
                                     enthalpy,
                                     age,
                                     bc_mask,
                                     bc_values,
                                     no_model_mask,
                                     no_model_ice_thickness,
                                     no_model_surface_elevation };
  // define
  for (const auto * vec : optional) {
    if (vec != nullptr) {
      for (const auto &var : vec->all_metadata()) {
        output.define_variable(var);
      }
    }
  }
  // write
  for (const auto *vec : optional) {
    if (vec != nullptr) {
      vec->write(output);
    }
  }
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
class VerticalVelocity {
public:
  VerticalVelocity(bool upstream_fd, double dx, double dy, const std::vector<double> &z)
    : m_upstream_fd(upstream_fd), m_Mz(z.size()), m_dx(dx), m_dy(dy), m_dz(m_Mz),
      m_work(m_Mz) {

    m_dz[0] = 0.0;
    for (unsigned int k = 1; k < m_Mz; ++k) {
      m_dz[k] = z[k] - z[k - 1];
    }
  }

  void compute(const stencils::Star<int> &cell_type, double basal_melt_rate, const double *u_w,
               const double *u_c, const double *u_e, const double *v_s, const double *v_c,
               const double *v_n, double *result) const {

    using mask::ice_free;
    using mask::icy;

    double west = 1.0, east = 1.0, south = 1.0, north = 1.0;
    double D_x = 0, // 1/(dx), 1/(2dx), or 0
        D_y    = 0; // 1/(dy), 1/(2dy), or 0

    // Switch between second-order centered differences in the interior and
    // first-order one-sided differences at ice margins.

    // x-derivative
    {
      // use basal velocity to determine FD direction ("upwind" when it's clear, centered when it's
      // not)
      if (m_upstream_fd) {
        const double uw = 0.5 * (u_w[0] + u_c[0]), ue = 0.5 * (u_c[0] + u_e[0]);

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

      if ((icy(cell_type.c) and ice_free(cell_type.e)) or
          (ice_free(cell_type.c) and icy(cell_type.e))) {
        east = 0;
      }
      if ((icy(cell_type.c) and ice_free(cell_type.w)) or
          (ice_free(cell_type.c) and icy(cell_type.w))) {
        west = 0;
      }

      if (east + west > 0) {
        D_x = 1.0 / (m_dx * (east + west));
      } else {
        D_x = 0.0;
      }
    }

    // y-derivative
    {
      // use basal velocity to determine FD direction ("upwind" when it's clear, centered when it's
      // not)
      if (m_upstream_fd) {
        const double vs = 0.5 * (v_s[0] + v_c[0]), vn = 0.5 * (v_c[0] + v_n[0]);

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

      if ((icy(cell_type.c) and ice_free(cell_type.n)) or
          (ice_free(cell_type.c) and icy(cell_type.n))) {
        north = 0;
      }
      if ((icy(cell_type.c) and ice_free(cell_type.s)) or
          (ice_free(cell_type.c) and icy(cell_type.s))) {
        south = 0;
      }

      if (north + south > 0) {
        D_y = 1.0 / (m_dy * (north + south));
      } else {
        D_y = 0.0;
      }
    }

    auto &u_x_plus_v_y = m_work;
    // compute u_x + v_y using a vectorizable loop
    for (unsigned int k = 0; k < m_Mz; ++k) {
      double u_x      = D_x * (west * (u_c[k] - u_w[k]) + east * (u_e[k] - u_c[k])),
             v_y      = D_y * (south * (v_c[k] - v_s[k]) + north * (v_n[k] - v_c[k]));
      u_x_plus_v_y[k] = u_x + v_y;
    }

    // at the base: include the basal melt rate
    result[0] = -basal_melt_rate;
    // within the ice and above:
    for (unsigned int k = 1; k < m_Mz; ++k) {
      result[k] = result[k - 1] - (0.5 * m_dz[k]) * (u_x_plus_v_y[k] + u_x_plus_v_y[k - 1]);
    }
  }

private:
  bool m_upstream_fd;
  unsigned int m_Mz;
  double m_dx;
  double m_dy;
  std::vector<double> m_dz;
  mutable std::vector<double> m_work;
};

void compute_vertical_velocity(const array::CellType1 &cell_type, const array::Array3D &u,
                               const array::Array3D &v, const array::Scalar *basal_melt_rate,
                               array::Array3D &result) {
  auto grid = u.grid();

  auto config = grid->ctx()->config();

  const bool use_upstream_fd =
      config->get_string("stress_balance.vertical_velocity_approximation") == "upstream";

  array::AccessScope list{&u, &v, &cell_type, &result};

  if (basal_melt_rate != nullptr) {
    list.add(*basal_melt_rate);
  }

  VerticalVelocity w(use_upstream_fd, grid->dx(), grid->dy(), u.levels());

  for (auto p : grid->points()) {
    const int i = p.i(), j = p.j();

    const double
      *u_w = u.get_column(i - 1, j),
      *u_c = u.get_column(i, j),
      *u_e = u.get_column(i + 1, j);
    const double
      *v_s = v.get_column(i, j - 1),
      *v_c = v.get_column(i, j),
      *v_n = v.get_column(i, j + 1);

    w.compute(cell_type.star_int(i, j),                                    //
              basal_melt_rate != nullptr ? (*basal_melt_rate)(i, j) : 0.0, //
              u_w, u_c, u_e,                                               //
              v_s, v_c, v_n,                                               //
              result.get_column(i, j));
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
 * compute_volumetric_strain_heating() for details.
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
 */
void compute_volumetric_strain_heating(const array::Array3D &u, const array::Array3D &v,
                                       const array::Array3D &enthalpy,
                                       const array::Scalar &thickness,
                                       const array::CellType1 &cell_type,
                                       const rheology::FlowLaw &flow_law,
                                       const EnthalpyConverter &EC,
                                       double enhancement_factor,
                                       array::Array3D &result) {

  using mask::icy;
  using mask::ice_free;

  assert(u.stencil_width() > 0);
  assert(v.stencil_width() > 0);

  auto grid = thickness.grid();

  double
    n = flow_law.exponent(),
    exponent = 0.5 * (1.0 / n + 1.0),
    e_to_a_power = pow(enhancement_factor, -1.0 / n);

  array::AccessScope list{&cell_type, &enthalpy, &result, &thickness, &u, &v};

  const std::vector<double> &z = grid->z();
  const unsigned int Mz = grid->Mz();
  std::vector<double> depth(Mz), pressure(Mz), hardness(Mz);

  ParallelSection loop(grid->com);
  try {
    for (auto p : grid->points()) {
      const int i = p.i(), j = p.j();

      auto M = cell_type.star_int(i, j);

      double H = thickness(i, j);
      int ks = grid->kBelowHeight(H);
      const double
        *u_c, *u_w, *u_n, *u_e, *u_s,
        *v_c, *v_w, *v_n, *v_e, *v_s;
      double *Sigma;
      const double *E_c;

      double west = 1, east = 1, south = 1, north = 1,
        D_x = 0,                // 1/(dx), 1/(2dx), or 0
        D_y = 0;                // 1/(dy), 1/(2dy), or 0

      // x-derivative
      {
        if ((icy(M.c) and ice_free(M.e)) or (ice_free(M.c) and icy(M.e))) {
          east = 0;
        }
        if ((icy(M.c) and ice_free(M.w)) or (ice_free(M.c) and icy(M.w))) {
          west = 0;
        }

        if (east + west > 0) {
          D_x = 1.0 / (grid->dx() * (east + west));
        } else {
          D_x = 0.0;
        }
      }

      // y-derivative
      {
        if ((icy(M.c) and ice_free(M.n)) or (ice_free(M.c) and icy(M.n))) {
          north = 0;
        }
        if ((icy(M.c) and ice_free(M.s)) or (ice_free(M.c) and icy(M.s))) {
          south = 0;
        }

        if (north + south > 0) {
          D_y = 1.0 / (grid->dy() * (north + south));
        } else {
          D_y = 0.0;
        }
      }

      u_c = u.get_column(i, j);
      u_w = u.get_column(i - 1, j);
      u_e = u.get_column(i + 1, j);
      u_s = u.get_column(i, j - 1);
      u_n = u.get_column(i, j + 1);

      v_c = v.get_column(i, j);
      v_w = v.get_column(i - 1, j);
      v_e = v.get_column(i + 1, j);
      v_s = v.get_column(i, j - 1);
      v_n = v.get_column(i, j + 1);

      E_c = enthalpy.get_column(i, j);
      Sigma = result.get_column(i, j);

      for (int k = 0; k <= ks; ++k) {
        depth[k] = H - z[k];
      }

      // pressure added by the ice (i.e. pressure difference between the
      // current level and the top of the column)
      EC.pressure(depth, ks, pressure); // FIXME issue #15

      flow_law.hardness_n(E_c, pressure.data(), ks + 1, hardness.data());

      for (int k = 0; k <= ks; ++k) {
        double dz;

        double u_z = 0.0, v_z = 0.0,
          u_x = D_x * (west  * (u_c[k] - u_w[k]) + east  * (u_e[k] - u_c[k])),
          u_y = D_y * (south * (u_c[k] - u_s[k]) + north * (u_n[k] - u_c[k])),
          v_x = D_x * (west  * (v_c[k] - v_w[k]) + east  * (v_e[k] - v_c[k])),
          v_y = D_y * (south * (v_c[k] - v_s[k]) + north * (v_n[k] - v_c[k]));

        if (k > 0) {
          dz = z[k+1] - z[k-1];
          u_z = (u_c[k+1] - u_c[k-1]) / dz;
          v_z = (v_c[k+1] - v_c[k-1]) / dz;
        } else {
          // use one-sided differences for u_z and v_z on the bottom level
          dz = z[1] - z[0];
          u_z = (u_c[1] - u_c[0]) / dz;
          v_z = (v_c[1] - v_c[0]) / dz;
        }

        Sigma[k] = 2.0 * e_to_a_power * hardness[k] * pow(D2(u_x, u_y, u_z, v_x, v_y, v_z), exponent);
      } // k-loop

      int remaining_levels = Mz - (ks + 1);
      if (remaining_levels > 0) {
        PetscErrorCode ierr = PetscArrayzero(&Sigma[ks+1], remaining_levels);
        PISM_CHK(ierr, "PetscArrayzero");
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
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
void compute_2D_principal_strain_rates(const array::Vector1 &V,
                                       const array::CellType1 &mask,
                                       array::Array2D<PrincipalStrainRates> &result) {

  using mask::ice_free;

  auto grid = result.grid();
  double dx = grid->dx();
  double dy = grid->dy();

  array::AccessScope list{&V, &mask, &result};

  for (auto p : grid->points()) {
    const int i = p.i(), j = p.j();

    if (mask.ice_free(i,j)) {
      result(i, j).eigen1 = 0.0;
      result(i, j).eigen2 = 0.0;
      continue;
    }

    auto m = mask.star(i,j);
    auto U = V.star(i,j);

    // strain in units s^-1
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
      u_x = 1.0 / (dx * (west + east)) * (west * (U.c.u - U[West].u) + east * (U[East].u - U.c.u));
      v_x = 1.0 / (dx * (west + east)) * (west * (U.c.v - U[West].v) + east * (U[East].v - U.c.v));
    }

    if (south + north > 0) {
      u_y = 1.0 / (dy * (south + north)) * (south * (U.c.u - U[South].u) + north * (U[North].u - U.c.u));
      v_y = 1.0 / (dy * (south + north)) * (south * (U.c.v - U[South].v) + north * (U[North].v - U.c.v));
    }

    const double A = 0.5 * (u_x + v_y),  // A = (1/2) trace(D)
      B   = 0.5 * (u_x - v_y),
      Dxy = 0.5 * (v_x + u_y),  // B^2 = A^2 - u_x v_y
      q   = sqrt(B*B + Dxy*Dxy);
    result(i, j).eigen1 = A + q;
    result(i, j).eigen2 = A - q; // q >= 0 so e1 >= e2

  }
}

//! @brief Compute 2D deviatoric stresses.
void compute_2D_stresses(const rheology::FlowLaw &flow_law,
                         const array::Vector1 &velocity,
                         const array::Scalar &hardness,
                         const array::CellType1 &cell_type,
                         array::Array2D<DeviatoricStresses> &result) {

  using mask::ice_free;

  auto grid = result.grid();

  const double
    dx = grid->dx(),
    dy = grid->dy();

  array::AccessScope list{&velocity, &hardness, &result, &cell_type};

  for (auto p : grid->points()) {
    const int i = p.i(), j = p.j();

    if (cell_type.ice_free(i, j)) {
      result(i,j).xx = 0.0;
      result(i,j).yy = 0.0;
      result(i,j).xy = 0.0;
      continue;
    }

    auto m = cell_type.star(i,j);
    auto U = velocity.star(i,j);

    // strain in units s^-1
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
    // Similarly in y-direction.
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
      u_x = 1.0 / (dx * (west + east)) * (west * (U.c.u - U[West].u) + east * (U[East].u - U.c.u));
      v_x = 1.0 / (dx * (west + east)) * (west * (U.c.v - U[West].v) + east * (U[East].v - U.c.v));
    }

    if (south + north > 0) {
      u_y = 1.0 / (dy * (south + north)) * (south * (U.c.u - U[South].u) + north * (U[North].u - U.c.u));
      v_y = 1.0 / (dy * (south + north)) * (south * (U.c.v - U[South].v) + north * (U[North].v - U.c.v));
    }

    double nu = 0.0;
    flow_law.effective_viscosity(hardness(i, j),
                                 secondInvariant_2D({u_x, v_x}, {u_y, v_y}),
                                 &nu, NULL);

    //get deviatoric stresses
    result(i,j).xx = 2.0*nu*u_x;
    result(i,j).yy = 2.0*nu*v_y;
    result(i,j).xy = nu*(u_y+v_x);
  }
}

//! \brief Compute the basal frictional heating.
/*!
  Ice shelves have zero basal friction heating.

  \param[in] V *basal* sliding velocity
  \param[in] tauc basal yield stress
  \param[in] mask (used to determine if floating or grounded)
  \param[out] result
 */
void compute_basal_frictional_heating(const IceBasalResistancePlasticLaw & sliding_law,
                                      const array::Vector &basal_velocity,
                                      const array::Scalar &tauc,
                                      const array::CellType &mask,
                                      array::Scalar &result) {

  auto grid = basal_velocity.grid();

  array::AccessScope list{ &basal_velocity, &result, &tauc, &mask };

  for (auto p : grid->points()) {
    const int i = p.i(), j = p.j();

    if (mask.ocean(i, j)) {
      result(i, j) = 0.0;
    } else {
      const auto &V  = basal_velocity(i, j);

      double C = sliding_law.drag(tauc(i, j), V.u, V.v);
      double basal_stress_x = -C * V.u, basal_stress_y = -C * V.v;

      result(i, j) = -basal_stress_x * V.u - basal_stress_y * V.v;
    }
  }
}

} // end of namespace stressbalance
} // end of namespace pism
