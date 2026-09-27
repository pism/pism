// Copyright (C) 2026 Andy Aschwanden and Constantine Khroulev
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

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "pism/coupler/debris/EnglacialTransport.hh"

#include "pism/coupler/debris/column_kernels.hh"
#include "pism/geometry/TransportScheme.hh"
#include "pism/util/Grid.hh"
#include "pism/util/array/CellType.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {
namespace debris {

EnglacialTransport::EnglacialTransport(std::shared_ptr<const Grid> grid,
                                       std::shared_ptr<TransportScheme3D> scheme,
                                       double solid_density, double cfl_ratio)
  : m_grid(grid),
    m_scheme(scheme),
    m_zi(column::interfaces(grid->z())),
    m_solid_density(solid_density),
    m_cfl_ratio(cfl_ratio),
    m_substeps(0),
    m_dt_vertical(std::numeric_limits<double>::infinity()),
    m_mass(std::make_shared<array::Array3D>(grid, "englacial_debris_mass", array::WITHOUT_GHOSTS,
                                            grid->z())),
    m_melt_out(grid, "debris_melt_out"),
    m_burial(grid, "debris_burial"),
    m_basal_loss(grid, "debris_basal_loss"),
    m_ice_free_loss(grid, "debris_ice_free_loss"),
    m_H_sub(grid, "debris_substep_ice_thickness") {

  m_mass->metadata(0)
      .long_name("englacial debris mass per unit area per vertical control volume")
      .units("kg m^-2");
  m_mass->set(0.0);

  m_melt_out.metadata(0).long_name("englacial debris released by surface melt").units("kg m^-2");
  m_burial.metadata(0).long_name("debris buried by accumulation").units("kg m^-2");
  m_basal_loss.metadata(0).long_name("englacial debris removed by basal melt").units("kg m^-2");
  m_ice_free_loss.metadata(0).long_name("englacial debris lost with the ice of columns that became ice-free").units("kg m^-2");

  if (not (cfl_ratio > 0.0)) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "the CFL ratio of the englacial debris advection has to be"
                                  " positive (got %f)", cfl_ratio);
  }
}

void EnglacialTransport::step(double dt, const array::Scalar1 &H_old,
                              const array::CellType1 &cell_type_old,
                              const array::Scalar &H_new, const array::CellType1 &cell_type_new,
                              const array::Scalar &top_smb, const array::Scalar &bottom_smb,
                              const array::Scalar &input_rate,
                              const array::Array3D &u, const array::Array3D &v,
                              const array::Array3D &w) {

  // number of sub-steps: each has to satisfy the vertical CFL condition in every column
  // (the horizontal CFL condition is the caller's responsibility, see
  // DebrisTransport::max_timestep_impl())
  {
    double dt_local = std::numeric_limits<double>::infinity();
    {
      array::AccessScope list{ &H_old, &cell_type_old, &w };
      for (auto p : m_grid->points()) {
        const int i = p.i(), j = p.j();
        if (cell_type_old.icy(i, j) and H_old(i, j) > 0.0) {
          dt_local = std::min(dt_local, column::vertical_dt_max(m_zi, H_old(i, j), w.get_column(i, j)));
        }
      }
    }
    m_dt_vertical = GlobalMin(m_grid->com, dt_local);

    m_substeps = 1;
    if (dt > 0.0 and std::isfinite(m_dt_vertical) and m_dt_vertical > 0.0) {
      const double n = std::ceil(dt / (m_cfl_ratio * m_dt_vertical));
      if (n > 1.0e9) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "englacial debris advection: the vertical CFL condition"
                                      " requires %.3g sub-steps in a time step of %.3g s"
                                      " (largest allowed sub-step: %.3g s); check the vertical velocity",
                                      n, dt, m_cfl_ratio * m_dt_vertical);
      }
      m_substeps = std::max(1U, (unsigned int)n);
    }
  }

  for (auto *v : { &m_melt_out, &m_burial, &m_basal_loss, &m_ice_free_loss }) {
    v->set(0.0);
  }

  const double f = 1.0 / m_substeps, dt_sub = dt * f;

  for (unsigned int s = 1; s <= m_substeps; ++s) {
    // 1. advection on the geometry at the start of the sub-step
    const array::Scalar1 *H_adv = &H_old;
    if (s > 1) {
      const double f_prev = (s - 1) * f;
      array::AccessScope list{ &H_old, &H_new, &m_H_sub };
      for (auto p : m_grid->points()) {
        const int i = p.i(), j = p.j();
        m_H_sub(i, j) = H_old(i, j) + f_prev * (H_new(i, j) - H_old(i, j));
      }
      m_H_sub.update_ghosts();
      H_adv = &m_H_sub;
    }
    m_scheme->update(dt_sub, *H_adv, cell_type_old, *m_mass, u, v, w);
    m_mass->copy_from(m_scheme->x());

    // 2-5. thickness changes of the sub-step
    column_step(s * f, f, dt_sub, H_old, H_new, cell_type_new, top_smb, bottom_smb, input_rate,
                s == m_substeps);
  }
}

void EnglacialTransport::column_step(double f_new, double f_sub, double dt_sub,
                                     const array::Scalar &H_old, const array::Scalar &H_new,
                                     const array::CellType1 &cell_type_new,
                                     const array::Scalar &top_smb, const array::Scalar &bottom_smb,
                                     const array::Scalar &input_rate, bool final) {
  const int Mz = (int)m_grid->Mz();

  array::AccessScope list{ &H_old, &H_new, &cell_type_new, &top_smb, &bottom_smb, &input_rate,
                           m_mass.get(), &m_melt_out, &m_burial, &m_basal_loss, &m_ice_free_loss };

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();

    double *m = m_mass->get_column(i, j);

    // thickness at the end of this sub-step and the mass balance applied during it
    const double H = H_old(i, j) + f_new * (H_new(i, j) - H_old(i, j)),
                 dH_top = f_sub * top_smb(i, j), dH_bottom = f_sub * bottom_smb(i, j);

    // 5. columns without ice give up their debris
    if (H <= 0.0 or (final and not cell_type_new.icy(i, j))) {
      m_ice_free_loss(i, j) += column::total_mass(Mz, m);
      for (int k = 0; k < Mz; ++k) {
        m[k] = 0.0;
      }
      continue;
    }

    // 2. thickness after the flow step, before mass fluxes were applied
    const double H_flow = std::max(H - dH_top - dH_bottom, 0.0);
    column::fold_above(m_zi, H_flow, m);

    // 3. basal mass balance: melt removes the bottom of the column; freeze-on adds clean
    // ice at the base (the debris keeps its height above the bed)
    if (dH_bottom < 0.0) {
      m_basal_loss(i, j) += column::remove_bottom(m_zi, H_flow, -dH_bottom, m);
    }
    const double H_base = std::max(H_flow + dH_bottom, 0.0);

    // 4. surface mass balance
    if (dH_top < 0.0) {
      m_melt_out(i, j) += column::remove_top(m_zi, H_base, H, m);
    } else if (dH_top > 0.0) {
      // accumulation buries the debris arriving at input cells (equation 16)
      const double M_add = std::max(input_rate(i, j), 0.0) * dt_sub * m_solid_density;
      column::add_top(m_zi, H_base, H, M_add, m);
      m_burial(i, j) += M_add;
    }

    // guard against mass left above the final surface by rounding
    column::fold_above(m_zi, H, m);
  }
}

array::Array3D &EnglacialTransport::mass() {
  return *m_mass;
}

const array::Array3D &EnglacialTransport::mass() const {
  return *m_mass;
}

const array::Scalar &EnglacialTransport::melt_out() const {
  return m_melt_out;
}

const array::Scalar &EnglacialTransport::burial() const {
  return m_burial;
}

const array::Scalar &EnglacialTransport::basal_loss() const {
  return m_basal_loss;
}

const array::Scalar &EnglacialTransport::ice_free_loss() const {
  return m_ice_free_loss;
}

void EnglacialTransport::concentration(const array::Scalar &H, array::Array3D &result) const {
  const int Mz = (int)m_grid->Mz();
  std::vector<double> C(Mz);

  array::AccessScope list{ &H, m_mass.get(), &result };

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();
    column::mass_to_concentration(m_zi, H(i, j), m_mass->get_column(i, j), C.data());
    result.set_column(i, j, C.data());
  }
}

void EnglacialTransport::set_concentration(const array::Array3D &C, const array::Scalar &H) {
  array::AccessScope list{ &H, &C, m_mass.get() };

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();
    column::concentration_to_mass(m_zi, H(i, j), C.get_column(i, j), m_mass->get_column(i, j));
  }
}

void EnglacialTransport::column_mass(array::Scalar &result) const {
  array::sum_columns(*m_mass, 0.0, 1.0, result);
}

unsigned int EnglacialTransport::substeps() const {
  return m_substeps;
}

double EnglacialTransport::vertical_dt_max() const {
  return m_dt_vertical;
}

double EnglacialTransport::local_mass() const {
  const int Mz = (int)m_grid->Mz();
  const double cell_area = m_grid->cell_area();

  array::AccessScope list{ m_mass.get() };

  double result = 0.0;
  for (auto p : m_grid->points()) {
    result += column::total_mass(Mz, m_mass->get_column(p.i(), p.j()));
  }
  return result * cell_area;
}

} // end of namespace debris
} // end of namespace pism
