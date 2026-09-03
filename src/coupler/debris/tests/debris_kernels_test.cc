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

/* Unit tests of the debris column and terminus kernels. Every case has a closed-form
   answer; no PISM, PETSc or MPI is needed. */

#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <vector>

#include "pism/coupler/debris/column_kernels.hh"
#include "pism/coupler/debris/terminus_kernels.hh"

using namespace pism::debris;

static int failures = 0;

static void check_close(const char *what, double actual, double expected, double tol) {
  if (std::fabs(actual - expected) > tol) {
    std::fprintf(stderr, "FAILED: %s: got %.12g, expected %.12g (tolerance %g)\n",
                 what, actual, expected, tol);
    failures += 1;
  }
}

static void check_true(const char *what, bool ok) {
  if (not ok) {
    std::fprintf(stderr, "FAILED: %s\n", what);
    failures += 1;
  }
}

//! PISM's "quadratic" vertical spacing: fine near the bed, coarse near the top.
static std::vector<double> quadratic_levels(double Lz, int Mz) {
  std::vector<double> z(Mz);
  for (int k = 0; k < Mz; ++k) {
    double s = (double)k / (Mz - 1);
    z[k] = Lz * s * s;
  }
  return z;
}

static void test_interfaces() {
  std::vector<double> z = { 0.0, 1.0, 3.0, 6.0 };
  auto zi = column::interfaces(z);

  check_true("interfaces: size", zi.size() == z.size() + 1);
  check_close("interfaces: zi[0]", zi[0], 0.0, 0.0);
  check_close("interfaces: zi[1]", zi[1], 0.5, 0.0);
  check_close("interfaces: zi[2]", zi[2], 2.0, 0.0);
  check_close("interfaces: zi[3]", zi[3], 4.5, 0.0);
  check_true("interfaces: top is unbounded", std::isinf(zi[4]));

  check_close("volume_thickness: full", column::volume_thickness(zi, 1, 10.0), 1.5, 1e-15);
  check_close("volume_thickness: truncated", column::volume_thickness(zi, 3, 5.0), 0.5, 1e-15);
  check_close("volume_thickness: above H", column::volume_thickness(zi, 3, 4.0), 0.0, 0.0);

  check_true("volume_index: bottom", column::volume_index(zi, 0.25) == 0);
  check_true("volume_index: interior", column::volume_index(zi, 3.0) == 2);
  check_true("volume_index: top", column::volume_index(zi, 100.0) == 3);
}

static void test_conversions() {
  const int Mz = 21;
  auto z  = quadratic_levels(100.0, Mz);
  auto zi = column::interfaces(z);
  const double H = 73.0;

  std::vector<double> C(Mz), m(Mz), C2(Mz);
  for (int k = 0; k < Mz; ++k) {
    C[k] = 1.0 + 0.1 * k;
  }

  column::concentration_to_mass(zi, H, C.data(), m.data());
  column::mass_to_concentration(zi, H, m.data(), C2.data());

  for (int k = 0; k < Mz; ++k) {
    if (column::volume_thickness(zi, k, H) > 0.0) {
      check_close("mass <-> concentration round trip", C2[k], C[k], 1e-12);
    } else {
      check_close("mass <-> concentration: empty volume", C2[k], 0.0, 0.0);
    }
  }

  // total mass = integral of the piecewise constant concentration over [0, H]
  double expected = 0.0;
  for (int k = 0; k < Mz; ++k) {
    expected += C[k] * column::volume_thickness(zi, k, H);
  }
  check_close("total_mass", column::total_mass(Mz, m.data()), expected, 1e-12);
}

static void test_remove_top() {
  const int Mz = 21;
  auto z  = quadratic_levels(100.0, Mz);
  auto zi = column::interfaces(z);

  // uniform concentration: removing dz of ice removes C dz of mass
  const double C0 = 2.5, H_old = 80.0, H_new = 61.3;
  std::vector<double> C(Mz, C0), m(Mz);
  column::concentration_to_mass(zi, H_old, C.data(), m.data());
  const double M_old = column::total_mass(Mz, m.data());

  double removed = column::remove_top(zi, H_old, H_new, m.data());

  check_close("remove_top: removed mass", removed, C0 * (H_old - H_new), 1e-10);
  check_close("remove_top: conservation", M_old - column::total_mass(Mz, m.data()), removed, 1e-10);

  // the remaining column still has uniform concentration C0
  std::vector<double> C2(Mz);
  column::mass_to_concentration(zi, H_new, m.data(), C2.data());
  for (int k = 0; k < Mz; ++k) {
    if (column::volume_thickness(zi, k, H_new) > 0.0) {
      check_close("remove_top: concentration unchanged", C2[k], C0, 1e-10);
    }
  }

  // removing everything: only volumes below the old surface are part of the column
  std::fill(m.begin(), m.end(), 1.0);
  int n_below = 0;
  for (int k = 0; k < Mz; ++k) {
    n_below += column::volume_thickness(zi, k, H_old) > 0.0 ? 1 : 0;
  }
  removed = column::remove_top(zi, H_old, 0.0, m.data());
  check_close("remove_top: remove all", removed, (double)n_below, 1e-12);
  check_close("remove_top: nothing left below H_old",
              column::total_mass(Mz, m.data()), (double)(Mz - n_below), 0.0);

  // no-op when the surface rises
  std::fill(m.begin(), m.end(), 1.0);
  removed = column::remove_top(zi, 50.0, 60.0, m.data());
  check_close("remove_top: rising surface is a no-op", removed, 0.0, 0.0);
  check_close("remove_top: rising surface keeps mass", column::total_mass(Mz, m.data()), (double)Mz, 0.0);
}

static void test_add_top() {
  const int Mz = 21;
  auto z  = quadratic_levels(100.0, Mz);
  auto zi = column::interfaces(z);

  const double M_add = 3.0;
  std::vector<double> m(Mz, 0.0);

  // Surfaces aligned with volume interfaces: the added ice fills whole volumes, so
  // removing it again returns exactly what was added.
  {
    const double H_old = zi[12], H_new = zi[15];
    column::add_top(zi, H_old, H_new, M_add, m.data());
    check_close("add_top: conservation", column::total_mass(Mz, m.data()), M_add, 1e-12);
    check_close("add_top: uniform concentration",
                m[12] / (zi[13] - zi[12]), M_add / (H_new - H_old), 1e-12);

    double removed = column::remove_top(zi, H_new, H_old, m.data());
    check_close("add_top then remove_top (aligned)", removed, M_add, 1e-10);
    check_close("add_top then remove_top (aligned): nothing left",
                column::total_mass(Mz, m.data()), 0.0, 1e-12);
  }

  // General surfaces: the mass added to the volume straddling H_old is mixed with that
  // whole volume (piecewise-constant representation), so a round trip conserves mass but
  // leaves part of it behind.
  {
    const double H_old = 40.0, H_new = 47.5;
    std::fill(m.begin(), m.end(), 0.0);
    column::add_top(zi, H_old, H_new, M_add, m.data());
    double removed = column::remove_top(zi, H_new, H_old, m.data());
    check_true("add_top then remove_top: removed <= added", removed <= M_add + 1e-12);
    check_close("add_top then remove_top: conservation",
                removed + column::total_mass(Mz, m.data()), M_add, 1e-12);
  }

  const double H_old = 40.0, H_new = 47.5;

  // mass goes into the surface volume when the thickness does not change
  std::fill(m.begin(), m.end(), 0.0);
  column::add_top(zi, H_old, H_old, M_add, m.data());
  check_close("add_top: no thickness change", m[column::volume_index(zi, H_old)], M_add, 0.0);

  // zero mass is a no-op
  column::add_top(zi, H_old, H_new, 0.0, m.data());
  check_close("add_top: zero mass", column::total_mass(Mz, m.data()), M_add, 0.0);
}

static void test_remove_bottom() {
  const int Mz = 21;
  auto z  = quadratic_levels(100.0, Mz);
  auto zi = column::interfaces(z);

  // uniform concentration: removing dH at the bottom removes C dH and leaves the
  // concentration unchanged
  const double C0 = 1.7, H = 65.0, dH = 4.25;
  std::vector<double> C(Mz, C0), m(Mz);
  column::concentration_to_mass(zi, H, C.data(), m.data());

  double removed = column::remove_bottom(zi, H, dH, m.data());
  check_close("remove_bottom: removed mass", removed, C0 * dH, 1e-10);
  check_close("remove_bottom: remaining mass", column::total_mass(Mz, m.data()), C0 * (H - dH), 1e-10);

  std::vector<double> C2(Mz);
  column::mass_to_concentration(zi, H - dH, m.data(), C2.data());
  for (int k = 0; k < Mz; ++k) {
    if (column::volume_thickness(zi, k, H - dH) > 0.0) {
      check_close("remove_bottom: concentration unchanged", C2[k], C0, 1e-10);
    }
  }

  // a debris layer at the bottom is melted out completely
  std::fill(m.begin(), m.end(), 0.0);
  m[0] = 5.0;                   // volume 0 spans [0, zi[1]]
  removed = column::remove_bottom(zi, H, zi[1], m.data());
  check_close("remove_bottom: bottom layer gone", removed, 5.0, 1e-12);
  check_close("remove_bottom: nothing left", column::total_mass(Mz, m.data()), 0.0, 1e-12);

  // a layer higher up moves down by dH: mass is conserved
  std::fill(m.begin(), m.end(), 0.0);
  m[10] = 2.0;
  removed = column::remove_bottom(zi, H, 1.0, m.data());
  check_close("remove_bottom: interior layer kept", removed, 0.0, 1e-12);
  check_close("remove_bottom: interior layer conserved", column::total_mass(Mz, m.data()), 2.0, 1e-12);

  // no-ops
  std::fill(m.begin(), m.end(), 1.0);
  check_close("remove_bottom: dH <= 0", column::remove_bottom(zi, H, -1.0, m.data()), 0.0, 0.0);
  check_close("remove_bottom: H <= 0", column::remove_bottom(zi, 0.0, 1.0, m.data()), 0.0, 0.0);
}

static void test_fold_above() {
  const int Mz = 21;
  auto z  = quadratic_levels(100.0, Mz);
  auto zi = column::interfaces(z);

  std::vector<double> m(Mz, 1.0);
  const double H = 30.0;
  const int k_top = column::volume_index(zi, H);

  double moved = column::fold_above(zi, H, m.data());
  check_close("fold_above: moved", moved, (double)(Mz - 1 - k_top), 0.0);
  check_close("fold_above: conservation", column::total_mass(Mz, m.data()), (double)Mz, 0.0);
  check_close("fold_above: top volume", m[k_top], 1.0 + moved, 0.0);
  for (int k = k_top + 1; k < Mz; ++k) {
    check_close("fold_above: emptied", m[k], 0.0, 0.0);
  }

  check_close("fold_above: no ice", column::fold_above(zi, 0.0, m.data()), 0.0, 0.0);
}

static void test_vertical_dt_max() {
  std::vector<double> z = { 0.0, 1.0, 3.0, 6.0, 10.0 };
  auto zi = column::interfaces(z);
  // volumes: [0,0.5] [0.5,2] [2,4.5] [4.5,8] [8,inf)

  std::vector<double> w = { 0.0, -2.0, 1.0, 4.0, 100.0 };
  // H = 7: thicknesses 0.5, 1.5, 2.5, 2.5, 0 -> dz/|w| = inf, 0.75, 2.5, 0.625
  check_close("vertical_dt_max", column::vertical_dt_max(zi, 7.0, w.data()), 0.625, 1e-15);
  // H = 3: 0.5, 1.5, 1.0, 0 -> inf, 0.75, 1.0
  check_close("vertical_dt_max: truncated", column::vertical_dt_max(zi, 3.0, w.data()), 0.75, 1e-15);

  std::vector<double> rest(5, 0.0);
  check_true("vertical_dt_max: at rest", std::isinf(column::vertical_dt_max(zi, 7.0, rest.data())));
}

static void test_upstream_chain() {
  // 7x7 surface rising toward the north-east with a ridge (maximum) at (5, 5); ice-free
  // outside [0, 6]
  auto surface = [](int i, int j) -> double {
    double x = i, y = j;
    return 100.0 - ((x - 5.0) * (x - 5.0) + (y - 5.0) * (y - 5.0));
  };
  auto icy = [](int i, int j) -> bool { return i >= 0 and i <= 6 and j >= 0 and j <= 6; };

  auto chain = terminus::upstream_chain(surface, icy, 1, 1, 3, 1.0, 1.0);
  check_true("upstream_chain: length", chain.size() == 3);
  if (chain.size() == 3) {
    check_true("upstream_chain: step 1", chain[0].di == 1 and chain[0].dj == 1);
    check_true("upstream_chain: step 2", chain[1].di == 2 and chain[1].dj == 2);
    check_true("upstream_chain: step 3", chain[2].di == 3 and chain[2].dj == 3);
  }

  // stops at the maximum
  chain = terminus::upstream_chain(surface, icy, 4, 4, 5, 1.0, 1.0);
  check_true("upstream_chain: stops at the maximum", chain.size() == 1);

  chain = terminus::upstream_chain(surface, icy, 5, 5, 5, 1.0, 1.0);
  check_true("upstream_chain: starting at the maximum", chain.empty());

  // does not step onto ice-free cells: from (0, 5) the steepest ascent would be (1, 5)
  // then (2, 5), ...; mark i = 1 ice-free and check that the walk goes around it
  auto icy2 = [](int i, int j) -> bool { return i >= 0 and i <= 6 and j >= 0 and j <= 6 and i != 1; };
  chain = terminus::upstream_chain(surface, icy2, 0, 5, 1, 1.0, 1.0);
  check_true("upstream_chain: avoids ice-free cells", chain.empty());

  check_true("upstream_chain: n = 0", terminus::upstream_chain(surface, icy, 1, 1, 0, 1.0, 1.0).empty());

  // a surface uniform in y rising in x: the steepest ascent is along x, not diagonal
  auto ramp = [](int i, int j) -> double { (void) j; return 10.0 * i; };
  chain = terminus::upstream_chain(ramp, icy, 2, 3, 2, 25.0, 25.0);
  check_true("upstream_chain: axis-aligned on a ramp",
             chain.size() == 2 and chain[0].di == 1 and chain[0].dj == 0 and
             chain[1].di == 2 and chain[1].dj == 0);
}

static void test_foreland_slope() {
  // ice for i <= 3, ice-free for i >= 4; surface decreasing in x with slope 0.1, flat in y
  auto surface  = [](int i, int j) -> double { (void) j; return 100.0 - 0.1 * 25.0 * i; };
  auto ice_free = [](int i, int j) -> bool { (void) j; return i >= 4; };

  check_close("foreland_slope: margin cell", terminus::foreland_slope(surface, ice_free, 3, 3, 25.0, 25.0), 0.1, 1e-14);
  check_close("foreland_slope: interior cell", terminus::foreland_slope(surface, ice_free, 2, 3, 25.0, 25.0), 0.0, 0.0);

  // a lower ice-free neighbor in both x and y: combine the two slopes
  auto surface2  = [](int i, int j) -> double { return 100.0 - 3.0 * i - 4.0 * j; };
  auto ice_free2 = [](int i, int j) -> bool { return i >= 4 or j >= 4; };
  check_close("foreland_slope: corner cell", terminus::foreland_slope(surface2, ice_free2, 3, 3, 1.0, 1.0), 5.0, 1e-14);

  // an ice-free neighbor that is *higher* does not count (a headwall)
  auto surface3  = [](int i, int j) -> double { (void) j; return 100.0 + 10.0 * i; };
  check_close("foreland_slope: headwall", terminus::foreland_slope(surface3, ice_free, 3, 3, 25.0, 25.0), 0.0, 0.0);
}

int main() {
  test_interfaces();
  test_conversions();
  test_remove_top();
  test_add_top();
  test_remove_bottom();
  test_fold_above();
  test_vertical_dt_max();
  test_upstream_chain();
  test_foreland_slope();

  if (failures > 0) {
    std::fprintf(stderr, "%d check(s) FAILED\n", failures);
    return 1;
  }

  std::printf("debris kernels: all checks passed\n");
  return 0;
}
