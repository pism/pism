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

#include <algorithm>            // std::min, std::max
#include <cmath>                // std::fabs
#include <limits>

#include "pism/coupler/debris/column_kernels.hh"

namespace pism {
namespace debris {
namespace column {

std::vector<double> interfaces(const std::vector<double> &z) {
  const size_t Mz = z.size();
  std::vector<double> zi(Mz + 1);

  zi[0] = 0.0;
  for (size_t k = 1; k < Mz; ++k) {
    zi[k] = 0.5 * (z[k - 1] + z[k]);
  }
  zi[Mz] = std::numeric_limits<double>::infinity();

  return zi;
}

double volume_thickness(const std::vector<double> &zi, int k, double H) {
  return std::max(0.0, std::min(zi[k + 1], H) - zi[k]);
}

int volume_index(const std::vector<double> &zi, double h) {
  const int Mz = (int)zi.size() - 1;
  for (int k = 0; k < Mz; ++k) {
    if (h < zi[k + 1]) {
      return k;
    }
  }
  return Mz - 1;
}

void mass_to_concentration(const std::vector<double> &zi, double H,
                           const double *m, double *C) {
  const int Mz = (int)zi.size() - 1;
  for (int k = 0; k < Mz; ++k) {
    double t = volume_thickness(zi, k, H);
    C[k] = t > 0.0 ? m[k] / t : 0.0;
  }
}

void concentration_to_mass(const std::vector<double> &zi, double H,
                           const double *C, double *m) {
  const int Mz = (int)zi.size() - 1;
  for (int k = 0; k < Mz; ++k) {
    m[k] = C[k] * volume_thickness(zi, k, H);
  }
}

double total_mass(int Mz, const double *m) {
  double result = 0.0;
  for (int k = 0; k < Mz; ++k) {
    result += m[k];
  }
  return result;
}

double remove_top(const std::vector<double> &zi, double H_old, double H_new, double *m) {
  if (H_new >= H_old) {
    return 0.0;
  }
  H_new = std::max(H_new, 0.0);

  const int Mz = (int)zi.size() - 1;

  double removed = 0.0;
  for (int k = Mz - 1; k >= 0; --k) {
    const double bottom = zi[k], top = std::min(zi[k + 1], H_old);

    if (top <= bottom) {
      // this volume is above the old surface: nothing here
      continue;
    }

    if (bottom >= H_new) {
      // entirely removed
      removed += m[k];
      m[k] = 0.0;
      continue;
    }

    // straddles the new surface: remove the fraction above it (uniform concentration)
    const double fraction = (top - H_new) / (top - bottom);
    removed += fraction * m[k];
    m[k] *= (1.0 - fraction);
    break;
  }

  return removed;
}

void add_top(const std::vector<double> &zi, double H_old, double H_new, double M_add,
             double *m) {
  if (M_add == 0.0) {
    return;
  }

  const int Mz = (int)zi.size() - 1;

  if (H_new <= H_old) {
    m[volume_index(zi, std::max(H_old, 0.0))] += M_add;
    return;
  }

  const double dH = H_new - H_old;
  for (int k = 0; k < Mz; ++k) {
    const double overlap = std::max(0.0, std::min(zi[k + 1], H_new) - std::max(zi[k], H_old));
    if (overlap > 0.0) {
      m[k] += M_add * (overlap / dH);
    }
  }
}

double remove_bottom(const std::vector<double> &zi, double H, double dH, double *m) {
  if (dH <= 0.0 or H <= 0.0) {
    return 0.0;
  }
  dH = std::min(dH, H);

  const int Mz = (int)zi.size() - 1;

  const double M_old = total_mass(Mz, m);

  // Concentration profile of the old column (uniform within each volume), re-sampled
  // after shifting the column down by dH: the ice that was at height h is now at h - dH.
  std::vector<double> C(Mz);
  mass_to_concentration(zi, H, m, C.data());

  const double H_new = H - dH;
  double M_new = 0.0;
  for (int k = 0; k < Mz; ++k) {
    // new volume k covers [zi[k], min(zi[k+1], H_new)] and contains the old ice from
    // [zi[k] + dH, min(zi[k+1], H_new) + dH]
    const double bottom = zi[k] + dH, top = std::min(zi[k + 1], H_new) + dH;
    double mass = 0.0;
    if (top > bottom) {
      for (int l = 0; l < Mz; ++l) {
        const double overlap = std::max(0.0, std::min(std::min(zi[l + 1], H), top) -
                                        std::max(zi[l], bottom));
        mass += C[l] * overlap;
      }
    }
    m[k] = mass;
    M_new += mass;
  }

  // everything that did not survive the shift was in the melted ice
  return std::max(0.0, M_old - M_new);
}

double fold_above(const std::vector<double> &zi, double H, double *m) {
  const int Mz = (int)zi.size() - 1;

  if (H <= 0.0) {
    // no ice: nothing to fold into; the caller decides what to do with the mass
    return 0.0;
  }

  const int k_top = volume_index(zi, H);

  double moved = 0.0;
  for (int k = k_top + 1; k < Mz; ++k) {
    moved += m[k];
    m[k] = 0.0;
  }
  m[k_top] += moved;

  return moved;
}

double vertical_dt_max(const std::vector<double> &zi, double H, const double *w) {
  const int Mz = (int)zi.size() - 1;

  double result = std::numeric_limits<double>::infinity();
  for (int k = 0; k < Mz; ++k) {
    const double dz = volume_thickness(zi, k, H);
    if (dz <= 0.0) {
      break;
    }
    const double speed = std::fabs(w[k]);
    if (speed > 0.0) {
      result = std::min(result, dz / speed);
    }
  }

  return result;
}

} // end of namespace column
} // end of namespace debris
} // end of namespace pism
