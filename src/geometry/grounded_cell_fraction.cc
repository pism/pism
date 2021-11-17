/* Copyright (C) 2018, 2020, 2021 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include <cassert>
#include <cmath>                // fabs, isnan

#include "grounded_cell_fraction.hh"

#include "pism/util/error_handling.hh"
#include "pism/util/pism_utilities.hh" // clip
#include "pism/util/iceModelVec.hh"

namespace pism {

struct Point {
  double x;
  double y;
};

/*!
 * Consider the right triangle A(0,0) - B(1,0) - C(0,1).
 *
 * Define a linear function z = a + (b - a) * x + (c - a) * y, where a, b, and c are its
 * values at points A, B, and C, respectively.
 *
 * Our goal is to find the fraction of the area of ABC in where z > 0.
 *
 * We note that z(x,y) is continuous, so unless a, b, and c have the same sign the line
 *
 * z = 0
 *
 * will intersect exactly two of the sides (possibly at a node of ABC).
 *
 * So, if the line (z = 0) does not intersect BC, for example, then it has to intersect AB
 * and AC.
 *
 * This method can be applied to arbitrary triangles. It does not even matter if values at
 * triangle nodes (a, b, c) are listed in the same order. (For any two triangles on a
 * plane there exists an affine map that takes one to the other. Also, affine maps
 * preserve ratios of areas of figures.)
 */

/*!
 * Compute the area of a triangle on a plane using the "shoelace formula."
 */
static inline double triangle_area(const Point &a, const Point &b, const Point &c) {
  // note: fabs should not be needed since we traverse all triangle nodes
  // counter-clockwise, but it is good to be safe
  return 0.5 * std::fabs((a.x - c.x) * (b.y - a.y) - (a.x - b.x) * (c.y - a.y));
}

/*!
 * Compute the coordinates of the intersection of (z = 0) with the side AB.
 */
Point intersect_ab(double a, double b) {
  if (a != b) {
    return {a / (a - b), 0.0};
  }
  return {-1.0, -1.0};        // no intersection
}

/*!
 * Compute the coordinates of the intersection of (z = 0) with the side BC.
 */
Point intersect_bc(double b, double c) {
  if (b != c) {
    return {c / (c - b), b / (b - c)};
  }
  return {-1.0, -1.0};        // no intersection
}

/*!
 * Compute the coordinates of the intersection of (z = 0) with the side AC.
 */
Point intersect_ac(double a, double c) {
  if (a != c) {
    return {0.0, a / (a - c)};
  }
  return {-1.0, -1.0};        // no intersection
}

/*!
 * Return true if a point p is not a valid point on a side of the triangle
 * (0,0)-(1,0)-(0,1).
 *
 * This is not a complete test (for example, it does not check if y = 1 - x for points on
 * the hypotenuse). The point of this is to exclude the kinds of invalid points we are
 * likely to see, not all of them.
 *
 * Note that we use (-1, -1) to indicate "invalid" points in the rest of the code and
 * these are easy to detect: they require only one comparison.
 */
bool invalid(const Point &p) {
  return (p.x < 0.0 or p.x > 1.0 or p.y < 0.0 or p.y > 1.0);
}

/*!
 * Return true if two points are the same.
 */
static bool same(const Point &a, const Point &b) {
  const double threshold = 1e-12;
  return std::fabs(a.x - b.x) < threshold and std::fabs(a.y - b.y) < threshold;
}

/*!
 * Consider the right triangle A(0,0) - B(1,0) - C(0,1).
 *
 * Define a linear function z = a + (b - a) * x + (c - a) * y, where a, b, and c are its
 * values at points A, B, and C, respectively.
 *
 * Our goal is to find the fraction of the triangle ABC in where z > 0.
 *
 * This corresponds to the grounded area fraction if z is the flotation criterion
 * function.
 */
double grounded_area_fraction(double a, double b, double c) {

  assert(not std::isnan(a));
  assert(not std::isnan(b));
  assert(not std::isnan(c));

  if (a > 0.0 and b > 0.0 and c > 0.0) {
    return 1.0;
  }

  if (a <= 0.0 and b <= 0.0 and c <= 0.0) {
    return 0.0;
  }

  // the area of the triangle (0,0)-(1,0)-(0,1)
  const double total_area = 0.5;

  const Point
    ab = intersect_ab(a, b),
    bc = intersect_bc(b, c),
    ac = intersect_ac(a, c);

  if (invalid(bc)) {
    assert(not (invalid(ab) or invalid(ac)));

    double ratio = triangle_area({0.0, 0.0}, ab, ac) / total_area;
    assert((ratio >= 0.0) and (ratio <= 1.0));

    if (a > 0.0) {
      return ratio;
    }
    return 1.0 - ratio;
  }

  if (invalid(ac)) {
    assert(not (invalid(ab) or invalid(bc)));

    double ratio = triangle_area({1.0, 0.0}, bc, ab) / total_area;
    assert((ratio >= 0.0) and (ratio <= 1.0));

    if (b > 0.0) {
      return ratio;
    }
    return 1.0 - ratio;
  }

  if (invalid(ab)) {
    assert(not (invalid(bc) or invalid(ac)));

    double ratio = triangle_area({0.0, 1.0}, ac, bc) / total_area;
    assert((ratio >= 0.0) and (ratio <= 1.0));

    if (c > 0.0) {
      return ratio;
    }
    return 1.0 - ratio;
  }

  // Note that we know that ab, bc, and ac are all valid.

  // the a == 0 case, the line F = 0 goes through A
  if (same(ab, ac)) {
    double ratio = triangle_area({1.0, 0.0}, bc, ab) / total_area;
    assert((ratio >= 0.0) and (ratio <= 1.0));

    if (b > 0.0) {
      return ratio;
    }
    return 1.0 - ratio;
  }

  // the b == 0 case and the c == 0 case
  if (same(ab, bc) or same(ac, bc)) {
    double ratio = triangle_area({0.0, 0.0}, ab, ac) / total_area;
    assert((ratio >= 0.0) and (ratio <= 1.0));

    if (a > 0.0) {
      return ratio;
    }
    return 1.0 - ratio;
  }

  // Note: the case of F=0 coinciding with a side of the triangle is covered by if clauses
  // above. For example, when F=0 coincides with AC, we have a = c = 0 and intersect_ac(a, c)
  // returns an invalid intersection point.

  throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                "the logic in grounded_area_fraction failed! Please submit a bug report.");
}

/*!
 * The flotation criterion.
 */
static double F(double SL, double B, double H, double alpha) {
  double
    water_depth = SL - B,
    shelf_depth = H * alpha;
  return shelf_depth - water_depth;
}


/*!
 * Compute the flotation criterion at all the points in the box stencil.
 */
typedef stencils::Box<double> Box;
static Box F(const Box &SL, const Box &B, const Box &H, double alpha) {
  return {F(SL.ij, B.ij, H.ij, alpha),
          F(SL.n,  B.n,  H.n,  alpha),
          F(SL.nw, B.nw, H.nw, alpha),
          F(SL.w,  B.w,  H.w,  alpha),
          F(SL.sw, B.sw, H.sw, alpha),
          F(SL.s,  B.s,  H.s,  alpha),
          F(SL.se, B.se, H.se, alpha),
          F(SL.e,  B.e,  H.e,  alpha),
          F(SL.ne, B.ne, H.ne, alpha)};
}

/*!
 * @param[in] ice_density ice density, kg/m3
 * @param[in] ocean_density ocean_density, kg/m3
 * @param[in] sea_level sea level (flotation) elevation, m
 * @param[in] ice_thickness ice thickness, m
 * @param[in] bed_topography bed elevation, m
 * @param[out] result grounded cell fraction, between 0 (floating) and 1 (grounded)
 */
void compute_grounded_cell_fraction(double ice_density,
                                    double ocean_density,
                                    const IceModelVec2S &sea_level,
                                    const IceModelVec2S &ice_thickness,
                                    const IceModelVec2S &bed_topography,
                                    IceModelVec2S &result) {
  IceGrid::ConstPtr grid = result.grid();
  double alpha = ice_density / ocean_density;

  IceModelVec::AccessList list{&sea_level, &ice_thickness, &bed_topography, &result};

  ParallelSection loop(grid->com);
  try {
    for (Points p(*grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      /*
        NW----------------N----------------NE
        |                 |                 |
        |                 |                 |
        |       nw--------n--------ne       |
        |        |        |        |        |
        |        |        |        |        |
        W--------w--------o--------e--------E
        |        |        |        |        |
        |        |        |        |        |
        |       sw--------s--------se       |
        |                 |                 |
        |                 |                 |
        SW----------------S----------------SE
      */

      // compute the floatation function at 8 points surrounding the current grid point
      stencils::Box<double> f;
      {
        auto S = sea_level.box(i, j);
        auto H = ice_thickness.box(i, j);
        auto B = bed_topography.box(i, j);

        auto x = F(S, B, H, alpha);

        f.ij = x.ij;
        f.sw = 0.25 * (x.sw + x.s + x.ij + x.w);
        f.se = 0.25 * (x.s + x.se + x.e + x.ij);
        f.ne = 0.25 * (x.ij + x.e + x.ne + x.n);
        f.nw = 0.25 * (x.w + x.ij + x.n + x.nw);

        f.s = 0.5 * (x.ij + x.s);
        f.e = 0.5 * (x.ij + x.e);
        f.n = 0.5 * (x.ij + x.n);
        f.w = 0.5 * (x.ij + x.w);
      }

      // compute the grounding fraction for the current cell by breaking it into 8
      // triangles
      double fraction = 0.125 * (grounded_area_fraction(f.ij, f.ne, f.n) +
                                 grounded_area_fraction(f.ij, f.n,  f.nw) +
                                 grounded_area_fraction(f.ij, f.nw, f.w) +
                                 grounded_area_fraction(f.ij, f.w,  f.sw) +
                                 grounded_area_fraction(f.ij, f.sw, f.s) +
                                 grounded_area_fraction(f.ij, f.s,  f.se) +
                                 grounded_area_fraction(f.ij, f.se, f.e) +
                                 grounded_area_fraction(f.ij, f.e,  f.ne));

      result(i, j) = clip(fraction, 0.0, 1.0);

    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

} // end of namespace pism
