/* Copyright (C) 2018 PISM Authors
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

#include <cmath>                // fabs

#include "grounded_cell_fraction_v2.hh"

#include "pism/util/error_handling.hh"
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
 * Our goal is to find the fraction of the triangle ABC in where z > 0.
 *
 * We note that z(x,y) is continuous, so unless a, b, and c have the same sign the line
 *
 * z = 0
 *
 * will intersect exactly two of the sides (possibly at a node of the triangle ABC).
 *
 * So, if the line (z = 0) does not intersect BC, for example, then it has to intersect AB
 * and AC.
 *
 * This method can be applied to arbitrary triangles as long as their nodes are listed in
 * the correct order. (For any two triangles on a plane there exists an affine map that
 * takes one to the other. Also, affine maps preserve ratios of areas of figures.)
 */

/*!
 * Compute the area of a triangle on a plane using the "shoelace formula."
 */
static inline double triangle_area(const Point &a, const Point &b, const Point &c) {
  // note: fabs should not be needed since we traverse all triangle nodes
  // counter-clockwise, but it is good to be safe
  return 0.5 * fabs((a.x - c.x) * (b.y - a.y) - (a.x - b.x) * (c.y - a.y));
}

/*!
 * Compute the coordinates of the intersection of (z = 0) with the side AB.
 */
Point intersect_ab(double a, double b) {
  if (a != b) {
    return {a / (a - b), 0.0};
  } else {
    return {-1.0, -1.0};        // no intersection
  }
}

/*!
 * Compute the coordinates of the intersection of (z = 0) with the side BC.
 */
Point intersect_bc(double b, double c) {
  if (b != c) {
    return {c / (c - b), b / (b - c)};
  } else {
    return {-1.0, -1.0};        // no intersection
  }
}

/*!
 * Compute the coordinates of the intersection of (z = 0) with the side AC.
 */
Point intersect_ac(double a, double c) {
  if (a != c) {
    return {0.0, a / (a - c)};
  } else {
    return {-1.0, -1.0};        // no intersection
  }
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
  if (p.x < 0.0 or p.x > 1.0 or p.y < 0.0 or p.y > 1.0) {
    return true;
  } else {
    return false;
  }
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

  if (a > 0.0 and b > 0.0 and c > 0.0) {
    return 1.0;
  }

  if (a < 0.0 and b < 0.0 and c < 0.0) {
    return 0.0;
  }

  // the area of the triangle (0,0)-(1,0)-(0,1)
  const double total_area = 0.5;

  Point
    ab = intersect_ab(a, b),
    bc = intersect_bc(b, c),
    ac = intersect_ac(a, c);

  if (invalid(bc)) {
    double ratio = triangle_area({0.0, 0.0}, ab, ac) / total_area;

    if (a > 0.0) {
      return ratio;
    } else {
      return 1.0 - ratio;
    }
  }

  if (invalid(ac)) {
    double ratio = triangle_area({1.0, 0.0}, bc, ab) / total_area;

    if (b > 0.0) {
      return ratio;
    } else {
      return 1.0 - ratio;
    }
  }

  if (invalid(ab)) {
    double ratio = triangle_area({0.0, 1.0}, ac, bc) / total_area;

    if (c > 0.0) {
      return ratio;
    } else {
      return 1.0 - ratio;
    }
  }

  throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                "the logic in grounded_area_fraction failed!");
}

/*!
 * Compute the grounded area fraction for a square cell by splitting it into 4 triangles.
 *
 * We split into 4 triangles instead of 2 to preserve symmetry.
 *
 * Triangles:
 *
 * EAB, EBC, ECD, EDA, where E is the center of the cell.
 *
 * The value at E is computed by averaging A, B, C, and D.
 */
double grounded_area_fraction_square(double A, double B, double C, double D) {

  double E = 0.25 * (A + B + C + D);

  return 0.25 * (grounded_area_fraction(E, A, B),
                 grounded_area_fraction(E, B, C),
                 grounded_area_fraction(E, C, D),
                 grounded_area_fraction(E, D, A));
}

// This structure extracts the box stencil information from an IceModelVec2S.
struct Box {
  double ij, n, nw, w, sw, s, se, e, ne;
  Box(const IceModelVec2S &X, int i, int j) {
    ij = X(i, j);
    n  = X(i, j + 1);
    nw = X(i - 1, j + 1);
    w  = X(i - 1, j);
    sw = X(i - 1, j - 1);
    s  = X(i, j - 1);
    se = X(i + 1, j - 1);
    e  = X(i + 1, j);
    ne = X(i + 1, j + 1);
  }
};

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
 * @param[in] ice_density ice density, kg/m3
 * @param[in] ocean_density ocean_density, kg/m3
 * @param[in] sea_level sea level (flotation) elevation, m
 * @param[in] ice_thickness ice thickness, m
 * @param[in] bed_topography bed elevation, m
 * @param[out] result grounded cell fraction, between 0 (floating) and 1 (grounded)
 */
void compute_grounded_cell_fraction_v2(double ice_density,
                                       double ocean_density,
                                       const IceModelVec2S &sea_level,
                                       const IceModelVec2S &ice_thickness,
                                       const IceModelVec2S &bed_topography,
                                       IceModelVec2S &result) {
  IceGrid::ConstPtr grid = result.grid();
  double alpha = ice_density / ocean_density;

  IceModelVec::AccessList list{&sea_level, &ice_thickness, &bed_topography};

  ParallelSection loop(grid->com);
  try {
    for (Points p(*grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      Box
        S(sea_level, i, j),
        H(ice_thickness, i, j),
        B(bed_topography, i, j);

      /*
     NW-----------------N-----------------NE
      |                 |                 |
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
      |                 |                 |
     SW-----------------S-----------------SE
      */

      double
        S_e = 0.5 * (S.ij + S.e),
        S_w = 0.5 * (S.ij + S.w),
        S_n = 0.5 *

      double
        F_ij = F(S.ij, B.ij, H.ij, alpha),
        F_e  = F(0.5 * (S.ij + S.e),
                 0.5 * (B.ij + B.e),
                 0.5)

    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

}

} // end of namespace pism
