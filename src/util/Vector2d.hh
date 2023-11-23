/* Copyright (C) 2015, 2016, 2020, 2021, 2022 PISM Authors
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

#ifndef PISM_VECTOR2D_HH
#define PISM_VECTOR2D_HH

#include <cmath>                // sqrt

namespace pism {

//! @brief This class represents a 2D vector field (such as ice
//! velocity) at a certain grid point.
class Vector2d {
public:
  Vector2d() : u(0), v(0) {}
  Vector2d(double a, double b) : u(a), v(b) {}
  Vector2d(const Vector2d &other) : u(other.u), v(other.v) {}

  //! Magnitude squared.
  inline double magnitude_squared() const {
    return u*u + v*v;
  }
  //! Magnitude.
  inline double magnitude() const {
    return sqrt(magnitude_squared());
  }

  inline Vector2d& operator=(const Vector2d &other) {
    // NOTE: we don't check for self-assignment because there is no memory
    // (de-)allocation here.
    u = other.u;
    v = other.v;
    return *this;
  }

  //! Set both components to the same number.
  inline Vector2d& operator=(const double &a) {
    u = a;
    v = a;
    return *this;
  }

  inline Vector2d& operator+=(const Vector2d &other) {
    u += other.u;
    v += other.v;
    return *this;
  }

  inline Vector2d& operator-=(const Vector2d &other) {
    u -= other.u;
    v -= other.v;
    return *this;
  }

  inline Vector2d& operator*=(const double &a) {
    u *= a;
    v *= a;
    return *this;
  }

  inline Vector2d& operator/=(const double &a) {
    u /= a;
    v /= a;
    return *this;
  }

  //! \brief Adds two vectors.
  inline Vector2d operator+(const Vector2d &other) const {
    return Vector2d(u + other.u, v + other.v);
  }

  //! \brief Substracts two vectors.
  inline Vector2d operator-(const Vector2d &other) const {
    return Vector2d(u - other.u, v - other.v);
  }

  //! \brief Scales a vector.
  inline Vector2d operator*(const double &a) const {
    return Vector2d(u * a, v * a);
  }

  //! \brief Scales a vector.
  inline Vector2d operator/(const double &a) const {
    return Vector2d(u / a, v / a);
  }

  double u, v;
};

// Multiplication of a vector by a constant is commutative.
inline Vector2d operator*(const double &a, const Vector2d &v) {
  return v * a;
}

} // end of namespace pism

#endif /* PISM_VECTOR2D_HH */
