/* Copyright (C) 2015, 2016, 2023 PISM Authors
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

#ifndef _VECTOR2_H_
#define _VECTOR2_H_

#include <cmath>                // sqrt

namespace pism {

//! @brief This class represents a 2D vector field (such as ice
//! velocity) at a certain grid point.
class Vector2 {
public:
  Vector2() : u(0), v(0) {}
  Vector2(double a, double b) : u(a), v(b) {}

  //! Magnitude squared.
  inline double magnitude_squared() const {
    return u*u + v*v;
  }
  //! Magnitude.
  inline double magnitude() const {
    return sqrt(magnitude_squared());
  }

  //! Set both components to the same number.
  inline Vector2& operator=(const double &a) {
    u = a;
    v = a;
    return *this;
  }

  inline Vector2& operator+=(const Vector2 &other) {
    u += other.u;
    v += other.v;
    return *this;
  }

  inline Vector2& operator-=(const Vector2 &other) {
    u -= other.u;
    v -= other.v;
    return *this;
  }

  inline Vector2& operator*=(const double &a) {
    u *= a;
    v *= a;
    return *this;
  }

  inline Vector2& operator/=(const double &a) {
    u /= a;
    v /= a;
    return *this;
  }

  //! \brief Adds two vectors.
  inline Vector2 operator+(const Vector2 &other) const {
    return Vector2(u + other.u, v + other.v);
  }

  //! \brief Substracts two vectors.
  inline Vector2 operator-(const Vector2 &other) const {
    return Vector2(u - other.u, v - other.v);
  }

  //! \brief Scales a vector.
  inline Vector2 operator*(const double &a) const {
    return Vector2(u * a, v * a);
  }

  //! \brief Scales a vector.
  inline Vector2 operator/(const double &a) const {
    return Vector2(u / a, v / a);
  }

  double u, v;
};

// Multiplication of a vector by a constant is commutative.
inline Vector2 operator*(const double &a, const Vector2 &v) {
  return v * a;
}

} // end of namespace pism

#endif /* _VECTOR2_H_ */
