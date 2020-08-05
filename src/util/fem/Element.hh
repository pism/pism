/* Copyright (C) 2020 PISM Authors
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
#ifndef PISM_ELEMENT_H
#define PISM_ELEMENT_H

#include <vector>
#include <cassert>              // assert

#include <petscdmda.h>          // DMDALocalInfo

#include "FEM.hh"
#include "pism/util/Vector2.hh"
#include "pism/util/petscwrappers/Mat.hh" // Mat, MatStencil


namespace pism {

class IceGrid;
class IceModelVec2Int;

struct Vector3 {
  double x, y, z;
};

namespace fem {

class Quadrature;

//! The mapping from global to local degrees of freedom.
/*! Computations of residual and Jacobian entries in the finite element method are done by
  iterating of the elements and adding the various contributions from each element. To do
  this, degrees of freedom from the global vector of unknowns must be remapped to
  element-local degrees of freedom for the purposes of local computation, (and the results
  added back again to the global residual and Jacobian arrays).

  An Element mediates the transfer between element-local and global degrees of freedom and
  provides values of shape functions at quadrature points. In this implementation, the
  global degrees of freedom are either scalars (double's) or vectors (Vector2's), one per
  node in the IceGrid, and the local degrees of freedom on the element are q1::n_chi or
  p1::n_chi (%i.e. four or three) scalars or vectors, one for each vertex of the element.

  In addition to this, the Element transfers locally computed contributions to residual
  and Jacobian matrices to their global counterparts.
*/
class Element {
public:
  ~Element();

  int n_chi() const {
    return m_n_chi;
  }

  /*!
   * `chi(q, k)` returns values and partial derivatives of the `k`-th shape function at a
   * quadrature point `q`.
   */
  const Germ& chi(unsigned int q, unsigned int k) const {
    assert(q < m_Nq);
    assert(k < m_n_chi);
    return m_germs[q * m_n_chi + k];
  }

  //! Number of quadrature points
  int n_pts() const {
    return m_Nq;
  }

  //! Weight of the quadrature point `q`
  double weight(unsigned int q) const {
    assert(q < m_Nq);
    return m_weights[q];
  }

  /*! @brief Given nodal values, compute the values at quadrature points.*/
  //! The output array `result` should have enough elements to hold values at all
  //! quadrature points.
  template <typename T>
  void evaluate(const T *x, T *result) const {
    for (unsigned int q = 0; q < m_Nq; q++) {
      result[q] = 0.0;
      for (unsigned int k = 0; k < m_n_chi; k++) {
        result[q] += m_germs[q * m_n_chi + k].val * x[k];
      }
    }
  }

protected:
  Element(const IceGrid &grid, int Nq, int n_chi, int block_size);
  Element(const DMDALocalInfo &grid_info, int Nq, int n_chi, int block_size);

  /*! @brief Add Jacobian contributions. */
  void add_contribution(const double *K, Mat J) const;

  void mark_row_invalid(int k);
  void mark_col_invalid(int k);

  DMDALocalInfo m_grid;

  // pointer to a shape function
  typedef Germ (*ShapeFunction)(unsigned int k, const QuadPoint &p);

  void initialize(const double J[3][3],
                  ShapeFunction f,
                  unsigned int n_chi,
                  const std::vector<QuadPoint>& points,
                  const std::vector<double>& W);

  //! grid offsets used to extract nodal values from the grid and add contributions from
  //! an element to the residual and the Jacobian.
  std::vector<int> m_i_offset;
  std::vector<int> m_j_offset;

  //! Number of nodes (and therefore the number of shape functions) in this particular
  //! type of elements.
  const unsigned int m_n_chi;

  //! Indices of the current element.
  int m_i, m_j;

  //! Number of quadrature points
  const unsigned int m_Nq;

  const int m_block_size;

  std::vector<Germ> m_germs;

  //! Stencils for updating entries of the Jacobian matrix.
  std::vector<MatStencil> m_row, m_col;

  //! Constant for marking invalid row/columns.
  //!
  //! Has to be negative because MatSetValuesBlockedStencil supposedly ignores negative indexes.
  //! Seems like it has to be negative and below the smallest allowed index, i.e. -2 and below with
  //! the stencil of width 1 (-1 *is* an allowed index). We use -2^30 and *don't* use PETSC_MIN_INT,
  //! because PETSC_MIN_INT depends on PETSc's configuration flags.
  static const int m_invalid_dof = -1073741824;

  //! Quadrature weights for a particular physical element
  std::vector<double> m_weights;
private:
  Element()
    : m_n_chi(0),
      m_Nq(0),
      m_block_size(0) {
    // empty
  }
};

class Element2 : public Element {
public:

  void reset(int i, int j);

  //! Convert a local degree of freedom index `k` to a global degree of freedom index (`i`,`j`).
  void local_to_global(unsigned int k, int &i, int &j) const {
    i = m_i + m_i_offset[k];
    j = m_j + m_j_offset[k];
  }

  Vector2 normal(unsigned int side) const {
    assert(side < m_n_chi);
    return m_normals[side];
  }

  unsigned int n_sides() const {
    return n_chi();
  }

  double side_length(unsigned int side) const {
    assert(side < m_n_chi);
    return m_side_lengths[side];
  }

  using Element::evaluate;
  /*! @brief Given nodal values, compute the values and partial derivatives at the
   *  quadrature points.*/
  //! Output arrays should have enough elements to hold values at all quadrature points.`
  template <typename T>
  void evaluate(const T *x, T *vals, T *dx, T *dy) {
    for (unsigned int q = 0; q < m_Nq; q++) {
      vals[q] = 0.0;
      dx[q]   = 0.0;
      dy[q]   = 0.0;
      for (unsigned int k = 0; k < m_n_chi; k++) {
        const Germ &psi = m_germs[q * m_n_chi + k];
        vals[q] += psi.val * x[k];
        dx[q]   += psi.dx  * x[k];
        dy[q]   += psi.dy  * x[k];
      }
    }
  }

  /*! @brief Get nodal values of an integer mask. */
  void nodal_values(const IceModelVec2Int &x_global, int *result) const;

  /*! @brief Extract nodal values for the element (`i`,`j`) from global array `x_global`
    into the element-local array `result`.
  */
  template<typename T>
  void nodal_values(T const* const* x_global, T* result) const {
    for (unsigned int k = 0; k < m_n_chi; ++k) {
      int i = 0, j = 0;
      local_to_global(k, i, j);
      result[k] = x_global[j][i];   // note the indexing order
    }
  }

  using Element::add_contribution;
  /*! @brief Add the values of element-local contributions `y` to the global vector `y_global`. */
  /*! The element-local array `local` should be an array of Nk values.
   *
   * Use this to add residual contributions.
   */
  template<typename T>
  void add_contribution(const T *local, T** y_global) const {
    for (unsigned int k = 0; k < m_n_chi; k++) {
      if (m_row[k].i == m_invalid_dof) {
        // skip rows marked as "invalid"
        continue;
      }
      const int
        i = m_row[k].i,
        j = m_row[k].j;
      y_global[j][i] += local[k];   // note the indexing order
    }
  }

  using Element::mark_row_invalid;
  using Element::mark_col_invalid;

protected:
  Element2(const IceGrid &grid, int Nq, int n_chi, int block_size);
  Element2(const DMDALocalInfo &grid_info, int Nq, int n_chi, int block_size);

  std::vector<Vector2> m_normals;

  std::vector<double> m_side_lengths;
};

//! Q1 element with sides parallel to X and Y axes
class Q1Element2 : public Element2 {
public:
  Q1Element2(const IceGrid &grid, const Quadrature &quadrature);
  Q1Element2(const DMDALocalInfo &grid_info,
             double dx, double dy,
             const Quadrature &quadrature);
};

//! P1 element embedded in Q1Element2
class P1Element2 : public Element2 {
public:
  P1Element2(const IceGrid &grid, const Quadrature &quadrature, int N);
};

class Element3 : public Element {
public:
  using Element::evaluate;
  /*! @brief Given nodal values, compute the values and partial derivatives at the
   *  quadrature points.*/
  //! Output arrays should have enough elements to hold values at all quadrature points.`
  template <typename T>
  void evaluate(const T *x, T *vals, T *dx, T *dy, T *dz) const {
    for (unsigned int q = 0; q < m_Nq; q++) {
      vals[q] = 0.0;
      dx[q]   = 0.0;
      dy[q]   = 0.0;
      dz[q]   = 0.0;
      for (unsigned int k = 0; k < m_n_chi; k++) {
        const Germ &psi = m_germs[q * m_n_chi + k];
        vals[q] += psi.val * x[k];
        dx[q]   += psi.dx  * x[k];
        dy[q]   += psi.dy  * x[k];
        dz[q]   += psi.dz  * x[k];
      }
    }
  }

  struct GlobalIndex {
    int i, j, k;
  };

  GlobalIndex local_to_global(int i, int j, int k, int n) const {
    return {i + m_i_offset[n],
            j + m_j_offset[n],
            k + m_k_offset[n]};
  }

  GlobalIndex local_to_global(int n) const {
    return {m_i + m_i_offset[n],
            m_j + m_j_offset[n],
            m_k + m_k_offset[n]};
  }

  /*! @brief Extract nodal values for the element (`i`,`j`,`k`) from global array `x_global`
    into the element-local array `result`.
  */
  template<typename T>
  void nodal_values(T const* const* const* x_global, T* result) const {
    for (unsigned int n = 0; n < m_n_chi; ++n) {
      auto I = local_to_global(n);
      result[n] = x_global[I.j][I.i][I.k]; // note the STORAGE_ORDER
    }
  }

  using Element::add_contribution;
  /*! @brief Add the values of element-local contributions `y` to the global vector `y_global`. */
  /*! The element-local array `local` should be an array of Nk values.
   *
   * Use this to add residual contributions.
   */
  template<typename T>
  void add_contribution(const T *local, T*** y_global) const {
    for (unsigned int n = 0; n < m_n_chi; n++) {
      if (m_row[n].i == m_invalid_dof) {
        // skip rows marked as "invalid"
        continue;
      }
      auto I = local_to_global(n);
      y_global[I.j][I.i][I.k] += local[n];   // note the STORAGE_ORDER
    }
  }
protected:
  Element3(const DMDALocalInfo &grid_info, int Nq, int n_chi, int block_size);
  Element3(const IceGrid &grid, int Nq, int n_chi, int block_size);

  std::vector<int> m_k_offset;

  int m_k;
};

//! @brief 3D Q1 element
/* Regular grid in the x and y directions, scaled in the z direction.
 *
 */
class Q1Element3 : public Element3 {
public:
  Q1Element3(const DMDALocalInfo &grid, double dx, double dy, const Quadrature &quadrature);
  Q1Element3(const IceGrid &grid, const Quadrature &quadrature);

  void reset(int i, int j, int k, const double *z);

  using Element::mark_row_invalid;
  using Element::mark_col_invalid;
private:
  double m_dx;
  double m_dy;

  // values of shape functions and their derivatives with respect to xi,eta,zeta at all
  // quadrature points
  std::vector<Germ> m_chi;

  // quadrature points (on the reference element)
  std::vector<QuadPoint> m_points;

  // quadrature weights (on the reference element)
  std::vector<double> m_w;
};


class Q1Element3Face {
public:
  Q1Element3Face(double dx, double dy, const Quadrature &quadrature);

  //! Number of quadrature points
  int n_pts() const {
    return m_Nq;
  }
  double weight(unsigned int q) const {
    assert(q < m_Nq);
    return m_weights[q];
  }

  const double& chi(unsigned int q, unsigned int k) const {
    assert(q < m_Nq);
    assert(k < m_n_chi);
    return m_chi[q * m_n_chi + k];
  }

  // outward pointing unit normal vector to an element face at the quadrature point q
  const Vector3& normal(unsigned int q) const {
    return m_normals[q];
  }

  // NB: here z contains nodal z coordinates for *all* nodes of the element
  void reset(int face, const double *z);

  template <typename T>
  void evaluate(const T *x, T *result) const {
    for (unsigned int q = 0; q < m_Nq; q++) {
      result[q] = 0.0;
      for (unsigned int k = 0; k < m_n_chi; k++) {
        result[q] += m_chi[q * m_n_chi + k] * x[k];
      }
    }
  }
protected:
  // grid spacing
  double m_dx;
  double m_dy;

  // values of shape functions at all quadrature points and their derivatives with respect
  // to xi, eta, zeta
  std::vector<double> m_chi;

  // quadrature points (on the reference element corresponding to a face)
  std::vector<QuadPoint> m_points;

  std::vector<Vector3> m_normals;

  // quadrature weights (on the reference element corresponding to a face)
  std::vector<double> m_w;

  // quadrature weights (on the face of a physical element)
  std::vector<double> m_weights;

  // number of quadrature points
  const unsigned int m_n_chi;
  const unsigned int m_Nq;
};

} // end of namespace fem
} // end of namespace pism

#endif /* PISM_ELEMENT_H */
