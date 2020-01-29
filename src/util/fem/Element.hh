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

#include "FEM.hh"
#include "pism/util/Vector2.hh"
#include "pism/util/petscwrappers/Mat.hh" // Mat, MatStencil

namespace pism {

class IceGrid;
class IceModelVec2Int;

namespace fem {

class Quadrature;

//! The mapping from global to local degrees of freedom.
/*! Computations of residual and Jacobian entries in the finite element method are done by iterating
  of the elements and adding the various contributions from each element. To do this, degrees of
  freedom from the global vector of unknowns must be remapped to element-local degrees of freedom
  for the purposes of local computation, (and the results added back again to the global residual
  and Jacobian arrays).

  An Element mediates the transfer between element-local and global degrees of freedom. In this
  very concrete implementation, the global degrees of freedom are either scalars (double's) or
  vectors (Vector2's), one per node in the IceGrid, and the local degrees of freedom on the element
  are q1::n_chi (%i.e. four) scalars or vectors, one for each vertex of the element.

  The Element is also (perhaps awkwardly) overloaded to mediate transfering locally computed
  contributions to residual and Jacobian matrices to their global counterparts.

  See also: \link FETools.hh FiniteElement/IceGrid background material\endlink.
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
  Germ chi(unsigned int q, unsigned int k) const {
    return m_germs[q][k];
  }

  //! Number of quadrature points
  int n_pts() const {
    return m_weights.size();
  }

  double weight(unsigned int q) const {
    return m_weights[q];
  }

  int n_sides() const {
    return n_chi();
  }

  Vector2 normal(int side) const {
    assert(side < m_n_chi);
    return m_normals[side];
  }

  double side_length(int side) const {
    return m_side_lengths[side];
  }

  /*! @brief Given nodal values, compute the values at quadrature points.*/
  //! The output array `result` should have enough elements to hold values at all
  //! quadrature points.
  template <typename T>
  void evaluate(const T *x, T *result) {
    for (unsigned int q = 0; q < m_Nq; q++) {
      result[q] = 0.0;
      for (unsigned int k = 0; k < q1::n_chi; k++) { // FIXME: q1::n_chi hard-wired
        result[q] += m_germs[q][k].val * x[k];
      }
    }
  }

  /*! @brief Given nodal values, compute the values and partial derivatives at the
   *  quadrature points.*/
  //! Output arrays should have enough elements to hold values at all quadrature points.`
  template <typename T>
  void evaluate(const T *x, T *vals, T *dx, T *dy) {
    for (unsigned int q = 0; q < m_Nq; q++) {
      vals[q] = 0.0;
      dx[q]   = 0.0;
      dy[q]   = 0.0;
      for (unsigned int k = 0; k < q1::n_chi; k++) { // FIXME: q1::n_chi hard-wired
        const Germ &psi = m_germs[q][k];
        vals[q] += psi.val * x[k];
        dx[q]   += psi.dx  * x[k];
        dy[q]   += psi.dy  * x[k];
      }
    }
  }

  /*! @brief Extract nodal values for the element (`i`,`j`) from global array `x_global`
    into the element-local array `result`.
  */
  template<typename T>
  void nodal_values(T const* const* x_global, T* result) const {
    for (int k = 0; k < m_n_chi; ++k) {
      int i = 0, j = 0;
      local_to_global(k, i, j);
      result[k] = x_global[j][i];   // note the indexing order
    }
  }

  /*! @brief Extract nodal values for the element (`i`,`j`) from global IceModelVec `x_global`
    into the element-local array `result`.
  */
  template<class C, typename T>
  void nodal_values(const C& x_global, T* result) const {
    for (int k = 0; k < m_n_chi; ++k) {
      int i = 0, j = 0;
      local_to_global(k, i, j);
      result[k] = x_global(i, j);   // note the indexing order
    }
  }

  /*! @brief Get nodal values of an integer mask. */
  void nodal_values(const IceModelVec2Int &x_global, int *result) const;

  /*! @brief Add the values of element-local contributions `y` to the global vector `y_global`. */
  /*! The element-local array `local` should be an array of Nk values.
   *
   * Use this to add residual contributions.
   */
  template<typename T>
  void add_contribution(const T *local, T** y_global) const {
    for (int k = 0; k < m_n_chi; k++) {
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

  template<class C, typename T>
  void add_contribution(const T *local, C& y_global) const {
    for (int k = 0; k < m_n_chi; k++) {
      if (m_row[k].i == m_invalid_dof) {
        // skip rows marked as "invalid"
        continue;
      }
      const int
        i = m_row[k].i,
        j = m_row[k].j;
      y_global(i, j) += local[k];   // note the indexing order
    }
  }

  void reset(int i, int j);

  /*! @brief Add Jacobian contributions. */
  void add_contribution(const double *K, Mat J) const;

  void mark_row_invalid(int k);
  void mark_col_invalid(int k);

  //! Convert a local degree of freedom index `k` to a global degree of freedom index (`i`,`j`).
  void local_to_global(int k, int &i, int &j) const {
    i = m_i + m_i_offset[k];
    j = m_j + m_j_offset[k];
  }

protected:
  Element(const IceGrid &grid, const Quadrature &q);

  // pointer to a 2D shape function
  typedef Germ (*ShapeFunction2)(unsigned int k, const QuadPoint &p);

  void initialize(ShapeFunction2 f,
                  unsigned int n_chi,
                  const std::vector<QuadPoint>& points,
                  const std::vector<double>& W);

  //! grid offsets used to extract nodal values from the grid and add contributions from
  //! an element to the residual and the Jacobian.
  std::vector<int> m_i_offset;
  std::vector<int> m_j_offset;

  //! Number of nodes (and therefore the number of shape functions) in this particular
  //! type of elements.
  int m_n_chi;

  std::vector<Vector2> m_normals;

  std::vector<double> m_side_lengths;

  //! Jacobian of the map from the reference element to a physical element
  double m_J[2][2];

private:
  Element(const IceGrid &grid) : m_Nq(0), m_grid(grid) {
    // empty
  }

  //! Constant for marking invalid row/columns.
  //!
  //! Has to be negative because MatSetValuesBlockedStencil supposedly ignores negative indexes.
  //! Seems like it has to be negative and below the smallest allowed index, i.e. -2 and below with
  //! the stencil of width 1 (-1 *is* an allowed index). We use -2^30 and *don't* use PETSC_MIN_INT,
  //! because PETSC_MIN_INT depends on PETSc's configuration flags.
  static const int m_invalid_dof = -1073741824;

  // maximum number of shape functions supported by this class
  static const int m_n_chi_max = q1::n_chi;

  //! Indices of the current element.
  int m_i, m_j;

  //! Number of quadrature points
  const unsigned int m_Nq;

  //! Stencils for updating entries of the Jacobian matrix.
  MatStencil m_row[m_n_chi_max], m_col[m_n_chi_max];

  // the grid (for marking ghost DOFs as "invalid")
  const IceGrid &m_grid;

  //! Quadrature weights
  std::vector<double> m_weights;

  // define Germs, which is an array of q1::n_chi "Germ"s
  typedef Germ Germs[q1::n_chi];

  Germs* m_germs;
};

//! Q1 element with sides parallel to X and Y axes
class Q1Element : public Element {
public:
  Q1Element(const IceGrid &grid, const Quadrature &quadrature);
};

//! P1 element embedded in Q1Element
class P1Element : public Element {
public:
  P1Element(const IceGrid &grid, const Quadrature &quadrature, int N);
};


} // end of namespace fem
} // end of namespace pism

#endif /* PISM_ELEMENT_H */
