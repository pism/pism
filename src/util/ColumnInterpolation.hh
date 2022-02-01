/* Copyright (C) 2014, 2015, 2022 PISM Authors
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

#ifndef _COLUMNINTERPOLATION_H_
#define _COLUMNINTERPOLATION_H_

#include <vector>

namespace pism {

class ColumnInterpolation {
public:
  ColumnInterpolation(const std::vector<double> &z_coarse,
                      const std::vector<double> &z_fine);

  void coarse_to_fine(const double *input, unsigned int ks, double *result) const;
  void fine_to_coarse(const double *input, double *result) const;

  // These two methods allocate fresh storage for the output.
  std::vector<double> coarse_to_fine(const std::vector<double> &input, unsigned int ks) const;
  std::vector<double> fine_to_coarse(const std::vector<double> &input) const;

  unsigned int Mz_coarse() const;
  const std::vector<double>& z_coarse() const;

  unsigned int Mz_fine() const;
  double dz_fine() const;
  const std::vector<double>& z_fine() const;
private:
  std::vector<double> m_z_fine, m_z_coarse;
  std::vector<double> m_constants;

  // Array m_coarse2fine contains indices of the ice coarse vertical grid
  // that are just below a level of the fine grid. I.e. m_coarse2fine[k] is
  // the coarse grid level just below fine-grid level k (zlevels_fine[k]).
  // Similarly for other arrays below.
  std::vector<unsigned int> m_coarse2fine, m_fine2coarse;
  bool m_use_linear_interpolation;
  bool m_identical_grids;

  void init_interpolation();
  void coarse_to_fine_linear(const double *input, unsigned int ks, double *result) const;
  void coarse_to_fine_quadratic(const double *input, unsigned int ks, double *result) const;
};

} // end of namespace pism

#endif /* _COLUMNINTERPOLATION_H_ */
