/* Copyright (C) 2018, 2020 PISM Authors
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

// Utilities for serial models using FFTW and extended computational grids.

#include <vector>
#include <complex>

#include <fftw3.h>

namespace pism {

namespace petsc {
class Vec;
} // end of namespace petsc

/*!
 * Template class for accessing the central part of an extended grid, i.e. PISM's grid
 * surrounded by "padding" necessary to reduce artifacts coming from interpreting model
 * inputs as periodic in space.
 *
 * Allows accessing a 1D array using 2D indexing.
 */
class FFTWArray {
public:
  // Access the central part of an array "a" of size My*Mx, using offsets i_offset and
  // j_offset specifying the corner of the part to be accessed.
  FFTWArray(fftw_complex *a, int Mx, int My, int i_offset, int j_offset);

  // Access the whole array "a" of size My*My.
  FFTWArray(fftw_complex *a, int Mx, int My);

  inline std::complex<double> &operator()(int i, int j) {
    return m_array[(m_j_offset + j) + m_My * (m_i_offset + i)];
  }

private:
  const int m_My, m_i_offset, m_j_offset;
  std::complex<double> *m_array;
};

std::vector<double> fftfreq(int M, double normalization);

//! Fill `input` with zeros.
void clear_fftw_array(fftw_complex *input, int Nx, int Ny);

//! Copy `source` to `destination`.
void copy_fftw_array(fftw_complex *source, fftw_complex *destination, int Nx, int Ny);

//! Set the real part of output to input. Input has the size of My*Mx, embedded in the
//! bigger (output) grid of size Ny*Nx. Offsets i0 and j0 specify the location of the
//! subset to set.
/*!
 * Sets the imaginary part to zero.
 */
void set_real_part(petsc::Vec &input,
                   double normalization,
                   int Mx, int My,
                   int Nx, int Ny,
                   int i0, int j0,
                   fftw_complex *output);

//! \brief Get the real part of input and put it in output.
/*!
 * See set_real_part for details.
 */
void get_real_part(fftw_complex *input,
                   double normalization,
                   int Mx, int My,
                   int Nx, int Ny,
                   int i0, int j0,
                   petsc::Vec &output);

} // end of namespace pism
