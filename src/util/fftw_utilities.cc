/* Copyright (C) 2018, 2019, 2020, 2021, 2023 PISM Authors
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

#include <cstring>              // memcpy

#include "pism/util/fftw_utilities.hh"

#include "pism/util/petscwrappers/Vec.hh"

namespace pism {

  // Access the central part of an array "a" of size My*Mx, using offsets i_offset and
  // j_offset specifying the corner of the part to be accessed.
FFTWArray::FFTWArray(fftw_complex *a, int Mx, int My, int i_offset, int j_offset)
  : m_My(My), m_i_offset(i_offset), m_j_offset(j_offset),
    m_array(reinterpret_cast<std::complex<double>*>(a)) {
  (void) Mx;
}

  // Access the whole array "a" of size My*My.
FFTWArray::FFTWArray(fftw_complex *a, int Mx, int My)
  : m_My(My), m_i_offset(0), m_j_offset(0),
    m_array(reinterpret_cast<std::complex<double>*>(a)) {
  (void) Mx;
}


/*!
 * Return the Discrete Fourier Transform sample frequencies.
 */
std::vector<double> fftfreq(int M, double normalization) {
  std::vector<double> result(M);

  int N = (M % 2) != 0 ? M / 2 : M / 2 - 1;

  for (int i = 0; i <= N; i++) {
    result[i] = i;
  }

  for (int i = N + 1; i < M; i++) {
    result[i] = i - M;
  }

  // normalize
  for (int i = 0; i < M; i++) {
    result[i] /= (M * normalization);
  }

  return result;
}

//! \brief Fill `input` with zeros.
void clear_fftw_array(fftw_complex *input, int Nx, int Ny) {
  FFTWArray fftw_in(input, Nx, Ny);
  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      fftw_in(i, j) = 0.0;
    }
  }
}

//! @brief Copy `source` to `destination`.
void copy_fftw_array(fftw_complex *source, fftw_complex *destination, int Nx, int Ny) {
  memcpy(destination, source, Nx * Ny * sizeof(fftw_complex));
}

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
                   fftw_complex *output) {
  petsc::VecArray2D in(input, Mx, My);
  FFTWArray out(output, Nx, Ny, i0, j0);

  for (int j = 0; j < My; ++j) {
    for (int i = 0; i < Mx; ++i) {
      out(i, j) = in(i, j) * normalization;
    }
  }
}


//! @brief Get the real part of input and put it in output.
/*!
 * See set_real_part for details.
 */
void get_real_part(fftw_complex *input,
                   double normalization,
                   int Mx, int My,
                   int Nx, int Ny,
                   int i0, int j0,
                   petsc::Vec &output) {
  petsc::VecArray2D out(output, Mx, My);
  FFTWArray in(input, Nx, Ny, i0, j0);
  for (int j = 0; j < My; ++j) {
    for (int i = 0; i < Mx; ++i) {
      out(i, j) = in(i, j).real() * normalization;
    }
  }
}

} // end of namespace pism
