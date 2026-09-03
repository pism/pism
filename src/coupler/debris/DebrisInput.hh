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

#ifndef PISM_DEBRIS_INPUT_HH
#define PISM_DEBRIS_INPUT_HH

#include <memory>

#include "pism/util/array/Scalar.hh"

namespace pism {

class Grid;

namespace array {
class Forcing;
}

namespace debris {

//! @brief Prescribed rate of debris input from the surrounding terrain.
/*!
  Reads `debris_input_rate` (thickness of solid debris per unit time) from
  `debris.transport.input.file`, possibly time-dependent. If no file is given the input is
  zero everywhere.

  Where the debris lands is decided by the caller: in the accumulation zone it is buried
  (equation 16 in Verhaegen and Huybrechts, 2026), in the ablation zone it is added to the
  supraglacial debris layer (equation 18).
*/
class DebrisInput {
public:
  DebrisInput(std::shared_ptr<const Grid> grid);

  void init();
  void update(double t, double dt);

  //! Input rate averaged over the last step (m of solid debris per second).
  const array::Scalar &rate() const;

  bool enabled() const;

private:
  std::shared_ptr<const Grid> m_grid;
  std::shared_ptr<array::Forcing> m_rate;
  array::Scalar m_zero;
  std::string m_filename;
  bool m_periodic;
};

} // end of namespace debris
} // end of namespace pism

#endif /* PISM_DEBRIS_INPUT_HH */
