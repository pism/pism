// Copyright (C) 2026 Constantine Khroulev and Andy Aschwanden
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

#ifndef PISM_DEBRIS_GIVEN_HH
#define PISM_DEBRIS_GIVEN_HH

#include "pism/coupler/DebrisModel.hh"

namespace pism {
namespace debris {

class Given : public DebrisModel {
public:
  Given(std::shared_ptr<const Grid>);
  virtual ~Given() = default;
protected:
  void init_impl(const Geometry &geometry);
  void update_impl(const Inputs &inputs, double t, double dt);

  const array::Scalar &debris_impl() const;

  std::set<VariableMetadata> state_impl() const;
  void write_state_impl(const OutputFile &output) const;

  std::shared_ptr<array::Forcing> m_debris_thickness;
};

} // end of namespace debris
} // end of namespace pism

#endif /* PISM_DEBRIS_GIVEN_HH */
