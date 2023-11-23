// Copyright (C) 2018, 2019, 2021, 2022, 2023 Constantine Khroulev and Andy Aschwanden
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

#ifndef __PISMFrontalMelt_hh
#define __PISMFrontalMelt_hh

#include <memory>

#include "pism/util/Component.hh"

namespace pism {

class FrontalMeltInputs {
public:
  FrontalMeltInputs();

  const Geometry *geometry;

  // used by DischargeRouting
  const array::Scalar *subglacial_water_flux;

};

//! @brief Frontal melt models and modifiers.
namespace frontalmelt {

//! A very rudimentary PISM frontal melt model.
class FrontalMelt : public Component {
public:
  // "modifier" constructor
  FrontalMelt(std::shared_ptr<const Grid> g, std::shared_ptr<FrontalMelt> input);
  // "model" constructor
  FrontalMelt(std::shared_ptr<const Grid> g);

  virtual ~FrontalMelt() = default;

  void init(const Geometry &geometry);

  void update(const FrontalMeltInputs &inputs, double t, double dt);

  const array::Scalar& frontal_melt_rate() const;

  const array::Scalar& retreat_rate() const;

protected:
  virtual void init_impl(const Geometry &geometry);

  // provides default (pass-through) implementations for "modifiers"
  virtual void update_impl(const FrontalMeltInputs &inputs, double t, double dt);
  virtual MaxTimestep max_timestep_impl(double t) const;
  virtual void define_model_state_impl(const File &output) const;
  virtual void write_model_state_impl(const File &output) const;

  virtual DiagnosticList diagnostics_impl() const;
  virtual TSDiagnosticList ts_diagnostics_impl() const;

  virtual const array::Scalar& frontal_melt_rate_impl() const = 0;

  void compute_retreat_rate(const Geometry &geometry, const array::Scalar &frontal_melt_rate,
                            array::Scalar &result) const;

  std::shared_ptr<FrontalMelt> m_input_model;

  bool apply(const array::CellType1 &M, int i, int j) const;

  array::Scalar m_retreat_rate;
  bool m_include_floating_ice;
};

} // end of namespace frontalmelt
} // end of namespace pism

#endif  // __PISMFrontalMelt_hh
