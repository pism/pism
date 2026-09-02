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

#ifndef PISM_DEBRIS_ICE_MELT_ENHANCEMENT_HH
#define PISM_DEBRIS_ICE_MELT_ENHANCEMENT_HH

#include <memory>

#include "pism/util/Component.hh"
#include "pism/util/array/Scalar.hh"

namespace pism {
namespace debris {

class DebrisModel;

//! @brief The effect of a supraglacial debris cover on the ice melt rate.
/*!
  Supplies a dimensionless factor multiplying the clean-ice melt rate: values
  greater than one mean that a thin, dispersed debris cover *enhances* melt,
  values less than one that a thicker layer insulates the ice.

  The parameterization is selected using `debris.ice_melt_enhancement.model`.
*/
class IceMeltEnhancement : public Component {
public:
  //! Available parameterizations.
  enum Model : int {
    //! No debris effect: the factor is 1 everywhere.
    NONE = 0,
    //! Read the factor from a file (variable `debris_melt_factor`).
    GIVEN = 1,
    //! Compute the factor from the debris thickness supplied by a DebrisModel
    //! using equation (14) of Verhaegen and Huybrechts (2026).
    VERHAEGEN = 2
  };

  /*!
   * `debris_model` is required by (and used only by) the "verhaegen" model,
   * which takes ownership of it: init() and update() are forwarded to it.
   */
  IceMeltEnhancement(std::shared_ptr<const Grid> grid);
  IceMeltEnhancement(std::shared_ptr<const Grid> grid,
                     std::shared_ptr<DebrisModel> debris_model);
  virtual ~IceMeltEnhancement() = default;

  void init(const Geometry &geometry);
  void update(const Geometry &geometry, double t, double dt);

  //! @brief The ratio of the sub-debris melt rate to the clean ice melt rate.
  const array::Scalar &ice_melt_enhancement() const;

  Model model() const;

  //! @brief Equation (14) of Verhaegen and Huybrechts (2026).
  static double verhaegen_melt_factor(double debris_thickness);

protected:
  MaxTimestep max_timestep_impl(double t, const CFLData *cfl_data) const;
  DiagnosticList spatial_diagnostics_impl() const;

  Model m_model;

  array::Scalar m_ice_melt_enhancement;

  //! forcing used by the "given" model
  std::shared_ptr<array::Forcing> m_melt_factor;

  //! debris thickness used by the "verhaegen" model
  std::shared_ptr<DebrisModel> m_debris_model;
};

} // end of namespace debris
} // end of namespace pism

#endif /* PISM_DEBRIS_ICE_MELT_ENHANCEMENT_HH */
