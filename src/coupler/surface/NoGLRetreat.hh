/* Copyright (C) 2021 PISM Authors
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

#ifndef PISM_NO_GL_RETREAT_HH
#define PISM_NO_GL_RETREAT_HH

#include "pism/coupler/SurfaceModel.hh"

namespace pism {
namespace surface {

class NoGLRetreat : public SurfaceModel {
public:
  NoGLRetreat(IceGrid::ConstPtr g, std::shared_ptr<SurfaceModel> input);
  virtual ~NoGLRetreat() = default;

  const IceModelVec2S& smb_adjustment() const;
protected:
  void init_impl(const Geometry &geometry);
  void update_impl(const Geometry &geometry, double t, double dt);

  const IceModelVec2S& mass_flux_impl() const;
  const IceModelVec2S& accumulation_impl() const;
  const IceModelVec2S& melt_impl() const;
  const IceModelVec2S& runoff_impl() const;

  DiagnosticList diagnostics_impl() const;
private:
  IceModelVec2S::Ptr m_mass_flux;
  IceModelVec2S m_smb_adjustment;
  IceModelVec2S m_min_ice_thickness;
};

} // end of namespace surface
} // end of namespace pism

#endif /* PISM_NO_GL_RETREAT_HH */
