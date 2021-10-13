/* Copyright (C) 2016, 2017, 2021 PISM Authors
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

#ifndef BTU_VERIFICATION_H
#define BTU_VERIFICATION_H

#include "pism/energy/BTU_Full.hh"

namespace pism {
namespace energy {

class BTU_Verification : public BTU_Full
{
public:
  BTU_Verification(IceGrid::ConstPtr g,
                   const BTUGrid &vertical_grid,
                   int test, bool bii);
  virtual ~BTU_Verification() = default;

protected:
  virtual void initialize_bottom_surface_flux();
  virtual void bootstrap(const IceModelVec2S &bedrock_top_temperature);
  int m_testname;
  bool m_bedrock_is_ice;
};

} // end of namespace energy
} // end of namespace pism


#endif /* BTU_VERIFICATION_H */
