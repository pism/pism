/* Copyright (C) 2016, 2017, 2018, 2019 PISM Authors
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

#ifndef HUYHURSTCALVING_H
#define HUYHURSTCALVING_H

#include "StressCalving.hh"

namespace pism {

class Geometry;

namespace calving {

class HuyhurstCalving : public StressCalving {
public:
  HuyhurstCalving(IceGrid::ConstPtr grid);
  virtual ~HuyhurstCalving();

  void init();

  void update(const IceModelVec2CellType &pism_mask, const IceModelVec2S &ice_thickness, 
              const IceModelVec2S &sealevel, const IceModelVec2S &surface);

protected:
  DiagnosticList diagnostics_impl() const;
  
protected:
  double m_B_tilde, m_exponent_r, m_sigma_threshold;

};

} // end of namespace calving
} // end of namespace pism

#endif /* HUYHURSTCALVING_H */
