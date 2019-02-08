/* Copyright (C) 2013, 2014, 2015, 2016, 2017, 2018, 2019 PISM Authors
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
#ifndef _PISMEIGENCALVING_H_
#define _PISMEIGENCALVING_H_

#include "StressCalving.hh"

namespace pism {

namespace calving {

class EigenCalving : public StressCalving {
public:
  EigenCalving(IceGrid::ConstPtr grid);
  virtual ~EigenCalving();

  void init();

  void update(double dt,
              const Geometry &geometry,
              const IceModelVec2Int &bc_mask,
              const IceModelVec2V &ice_velocity,
              IceModelVec2CellType &cell_type,
              IceModelVec2S &Href,
              IceModelVec2S &ice_thickness);

protected:
  virtual DiagnosticList diagnostics_impl() const;

  void compute_retreat_rate(const IceModelVec2CellType &cell_type,
                            const IceModelVec2V &ice_velocity,
                            IceModelVec2S &result) const;

  void compute_retreat_rate(const FrontRetreatInputs &inputs,
                            IceModelVec2S &result) const;

protected:
  double m_K;
};

} // end of namespace calving
} // end of namespace pism

#endif /* _PISMEIGENCALVING_H_ */
