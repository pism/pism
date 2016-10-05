/* Copyright (C) 2013, 2014, 2015, 2016 PISM Authors
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

#ifndef _PO_DELTA_MBP_H_
#define _PO_DELTA_MBP_H_

#include "coupler/util/PScalarForcing.hh"
#include "POModifier.hh"

namespace pism {
namespace ocean {

/**
 * Scalar melange back-pressure fraction forcing.
 * 
 */
class Delta_MBP : public PScalarForcing<OceanModel,OceanModifier>
{
public:
  Delta_MBP(IceGrid::ConstPtr g, OceanModel* in);
  virtual ~Delta_MBP();

protected:
  virtual MaxTimestep max_timestep_impl(double t);
  virtual void init_impl();
  virtual void melange_back_pressure_fraction_impl(IceModelVec2S &result) const;
};

} // end of namespace ocean
} // end of namespace pism

#endif /* _PO_DELTA_MBP_H_ */
