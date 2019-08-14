/* Copyright (C) 2015, 2017, 2018, 2019 PISM Authors
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

#ifndef _REGIONAL_YIELD_STRESS_H_
#define _REGIONAL_YIELD_STRESS_H_

#include "pism/basalstrength/YieldStress.hh"

namespace pism {

class RegionalYieldStress : public YieldStress {
public:
  RegionalYieldStress(std::shared_ptr<YieldStress> input);
  virtual ~RegionalYieldStress();
protected:
  virtual void init_impl(const YieldStressInputs &inputs);
  virtual void update_impl(const YieldStressInputs &inputs, double t, double dt);
private:
  std::shared_ptr<YieldStress> m_input;
};

} // end of namespace pism

#endif /* _REGIONAL_YIELD_STRESS_H_ */
