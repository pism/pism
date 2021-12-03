/* Copyright (C) 2015, 2017, 2018, 2019, 2021 PISM Authors
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

/*!
 * Regional version of yield stress models. Sets high tauc in "no model" areas.
 *
 * Note: this class has to implement all the virtual methods of `Component` because it has
 * to forward these calls to the model provided to its constructor.
 */
class RegionalYieldStress : public YieldStress {
public:
  RegionalYieldStress(std::shared_ptr<YieldStress> input);
  virtual ~RegionalYieldStress() = default;
private:
  void restart_impl(const File &input_file, int record);

  void bootstrap_impl(const File &input_file, const YieldStressInputs &inputs);

  void init_impl(const YieldStressInputs &inputs);

  void update_impl(const YieldStressInputs &inputs, double t, double dt);

  MaxTimestep max_timestep_impl(double t) const;

  void define_model_state_impl(const File &output) const;

  void write_model_state_impl(const File &output) const;

  DiagnosticList diagnostics_impl() const;

  TSDiagnosticList ts_diagnostics_impl() const;

  std::shared_ptr<YieldStress> m_input;

  double m_high_tauc;
};

} // end of namespace pism

#endif /* _REGIONAL_YIELD_STRESS_H_ */
