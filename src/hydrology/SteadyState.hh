/* Copyright (C) 2019 PISM Authors
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

#ifndef STEADY_STATE_H
#define STEADY_STATE_H

#include "NullTransport.hh"

namespace pism {
namespace hydrology {

class EmptyingProblem;

/*!
 * A version of the "null-transport" hydrology model that adds the steady state water
 * flux.
 */
class SteadyState : public NullTransport {
public:
  SteadyState(std::shared_ptr<const Grid> g);
  virtual ~SteadyState() = default;

protected:
  void initialization_message() const;

  void init_time(const std::string &input_file);

  void init_impl(const array::Scalar &W_till, const array::Scalar &W, const array::Scalar &P);

  void bootstrap_impl(const File &input_file,
                      const array::Scalar &ice_thickness);
  void restart_impl(const File &input_file, int record);

  void update_impl(double t, double dt, const Inputs& inputs);

  std::map<std::string, Diagnostic::Ptr> diagnostics_impl() const;

  MaxTimestep max_timestep_impl(double t) const;
  void define_model_state_impl(const File &output) const;
  void write_model_state_impl(const File &output) const;

  std::shared_ptr<EmptyingProblem> m_emptying_problem;

  //! time of the last water flux update
  double m_t_last;
  //! Update interval in seconds
  double m_update_interval;
  //! Temporal resolution to use when checking whether it's time to update
  double m_t_eps;
  //! Name of the variable used to store the last update time.
  std::string m_time_name;

  //! Times corresponding to records in the input file
  std::vector<double> m_time;
  //! Time bounds corresponding to records in the input file
  std::vector<double> m_time_bounds;

  //! Set to true in bootstrap_impl() if update_impl() has to bootstrap m_Q.
  bool m_bootstrap;
};

} // end of namespace hydrology
} // end of namespace pism

#endif /* STEADY_STATE_H */
