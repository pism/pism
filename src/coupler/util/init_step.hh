/* Copyright (C) 2017, 2018 PISM Authors
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

#include "pism/util/Time.hh"
#include "pism/util/MaxTimestep.hh"

namespace pism {

/*!
 * Update a `model` by asking it to perform time-stepping from the current time to one year in the
 * future (or as far as the time step restriction allows).
 *
 * This is sometimes necessary during initialization, but should be avoided if possible.
 */
template<class M>
void init_step(M *model, const Geometry &geometry, const Time& time) {
  const double
    now               = time.current(),
    one_year_from_now = time.increment_date(now, 1.0);

  // Take a one year long step if we can.
  MaxTimestep max_dt(one_year_from_now - now);

  max_dt = std::min(max_dt, model->max_timestep(now));

  // Do not take time-steps shorter than 1 second
  if (max_dt.value() < 1.0) {
    max_dt = MaxTimestep(1.0);
  }

  assert(max_dt.finite());

  model->update(geometry, now, max_dt.value());
}

} // end of namespace pism
