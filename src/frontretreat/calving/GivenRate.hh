/* Copyright (C) 2013, 2014, 2015, 2016, 2017, 2018, 2021, 2022 PISM Authors
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

#ifndef _PISMGIVENRATE_H_
#define _PISMGIVENRATE_H_

#include "pism/util/Component.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/IceModelVec2CellType.hh"

namespace pism {

class IceModelVec2T;

namespace calving {

/*! \brief Calving mechanism removing the ice at the shelf front for
    a given rate (2D map). */
class GivenRate : public Component
{
public:
  GivenRate(IceGrid::ConstPtr g);
  virtual ~GivenRate() = default;

  void init();
  void update(double t,
              double dt);

  const IceModelVec2S& calvingrate() const;

protected:
  DiagnosticList diagnostics_impl() const;
  std::shared_ptr<IceModelVec2T> m_calving_rate;
};

} // end of namespace calving
} // end of namespace pism

#endif /* _PISMGIVENRATE_H_ */
