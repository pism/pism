// Copyright (C)  2009-2018 Ricarda Winkelmann, Torsten Albrecht, Constantine Khrulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef __PAPIK_hh
#define __PAPIK_hh

#include "YearlyCycle.hh"

namespace pism {
namespace atmosphere {

class PIK : public YearlyCycle {
public:
  PIK(IceGrid::ConstPtr g);
  virtual ~PIK();

private:
  void init_impl(const Geometry &geometry);

  MaxTimestep max_timestep_impl(double t) const;
  void update_impl(const Geometry &geometry, double t, double dt);

  enum Parameterization {MARTIN, HUYBRECHTS_DEWOLDE, MARTIN_HUYBRECHTS_DEWOLDE,
                         ERA_INTERIM, ERA_INTERIM_SIN, ERA_INTERIM_LON};

  Parameterization m_parameterization;
};

} // end of namespace atmosphere
} // end of namespace pism

#endif // __PAPIK_hh
