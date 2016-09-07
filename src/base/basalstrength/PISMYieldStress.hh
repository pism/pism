// Copyright (C) 2004--2012, 2014, 2015, 2016 Jed Brown, Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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

#ifndef _PISMYIELDSTRESS_H_
#define _PISMYIELDSTRESS_H_

#include "base/util/PISMComponent.hh"
#include "base/util/iceModelVec.hh"

namespace pism {

//! \brief The PISM basal yield stress model interface (virtual base class)
class YieldStress : public Component {
public:
  YieldStress(IceGrid::ConstPtr g);
  virtual ~YieldStress();

  void init();
  const IceModelVec2S& basal_material_yield_stress();
  void update();
protected:
  virtual void init_impl() = 0;
  virtual void update_impl() = 0;
  IceModelVec2S m_basal_yield_stress;
};

} // end of namespace pism

#endif /* _PISMYIELDSTRESS_H_ */
