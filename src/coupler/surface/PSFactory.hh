// Copyright (C) 2011, 2014 PISM Authors
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

#ifndef _PSFACTORY_H_
#define _PSFACTORY_H_

#include "PSModifier.hh"

class PSFactory : public PCFactory<PISMSurfaceModel,PSModifier> {
public:
  PSFactory(IceGrid& g, const PISMConfig& conf)
    : PCFactory<PISMSurfaceModel,PSModifier>(g, conf)
  {
    add_standard_types();
    option = "surface";
  }

  virtual ~PSFactory() {}
  virtual void add_standard_types();
};

#endif /* _PSFACTORY_H_ */
