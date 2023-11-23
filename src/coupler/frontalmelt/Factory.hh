// Copyright (C) 2011, 2012, 2014, 2015, 2017, 2018, 2021 PISM Authors
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

#ifndef _PFMFACTORY_H_
#define _PFMFACTORY_H_

#include "pism/coupler/util/PCFactory.hh"
#include "pism/coupler/FrontalMelt.hh"

namespace pism {
namespace frontalmelt {
class Factory : public PCFactory<frontalmelt::FrontalMelt> {
public:
  Factory(std::shared_ptr<const Grid> g);
  ~Factory() = default;
};
} // end of namespace frontalmelt
} // end of namespace pism

#endif /* _PFMFACTORY_H_ */
