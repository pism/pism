// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2021 PISM Authors
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

#ifndef _PISMCONSTANTYIELDSTRESS_H_
#define _PISMCONSTANTYIELDSTRESS_H_

#include "YieldStress.hh"
#include "pism/util/iceModelVec.hh"

namespace pism {

class IceGrid;

class ConstantYieldStress : public YieldStress {
public:
  ConstantYieldStress(IceGrid::ConstPtr g);
  virtual ~ConstantYieldStress() = default;
private:
  void restart_impl(const File &input_file, int record);

  void bootstrap_impl(const File &input_file, const YieldStressInputs &inputs);

  void init_impl(const YieldStressInputs &inputs);

  void update_impl(const YieldStressInputs &inputs, double t, double dt);

  MaxTimestep max_timestep_impl(double t) const;
};

} // end of namespace pism

#endif /* _PISMCONSTANTYIELDSTRESS_H_ */
