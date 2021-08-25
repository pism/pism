// Copyright (C) 2009--2021 Constantine Khroulev
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

#ifndef __IceModelVec2T_hh
#define __IceModelVec2T_hh

#include "iceModelVec.hh"
#include "MaxTimestep.hh"

namespace pism {

//! A class for storing and accessing 2D time-series (for climate forcing)
/*! This class was created to read time-dependent and spatially-varying climate
  forcing data, in particular snow temperatures and precipitation.

  If requests (calls to update()) go in sequence, every records should be read
  only once.

  Note that this class is optimized for use with a PDD scheme -- it stores
  records so that data corresponding to a grid point are stored in adjacent
  memory locations.

  IceModelVec2T is always global (%i.e. has no ghosts).

  Both versions of interp() use piecewise-constant interpolation and
  extrapolate (by a constant) outside the available range.
*/
class IceModelVec2T : public IceModelVec2S {
public:
  static std::shared_ptr<IceModelVec2T>
  ForcingField(IceGrid::ConstPtr grid,
               const File &file,
               const std::string &short_name,
               const std::string &standard_name,
               int max_buffer_size,
               int evaluations_per_year,
               bool periodic,
               InterpolationType interpolation_type = PIECEWISE_CONSTANT);

  static std::shared_ptr<IceModelVec2T> Constant(IceGrid::ConstPtr grid,
                                                 const std::string &short_name,
                                                 double value);

  IceModelVec2T(IceGrid::ConstPtr grid,
                const std::string &short_name,
                unsigned int buffer_size,
                unsigned int n_evaluations_per_year,
                InterpolationType interpolation_type = PIECEWISE_CONSTANT);
  virtual ~IceModelVec2T();

  unsigned int buffer_size();

  void init(const std::string &filename, bool periodic);

  void update(double t, double dt);
  MaxTimestep max_timestep(double t) const;

  void interp(double t);

  void interp(int i, int j, std::vector<double> &results);

  void average(double t, double dt);

  void begin_access() const;
  void end_access() const;
  void init_interpolation(const std::vector<double> &ts);

private:
  struct Data;

  Data *m_data;

  double*** array3();
  void update(unsigned int start);
  void discard(int N);
  double average(int i, int j);
  void set_record(int n);
  void get_record(int n);
};


} // end of namespace pism

#endif // __IceModelVec2T_hh
