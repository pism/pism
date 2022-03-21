// Copyright (C) 2009--2022 Constantine Khrulev
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

#ifndef PISM_ARRAY_FORCING
#define PISM_ARRAY_FORCING

#include "array/Scalar.hh"
#include "MaxTimestep.hh"
#include "interpolation.hh"     // InterpolationType

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
class IceModelVec2T : public array::Scalar {
public:

  IceModelVec2T(IceGrid::ConstPtr grid,
                const File &file,
                const std::string &short_name,
                const std::string &standard_name,
                int max_buffer_size,
                bool periodic,
                InterpolationType interpolation_type = PIECEWISE_CONSTANT);

  virtual ~IceModelVec2T();

  static std::shared_ptr<IceModelVec2T>
  Constant(IceGrid::ConstPtr grid, const std::string &short_name, double value);

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

  IceModelVec2T(IceGrid::ConstPtr grid,
                const std::string &short_name,
                unsigned int buffer_size,
                bool dummy,
                InterpolationType interpolation_type);
  void allocate(IceGrid::ConstPtr grid,
                const std::string &short_name,
                unsigned int buffer_size,
                InterpolationType interpolation_type);

  double*** array3();
  void update(unsigned int start);
  void discard(int N);
  void set_record(int n);
  void init_periodic_data(const File &file);
};


} // end of namespace pism

#endif // PISM_ARRAY_FORCING
