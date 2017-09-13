// Copyright (C) 2009--2017 Constantine Khroulev
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

  It allocates a number (given as an argument to the set_n_records() method) of
  records and reads them from a file if necessary.

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
  IceModelVec2T();
  virtual ~IceModelVec2T();

  virtual void set_n_records(unsigned int N);
  virtual void set_n_evaluations_per_year(unsigned int N);
  virtual unsigned int get_n_records();
  void create(IceGrid::ConstPtr mygrid, const std::string &my_short_name);
  virtual void init(const std::string &filename, unsigned int period,
                              double reference_time);
  virtual void init_constant(double value);
  virtual void update(double my_t, double my_dt);
  virtual void set_record(int n);
  virtual void get_record(int n);
  MaxTimestep max_timestep(double my_t) const;

  virtual void interp(double my_t);

  virtual void interp(int i, int j, std::vector<double> &results);

  virtual void average(double my_t, double my_dt);
  virtual double average(int i, int j);

  virtual void begin_access() const;
  virtual void end_access() const;
  virtual void init_interpolation(const std::vector<double> &ts);

protected:
  std::vector<double> m_time,             //!< all the times available in filename
    m_time_bounds;                //!< time bounds
  std::string m_filename;         //!< file to read (regrid) from
  petsc::DM::Ptr m_da3;
  petsc::Vec m_v3;                       //!< a 3D Vec used to store records
  mutable void ***m_array3;
  unsigned int m_n_records, //!< maximum number of records to store in memory
    m_N,                    //!< number of records kept in memory
    m_n_evaluations_per_year;     //!< number of evaluations per year
  //!< used to compute temporal averages
  int m_first; //!< in-file index of the first record stored in memory
  //!< ("int" to allow first==-1 as an "invalid" first value)

  std::vector<unsigned int> m_interp_indices;
  unsigned int m_period;        // in years
  double m_reference_time;      // in seconds

  double*** get_array3();
  virtual void update(unsigned int start);
  virtual void discard(int N);
};


} // end of namespace pism

#endif // __IceModelVec2T_hh
