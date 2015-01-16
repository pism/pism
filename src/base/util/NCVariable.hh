// Copyright (C) 2009--2015 Constantine Khroulev
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

#ifndef __NCVariable_hh
#define __NCVariable_hh

#include <set>
#include <map>
#include <vector>
#include <string>

#include "PISMUnits.hh"

// We use PIO and IO_Type here. (I should move methods using this out
// of NCSpatialVariable. -- CK)
#include "PIO.hh"

namespace pism {

class Time;

//! @brief A class for handling variable metadata, reading, writing and converting
//! from input units and to output units.
/*! A NetCDF variable can have any number of attributes, but some of them get
  special treatment:

  - units: specifies internal units. When read, a variable is
  converted to these units. When written, it is converted from these
  to glaciological_units if write_in_glaciological_units is true.

  - glaciological_units: is never written to a file; replaces 'units'
  in the output if write_in_glaciological_units is true.

  - valid_min, valid_max: specify the valid range of a variable. Are
  read from an input file *only* if not specified previously. If
  both are set, then valid_range is used in the output instead.

  Also:

  - empty string attributes are ignored (they are not written to the
  output; file and has_attribute("foo") returns false if "foo" is
  absent or equal to an empty string).

  Typical attributes stored here:

  - long_name
  - standard_name
  - pism_intent
  - units
  - glaciological_units (saved to files as "units")

  Use the `name` of "PISM_GLOBAL" to read and write global attributes.
  (See also PIO.)

*/

class NCVariable {
public:
  NCVariable(const std::string &name, const UnitSystem &system, unsigned int ndims = 0);
  virtual ~NCVariable();

  // setters
  void set_units(const std::string &unit_spec);
  void set_glaciological_units(const std::string &unit_spec);

  void set_double(const std::string &name, double value);
  void set_doubles(const std::string &name, const std::vector<double> &values);
  void set_name(const std::string &name);
  void set_string(const std::string &name, const std::string &value);

  void clear_all_doubles();
  void clear_all_strings();

  // getters
  Unit get_units() const;
  Unit get_glaciological_units() const;

  double get_double(const std::string &name) const;
  std::vector<double> get_doubles(const std::string &name) const;
  std::string get_name() const;
  std::string get_string(const std::string &name) const;

  unsigned int get_n_spatial_dimensions() const;

  bool has_attribute(const std::string &name) const;

  typedef std::map<std::string,std::string> StringAttrs;
  const StringAttrs& get_all_strings() const;

  typedef std::map<std::string,std::vector<double> > DoubleAttrs;
  const DoubleAttrs& get_all_doubles() const;

  void report_to_stdout(MPI_Comm com, int verbosity_threshold) const;

protected:
  unsigned int m_n_spatial_dims;

private:
  //! internal (PISM) units
  Unit m_units;

  //! @brief for diagnostic variables: units to use when writing to a
  //! NetCDF file and for standard out reports
  Unit m_glaciological_units;

  //! string and boolean attributes
  std::map<std::string, std::string> m_strings;

  //! scalar and array attributes
  std::map<std::string, std::vector<double> > m_doubles;
  std::string m_short_name;
};

class LocalInterpCtx;
class IceGrid;

enum RegriddingFlag {OPTIONAL, OPTIONAL_FILL_MISSING, CRITICAL, CRITICAL_FILL_MISSING};

//! Spatial NetCDF variable (corresponding to a 2D or 3D scalar field).
class NCSpatialVariable : public NCVariable {
public:
  NCSpatialVariable(const UnitSystem &system, const std::string &name,
                    const IceGrid &g);
  NCSpatialVariable(const UnitSystem &system, const std::string &name,
                    const IceGrid &g, const std::vector<double> &zlevels);
  NCSpatialVariable(const NCSpatialVariable &other);
  virtual ~NCSpatialVariable();
  void set_levels(const std::vector<double> &levels);

  void set_time_independent(bool flag);

  void read(const PIO &file, unsigned int time, double *output);
  void write(const PIO &file, IO_Type nctype,
             bool write_in_glaciological_units, const double *input) const;

  void regrid(const PIO &file,
              RegriddingFlag flag,
              bool report_range,
              double default_value, double *output);
  void regrid(const PIO &file,
              unsigned int t_start,
              RegriddingFlag flag,
              bool report_range,
              double default_value, double *output);

  void define(const PIO &nc, IO_Type nctype,
              bool write_in_glaciological_units) const;

  NCVariable& get_x();
  NCVariable& get_y();
  NCVariable& get_z();

  const NCVariable& get_x() const;
  const NCVariable& get_y() const;
  const NCVariable& get_z() const;

private:
  MPI_Comm m_com;
  std::string m_variable_order;        //!< variable order in output files;
  std::string m_time_dimension_name;
  NCVariable m_x, m_y, m_z;
  std::vector<double> m_zlevels;
  const IceGrid *m_grid;
  void report_range(double min, double max, bool found_by_standard_name);
  void check_range(const std::string &filename, double min, double max);
  void define_dimensions(const PIO &nc) const;

  void init_internal(const std::string &name, const IceGrid &g,
                     const std::vector<double> &z_levels);
};

//! An internal class for reading, writing and converting time-series.
class NCTimeseries : public NCVariable {
public:
  NCTimeseries(const std::string &name, const std::string &dimension_name,
               const UnitSystem &system);
  virtual ~NCTimeseries();

  std::string get_dimension_name() const;

  virtual void define(const PIO &nc, IO_Type nctype, bool) const;
private:
  std::string m_dimension_name;        //!< the name of the NetCDF dimension this timeseries depends on
};

class NCTimeBounds : public NCTimeseries
{
public:
  NCTimeBounds(const std::string &name, const std::string &dimension_name,
               const UnitSystem &system);
  virtual ~NCTimeBounds();
  virtual void define(const PIO &nc, IO_Type nctype, bool) const;
private:
  std::string m_bounds_name;
};

} // end of namespace pism

#endif  // __NCVariable_hh
