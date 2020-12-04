// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2019 PISM Authors
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

#ifndef _PISM_FILE_ACCESS_H_
#define _PISM_FILE_ACCESS_H_

#include <vector>
#include <string>
#include <mpi.h>
#include <map>
#include <set>

#include "pism/util/Units.hh"
#include "pism/util/io/IO_Flags.hh"

namespace pism {

enum AxisType {X_AXIS, Y_AXIS, Z_AXIS, T_AXIS, UNKNOWN_AXIS};

class IceGrid;

/*!
 * Convert a string to PISM's backend type.
 */
IO_Backend string_to_backend(const std::string &backend);

struct VariableLookupData {
  bool exists;
  bool found_using_standard_name;
  std::string name;
};

//! \brief High-level PISM I/O class.
/*!
 * Hides the low-level NetCDF wrapper.
 */
class File
{
public:
  File(MPI_Comm com, const std::string &filename, IO_Backend backend, IO_Mode mode,
       int iosysid = -1, const std::map<std::string, int> &varsi = std::map<std::string, int>(),
       const std::vector<int>& gridIDs = std::vector<int>(), int FileID = -1);
  ~File();

  IO_Backend backend() const;

  MPI_Comm com() const;

  void set_dimatt() const;

  std::string get_dimatt_value(std::string &dim_name) const; 

  void set_dimatt_value(std::string &dim_name, std::string &dimatt_value) const; 

  void close();

  void redef() const;

  void enddef() const;

  void sync() const;

  std::string filename() const;

  unsigned int nrecords() const;

  unsigned int nrecords(const std::string &name, const std::string &std_name,
                        units::System::Ptr unit_system) const;

  unsigned int nvariables() const;

  unsigned int nattributes(const std::string &var_name) const;

  // dimensions

  void define_dimension(const std::string &name, size_t length) const;

  unsigned int dimension_length(const std::string &name) const;

  std::vector<std::string> dimensions(const std::string &variable_name) const;

  bool find_dimension(const std::string &name) const;

  AxisType dimension_type(const std::string &name,
                          units::System::Ptr unit_system) const;

  std::vector<double> read_dimension(const std::string &name) const;

  // variables

  std::string variable_name(unsigned int id) const;

  void define_variable(const std::string &name, IO_Type nctype,
                       const std::vector<std::string> &dims) const;

  VariableLookupData find_variable(const std::string &short_name, const std::string &std_name) const;

  bool find_variable(const std::string &short_name) const;

  void read_variable(const std::string &variable_name,
                       const std::vector<unsigned int> &start,
                       const std::vector<unsigned int> &count,
                       double *ip) const;

  void read_variable_transposed(const std::string &variable_name,
                                const std::vector<unsigned int> &start,
                                const std::vector<unsigned int> &count,
                                const std::vector<unsigned int> &imap, double *ip) const;

  void write_variable(const std::string &variable_name,
                      const std::vector<unsigned int> &start,
                      const std::vector<unsigned int> &count,
                      const double *op) const;

  void write_distributed_array(const std::string &variable_name,
                               const IceGrid &grid,
                               unsigned int z_count,
                               const double *input) const;

  // attributes

  void remove_attribute(const std::string &variable_name, const std::string &att_name) const;

  std::string attribute_name(const std::string &var_name, unsigned int n) const;

  IO_Type attribute_type(const std::string &var_name, const std::string &att_name) const;

  void write_attribute(const std::string &var_name, const std::string &att_name,
                      IO_Type nctype, const std::vector<double> &values) const;

  void write_attribute(const std::string &var_name, const std::string &att_name,
                       const std::string &value) const;

  std::vector<double> read_double_attribute(const std::string &var_name,
                                            const std::string &att_name) const;

  std::string read_text_attribute(const std::string &var_name, const std::string &att_name) const;

  void append_history(const std::string &history) const;

  //new functions because of CDI class
  void new_grid(int lengthx, int lengthy) const;
  void new_timestep(int tsID) const;
  void expose_windows() const;
  void reference_date(double time) const;
  std::map<std::string, int> get_variables_map() const;
  void define_vlist() const;
  void send_diagnostics(const std::set<std::string> &variables) const;
  void set_beforediag(bool value) const;
  void set_gridIDs(const std::vector<int>& gridIDs = std::vector<int>()) const;
  std::vector<int> get_gridIDs() const;
  int get_streamID() const;
  int get_vlistID() const;

private:
  struct Impl;
  Impl *m_impl;
  mutable std::map<std::string, std::string> m_dimatt;
  void open(const std::string &filename, IO_Mode mode, const std::map<std::string, int> &varsi = std::map<std::string, int>(), int FileID = -1);

  // disable copying and assignments
  File(const File &other);
  File & operator=(const File &);
};

} // end of namespace pism

#endif /* _PISM_FILE_ACCESS_H_ */
