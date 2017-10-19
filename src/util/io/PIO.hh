// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017 PISM Authors
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

#ifndef _PIO_H_
#define _PIO_H_

#include <vector>
#include <string>
#include <mpi.h>

#include "pism/util/Units.hh"
#include "pism/util/io/IO_Flags.hh"

namespace pism {

enum AxisType {X_AXIS, Y_AXIS, Z_AXIS, T_AXIS, UNKNOWN_AXIS};

//! \brief High-level PISM I/O class.
/*!
 * Hides the low-level NetCDF wrapper.
 */
class PIO
{
public:
  PIO(MPI_Comm com, const std::string &backend, const std::string &filename, IO_Mode mode);
  ~PIO();

  MPI_Comm com() const;

  void close();

  void redef() const;

  void enddef() const;

  std::string inq_filename() const;

  unsigned int inq_nrecords() const;

  unsigned int inq_nrecords(const std::string &name, const std::string &std_name,
                            units::System::Ptr unit_system) const;

  void inq_var(const std::string &short_name, const std::string &std_name, bool &exists,
               std::string &result, bool &found_by_standard_name) const;

  bool inq_var(const std::string &short_name) const;

  std::vector<std::string> inq_vardims(const std::string &name) const;

  bool inq_dim(const std::string &name) const;

  unsigned int inq_dimlen(const std::string &name) const;

  AxisType inq_dimtype(const std::string &name,
                       units::System::Ptr unit_system) const;

  void inq_dim_limits(const std::string &name, double *min, double *max) const;

  void def_dim(const std::string &name, size_t length) const;

  void def_var(const std::string &name, IO_Type nctype,
               const std::vector<std::string> &dims) const;

  void get_dim(const std::string &name, std::vector<double> &result) const;

  void get_1d_var(const std::string &name, unsigned int start, unsigned int count,
                  std::vector<double> &result) const;

  void put_1d_var(const std::string &name, unsigned int start, unsigned int count,
                  const std::vector<double> &data) const;

  void append_history(const std::string &history) const;

  unsigned int inq_nattrs(const std::string &var_name) const;

  std::string inq_attname(const std::string &var_name, unsigned int n) const;

  IO_Type inq_atttype(const std::string &var_name, const std::string &att_name) const;

  void put_att_double(const std::string &var_name, const std::string &att_name,
                      IO_Type nctype, const std::vector<double> &values) const;

  void put_att_double(const std::string &var_name, const std::string &att_name,
                      IO_Type nctype, double value) const;

  void put_att_text(const std::string &var_name, const std::string &att_name,
                    const std::string &value) const;

  std::vector<double> get_att_double(const std::string &var_name,
                                     const std::string &att_name) const;

  std::string get_att_text(const std::string &var_name, const std::string &att_name) const;

  void get_vara_double(const std::string &variable_name,
                       const std::vector<unsigned int> &start,
                       const std::vector<unsigned int> &count,
                       double *ip) const;

  void put_vara_double(const std::string &variable_name,
                       const std::vector<unsigned int> &start,
                       const std::vector<unsigned int> &count,
                       const double *op) const;

  void get_varm_double(const std::string &variable_name,
                       const std::vector<unsigned int> &start,
                       const std::vector<unsigned int> &count,
                       const std::vector<unsigned int> &imap, double *ip) const;

  void put_varm_double(const std::string &variable_name,
                       const std::vector<unsigned int> &start,
                       const std::vector<unsigned int> &count,
                       const std::vector<unsigned int> &imap,
                       const double *op) const;

  std::string backend_type() const;
private:
  struct Impl;
  Impl *m_impl;

  void open(const std::string &filename, IO_Mode mode);

  void detect_mode(const std::string &filename);

  // disable copying and assignments
  PIO(const PIO &other);
  PIO & operator=(const PIO &);
};

} // end of namespace pism

#endif /* _PIO_H_ */
