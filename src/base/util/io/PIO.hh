// Copyright (C) 2012, 2013, 2014, 2015 PISM Authors
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

#include "PISMUnits.hh"

#include <petscvec.h>
#include <map>
#include <vector>
#include <string>

#include "PISMNCFile.hh"

namespace pism {

enum AxisType {X_AXIS, Y_AXIS, Z_AXIS, T_AXIS, UNKNOWN_AXIS};

class IceGrid;
class LocalInterpCtx;
class NCVariable;
class NCTimeseries;
class NCTimeBounds;
class Time;

//! \brief High-level PISM I/O class.
/*!
 * Hides the low-level NetCDF wrapper.
 */
class PIO
{
public:
  PIO(MPI_Comm com, const std::string &mode, const UnitSystem &units_system);
  PIO(const IceGrid &g, const std::string &mode);
  PIO(const PIO &other);
  ~PIO();

  static bool check_if_exists(MPI_Comm com, const std::string &filename);

  void open(const std::string &filename, IO_Mode mode);

  void close();

  void redef() const;

  void enddef() const;

  std::string inq_filename() const;

  unsigned int inq_nrecords() const;

  unsigned int inq_nrecords(const std::string &name, const std::string &std_name) const;

  void inq_var(const std::string &short_name, const std::string &std_name, bool &exists,
               std::string &result, bool &found_by_standard_name) const;

  bool inq_var(const std::string &short_name) const;

  std::vector<std::string> inq_vardims(const std::string &name) const;

  bool inq_dim(const std::string &name) const;

  unsigned int inq_dimlen(const std::string &name) const;

  AxisType inq_dimtype(const std::string &name) const;

  void inq_dim_limits(const std::string &name, double *min, double *max) const;

  void inq_units(const std::string &name, bool &has_units, Unit &units) const;

  void def_dim(unsigned long int length, const NCVariable &metadata) const;

  void def_var(const std::string &name, IO_Type nctype, const std::vector<std::string> &dims) const;

  void get_dim(const std::string &name, std::vector<double> &result) const;

  void get_1d_var(const std::string &name, unsigned int start, unsigned int count,
                  std::vector<double> &result) const;

  void put_1d_var(const std::string &name, unsigned int start, unsigned int count,
                  const std::vector<double> &data) const;

  void put_dim(const std::string &name, const std::vector<double> &data) const;

  void append_time(const std::string &var_name, double value) const;

  void def_time(const std::string &name, const std::string &calendar, const std::string &units) const;

  void append_history(const std::string &history) const;

  unsigned int inq_nattrs(const std::string &var_name) const;

  std::string inq_attname(const std::string &var_name, unsigned int n) const;

  IO_Type inq_atttype(const std::string &var_name, const std::string &att_name) const;

  void put_att_double(const std::string &var_name, const std::string &att_name, IO_Type nctype,
                      const std::vector<double> &values) const;

  void put_att_double(const std::string &var_name, const std::string &att_name, IO_Type nctype,
                      double value) const;

  void put_att_text(const std::string &var_name, const std::string &att_name, const std::string &value) const;

  std::vector<double> get_att_double(const std::string &var_name, const std::string &att_name) const;

  std::string get_att_text(const std::string &var_name, const std::string &att_name) const;

  void get_vec(const IceGrid *grid, const std::string &var_name, unsigned int z_count,
               unsigned int t, Vec g) const;

  void put_vec(const IceGrid *grid, const std::string &var_name,
               unsigned int z_count, Vec g) const;

  void regrid_vec(const IceGrid *grid, const std::string &var_name,
                  const std::vector<double> &zlevels_out,
                  unsigned int t_start, Vec g) const;

  void regrid_vec_fill_missing(const IceGrid *grid, const std::string &var_name,
                               const std::vector<double> &zlevels_out,
                               unsigned int t_start,
                               double default_value,
                               Vec g) const ;

  void get_vara_double(const std::string &variable_name,
                       const std::vector<unsigned int> &start,
                       const std::vector<unsigned int> &count,
                       double *ip) const;

  void put_vara_double(const std::string &variable_name,
                       const std::vector<unsigned int> &start,
                       const std::vector<unsigned int> &count,
                       double *op) const;

  void get_varm_double(const std::string &variable_name,
                       const std::vector<unsigned int> &start,
                       const std::vector<unsigned int> &count,
                       const std::vector<unsigned int> &imap, double *ip) const;

  void put_varm_double(const std::string &variable_name,
                       const std::vector<unsigned int> &start,
                       const std::vector<unsigned int> &count,
                       const std::vector<unsigned int> &imap, double *op) const;

  void set_local_extent(unsigned int xs, unsigned int xm,
                        unsigned int ys, unsigned int ym);

  void read_timeseries(const NCTimeseries &metadata,
                       Time *time,
                       std::vector<double> &data) const;


  void write_timeseries(const NCTimeseries &metadata, size_t t_start,
                        double data,
                        IO_Type nctype = PISM_DOUBLE) const;
  void write_timeseries(const NCTimeseries &metadata, size_t t_start,
                        std::vector<double> &data,
                        IO_Type nctype = PISM_DOUBLE) const;

  void read_time_bounds(const NCTimeBounds &metadata,
                        Time *time,
                        std::vector<double> &data) const;

  void write_time_bounds(const NCTimeBounds &metadata, size_t t_start,
                         std::vector<double> &data,
                         IO_Type nctype = PISM_DOUBLE) const;

  void read_attributes(const std::string &name, NCVariable &variable) const;
  void write_attributes(const NCVariable &var, IO_Type nctype,
                        bool write_in_glaciological_units) const;

  void write_global_attributes(const NCVariable &var) const;

  void read_valid_range(const std::string &name, NCVariable &variable) const;

  std::string backend_type() const;

private:
  MPI_Comm m_com;
  std::string m_backend_type;
  NCFile::Ptr m_nc;
  int m_xs, m_xm, m_ys, m_ym;
  UnitSystem m_unit_system;

  void use_mapped_io(std::string var_name, bool &result) const;

  LocalInterpCtx* get_interp_context(const std::string &name,
                                     const IceGrid &grid,
                                     const std::vector<double> &zlevels) const;

  void compute_start_and_count(const std::string &name,
                               unsigned int t_start, unsigned int t_count,
                               unsigned int x_start, unsigned int x_count,
                               unsigned int y_start, unsigned int y_count,
                               unsigned int z_start, unsigned int z_count,
                               std::vector<unsigned int> &start,
                               std::vector<unsigned int> &count,
                               std::vector<unsigned int> &imap) const;

  int k_below(double z, const std::vector<double> &zlevels) const;

  void regrid(const IceGrid *grid, const std::vector<double> &zlevels_out,
              LocalInterpCtx *lic, double *output_array) const;

  void detect_mode(const std::string &filename);

  void constructor(MPI_Comm com, const std::string &mode);
};

} // end of namespace pism

#endif /* _PIO_H_ */
