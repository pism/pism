// Copyright (C) 2012, 2013, 2014 PISM Authors
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

#include "IceGrid.hh"           // Needed for Periodicity enum declaration.
#include "PISMNCFile.hh"

namespace pism {

enum AxisType {X_AXIS, Y_AXIS, Z_AXIS, T_AXIS, UNKNOWN_AXIS};

class grid_info;
class LocalInterpCtx;
class NCVariable;
class NCTimeseries;
class NCTimeBounds;

//! \brief High-level PISM I/O class.
/*!
 * Hides the low-level NetCDF wrapper.
 */
class PIO
{
public:
  PIO(MPI_Comm com, const std::string &mode, const UnitSystem &units_system);
  PIO(IceGrid &g, const std::string &mode);
  PIO(const PIO &other);
  ~PIO();

  void check_if_exists(const std::string &filename, bool &result);

  void open(const std::string &filename, IO_Mode mode);

  void close();

  void redef() const;

  void enddef() const;

  std::string inq_filename() const;

  void inq_nrecords(unsigned int &result) const;

  void inq_nrecords(const std::string &name, const std::string &std_name,
                    unsigned int &result) const;

  void inq_var(const std::string &short_name, const std::string &std_name, bool &exists,
               std::string &result, bool &found_by_standard_name) const;

  void inq_var(const std::string &short_name, bool &exists) const;

  void inq_vardims(const std::string &name, std::vector<std::string> &result) const;

  void inq_dim(const std::string &name, bool &exists) const;

  void inq_dimlen(const std::string &name, unsigned int &result) const;

  void inq_dimtype(const std::string &name, AxisType &result) const;

  void inq_dim_limits(const std::string &name, double *min, double *max) const;

  PetscErrorCode inq_grid(const std::string &var_name, IceGrid *grid, Periodicity periodicity) const;

  void inq_units(const std::string &name, bool &has_units, Unit &units) const;

  void inq_grid_info(const std::string &name, Periodicity p,
                     grid_info &g) const;

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

  void inq_nattrs(const std::string &var_name, int &result) const;

  void inq_attname(const std::string &var_name, unsigned int n, std::string &result) const;

  void inq_atttype(const std::string &var_name, const std::string &att_name, IO_Type &result) const;

  void put_att_double(const std::string &var_name, const std::string &att_name, IO_Type nctype,
                      const std::vector<double> &values) const;

  void put_att_double(const std::string &var_name, const std::string &att_name, IO_Type nctype,
                      double value) const;

  void put_att_text(const std::string &var_name, const std::string &att_name, const std::string &value) const;

  void get_att_double(const std::string &var_name, const std::string &att_name,
                      std::vector<double> &result) const;

  void get_att_text(const std::string &var_name, const std::string &att_name, std::string &result) const;

  PetscErrorCode get_vec(IceGrid *grid, const std::string &var_name, unsigned int z_count, unsigned int t, Vec g) const;

  PetscErrorCode put_vec(IceGrid *grid, const std::string &var_name, unsigned int z_count, Vec g) const;

  void regrid_vec(IceGrid *grid, const std::string &var_name,
                  const std::vector<double> &zlevels_out,
                  unsigned int t_start, Vec g) const;

  void regrid_vec_fill_missing(IceGrid *grid, const std::string &var_name,
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

private:
  MPI_Comm m_com;
  std::string m_mode;
  NCFile::Ptr m_nc;
  int m_xs, m_xm, m_ys, m_ym;
  UnitSystem m_unit_system;

  void use_mapped_io(std::string var_name, bool &result) const;

  void get_interp_context(const std::string &name,
                          const IceGrid &grid,
                          const std::vector<double> &zlevels,
                          LocalInterpCtx* &lic) const;

  void compute_start_and_count(const std::string &name,
                               unsigned int t_start, unsigned int t_count,
                               unsigned int x_start, unsigned int x_count,
                               unsigned int y_start, unsigned int y_count,
                               unsigned int z_start, unsigned int z_count,
                               std::vector<unsigned int> &start,
                               std::vector<unsigned int> &count,
                               std::vector<unsigned int> &imap) const;

  int k_below(double z, const std::vector<double> &zlevels) const;

  PetscErrorCode regrid(IceGrid *grid, const std::vector<double> &zlevels_out,
                        LocalInterpCtx *lic, Vec g) const;

  void detect_mode(const std::string &filename);

  void constructor(MPI_Comm com, const std::string &mode);
};

//! \brief Contains parameters of an input file grid.
class grid_info {
public:
  grid_info();
  // dimension lengths
  unsigned int t_len, x_len, y_len, z_len;
  double time,                  //!< current time (seconds)
    x0,                         //!< x-coordinate of the domain center
    y0,                         //!< y-coordinate of the domain center
    Lx,                         //!< domain half-width
    Ly,                         //!< domain half-height
    z_min,                      //!< minimal value of the z dimension
    z_max;                      //!< maximal value of the z dimension
  std::vector<double> x, y, z;       //!< coordinates
};

} // end of namespace pism

#endif /* _PIO_H_ */
