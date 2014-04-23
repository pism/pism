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
  PIO(MPI_Comm com, std::string mode, PISMUnitSystem units_system);
  PIO(IceGrid &g, std::string mode);
  PIO(const PIO &other);
  ~PIO();

  PetscErrorCode check_if_exists(std::string filename, bool &result);

  PetscErrorCode open(std::string filename, PISM_IO_Mode mode);

  PetscErrorCode close();

  PetscErrorCode redef() const;

  PetscErrorCode enddef() const;

  std::string inq_filename() const;

  PetscErrorCode inq_nrecords(unsigned int &result) const;

  PetscErrorCode inq_nrecords(std::string name, std::string std_name, unsigned int &result) const;

  PetscErrorCode inq_var(std::string short_name, std::string std_name, bool &exists,
                         std::string &result, bool &found_by_standard_name) const;

  PetscErrorCode inq_var(std::string short_name, bool &exists) const;

  PetscErrorCode inq_vardims(std::string name, std::vector<std::string> &result) const;

  PetscErrorCode inq_dim(std::string name, bool &exists) const;

  PetscErrorCode inq_dimlen(std::string name, unsigned int &result) const;

  PetscErrorCode inq_dimtype(std::string name, AxisType &result) const;

  PetscErrorCode inq_dim_limits(std::string name, double *min, double *max) const;

  PetscErrorCode inq_grid(std::string var_name, IceGrid *grid, Periodicity periodicity) const;

  PetscErrorCode inq_units(std::string name, bool &has_units, PISMUnit &units) const;

  PetscErrorCode inq_grid_info(std::string name, grid_info &g) const;

  PetscErrorCode def_dim(unsigned long int length, const NCVariable &metadata) const;

  PetscErrorCode def_var(std::string name, PISM_IO_Type nctype, std::vector<std::string> dims) const;

  PetscErrorCode get_dim(std::string name, std::vector<double> &result) const;

  PetscErrorCode get_1d_var(std::string name, unsigned int start, unsigned int count,
                            std::vector<double> &result) const;

  PetscErrorCode put_1d_var(std::string name, unsigned int start, unsigned int count,
                            const std::vector<double> &data) const;

  PetscErrorCode put_dim(std::string name, const std::vector<double> &data) const;

  PetscErrorCode append_time(std::string var_name, double value) const;

  PetscErrorCode def_time(std::string name, std::string calendar, std::string units) const;

  PetscErrorCode append_history(std::string history) const;

  PetscErrorCode inq_nattrs(std::string var_name, int &result) const;

  PetscErrorCode inq_attname(std::string var_name, unsigned int n, std::string &result) const;

  PetscErrorCode inq_atttype(std::string var_name, std::string att_name, PISM_IO_Type &result) const;

  PetscErrorCode put_att_double(std::string var_name, std::string att_name, PISM_IO_Type nctype,
                                const std::vector<double> values) const;

  PetscErrorCode put_att_double(std::string var_name, std::string att_name, PISM_IO_Type nctype,
                                double value) const;

  PetscErrorCode put_att_text(std::string var_name, std::string att_name, std::string value) const;

  PetscErrorCode get_att_double(std::string var_name, std::string att_name,
                                std::vector<double> &result) const;

  PetscErrorCode get_att_text(std::string var_name, std::string att_name, std::string &result) const;

  PetscErrorCode get_vec(IceGrid *grid, std::string var_name, unsigned int z_count, unsigned int t, Vec g) const;

  PetscErrorCode put_vec(IceGrid *grid, std::string var_name, unsigned int z_count, Vec g) const;

  PetscErrorCode regrid_vec(IceGrid *grid, std::string var_name,
                            const std::vector<double> &zlevels_out,
                            unsigned int t_start, Vec g) const;

  PetscErrorCode get_vara_double(std::string variable_name,
                                 std::vector<unsigned int> start,
                                 std::vector<unsigned int> count,
                                 double *ip) const;

  PetscErrorCode put_vara_double(std::string variable_name,
                                 std::vector<unsigned int> start,
                                 std::vector<unsigned int> count,
                                 double *op) const;

  PetscErrorCode get_varm_double(std::string variable_name,
                                 std::vector<unsigned int> start,
                                 std::vector<unsigned int> count,
                                 std::vector<unsigned int> imap, double *ip) const;

  PetscErrorCode put_varm_double(std::string variable_name,
                                 std::vector<unsigned int> start,
                                 std::vector<unsigned int> count,
                                 std::vector<unsigned int> imap, double *op) const;

  void set_local_extent(unsigned int xs, unsigned int xm,
                        unsigned int ys, unsigned int ym);

  PetscErrorCode read_timeseries(const NCTimeseries &metadata,
                                 PISMTime *time,
                                 std::vector<double> &data) const;


  PetscErrorCode write_timeseries(const NCTimeseries &metadata, size_t t_start,
                                  double data,
                                  PISM_IO_Type nctype = PISM_DOUBLE) const;
  PetscErrorCode write_timeseries(const NCTimeseries &metadata, size_t t_start,
                                  std::vector<double> &data,
                                  PISM_IO_Type nctype = PISM_DOUBLE) const;

  PetscErrorCode read_time_bounds(const NCTimeBounds &metadata,
                                  PISMTime *time,
                                  std::vector<double> &data) const;

  PetscErrorCode write_time_bounds(const NCTimeBounds &metadata, size_t t_start,
                                   std::vector<double> &data,
                                   PISM_IO_Type nctype = PISM_DOUBLE) const;

  PetscErrorCode read_attributes(std::string name, NCVariable &variable) const;
  PetscErrorCode write_attributes(const NCVariable &var, PISM_IO_Type nctype,
                                  bool write_in_glaciological_units) const;

  PetscErrorCode write_global_attributes(const NCVariable &var) const;

  PetscErrorCode read_valid_range(std::string name, NCVariable &variable) const;

private:
  MPI_Comm m_com;
  std::string m_mode;
  bool shallow_copy;
  PISMNCFile *nc;
  int m_xs, m_xm, m_ys, m_ym;
  PISMUnitSystem m_unit_system;

  PetscErrorCode get_interp_context(std::string name,
                                    const IceGrid &grid,
                                    const std::vector<double> &zlevels,
                                    LocalInterpCtx* &lic) const;

  PetscErrorCode compute_start_and_count(std::string name, unsigned int t_start,
                                         unsigned int x_start, unsigned int x_count,
                                         unsigned int y_start, unsigned int y_count,
                                         unsigned int z_start, unsigned int z_count,
                                         std::vector<unsigned int> &start,
                                         std::vector<unsigned int> &count,
                                         std::vector<unsigned int> &imap) const;

  PetscErrorCode k_below(double z, const std::vector<double> &zlevels) const;

  PetscErrorCode regrid(IceGrid *grid, const std::vector<double> &zlevels_out, LocalInterpCtx *lic, Vec g) const;

  PetscErrorCode detect_mode(std::string filename);

  void constructor(MPI_Comm com, std::string mode);
};

//! \brief Contains parameters of an input file grid.
class grid_info {
public:
  grid_info();
  // dimension lengths
  unsigned int t_len, x_len, y_len, z_len;
  double time,                  //!< current time (seconds)
    x_min,                      //!< [x_min, x_max] is the X extent of the grid
    x_max,                      //!< [x_min, x_max] is the X extent of the grid
    y_min,                      //!< [y_min, y_max] is the Y extent of the grid
    y_max,                      //!< [y_min, y_max] is the Y extent of the grid
    z_min,                      //!< minimal value of the z dimension
    z_max;                      //!< maximal value of the z dimension
  std::vector<double> x, y, z;       //!< coordinates
};

} // end of namespace pism

#endif /* _PIO_H_ */
