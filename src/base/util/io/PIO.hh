// Copyright (C) 2012, 2013 PISM Authors
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

enum AxisType {X_AXIS, Y_AXIS, Z_AXIS, T_AXIS, UNKNOWN_AXIS};

class grid_info;
class LocalInterpCtx;

//! \brief High-level PISM I/O class.
/*!
 * Hides the low-level NetCDF wrapper.
 */
class PIO
{
public:
  PIO(MPI_Comm com, int rank, std::string mode, PISMUnitSystem units_system);
  PIO(IceGrid &g, std::string mode);
  PIO(const PIO &other);
  virtual ~PIO();

  virtual PetscErrorCode check_if_exists(std::string filename, bool &result);

  virtual PetscErrorCode open(std::string filename, int mode, bool append = false);

  virtual PetscErrorCode close();

  virtual PetscErrorCode redef() const;

  virtual PetscErrorCode enddef() const;

  virtual std::string inq_filename() const;

  virtual PetscErrorCode inq_nrecords(unsigned int &result) const;

  virtual PetscErrorCode inq_nrecords(std::string name, std::string std_name, unsigned int &result) const;

  virtual PetscErrorCode inq_var(std::string short_name, std::string std_name, bool &exists,
                                 std::string &result, bool &found_by_standard_name) const;

  virtual PetscErrorCode inq_var(std::string short_name, bool &exists) const;

  virtual PetscErrorCode inq_vardims(std::string name, std::vector<std::string> &result) const;

  virtual PetscErrorCode inq_dim(std::string name, bool &exists) const;

  virtual PetscErrorCode inq_dimlen(std::string name, unsigned int &result) const;

  virtual PetscErrorCode inq_dimtype(std::string name, AxisType &result) const;

  virtual PetscErrorCode inq_dim_limits(std::string name, double *min, double *max) const;

  virtual PetscErrorCode inq_grid(std::string var_name, IceGrid *grid, Periodicity periodicity) const;

  virtual PetscErrorCode inq_units(std::string name, bool &has_units, PISMUnit &units) const;

  virtual PetscErrorCode inq_grid_info(std::string name, grid_info &g) const;

  virtual PetscErrorCode def_dim(std::string name, long int length, std::map<std::string,std::string> attrs) const;

  virtual PetscErrorCode def_var(std::string name, PISM_IO_Type nctype, std::vector<std::string> dims) const;

  virtual PetscErrorCode get_dim(std::string name, std::vector<double> &result) const;

  virtual PetscErrorCode get_1d_var(std::string name, unsigned int start, unsigned int count,
                                    std::vector<double> &result) const;

  virtual PetscErrorCode put_1d_var(std::string name, unsigned int start, unsigned int count,
                                    const std::vector<double> &data) const;

  virtual PetscErrorCode put_dim(std::string name, const std::vector<double> &data) const;

  virtual PetscErrorCode append_time(std::string var_name, double value) const;

  virtual PetscErrorCode def_time(std::string name, std::string calendar, std::string units) const;

  virtual PetscErrorCode append_history(std::string history) const;

  virtual PetscErrorCode inq_nattrs(std::string var_name, int &result) const;

  virtual PetscErrorCode inq_attname(std::string var_name, unsigned int n, std::string &result) const;

  virtual PetscErrorCode inq_atttype(std::string var_name, std::string att_name, PISM_IO_Type &result) const;

  virtual PetscErrorCode put_att_double(std::string var_name, std::string att_name, PISM_IO_Type nctype,
                                        const std::vector<double> values) const;

  virtual PetscErrorCode put_att_double(std::string var_name, std::string att_name, PISM_IO_Type nctype,
                                        double value) const;

  virtual PetscErrorCode put_att_text(std::string var_name, std::string att_name, std::string value) const;

  virtual PetscErrorCode get_att_double(std::string var_name, std::string att_name,
                                        std::vector<double> &result) const;

  virtual PetscErrorCode get_att_text(std::string var_name, std::string att_name, std::string &result) const;

  virtual PetscErrorCode get_vec(IceGrid *grid, std::string var_name, unsigned int z_count, int t, Vec g) const;

  virtual PetscErrorCode put_vec(IceGrid *grid, std::string var_name, unsigned int z_count, Vec g) const;

  virtual PetscErrorCode regrid_vec(IceGrid *grid, std::string var_name,
                                    const std::vector<double> &zlevels_out, LocalInterpCtx *lic, Vec g) const;

  virtual PetscErrorCode get_vara_double(std::string variable_name,
                                         std::vector<unsigned int> start,
                                         std::vector<unsigned int> count,
                                         double *ip) const;

  virtual PetscErrorCode put_vara_double(std::string variable_name,
                                         std::vector<unsigned int> start,
                                         std::vector<unsigned int> count,
                                         double *op) const;

  virtual PetscErrorCode get_varm_double(std::string variable_name,
                                         std::vector<unsigned int> start,
                                         std::vector<unsigned int> count,
                                         std::vector<unsigned int> imap, double *ip) const;

  virtual PetscErrorCode put_varm_double(std::string variable_name,
                                         std::vector<unsigned int> start,
                                         std::vector<unsigned int> count,
                                         std::vector<unsigned int> imap, double *op) const;

  void set_local_extent(unsigned int xs, unsigned int xm,
                        unsigned int ys, unsigned int ym);
protected:
  MPI_Comm com;
  int rank;
  std::string m_mode;
  bool shallow_copy;
  PISMNCFile *nc;
  int m_xs, m_xm, m_ys, m_ym;
  PISMUnitSystem m_unit_system;

  PetscErrorCode compute_start_and_count(std::string name, int t_start,
                                         int x_start, int x_count,
                                         int y_start, int y_count,
                                         int z_start, int z_count,
                                         std::vector<unsigned int> &start,
                                         std::vector<unsigned int> &count,
                                         std::vector<unsigned int> &imap) const;

  PetscErrorCode k_below(double z, const std::vector<double> &zlevels) const;

  PetscErrorCode regrid(IceGrid *grid, const std::vector<double> &zlevels_out, LocalInterpCtx *lic, Vec g) const;

  PetscErrorCode detect_mode(std::string filename);
private:
  void constructor(MPI_Comm com, int rank, std::string mode);
};

#endif /* _PIO_H_ */
