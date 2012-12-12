// Copyright (C) 2012 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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

#include "udunits.h"

#include <petscvec.h>
#include <map>
#include <vector>
#include <string>

#include "IceGrid.hh"           // Needed for Periodicity enum declaration.
#include "PISMNCFile.hh"

// use namespace std BUT remove trivial namespace browser from doxygen-erated HTML source browser
/// @cond NAMESPACE_BROWSER
using namespace std;
/// @endcond

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
  PIO(MPI_Comm com, int rank, string mode);
  PIO(IceGrid &g, string mode);
  PIO(const PIO &other);
  virtual ~PIO();

  virtual PetscErrorCode open(string filename, int mode, bool append = false);

  virtual PetscErrorCode close();

  virtual PetscErrorCode redef() const;

  virtual PetscErrorCode enddef() const;

  virtual string inq_filename() const;

  virtual PetscErrorCode inq_nrecords(unsigned int &result) const;

  virtual PetscErrorCode inq_nrecords(string name, string std_name, unsigned int &result) const;

  virtual PetscErrorCode inq_var(string short_name, string std_name, bool &exists,
                                 string &result, bool &found_by_standard_name) const;

  virtual PetscErrorCode inq_var(string short_name, bool &exists) const;

  virtual PetscErrorCode inq_vardims(string name, vector<string> &result) const;

  virtual PetscErrorCode inq_dim(string name, bool &exists) const;

  virtual PetscErrorCode inq_dimlen(string name, unsigned int &result) const;

  virtual PetscErrorCode inq_dimtype(string name, AxisType &result) const;

  virtual PetscErrorCode inq_dim_limits(string name, double *min, double *max) const;

  virtual PetscErrorCode inq_grid(string var_name, IceGrid *grid, Periodicity periodicity) const;

  virtual PetscErrorCode inq_units(string name, bool &has_units, utUnit &units,
                                   bool use_reference_date = false) const;

  virtual PetscErrorCode inq_grid_info(string name, grid_info &g) const;

  virtual PetscErrorCode def_dim(string name, long int length, map<string,string> attrs) const;

  virtual PetscErrorCode def_var(string name, PISM_IO_Type nctype, vector<string> dims) const;

  virtual PetscErrorCode get_dim(string name, vector<double> &result) const;

  virtual PetscErrorCode get_1d_var(string name, unsigned int start, unsigned int count,
                                    vector<double> &result) const;

  virtual PetscErrorCode put_1d_var(string name, unsigned int start, unsigned int count,
                                    const vector<double> &data) const;

  virtual PetscErrorCode put_dim(string name, const vector<double> &data) const;

  virtual PetscErrorCode append_time(string var_name, double value) const;

  virtual PetscErrorCode def_time(string name, string calendar, string units) const;

  virtual PetscErrorCode append_history(string history) const;

  virtual PetscErrorCode inq_nattrs(string var_name, int &result) const;

  virtual PetscErrorCode inq_attname(string var_name, unsigned int n, string &result) const;

  virtual PetscErrorCode inq_atttype(string var_name, string att_name, PISM_IO_Type &result) const;

  virtual PetscErrorCode put_att_double(string var_name, string att_name, PISM_IO_Type nctype,
                                        const vector<double> values) const;

  virtual PetscErrorCode put_att_double(string var_name, string att_name, PISM_IO_Type nctype,
                                        double value) const;

  virtual PetscErrorCode put_att_text(string var_name, string att_name, string value) const;

  virtual PetscErrorCode get_att_double(string var_name, string att_name,
                                        vector<double> &result) const;

  virtual PetscErrorCode get_att_text(string var_name, string att_name, string &result) const;

  virtual PetscErrorCode get_vec(IceGrid *grid, string var_name, unsigned int z_count, int t, Vec g) const;

  virtual PetscErrorCode put_vec(IceGrid *grid, string var_name, unsigned int z_count, Vec g) const;

  virtual PetscErrorCode regrid_vec(IceGrid *grid, string var_name,
                                    const vector<double> &zlevels_out, LocalInterpCtx *lic, Vec g) const;

  virtual PetscErrorCode get_vara_double(string variable_name,
                                         vector<unsigned int> start,
                                         vector<unsigned int> count,
                                         double *ip) const;

  virtual PetscErrorCode put_vara_double(string variable_name,
                                         vector<unsigned int> start,
                                         vector<unsigned int> count,
                                         double *op) const;

  virtual PetscErrorCode get_varm_double(string variable_name,
                                         vector<unsigned int> start,
                                         vector<unsigned int> count,
                                         vector<unsigned int> imap, double *ip) const;

  virtual PetscErrorCode put_varm_double(string variable_name,
                                         vector<unsigned int> start,
                                         vector<unsigned int> count,
                                         vector<unsigned int> imap, double *op) const;

  void set_local_extent(unsigned int xs, unsigned int xm,
                        unsigned int ys, unsigned int ym);
protected:
  MPI_Comm com;
  int rank;
  string m_mode;
  bool shallow_copy;
  PISMNCFile *nc;
  int m_xs, m_xm, m_ys, m_ym;

  PetscErrorCode compute_start_and_count(string name, int t_start,
                                         int x_start, int x_count,
                                         int y_start, int y_count,
                                         int z_start, int z_count,
                                         vector<unsigned int> &start,
                                         vector<unsigned int> &count,
                                         vector<unsigned int> &imap) const;

  PetscErrorCode k_below(double z, const vector<double> &zlevels) const;

  PetscErrorCode regrid(IceGrid *grid, const vector<double> &zlevels_out, LocalInterpCtx *lic, Vec g) const;

  PetscErrorCode detect_mode(string filename);
private:
  void constructor(MPI_Comm com, int rank, string mode);
};

#endif /* _PIO_H_ */
