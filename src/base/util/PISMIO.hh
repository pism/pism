// Copyright (C) 2006--2011 Jed Brown, Ed Bueler and Constantine Khroulev
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

#ifndef __PISMIO
#define __PISMIO

#include <petscvec.h>
#include "NCTool.hh"

class IceGrid;
class LocalInterpCtx;
struct grid_info;

//! A class containing IO functions used to read and write spatial variables.
/*!
  \section pismio_overview Reading and writing gridded data (PISMIO)

  PISMIO is used to read gridded data from NetCDF files. Most uses of PISMIO
  are hidden inside NCSpatialVariable::read(), NCSpatialVariable::write() and
  NCSpatialVariable::regrid(), which are, in turn, used by IceModelVec.

  You need to use this class to "prepare" a NetCDF file for writing:
  \code
  PISMIO nc(&grid);

  ierr = nc.open_for_writing(filename, false, true); CHKERRQ(ierr);
  ierr = nc.append_time(grid.year); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);
  \endcode

  No action is needed to be able to write to an output ("-o") file, a snapshot
  file or the like, though; IceModel has already prepared it.
 */
class PISMIO : public NCTool {
public:
  PISMIO(IceGrid *my_grid);
  virtual ~PISMIO() {}

  using NCTool::open_for_writing;
  PetscErrorCode open_for_writing(string filename, bool append,
				  bool check_dims = false);
  PetscErrorCode get_grid(string filename, string var_name);
  PetscErrorCode create_dimensions() const;

  PetscErrorCode get_grid_info(string name, grid_info &g) const;

  PetscErrorCode get_var(int varid, Vec g, int z_count, int t) const;
  PetscErrorCode put_var(int varid, Vec g, int z_count) const;

  PetscErrorCode regrid_var(int varid, const vector<double> &zlevels_out, LocalInterpCtx *lic, Vec g) const;
private:
  int compute_block_size(vector<int> count) const;
  PetscErrorCode compute_start_and_count(int varid, int t_start,
                                         int x_start, int x_count,
                                         int y_start, int y_count,
                                         int z_start, int z_count,
                                         vector<int> &start,
                                         vector<int> &count, 
                                         vector<int> &imap) const;

  bool check_dimensions() const;
  int k_below(double z, const vector<double> &zlevels) const;
  IceGrid* grid;

  int event_write,
    event_write_proc0,
    event_write_send_and_receive;
};

#endif	// __PISMIO
