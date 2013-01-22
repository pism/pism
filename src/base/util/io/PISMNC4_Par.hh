// Copyright (C) 2012, 2013 PISM Authors
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

#ifndef _PISMNC4_PAR_H_
#define _PISMNC4_PAR_H_

#include "PISMNC4File.hh"

class PISMNC4_Par : public PISMNC4File
{
public:
  PISMNC4_Par(MPI_Comm c, int r)
    : PISMNC4File(c, r, 0) {}
  virtual ~PISMNC4_Par() {}

  // open/create/close
  virtual int open(string filename, int mode);

  virtual int create(string filename);

protected:
  virtual int set_access_mode(int varid, bool mapped) const;
};


#endif /* _PISMNC4_PAR_H_ */
