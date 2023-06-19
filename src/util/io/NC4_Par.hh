// Copyright (C) 2012, 2013, 2014, 2015, 2017, 2019, 2020, 2023 PISM Authors
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

#ifndef _PISMNC4_PAR_H_
#define _PISMNC4_PAR_H_

#include "pism/util/io/NC4File.hh"

namespace pism {
namespace io {

class NC4_Par : public NC4File
{
public:
  NC4_Par(MPI_Comm c)
    : NC4File(c, 0) {}
  virtual ~NC4_Par() {}
protected:
  // open/create/close
  virtual void open_impl(const std::string &filename, io::Mode mode);

  virtual void create_impl(const std::string &filename);

  virtual void set_access_mode(int varid, bool mapped) const;

  virtual void set_compression_level_impl(int level) const;
};


} // end of namespace io
} // end of namespace pism

#endif /* _PISMNC4_PAR_H_ */
