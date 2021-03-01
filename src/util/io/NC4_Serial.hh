// Copyright (C) 2012, 2013, 2014, 2015, 2019, 2020, 2021 PISM Authors
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

#ifndef _NC4_Serial_H_
#define _NC4_Serial_H_

#include "NC_Serial.hh"

namespace pism {
namespace io {

class NC4_Serial : public NC_Serial
{
public:
  NC4_Serial(MPI_Comm c);
  virtual ~NC4_Serial();

protected:
  void set_compression_level_impl(int level) const;

  void create_impl(const std::string &filename, int FileID = -1, int filetype);
};

  void def_var_impl(const std::string &name, IO_Type nctype,
                    const std::vector<std::string> &dims) const;

  mutable int m_compression_level;
};

} // end of namespace io
} // end of namespace pism

#endif /* _NC4_Serial_H_ */
