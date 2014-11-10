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

#ifndef _PISMNC4FILE_1FPP_H_
#define _PISMNC4FILE_1FPP_H_

#include "PISMNC4File.hh"

namespace pism {

class NC4_Quilt : public NC4File
{
public:
  NC4_Quilt(MPI_Comm c, unsigned int compression_level)
    : NC4File(c, compression_level), suffix("_patch")
  {
  }
  virtual ~NC4_Quilt() {}

protected:
  // implementations:

  // open/create/close
  int open_impl(const std::string &filename, IO_Mode mode);

  int create_impl(const std::string &filename);

  int close_impl();

  // dim
  int def_dim_impl(const std::string &name, size_t length) const;

  // var
  int def_var_impl(const std::string &name, IO_Type nctype, const std::vector<std::string> &dims) const;
  // att
  using NCFile::put_att_double_impl;
  int put_att_double_impl(const std::string &variable_name, const std::string &att_name, IO_Type xtype, const std::vector<double> &data) const;

  int put_att_text_impl(const std::string &variable_name, const std::string &att_name, const std::string &value) const;

  int move_if_exists_impl(const std::string &filename, int rank_to_use = 0);

private:
  int integer_open_mode(IO_Mode input) const;

  int get_put_var_double(const std::string &variable_name,
                         const std::vector<unsigned int> &start,
                         const std::vector<unsigned int> &count,
                         const std::vector<unsigned int> &imap, double *ip,
                         bool get,
                         bool mapped) const;
  
  void correct_start_and_count(const std::string &name,
                               std::vector<unsigned int> &start,
                               std::vector<unsigned int> &count) const;
  int global_stat(int stat) const;

  const std::string suffix;
};

} // end of namespace pism

#endif /* _PISMNC4FILE_1FPP_H_ */
