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

#ifndef _PISMNC4FILE_1FPP_H_
#define _PISMNC4FILE_1FPP_H_

#include "PISMNC4File.hh"

class PISMNC4_Quilt : public PISMNC4File
{
public:
  PISMNC4_Quilt(MPI_Comm c, int r, unsigned int compression_level)
    : PISMNC4File(c, r, compression_level), suffix("_patch")
  {
  }
  virtual ~PISMNC4_Quilt() {}

  // open/create/close
  virtual int open(string filename, int mode);

  virtual int create(string filename);

  virtual int close();

  // dim
  virtual int def_dim(string name, size_t length) const;

  // var
  virtual int def_var(string name, PISM_IO_Type nctype, vector<string> dims) const;
  // att
  using PISMNCFile::put_att_double;
  virtual int put_att_double(string variable_name, string att_name, PISM_IO_Type xtype, vector<double> &data) const;

  virtual int put_att_text(string variable_name, string att_name, string value) const;

  virtual int move_if_exists(string filename, int rank_to_use = 0);
private:
  virtual int get_put_var_double(string variable_name,
                                 vector<unsigned int> start,
                                 vector<unsigned int> count,
                                 vector<unsigned int> imap, double *ip,
                                 bool get,
                                 bool mapped) const;
  
  void correct_start_and_count(string name,
                               vector<unsigned int> &start,
                               vector<unsigned int> &count) const;
  int global_stat(int stat) const;

  const string suffix;
};

#endif /* _PISMNC4FILE_1FPP_H_ */
