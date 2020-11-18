// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2019 PISM Authors
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

#include "cdi.h"
#include "cdipio.h"
#include <yaxt.h>

#include <sstream>
#include <string.h>
#include <map>

#include "CDI.hh"

#include "pism/util/io/pism_cdi_type_conversion.hh"

#include "pism/util/error_handling.hh"

namespace pism {
namespace io {

CDI::CDI(MPI_Comm c) : NCFile(c) {
	m_gridID = -1;
	m_vlistID = -1;
	m_zID = -1;
	m_zbID = -1;
	m_zsID = -1;
	m_tID = -1;
}

CDI::~CDI() {
}

void CDI::open_impl(const std::string &fname, IO_Mode mode) {
	if (mode == PISM_READONLY) {
        m_file_id = streamOpenRead(fname.c_str());
	} else {
        m_file_id = streamOpenWrite(fname.c_str(), CDI_FILETYPE_NC4);
	}
}

void CDI::create_impl(const std::string &filename) {
	int mode = CDI_FILETYPE_NC4;
	m_file_id = streamOpenWrite(filename.c_str(), mode);
}

void CDI::sync_impl() const {
}

void CDI::close_impl() {
	streamClose(m_file_id);
	m_file_id = -1;
	destroy_objs();
}

void CDI::destroy_objs() {
	vlistDestroy(m_vlistID);
	taxisDestroy(m_tID);
	zaxisDestroy(m_zID);
	zaxisDestroy(m_zbID);
	zaxisDestroy(m_zsID);
	gridDestroy(m_gridID);
	m_varsID.clear();
}

void CDI::enddef_impl() const {
}

void CDI::redef_impl() const {
}

void CDI::def_dim_impl(const std::string &name, size_t length) const {
	if (m_gridID != -1) {
		if (strcmp(name.c_str(),"x")==0) {
			gridDefXsize(m_gridID, (int)length);
			gridDefXname(m_gridID, name.c_str());
		} else if (strcmp(name.c_str(),"y")==0) {
			gridDefYsize(m_gridID, (int)length);
			gridDefYname(m_gridID, name.c_str());
		} else if (strcmp(name.c_str(),"z")==0 && m_zID==-1) {
			m_zID = zaxisCreate(ZAXIS_GENERIC, (int)length);
			zaxisDefName(m_zID,name.c_str());
		} else if (strcmp(name.c_str(),"zb")==0 && m_zbID==-1) {
			m_zbID = zaxisCreate(ZAXIS_GENERIC, (int)length);
			zaxisDefName(m_zbID,name.c_str());
		} else if (strcmp(name.c_str(),"time")==0 && m_tID==-1) {
			m_tID = taxisCreate(TAXIS_ABSOLUTE);
		}
		if (m_zsID == -1) {
			m_zsID = zaxisCreate(ZAXIS_SURFACE, 1);
		}
	}
}

void CDI::inq_dimid_impl(const std::string &dimension_name, bool &exists) const {
	(void) dimension_name;
	(void) exists;
}

void CDI::inq_dimlen_impl(const std::string &dimension_name, unsigned int &result) const {
	if (strcmp(dimension_name.c_str(),"x")==0) {
		result = gridInqXsize(m_gridID);
	} else if (strcmp(dimension_name.c_str(),"y")==0) {
		result = gridInqYsize(m_gridID);
	} else if (strcmp(dimension_name.c_str(),"z")==0) {
		result = zaxisInqSize(m_zID);
	} else if (strcmp(dimension_name.c_str(),"zb")==0) {
		result = zaxisInqSize(m_zbID);
	}
}

void CDI::inq_unlimdim_impl(std::string &result) const {
	(void) result;
}

void CDI::def_var_impl(const std::string &name, IO_Type nctype, const std::vector<std::string> &dims) const {
	if (m_vlistID == -1) {
		m_vlistID = vlistCreate();
	}
	int zaxisID = -1;
	int tsteptype = -1;

	for (auto d : dims) {
		if (strcmp(d.c_str(),"z")==0) {
			zaxisID = m_zID;
		} else if (strcmp(d.c_str(),"zb")==0) {
			zaxisID = m_zbID;
		}
		if (strcmp(d.c_str(),"time")==0) {
			tsteptype = TIME_VARYING;
		} else {
			tsteptype = TIME_CONSTANT;
		}
    }
    if (zaxisID == -1) {
    	zaxisID = m_zsID;
    }

	int varID = vlistDefVar(m_vlistID, m_gridID, zaxisID, tsteptype);
	vlistDefVarName(m_vlistID, varID, name.c_str());
	vlistDefVarDatatype(m_vlistID, varID, pism_type_to_cdi_type(nctype));
	m_varsID[name] = varID;
}

void CDI::get_vara_double_impl(const std::string &variable_name,
                                 const std::vector<unsigned int> &start,
                                 const std::vector<unsigned int> &count,
                                 double *ip) const {
	(void) variable_name;
	(void) start;
	(void) count;
	(void) ip;
}

void CDI::put_vara_double_impl(const std::string &variable_name,
                                  const std::vector<unsigned int> &start,
                                  const std::vector<unsigned int> &count,
                                  const double *op) const {
	(void) variable_name;
	(void) start;
	(void) count;
	(void) op;
}

void CDI::get_varm_double_impl(const std::string &variable_name,
                                 const std::vector<unsigned int> &start,
                                 const std::vector<unsigned int> &count,
                                 const std::vector<unsigned int> &imap,
                                 double *ip) const {
	(void) variable_name;
	(void) start;
	(void) count;
	(void) imap;
	(void) ip;
}

void CDI::inq_nvars_impl(int &result) const {
	result = m_varsID.size();
}

void CDI::inq_vardimid_impl(const std::string &variable_name, std::vector<std::string> &result) const {
	(void) variable_name;
	(void) result;
}

void CDI::inq_varnatts_impl(const std::string &variable_name, int &result) const {
	cdiInqNatts(m_vlistID, m_varsID[variable_name], &result);
}

void CDI::inq_varid_impl(const std::string &variable_name, bool &exists) const {
	exists = m_varsID.count(variable_name) > 0;
}

void CDI::inq_varname_impl(unsigned int j, std::string &result) const {
	for (auto &i : m_varsID) {
      if (i.second == j) {
         result = i.first;
         break; // to stop searching
      }
   }
}

void CDI::get_att_double_impl(const std::string &variable_name,
                                  const std::string &att_name,
                                  std::vector<double> &result) const {
	(void) variable_name;
	(void) att_name;
	(void) result;
}

void CDI::get_att_text_impl(const std::string &variable_name, const std::string &att_name, std::string &result) const {
	(void) variable_name;
	(void) att_name;
	(void) result;
}

void CDI::del_att_impl(const std::string &variable_name, const std::string &att_name) const {
	(void) variable_name;
	(void) att_name;
}

void CDI::put_att_double_impl(const std::string &variable_name, const std::string &att_name, IO_Type nctype, const std::vector<double> &data) const {
	cdiDefAttFlt(m_vlistID, m_varsID[variable_name], att_name.c_str(), CDI_DATATYPE_FLT64, data.size(), &data[0]);
}

void CDI::put_att_text_impl(const std::string &variable_name,
                                const std::string &att_name,
                                const std::string &value) const {
	cdiDefAttTxt(m_vlistID, m_varsID[variable_name], att_name.c_str(), value.size(), value.c_str());
}

void CDI::inq_attname_impl(const std::string &variable_name,
                               unsigned int n,
                               std::string &result) const {
	(void) variable_name;
	(void) n;
	(void) result;
}

void CDI::inq_atttype_impl(const std::string &variable_name,
                               const std::string &att_name,
                               IO_Type &result) const {
	(void) variable_name;
	(void) att_name;
	(void) result;
}

void CDI::set_fill_impl(int fillmode, int &old_modep) const {
	(void) fillmode;
	(void) old_modep;
}

void CDI::create_grid_impl(int lengthx, int lengthy) const {
	m_gridID = gridCreate(GRID_GENERIC, lengthx*lengthy);
}

void CDI::define_timestep_impl(int tsID) const {
	streamDefTimestep(m_file_id, tsID);
}

void CDI::write_darray_impl(const std::string &variable_name,
                                   const IceGrid &grid,
                                   unsigned int z_count,
                                   unsigned int record,
                                   const double *input) {
  
  // ATTENTION: need to transpose the input before writing: not done yet
  int varid = m_varsID[variable_name];
  Xt_idxlist decomp = grid.yaxt_decomposition(z_count);
  size_t nmiss = 0;
  
  streamWriteVarPart(m_file_id, varid, input, nmiss, decomp);
}


} // end of namespace io
} // end of namespace pism
