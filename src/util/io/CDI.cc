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
#include <mpi.h>
#include <sstream>
#include <string.h>
#include <algorithm>
#include <map>
 #include <iostream>

#include "CDI.hh"

#include "pism/util/io/pism_cdi_type_conversion.hh"
#include "pism/util/IceGrid.hh"

#include "pism/util/error_handling.hh"

extern "C"{
#include "cdipio.h"
#include "cdi.h"
#include "yaxt.h"
}

namespace pism {
namespace io {

CDI::CDI(MPI_Comm c) : NCFile(c) {
//	m_gridID = -1;
//	m_gridsID = -1;
	m_vlistID = -1;
//	m_zID = -1;
//	m_zbID = -1;
//	m_zsID = -1;
//	m_tID = -1;
	m_varsID.clear();
	m_beforediag = true;
	m_gridexist = false;
	m_istimedef = false;
}

CDI::~CDI() {
}

void CDI::open_impl(const std::string &fname, IO_Mode mode, const std::map<std::string, int> &varsi) {
	if (mode == PISM_READONLY) {
        m_file_id = streamOpenRead(fname.c_str());
	} else {
        m_file_id = streamOpenAppend(fname.c_str());
        m_vlistID = streamInqVlist(m_file_id);
	m_tID = vlistInqTaxis(m_vlistID);
        m_varsID = varsi;
	}
	m_firststep = false;
}

void CDI::set_ncgridIDs_impl(const std::vector<int>& gridIDs) const {
	if (gridIDs.empty()) {
		m_gridID = -1;
		m_gridsID = -1;
		m_tID = -1;
		m_zID = -1;
		m_zbID = -1;
		m_zsID = -1;		
	} else {
		m_gridID = gridIDs[0];
		m_gridsID = gridIDs[1];
		m_tID = gridIDs[2];
		m_zID = gridIDs[3];
		m_zbID = gridIDs[4];
		m_zsID = gridIDs[5];
		m_gridexist = true;
	}
}

std::vector<int> CDI::get_ncgridIDs_impl() const {
	std::vector<int> gridIDs{m_gridID, m_gridsID, m_tID, m_zID, m_zbID, m_zsID};
	return gridIDs;
}


void CDI::create_impl(const std::string &filename) {
	int mode = CDI_FILETYPE_NC4;
	m_file_id = streamOpenWrite(filename.c_str(), mode);
	m_firststep = true;
}

void CDI::sync_impl() const {
}

void CDI::close_impl() {
	//double time_start = MPI_Wtime();
	//streamClose(m_file_id);
	//double time_end = MPI_Wtime();
	//std::cout << "streamClose() time = " << time_end-time_start;
	m_file_id = -1;
	destroy_objs();
}

void CDI::destroy_objs() {
	m_varsID.clear();
}

void CDI::enddef_impl() const {
//	pioEndDef(); concept of STAGES does not exist anymore
}

void CDI::redef_impl() const {
}

void CDI::def_dim_impl(const std::string &name, size_t length) const {
	if (m_vlistID == -1) {
        	m_vlistID = vlistCreate();
    	}

	if (!m_gridexist) {
	if (strcmp(name.c_str(),"x")==0) {
		gridDefXsize(m_gridID, (int)length);
		gridDefXname(m_gridID, name.c_str());
                m_dims_name.push_back(name);
	} else if (strcmp(name.c_str(),"y")==0) {
		gridDefYsize(m_gridID, (int)length);
		gridDefYname(m_gridID, name.c_str());
		m_dims_name.push_back(name);
	} else if (strcmp(name.c_str(),"z")==0 && m_zID==-1) {
		m_zID = zaxisCreate(ZAXIS_GENERIC, (int)length);
		zaxisDefName(m_zID,name.c_str());
		m_dims_name.push_back(name);
	} else if (strcmp(name.c_str(),"zb")==0 && m_zbID==-1) {
		m_zbID = zaxisCreate(ZAXIS_GENERIC, (int)length);
		zaxisDefName(m_zbID,name.c_str());
		m_dims_name.push_back(name);
	} else if (strcmp(name.c_str(),"time")==0 && m_tID==-1) {
		m_tID = taxisCreate(TAXIS_ABSOLUTE);
		vlistDefTaxis(m_vlistID, m_tID);
		m_dims_name.push_back(name);
	}
	if (m_zsID == -1) {
		m_zsID = zaxisCreate(ZAXIS_SURFACE, 1);
	}
	} else {
		if (!m_istimedef) {
			m_dims_name.push_back("x");
			m_dims_name.push_back("y");
			m_dims_name.push_back("z");
			m_dims_name.push_back("zb");
			m_dims_name.push_back("time");
			vlistDefTaxis(m_vlistID, m_tID);
			m_istimedef = true;
		}
	}
}

void CDI::def_ref_date_impl(double time) const {
	// Taking into account only standard calendar right now
	double nyearsf = - time / 365 / 24 / 60 / 60;
	taxisDefVdate(m_tID, -(long int)nyearsf);
	long int seconds = (nyearsf - (long int)nyearsf) * 86400;
	long int minutes, hours;
	hours = seconds / 3600;
	minutes = (seconds - (3600*hours)) / 60;
	seconds = (seconds - (3600*hours) - minutes*60);
	long int daytime = hours * 10000 + minutes * 100 + seconds;
	taxisDefVtime(m_tID, daytime);
}

void CDI::inq_dimid_impl(const std::string &dimension_name, bool &exists) const {
	if (strcmp(dimension_name.c_str(),"x")==0 || strcmp(dimension_name.c_str(),"y")==0) {
		exists = m_gridID != -1 ? true : false;
	} else if (strcmp(dimension_name.c_str(),"z")==0) {
		exists = m_zID != -1 ? true : false;
	} else if (strcmp(dimension_name.c_str(),"zb")==0) {
		exists = m_zbID != -1 ? true : false;
	} else if (strcmp(dimension_name.c_str(),"time")==0) {
		exists = m_tID != -1 ? true : false;
	}
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
	} else if (strcmp(dimension_name.c_str(),"time")==0) {
		result = inq_current_timestep() + 1;
	}
}

int CDI::inq_current_timestep() const {
	int timesID = -1, nrec = -1;
	if (m_firststep) {
		timesID = 0;
	} else { 
	while (nrec != 0) {
		timesID++;
		nrec = streamInqTimestep(m_file_id, timesID);
	}
	}
	return timesID;
}

void CDI::inq_unlimdim_impl(std::string &result) const {
	result = "time"; // limitation of CDI: cannot set time variable name
}

void CDI::def_var_impl(const std::string &name, IO_Type nctype, const std::vector<std::string> &dims) const {
    // No need to define the dimensions as variables
    if (std::find(m_dims_name.begin(), m_dims_name.end(), name) != m_dims_name.end())
    {
        return;
    }
    
    // Define variables
    if (m_vlistID == -1) {
	m_vlistID = vlistCreate();
    }

    if (dims.empty() || dims.size()<2) { // scalar variable
		def_var_scalar_impl(name, nctype, dims);
    } else {         // multi-dimensional variable
		def_var_multi_impl(name, nctype, dims);
    }
}

void CDI::def_var_scalar_impl(const std::string &name, IO_Type nctype,  const std::vector<std::string> &dims) const {
    if (m_gridsID == -1) {
        m_gridsID = gridCreate(GRID_LONLAT, 1);
        gridDefXsize(m_gridsID, 1);
        gridDefYsize(m_gridsID, 0);
	gridDefXname(m_gridsID, "xd");
    }
    if (m_zsID == -1) {
        m_zsID = zaxisCreate(ZAXIS_SURFACE, 1);
    }

    int tsteptype;
    if (dims.empty()) {
        tsteptype = TIME_CONSTANT;
    } else {
        for (auto d : dims) {
            if (strcmp(d.c_str(),"time")==0) tsteptype = TIME_VARYING;
        }
    }

    int varID = vlistDefVar(m_vlistID, m_gridsID, m_zsID, tsteptype);
    vlistDefVarName(m_vlistID, varID, name.c_str());
    int type = pism_type_to_cdi_type(nctype);
    vlistDefVarDatatype(m_vlistID, varID, type);
    m_varsID[name] = varID;
}

void CDI::def_var_multi_impl(const std::string &name, IO_Type nctype, const std::vector<std::string> &dims) const {
    int zaxisID = -1;
    int tsteptype = TIME_CONSTANT;

    for (auto d : dims) {
        if (strcmp(d.c_str(),"z")==0) {
            zaxisID = m_zID;
        } else if (strcmp(d.c_str(),"zb")==0) {
            zaxisID = m_zbID;
        } else if (strcmp(d.c_str(),"time")==0) {
            tsteptype = TIME_VARYING;
    	}
    }
    if (zaxisID == -1) {
        zaxisID = m_zsID;
    }

    int varID = vlistDefVar(m_vlistID, m_gridID, zaxisID, tsteptype);
    vlistDefVarName(m_vlistID, varID, name.c_str());
    int type = pism_type_to_cdi_type(nctype);
    vlistDefVarDatatype(m_vlistID, varID, type);
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
	// write dimensions values and scalar variables
	if (!m_gridexist) {	
	if (strcmp(variable_name.c_str(),"x")==0) {
		gridDefXvals(m_gridID, op);
	} else if (strcmp(variable_name.c_str(),"y")==0) {
		gridDefYvals(m_gridID, op);
	} else if (strcmp(variable_name.c_str(),"z")==0) {
		zaxisDefLevels(m_zID, op);
	} else if (strcmp(variable_name.c_str(),"zb")==0) {
		zaxisDefLevels(m_zbID, op);
	}
	}
	// need to be generalised to write scalars
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
	int varID = m_varsID[variable_name];
	int type = 0;
	int current_grid = vlistInqVarGrid(m_vlistID, varID);
	if (current_grid == m_gridID) {
		int current_z = vlistInqVarZaxis(m_vlistID, varID);
		if (current_z == m_zID) {
			type = type + 2;
		} else if (current_z == m_zbID) {
			type = type + 4;
		}
		int current_t = vlistInqVarTimetype(m_vlistID, varID);
		if (current_t == TIME_VARYING) {
			type++;
		}
	} else {
		result.clear();
		return;
	}
    // todo: inquire the name - don't write it directly
	switch (type) {
		case 0:
			result.resize(2); result[0] = "y"; result[1] = "x";
			break;
		case 1:
			result.resize(3); result[0] = "time"; result[1] = "y"; result[2] = "x";
			break;
		case 2:
			result.resize(3); result[0] = "z"; result[1] = "y"; result[2] = "x";
			break;
		case 3:
			result.resize(4); result[0] = "time"; result[1] = "z"; result[2] = "y"; result[3] = "x";
			break;
		case 4:
			result.resize(3); result[0] = "zb"; result[1] = "y"; result[2] = "x";
			break;
		case 5:
			result.resize(4); result[0] = "time"; result[1] = "zb"; result[2] = "y"; result[3] = "x";
			break;
	}
}

void CDI::inq_varnatts_impl(const std::string &variable_name, int &result) const {
	int varID = -1;
	if (variable_name == "PISM_GLOBAL") {
		varID = CDI_GLOBAL;
	} else {
		varID = m_varsID[variable_name];
	}
	cdiInqNatts(m_vlistID, varID, &result);
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
	int varID = -1;
	if (variable_name == "PISM_GLOBAL") {
		varID = CDI_GLOBAL;
	} else {
		varID = m_varsID[variable_name];
	}
	int cdilen = cdiInqAttLen(m_vlistID, varID, att_name.c_str());

	if (cdilen == 0) {
    	result.clear();
    	return;
  	}
  	result.resize(cdilen);

  	cdiInqAttFlt(m_vlistID, varID, att_name.c_str(), cdilen, result.data());
}

void CDI::get_att_text_impl(const std::string &variable_name, const std::string &att_name, std::string &result) const {
	int varID = -1;
	if (variable_name == "PISM_GLOBAL") {
		varID = CDI_GLOBAL;
	} else {
		varID = m_varsID[variable_name];
	}
	int cdilen = cdiInqAttLen(m_vlistID, varID, att_name.c_str());
	if (cdilen == -1 || cdilen==0) {
		result.clear();
    		return;
	}
	result.resize(cdilen);
	std::vector<char> str(cdilen + 1, 0);
	cdiInqAttTxt(m_vlistID, varID, att_name.c_str(), cdilen, str.data());
	result = str.data();
}

void CDI::del_att_impl(const std::string &variable_name, const std::string &att_name) const {
	int varID = -1;
	if (variable_name == "PISM_GLOBAL") {
		varID = CDI_GLOBAL;
	} else {
		varID = m_varsID[variable_name];
	}
	cdiDelAtt(m_vlistID, varID, att_name.c_str());
}

void CDI::put_att_double_impl(const std::string &variable_name, const std::string &att_name, IO_Type nctype, const std::vector<double> &data) const {
	if (std::find(m_dims_name.begin(), m_dims_name.end(), variable_name) != m_dims_name.end())
    	{
        	return;
    	}
	int varID = -1;
	if (variable_name == "PISM_GLOBAL") {
		varID = CDI_GLOBAL;
	} else {
		varID = m_varsID[variable_name];
	}
	int type = vlistInqVarDatatype(m_vlistID,varID);
	if (type != CDI_DATATYPE_FLT32 && type != CDI_DATATYPE_FLT64) {
		type = CDI_DATATYPE_FLT32;
	}
	cdiDefAttFlt(m_vlistID, varID, att_name.c_str(), type, data.size(), &data[0]);
}

void CDI::put_att_text_impl(const std::string &variable_name,
                                const std::string &att_name,
                                const std::string &value) const {
        if (std::find(m_dims_name.begin(), m_dims_name.end(), variable_name) != m_dims_name.end())
        {
                return;
        }
	int varID = -1;
	if (variable_name == "PISM_GLOBAL") {
		varID = CDI_GLOBAL;
	} else {
		varID = m_varsID[variable_name];
	}
	cdiDefAttTxt(m_vlistID, varID, att_name.c_str(), value.size(), value.c_str());
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
	int varID = -1;
	if (variable_name == "PISM_GLOBAL") {
		varID = CDI_GLOBAL;
	} else {
		varID = m_varsID[variable_name];
	}
	int cditype = cdiInqAttType(m_vlistID, varID, att_name.c_str());
	result = cdi_type_to_pism_type(cditype);
}

void CDI::set_fill_impl(int fillmode, int &old_modep) const {
	(void) fillmode;
	(void) old_modep;
}

void CDI::create_grid_impl(int lengthx, int lengthy) const {
	if (m_gridID ==-1)
		m_gridID = gridCreate(GRID_LONLAT, lengthx*lengthy);
}

void CDI::define_timestep_impl(int tsID) const {
	streamDefTimestep(m_file_id, tsID);
	if (tsID==0) m_firststep = false;
}

void CDI::write_timestep_impl() const {
	pioWriteTimestep();
}

void CDI::write_darray_impl(const std::string &variable_name,
                                   const IceGrid &grid,
                                   unsigned int z_count,
                                   unsigned int record,
                                   const double *input) {
	// transpose input data
	int dim = grid.local_length((int)z_count);
	double *inputIO;
	inputIO = (double*) malloc(dim * sizeof(double));
	grid.io_transpose(input, inputIO, (int)z_count);

	int varid = m_varsID[variable_name];
	
	// create decomposition if new
	Xt_idxlist decompid = grid.yaxt_decomposition((int)z_count);
	size_t nmiss = 0;
	if (!m_beforediag || m_diagvars.count(variable_name)==0) {
		streamWriteVarPart(m_file_id, varid, inputIO, nmiss, decompid);
	}
}

std::map<std::string, int> CDI::get_var_map_impl() {
	return m_varsID;
}

void CDI::def_vlist_impl() const {
	streamDefVlist(m_file_id, m_vlistID);
}

void CDI::set_diagvars_impl(const std::set<std::string> &variables) const {
	m_diagvars = variables;
}

void CDI::set_bdiag_impl(bool value) const {
	m_beforediag = value;
}

int CDI::get_ncstreamID_impl() const {
	return m_file_id;
}

int CDI::get_ncvlistID_impl() const {
	return m_vlistID;
}

} // end of namespace io
} // end of namespace pism
