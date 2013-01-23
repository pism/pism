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

// Defining H5_NO_DEPRECATED_SYMBOLS is supposed to be enough
// to select the v1.8 HDF5 API even if the library was built
// with v1.6 API set as the default...
#define H5_NO_DEPRECATED_SYMBOLS 1

// but we need the following seven lines because of a bug in H5version.h. :-(
#define H5Dopen_vers 2
#define H5Dcreate_vers 2
#define H5Eprint_vers 2
#define H5Eclear_vers 2
#define H5Eget_auto_vers 2
#define H5Eset_auto_vers 2
#define H5E_auto_t_vers 2

#include <hdf5.h>
#include <hdf5_hl.h>


#include <assert.h>
#include <cmath>
#include <stdlib.h>

#define TEMPORARY_STRING_LENGTH 32768

#include "PISMNC4_HDF5.hh"

static hid_t pism_type_to_hdf5_type(PISM_IO_Type xtype) {
  switch(xtype) {
  case PISM_BYTE:
  case PISM_CHAR:
    return H5T_NATIVE_CHAR;
  case PISM_FLOAT:
    return H5T_NATIVE_FLOAT;
  case PISM_DOUBLE:
    return H5T_NATIVE_DOUBLE;
  case PISM_SHORT:
    return H5T_NATIVE_SHORT;
  default:
  case PISM_INT:
    return H5T_NATIVE_INT;
  }
}

// Callback functions

// This functions finds an inlimited dimension. (We assume that there is only one.)
static herr_t find_unlimdim(hid_t loc_id, const char *name, const H5L_info_t *linfo, void *opdata) {
  string *dim_name = (string*)opdata;
  int return_value = 0;
  herr_t stat;
  H5O_info_t object_info;

  (void)linfo;

  stat = H5Oget_info_by_name(loc_id, name, &object_info, H5P_DEFAULT);
  assert(stat >= 0);

  if (object_info.type != H5O_TYPE_DATASET)
    return 0;

  hid_t did = H5Dopen(loc_id, name, H5P_DEFAULT);
  assert(did >= 0);

  int is_scale = H5DSis_scale(did);
  assert(is_scale >= 0);

  if (is_scale > 0) {
    hid_t space_id = H5Dget_space(did);
    assert(space_id >= 0);

    int rank = H5Sget_simple_extent_ndims(space_id);
    assert(rank >= 0);
    if (rank == 1) {
      hsize_t max_dim;

      stat = H5Sget_simple_extent_dims(space_id, NULL, &max_dim);

      if (max_dim == H5S_UNLIMITED) {
        (*dim_name) = name;
        return_value = 1;
      }
    }
    H5Sclose(space_id);
  }

  H5Dclose(did);

  return return_value;
}

// This function builds the list of all the dimensions in a file.
static herr_t get_dimensions(hid_t loc_id, const char *name, const H5L_info_t *linfo, void *opdata)
{
  vector<string> *dim_names = (vector<string>*)opdata;
  herr_t stat;
  H5O_info_t object_info;

  (void)linfo;

  stat = H5Oget_info_by_name(loc_id, name, &object_info, H5P_DEFAULT);
  assert(stat >= 0);

  if (object_info.type != H5O_TYPE_DATASET)
    return 0;

  hid_t did = H5Dopen(loc_id, name, H5P_DEFAULT);
  assert(did >= 0);

  int is_scale = H5DSis_scale(did);
  assert(is_scale >= 0);

  if (is_scale > 0)
    dim_names->push_back(name);

  H5Dclose(did);

  return 0;
}

/* attribute type of a DS dataset */
typedef struct ds_list_t {
  hobj_ref_t ref;     /* object reference  */
  int        dim_idx; /* dimension index of the dataset */
} ds_list_t;

//! \brief Creates the data type for the REFERENCE_LIST attribute.
/*!
 * This function was copied from the HDF5 library sources.
 */
static hid_t H5DS_get_REFLIST_type(void)
{
  hid_t ntid_t = -1;

  /* Build native type that corresponds to compound datatype
     used to store ds_list_t structure in the REFERENCE_LIST
     attribute */

  if((ntid_t = H5Tcreate(H5T_COMPOUND, sizeof(ds_list_t))) < 0)
    goto out;

  if(H5Tinsert(ntid_t, "dataset", HOFFSET(ds_list_t,ref), H5T_STD_REF_OBJ) < 0)
    goto out;

  if(H5Tinsert(ntid_t, "dimension", HOFFSET(ds_list_t, dim_idx), H5T_NATIVE_INT) < 0)
    goto out;

  return ntid_t;
 out:
  H5E_BEGIN_TRY {
    H5Tclose(ntid_t);
  } H5E_END_TRY;
  return -1;
}

/*!
 * All tuning parameters are here and in create_file_access_plist().
 */
static unsigned int compute_stripe_size(MPI_Comm com, int xm, int ym) {
  unsigned int stripe_size;

  assert(xm >= 0 && ym >= 0);

  int max_xm, max_ym;
  MPI_Allreduce(&xm, &max_xm, 1, MPI_INT, MPI_MAX, com);
  MPI_Allreduce(&ym, &max_ym, 1, MPI_INT, MPI_MAX, com);

  // round the chunk size up to the nearest megabyte
  stripe_size = (unsigned int)ceil(static_cast<double>(sizeof(double) * max_xm * max_ym) / (1024.0 * 1024.0));

  // don't use stipes of more that 32 Mbytes
  if (stripe_size > 32)
    stripe_size = 32;

  // Convert to bytes
  stripe_size *= 1024 * 1024;

  return stripe_size;
}

/*!
 * All tuning parameters are here and in compute_stripe_size().
 */
static hid_t create_file_access_plist(MPI_Comm com, MPI_Info info,
                                      int xm, int ym) {
  hid_t plist_id;

  unsigned int stripe_size = compute_stripe_size(com, xm, ym);

  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  assert(plist_id >= 0);

  // Use MPI I/O
  H5Pset_fapl_mpio(plist_id, com, info);

  // These settings should help large parallel runs. They also make files
  // bigger, sometimes significantly.
  {
    // Allocate a 1 Mb block for metadata (I bet this is overkill.)
    H5Pset_meta_block_size(plist_id, 1024*1024);

    // Align anything larger than 16Kb to stripe boundaries
    H5Pset_alignment(plist_id, 16*1024, stripe_size);
  }

  // Disable metadata cache evictions. This may not be necessary.
  if (true) {
    H5AC_cache_config_t mdc_config;

    mdc_config.version = H5AC__CURR_CACHE_CONFIG_VERSION;
    H5Pget_mdc_config(plist_id, &mdc_config);

    mdc_config.evictions_enabled = false;
    mdc_config.incr_mode = H5C_incr__off;
    mdc_config.decr_mode = H5C_decr__off;
    mdc_config.flash_incr_mode = H5C_flash_incr__off;

    H5Pset_mdc_config(plist_id, &mdc_config);
  }

  // Set MPI hints
  {
    char hint_value[100];
    snprintf(hint_value, 100, "%d", stripe_size);
    MPI_Info_set(info, const_cast<char*>("striping_unit"), hint_value);

    MPI_Info_set(info, const_cast<char*>("striping_factor"), const_cast<char*>("-1"));
  }

  return plist_id;
}

//! \brief Get the list of dimensions a variable depends on.
static herr_t inq_dimensions(hid_t dsid, vector<string> &dims) {
  hid_t space_id, attr_id, tid;
  int rank;
  hvl_t *buf = NULL;
  int stat;

  dims.clear();

  int is_scale = H5DSis_scale(dsid);
  if (is_scale > 0) {
    char name[TEMPORARY_STRING_LENGTH];
    stat = H5DSget_scale_name(dsid, name, TEMPORARY_STRING_LENGTH);

    dims.push_back(name);

    return 0;
  }

  space_id = H5Dget_space(dsid);
  assert(space_id >= 0);

  // get the number of dimensions
  rank = H5Sget_simple_extent_ndims(space_id);
  assert(rank >= 0);

  stat = H5Sclose(space_id);
  assert(stat >= 0);

  attr_id = H5Aopen(dsid, DIMENSION_LIST, H5P_DEFAULT);
  assert(attr_id >= 0);

  tid = H5Aget_type(attr_id);
  assert(tid >= 0);

  // space_id is needed by the H5Dvlen_reclaim() call below.
  space_id = H5Aget_space(attr_id);
  assert(space_id >= 0);

  buf = (hvl_t *)malloc((size_t)rank * sizeof(hvl_t));
  if(buf == NULL)
    goto out;

  stat = H5Aread(attr_id, tid, buf);
  assert(stat >= 0);

  for (int ii = 0; ii < rank; ++ii) {
    assert(buf[ii].len == 1);

    size_t len = H5Rget_name(dsid, H5R_OBJECT, buf[ii].p, NULL, 0);

    string name;
    name.resize(len + 2);

    H5Rget_name(dsid, H5R_OBJECT, buf[ii].p, &name[0], len + 1);

    // remove the leading "/" and save the name
    dims.push_back(&name[1]);
  }

  stat = H5Dvlen_reclaim(tid, space_id, H5P_DEFAULT, buf);
  assert(stat >= 0);

  free(buf);

 out:
  H5Tclose(tid);
  H5Sclose(space_id);
  H5Aclose(attr_id);

  return 0;
}


//! \brief Extend a variable along an unlimited dimension.
/*!
 * Does *not* work if "dataset" is a coordinate variable.
 */
herr_t extend_dataset(hid_t dataset, string dimension, int increment) {
  hid_t space_id;
  int rank, stat;

  space_id = H5Dget_space(dataset);
  assert(space_id >= 0);

  rank = H5Sget_simple_extent_ndims(space_id);
  assert(rank >= 0);

  vector<hsize_t> dims(rank);
  stat = H5Sget_simple_extent_dims(space_id, &dims[0], NULL);
  assert(stat == rank);

  vector<string> dim_names;
  stat = inq_dimensions(dataset, dim_names);
  assert(stat >= 0);

  for (unsigned int j = 0; j < dim_names.size(); ++j) {
    if (dim_names[j] == dimension)
      dims[j] += increment;
  }

  stat = H5Dset_extent(dataset, &dims[0]);
  assert(stat >= 0);

  stat = H5Sclose(space_id);
  assert(stat >= 0);

  return 0;
}

//! \brief Extends a "dimension scale" (coordinate variable), including all the
//! variables that use it.
static herr_t extend_dimension(hid_t dim_id, int increment) {

  hid_t reflist_t, attr_id, attr_space_id;
  int has_reflist, nelmts;
  ds_list_t  *dsbuf = NULL;
  int stat;
  char name[TEMPORARY_STRING_LENGTH];

  // extend the coordinate variable
  {
    hsize_t dim, max_dim;
    hid_t dim_space_id = H5Dget_space(dim_id);

    stat = H5Sget_simple_extent_ndims(dim_space_id);
    assert(stat == 1);

    stat = H5Sget_simple_extent_dims(dim_space_id, &dim, &max_dim);
    assert(stat == 1);

    dim += increment;

    stat = H5Dset_extent(dim_id, &dim);
    assert(stat >= 0);

    H5Sclose(dim_space_id);
  }

  /* try to find the attribute "REFERENCE_LIST" on the dataset */
  if((has_reflist = H5LTfind_attribute(dim_id, REFERENCE_LIST)) < 0)
    return has_reflist;

  if (has_reflist == 0)
    return 0;

  stat = H5DSget_scale_name(dim_id, name, TEMPORARY_STRING_LENGTH);
  assert(stat >= 0);

  attr_id = H5Aopen(dim_id, REFERENCE_LIST, H5P_DEFAULT);
  assert(attr_id >= 0);

  if((attr_space_id = H5Aget_space(attr_id)) < 0)
    return -1;

  if((nelmts = H5Sget_simple_extent_npoints(attr_space_id)) < 0)
    return -1;

  dsbuf = (ds_list_t*) malloc((size_t)nelmts * sizeof(ds_list_t));
  if(dsbuf == NULL)
    return -1;

  reflist_t = H5DS_get_REFLIST_type();
  assert(reflist_t > 0);

  if(H5Aread(attr_id, reflist_t, dsbuf) < 0)
    goto out;

  for(int ii=0; ii<nelmts; ii++) {

    hid_t dim_id_i;

    /* get the dataset id */
    dim_id_i = H5Rdereference(dim_id, H5R_OBJECT, &dsbuf[ii].ref);

    // extend this dataset
    stat = extend_dataset(dim_id_i, name, increment);
    assert(stat >= 0);

    /* close the dereferenced dataset */
    H5Dclose(dim_id_i);
  } /* ii */

  /* close space and attribute */
  if (H5Sclose(attr_space_id) < 0)
    goto out;
  if (H5Aclose(attr_id) < 0)
    goto out;

 out:
  free(dsbuf);
  H5Tclose(reflist_t);
  return 0;
}

PISMNC4_HDF5::PISMNC4_HDF5(MPI_Comm c, int r)
  : PISMNCFile(c, r) {
  file_id = -1;
}

PISMNC4_HDF5::~PISMNC4_HDF5() {
}


// open/create/close
int PISMNC4_HDF5::open(string filename, int mode) {
  herr_t stat;
  hid_t plist_id;

  m_filename = filename;

  plist_id = create_file_access_plist(com, MPI_INFO_NULL, m_xm, m_ym);

  stat = H5Fis_hdf5(filename.c_str()); check(stat);
  if (stat == 0) {
    if (rank == 0)
      fprintf(stderr, "ERROR: %s does not exist or is not a HDF5 file.\n",
              filename.c_str());
  }

  file_id = H5Fopen(filename.c_str(), mode == PISM_NOWRITE ? H5F_ACC_RDONLY : H5F_ACC_RDWR,
                    plist_id);

  H5Pclose(plist_id);

  if (file_id >= 0)
    return 0;
  else
    return 1;
}


int PISMNC4_HDF5::create(string filename) {
  hid_t plist_id;
  MPI_Info info;


  MPI_Info_create(&info);
  plist_id = create_file_access_plist(com, info, m_xm, m_ym); check(plist_id);

  file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

  H5Pclose(plist_id);
  MPI_Info_free(&info);

  if (file_id >= 0)
    return 0;
  else
    return 1;
}


int PISMNC4_HDF5::close() {

  int stat = H5Fclose(file_id);

  file_id = -1;
  m_filename.clear();

  if (stat >= 0)
    return 0;
  else
    return 1;
}


// redef/enddef
int PISMNC4_HDF5::enddef() const {
  // A no-op.
  return 0;
}


int PISMNC4_HDF5::redef() const {
  // A no-op.
  return 0;
}


// dim
//! \brief Defines a dimension and the associated coordinate variable.
int PISMNC4_HDF5::def_dim(string name, size_t length) const {

  herr_t stat;
  hid_t dataspace_id, dim_id;

  stat = H5LTfind_dataset(file_id, name.c_str()); check(stat);

  // Check if this variable already exists and return if it does.
  if (stat > 0)
    return 0;

  if (length == PISM_UNLIMITED) {
    hsize_t zero = 0, one = 1, unlimited = H5S_UNLIMITED;
    dataspace_id = H5Screate_simple(1, &zero, &unlimited); check(dataspace_id);

    hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE); check(plist_id);
    stat = H5Pset_chunk(plist_id, 1, &one); check(stat);

    dim_id = H5Dcreate(file_id, name.c_str(), H5T_NATIVE_DOUBLE, dataspace_id,
                       H5P_DEFAULT, plist_id, H5P_DEFAULT);

    stat = H5Pclose(plist_id); check(stat);
  } else {
    hsize_t hdf5_length = length;
    dataspace_id = H5Screate_simple(1, &hdf5_length, NULL); check(dataspace_id);

    dim_id = H5Dcreate(file_id, name.c_str(), H5T_NATIVE_DOUBLE, dataspace_id,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  }
  check(dim_id);

  stat = H5DSset_scale(dim_id, name.c_str()); check(stat);

  stat = H5Sclose(dataspace_id); check(stat);
  stat = H5Dclose(dim_id); check(stat);

  if (stat >= 0)
    return 0;
  else
    return 1;
}


/*!
 * We assume that every dimension has an associated coordinate variable, so
 * checking if a dimension exist is the same as checking if the corresponding
 * variable exists.
 */
int PISMNC4_HDF5::inq_dimid(string dimension_name, bool &exists) const {
  int stat = this->inq_varid(dimension_name, exists);
  return stat;
}

/*!
 * Get the length of a dimension.
 *
 * We assume that all dimensions are one-dimensional.
 */
int PISMNC4_HDF5::inq_dimlen(string dimension_name, unsigned int &result) const {
  hid_t dim_id, space_id;
  herr_t stat;

  dim_id = H5Dopen(file_id, dimension_name.c_str(), H5P_DEFAULT); check(dim_id);

  space_id = H5Dget_space(dim_id); check(space_id);

  int variable_rank = H5Sget_simple_extent_ndims(space_id); check(variable_rank);

  if (variable_rank == 1) {
    hsize_t dim;
    stat = H5Sget_simple_extent_dims(space_id, &dim, NULL); check(stat);
    result = dim;
  } else {
    result = 0;
  }

  stat = H5Sclose(space_id); check(stat);
  stat = H5Dclose(dim_id); check(stat);

  return 0;
}

//! \brief Finds the unlimited dimension by iterating over all dimensions.
int PISMNC4_HDF5::inq_unlimdim(string &result) const {
  hsize_t idx;
  herr_t stat = H5Literate_by_name(file_id, "/", H5_INDEX_NAME, H5_ITER_INC,
                                   &idx, find_unlimdim, &result, H5P_DEFAULT);

  if (stat <= 0)
    result.clear();

  return 0;
}

/*!
 * Gets the name of the j-th dimension.
 *
 * We want dimension indices to vary from zero to N-1, where N is the number of
 * dimensions, to allow for iterating over dimensions.
 *
 * HDF5 dimension scales are *variables*, so it is possible that the first
 * dimension has the index of 0, but the second one has the index of, say 2
 * (index 1 belonging to a regular variable, that is).
 *
 * To work around this we build a list of dimensions and pick from there.
 */
int PISMNC4_HDF5::inq_dimname(int j, string &result) const {
  herr_t stat;
  vector<string> dim_names;
  hsize_t idx;

  stat = H5Literate_by_name(file_id, "/", H5_INDEX_NAME, H5_ITER_INC, &idx, get_dimensions,
                            &dim_names, H5P_DEFAULT);
  check(stat);

  if ((size_t)j < dim_names.size())
    result = dim_names[j];
  else
    result.clear();

  return 0;
}

/*!
 * Gets the total number of dimensions. See the comment documenting inq_dimname().
 */
int PISMNC4_HDF5::inq_ndims(int &result) const {
  // create a list of dimensions, return the length
  herr_t stat;
  vector<string> dim_names;
  hsize_t idx;

  stat = H5Literate_by_name(file_id, "/", H5_INDEX_NAME, H5_ITER_INC, &idx, get_dimensions,
                            &dim_names, H5P_DEFAULT);
  check(stat);

  result = dim_names.size();

  return 0;
}


// var
/*!
 * FIXME: I need to re-think chunking for 3D variables (it should match the
 * in-memory storage order).
 */
int PISMNC4_HDF5::def_var(string name, PISM_IO_Type xtype, vector<string> dims) const {
  herr_t stat;

  stat = H5LTfind_dataset(file_id, name.c_str()); check(stat);

  // Check if this variable already exists and return if it does.
  if (stat > 0)
    return 0;

  int max_xm, max_ym;
  MPI_Allreduce(&m_xm, &max_xm, 1, MPI_INT, MPI_MAX, com);
  MPI_Allreduce(&m_ym, &max_ym, 1, MPI_INT, MPI_MAX, com);

  vector<hsize_t> extent, max_extent, chunk;

  vector<string>::iterator j;
  for (j = dims.begin(); j != dims.end(); ++j) {
    hid_t dim_id = H5Dopen(file_id, j->c_str(), H5P_DEFAULT);
    hid_t ds_id = H5Dget_space(dim_id);

    int variable_rank = H5Sget_simple_extent_ndims(ds_id);
    assert(variable_rank == 1);

    hsize_t dim_extent, dim_maxextent;
    stat = H5Sget_simple_extent_dims(ds_id, &dim_extent, &dim_maxextent);
    assert(stat == 1);

    extent.push_back(dim_extent);
    max_extent.push_back(dim_maxextent);

    if (*j == "x")
      chunk.push_back(max_xm);
    else if (*j == "y")
      chunk.push_back(max_ym);
    else
      chunk.push_back(1);

    H5Sclose(ds_id);
    H5Dclose(dim_id);
  }

  hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE); check(plist_id);

  hid_t dataspace;
  if (extent.size() > 0) {
    // default case
    dataspace = H5Screate_simple(extent.size(), &extent[0], &max_extent[0]); check(dataspace);
    stat = H5Pset_chunk(plist_id, chunk.size(), &chunk[0]); check(stat);
  } else {
    // scalar variables (such as "mapping")
    dataspace = H5Screate(H5S_SCALAR); check(dataspace);
  }

  hid_t file_type = pism_type_to_hdf5_type(xtype);

  hid_t dset_id = H5Dcreate(file_id, name.c_str(), file_type, dataspace,
                            H5P_DEFAULT, plist_id, H5P_DEFAULT);
  check(dset_id);

  H5Sclose(dataspace);
  H5Pclose(plist_id);

  for (unsigned int k = 0; k < dims.size(); ++k) {
    hid_t dim_id = H5Dopen(file_id, dims[k].c_str(), H5P_DEFAULT);

    H5DSattach_scale(dset_id, dim_id, k);

    H5Dclose(dim_id);
  }

  H5Dclose(dset_id);

  return 0;
}


int PISMNC4_HDF5::get_vara_double(string variable_name,
                                  vector<unsigned int> start,
                                  vector<unsigned int> count,
                                  double *ip) const {
  herr_t stat;
  hid_t plist_id;

  assert(start.size() == count.size());

  hid_t dset_id = H5Dopen(file_id, variable_name.c_str(), H5P_DEFAULT); check(dset_id);

  // Enable collective I/O.
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE); // Tuning: try changing this.

  vector<hsize_t> hdf5_start(count.size()),
    hdf5_count(count.size());

  for (unsigned int k = 0; k < count.size(); ++k) {
    hdf5_start[k] = start[k];
    hdf5_count[k] = count[k];
  }

  hid_t mem_space = H5Screate_simple(hdf5_count.size(), &hdf5_count[0], NULL);
  check(mem_space);

  hid_t file_space = H5Dget_space(dset_id);
  check(file_space);

  stat = H5Sselect_hyperslab(file_space, H5S_SELECT_SET,
                             &hdf5_start[0], NULL, &hdf5_count[0], NULL);
  check(stat);

  stat = H5Dread(dset_id, H5T_NATIVE_DOUBLE, mem_space, file_space,
                 plist_id, ip);
  check(stat);

  stat = H5Pclose(plist_id); check(stat);
  stat = H5Sclose(file_space); check(stat);
  stat = H5Sclose(mem_space); check(stat);
  stat = H5Dclose(dset_id); check(stat);

  if (stat >= 0)
    return 0;
  else
    return 1;
}



int PISMNC4_HDF5::put_vara_double(string variable_name,
                                  vector<unsigned int> start,
                                  vector<unsigned int> count,
                                  double *op) const {
  herr_t stat;
  hid_t plist_id;

  assert(start.size() == count.size());

  hid_t dset_id = H5Dopen(file_id, variable_name.c_str(), H5P_DEFAULT); check(dset_id);

  // Check if dset_id points to a dimension scale (we may need to extend it).
  int is_scale = H5DSis_scale(dset_id);
  check(is_scale);

  if (is_scale > 0) {
    // dset_id points to a scale
    assert(start.size() == 1);

    hid_t space_id = H5Dget_space(dset_id);

    stat = H5Sget_simple_extent_ndims(space_id);
    assert(stat == 1);

    hsize_t dim, max_dim;
    stat = H5Sget_simple_extent_dims(space_id, &dim, &max_dim);
    assert(stat == 1);

    int increment = start[0] + count[0] - dim;

    if (increment > 0) {
      assert(max_dim == H5S_UNLIMITED);

      stat = extend_dimension(dset_id, increment);
      check(stat);
    }

    H5Sclose(space_id);
  } // end of "if(is_scale > 0)"

  // Enable collective I/O.
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE); // Tuning: try changing this.

  vector<hsize_t> hdf5_start(count.size()),
    hdf5_count(count.size());

  for (unsigned int k = 0; k < count.size(); ++k) {
    hdf5_start[k] = start[k];
    hdf5_count[k] = count[k];
  }

  hid_t mem_space = H5Screate_simple(hdf5_count.size(), &hdf5_count[0], NULL); check(mem_space);

  hid_t file_space = H5Dget_space(dset_id); check(file_space);

  stat = H5Sselect_hyperslab(file_space, H5S_SELECT_SET,
                             &hdf5_start[0], NULL, &hdf5_count[0], NULL); check(stat);

  stat = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, mem_space, file_space,
                  plist_id, op); check(stat);

  stat = H5Pclose(plist_id); check(stat);
  stat = H5Sclose(file_space); check(stat);
  stat = H5Sclose(mem_space); check(stat);
  stat = H5Dclose(dset_id); check(stat);

  if (stat >= 0)
    return 0;
  else
    return 1;
}


int PISMNC4_HDF5::get_varm_double(string /*variable_name*/,
                                  vector<unsigned int> /*start*/,
                                  vector<unsigned int> /*count*/,
                                  vector<unsigned int> /*imap*/, double */*ip*/) const {
  // Not implemented.
  return -1;
}


int PISMNC4_HDF5::put_varm_double(string /*variable_name*/,
                                  vector<unsigned int> /*start*/,
                                  vector<unsigned int> /*count*/,
                                  vector<unsigned int> /*imap*/, double */*op*/) const {
  // Not implemented.
  return -1;
}


int PISMNC4_HDF5::inq_nvars(int &result) const {
  // For now assume that the number of variables is the number of links in the root group.
  herr_t stat;
  H5G_info_t group_info;
  stat = H5Gget_info_by_name(file_id, "/", &group_info, H5P_DEFAULT); check(stat);

  result = group_info.nlinks;

  return 0;
}


int PISMNC4_HDF5::inq_vardimid(string variable_name, vector<string> &result) const {
  herr_t stat;

  hid_t var_id = H5Dopen(file_id, variable_name.c_str(), H5P_DEFAULT); check(var_id);

  stat = inq_dimensions(var_id, result); check(stat);

  stat = H5Dclose(var_id); check(stat);

  if (stat >= 0)
    return 0;
  else
    return 1;
}


int PISMNC4_HDF5::inq_varnatts(string variable_name, int &result) const {
  H5O_info_t info;
  herr_t stat = H5Oget_info_by_name(file_id, variable_name.c_str(), &info, H5P_DEFAULT); check(stat);

  result = (int)info.num_attrs;

  return 0;
}


int PISMNC4_HDF5::inq_varid(string variable_name, bool &exists) const {
  herr_t stat = H5LTfind_dataset(file_id, variable_name.c_str()); check(stat);

  if (stat > 0)
    exists = true;
  else
    exists = false;

  return 0;
}


int PISMNC4_HDF5::inq_varname(unsigned int j, string &result) const {
  herr_t stat;

  size_t len = H5Lget_name_by_idx(file_id, "/", H5_INDEX_NAME, H5_ITER_INC, j, NULL, 1, H5P_DEFAULT);
  result.resize(len + 2);

  stat = H5Lget_name_by_idx(file_id, "/", H5_INDEX_NAME, H5_ITER_INC, j, &result[0], len+1, H5P_DEFAULT); check(stat);

  return 0;
}


int PISMNC4_HDF5::inq_vartype(string variable_name, PISM_IO_Type &result) const {
  hid_t var_id = H5Dopen(file_id, variable_name.c_str(), H5P_DEFAULT); check(var_id);

  hid_t type_id = H5Dget_type(var_id); check(type_id);

  hid_t native_type_id = H5Tget_native_type(type_id, H5T_DIR_ASCEND);
  check(native_type_id);

  // We use byte, float, and double. Others are lumped with double.

  if (H5Tequal(native_type_id, H5T_NATIVE_CHAR)) {
    result = PISM_CHAR;
  } else if (H5Tequal(native_type_id, H5T_NATIVE_FLOAT)) {
    result = PISM_FLOAT;
  } else {
    result = PISM_DOUBLE;
  }

  H5Tclose(native_type_id);
  H5Tclose(type_id);
  H5Dclose(var_id);

  return 0;
}


// att
int PISMNC4_HDF5::get_att_double(string variable_name, string att_name, vector<double> &result) const {

  if (variable_name == "PISM_GLOBAL")
    variable_name = "/";

  hid_t attr_id = H5Aopen_by_name(file_id, variable_name.c_str(), att_name.c_str(),
                                  H5P_DEFAULT, H5P_DEFAULT); check(attr_id);
  hid_t space_id = H5Aget_space(attr_id); check(space_id);

  size_t len = H5Sget_simple_extent_npoints(space_id); check(len);

  result.resize(len);

  herr_t stat = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &result[0]); check(stat);

  stat = H5Sclose(space_id); check(stat);
  stat = H5Aclose(attr_id); check(stat);

  return 0;
}


//! \brief Get a string (text) attibute.
/*!
 * String attributes are weird: they are considered "scalar" datasets using a
 * datatype of which has strlen(str) characters.
 */
int PISMNC4_HDF5::get_att_text(string variable_name, string att_name, string &result) const {

  hid_t attr_id;
  herr_t stat;

  if (variable_name == "PISM_GLOBAL")
    variable_name = "/";

  stat = H5Aexists_by_name(file_id, variable_name.c_str(), att_name.c_str(),
                           H5P_DEFAULT); check(stat);
  if (stat == 0) {
    result.clear();
    return 0;
  }

  attr_id = H5Aopen_by_name(file_id, variable_name.c_str(), att_name.c_str(),
                            H5P_DEFAULT, H5P_DEFAULT); check(attr_id);

  hid_t type_id = H5Aget_type(attr_id); check(type_id);
  size_t len = H5Tget_size(type_id); check(len);

  // len includes the terminating zero
  result.resize(len);

  hid_t mem_type_id = H5Tget_native_type(type_id, H5T_DIR_ASCEND);
  stat = H5Aread(attr_id, mem_type_id, &result[0]); check(stat);

  stat = H5Tclose(mem_type_id); check(stat);
  stat = H5Tclose(type_id); check(stat);
  stat = H5Aclose(attr_id); check(stat);

  return 0;
}

int PISMNC4_HDF5::put_att_double(string variable_name, string att_name, PISM_IO_Type xtype, vector<double> &data) const {

  herr_t stat;

  if (variable_name == "PISM_GLOBAL")
    variable_name = "/";

  // Remove the attribute if it already exists
  stat = H5Aexists_by_name(file_id, variable_name.c_str(), att_name.c_str(), H5P_DEFAULT);
  check(stat);

  if (stat > 0) {
    stat = H5Adelete_by_name(file_id, variable_name.c_str(), att_name.c_str(), H5P_DEFAULT);
    check(stat);
  }

  hid_t file_type = pism_type_to_hdf5_type(xtype);

  hsize_t len = data.size();
  hid_t att_space = H5Screate_simple(1, &len, NULL);

  hid_t att_id = H5Acreate_by_name(file_id, variable_name.c_str(), att_name.c_str(),
                                   file_type, att_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  H5Sclose(att_space);

  stat = H5Awrite(att_id, H5T_NATIVE_DOUBLE, &data[0]); check(stat);

  H5Aclose(att_id);

  return 0;
}


int PISMNC4_HDF5::put_att_text(string variable_name, string att_name, string value) const {

  if (variable_name == "PISM_GLOBAL")
    variable_name = "/";

  herr_t stat = H5LTset_attribute_string(file_id, variable_name.c_str(),
                                         att_name.c_str(), value.c_str()); check(stat);

  return 0;
}


int PISMNC4_HDF5::inq_attname(string variable_name, unsigned int n, string &result) const {
  size_t len = H5Aget_name_by_idx(file_id, variable_name.c_str(),
                                  H5_INDEX_NAME, H5_ITER_INC,
                                  n, NULL, 0, H5P_DEFAULT);
  result.resize(len + 2);

  herr_t stat = H5Aget_name_by_idx(file_id, variable_name.c_str(),
                                   H5_INDEX_NAME, H5_ITER_INC,
                                   n, &result[0], len + 1, H5P_DEFAULT); check(stat);

  return 0;
}


int PISMNC4_HDF5::inq_atttype(string variable_name, string att_name, PISM_IO_Type &result) const {

  herr_t stat;
  hid_t att_id = H5Aopen_by_name(file_id,
                                 variable_name.c_str(),
                                 att_name.c_str(), H5P_DEFAULT, H5P_DEFAULT); check(att_id);

  hid_t att_type = H5Aget_type(att_id); check(att_type);

  H5T_class_t type_class = H5Tget_class(att_type);

  if (type_class == H5T_STRING) {
    result = PISM_CHAR;
  } else {
    hid_t native_type_id = H5Tget_native_type(att_type, H5T_DIR_ASCEND); check(native_type_id);

    if (H5Tequal(native_type_id, H5T_NATIVE_DOUBLE) > 0)
      result = PISM_DOUBLE;
    else if (H5Tequal(native_type_id, H5T_NATIVE_FLOAT) > 0)
      result = PISM_FLOAT;
    else if (H5Tequal(native_type_id, H5T_NATIVE_CHAR) > 0)
      result = PISM_BYTE;
    else if (H5Tequal(native_type_id, H5T_NATIVE_SHORT) > 0)
      result = PISM_SHORT;
    else
      result = PISM_INT;

    stat = H5Tclose(native_type_id); check(stat);
  }

  stat = H5Tclose(att_type); check(stat);

  stat = H5Aclose(att_id); check(stat);

  return 0;
}


// misc
int PISMNC4_HDF5::set_fill(int /*fillmode*/, int &/*old_modep*/) const {
  // A no-op. Pre-filling is disabled by default and there is no reason to
  // enable it.
  return 0;
}


void PISMNC4_HDF5::check(int return_code) const {
  if (return_code < 0) {
    H5Eprint(H5E_DEFAULT, stderr);
    H5Eclear(H5E_DEFAULT);
  }
}

string PISMNC4_HDF5::get_format() const {
  return "netcdf4";
}
