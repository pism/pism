/* Copyright (C) 2013 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef _ICEMODELVEC3CUSTOM_H_
#define _ICEMODELVEC3CUSTOM_H_

#include "iceModelVec.hh"

/**
 * This class allows storing and saving 3D (or 2D with dof>1) data
 * using *one* NetCDF variable. (Using a "custom" vertical grid.) It
 * is used to store 3D velocity components on the scaled ("sigma")
 * vertical grid by the BlatterStressBalance class and for latitude
 * and longitude bounds by IceModel.
 *
 * \note DOF>1 2D data that should be stored in and read from several
 * variables can be stored using IceModelVec2.
 */
class IceModelVec3Custom : public IceModelVec3D {
public:
  IceModelVec3Custom();
  virtual ~IceModelVec3Custom();

  virtual PetscErrorCode create(IceGrid &mygrid,
                                string short_name,
                                string z_name,
                                vector<double> my_zlevels,
                                map<string, string> z_attrs);
};

#endif /* _ICEMODELVEC3CUSTOM_H_ */
