/* Copyright (C) 2015 Florian Ziemen and the PISM Authors
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

cmake_minimum_required(VERSION 2.8)

PROJECT(WriteVTU)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})


add_executable(WriteVTU MACOSX_BUNDLE WriteVTU)

if(VTK_LIBRARIES)
  target_link_libraries(WriteVTU ${VTK_LIBRARIES})
else()
  target_link_libraries(WriteVTU vtkHybrid vtkWidgets)
endif()

target_link_libraries (WriteVTU netcdf_c++4 boost_program_options)
