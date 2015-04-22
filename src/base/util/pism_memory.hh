/* Copyright (C) 2015 PISM Authors
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

//! @file pism_memory.hh
//!
//! We're trying to support both older compilers that provide C++ shared_ptr and weak_ptr in
//! the std::tr1 namespace and newer, C++11-standard-compliant ones providing these smart pointers
//! in std. This header should be removed once we decide to require a C++11-supporting compiler.

#ifndef _PISM_MEMORY_H_
#define _PISM_MEMORY_H_

#ifdef PISM_USE_TR1
#include <tr1/memory>
#define PISM_SHARED_PTR_NSPACE std::tr1
#else
#include <memory>
#define PISM_SHARED_PTR_NSPACE std
#endif

#define PISM_SHARED_PTR(TYPE) PISM_SHARED_PTR_NSPACE::shared_ptr<TYPE>
#define PISM_WEAK_PTR(TYPE) PISM_SHARED_PTR_NSPACE::weak_ptr<TYPE>

#endif /* _PISM_MEMORY_H_ */
