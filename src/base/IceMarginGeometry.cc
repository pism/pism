// Copyright (C) 2004-2011 Jed Brown, Ed Bueler, Constantine Khroulev
//                         and David Maxwell
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

#include "IceMarginGeometry.hh"

IceMarginGeometry::IceMarginGeometry( PetscReal seaLevel, 
                                      const IceFlowLaw &ice, 
                                      const NCConfigVariable &config)
{
  PetscReal ice_rho = ice.rho,
    ocean_rho = config.get("sea_water_density");
  m_seaLevel = seaLevel;
  m_alpha = 1.0 - ice_rho/ocean_rho;
}
