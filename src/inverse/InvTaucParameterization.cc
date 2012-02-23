// Copyright (C) 2011 David Maxwell
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

#include "InvTaucParameterization.hh"

InvTaucParamIdent g_InvTaucParamIdent;
InvTaucParamSquare g_InvTaucParamSquare;
InvTaucParamExp g_InvTaucParamExp;

PetscErrorCode InvTaucParamIdent::toTauc( PetscReal p, PetscReal *value, 
                                          PetscReal *derivative)
{
  if(value != NULL) *value = p;
  if(derivative != NULL) *derivative = 1.;
  return 0;
}

PetscErrorCode InvTaucParamIdent::fromTauc( PetscReal tauc, PetscReal *OUTPUT)
{
  *OUTPUT = tauc;
  return 0;
}


PetscErrorCode InvTaucParamSquare::toTauc( PetscReal p, PetscReal *value, 
                                           PetscReal *derivative)
{
  if(value != NULL) *value = p*p;
  if(derivative != NULL) *derivative = 2*p;
  return 0;
}

PetscErrorCode InvTaucParamSquare::fromTauc( PetscReal tauc, PetscReal *OUTPUT)
{
  *OUTPUT = sqrt(tauc);
  // if(tauc>=0) {
  //   *OUTPUT = sqrt(tauc); 
  // } else {
  //   *OUTPUT = NaN;
  // }
  return 0;
}

PetscReal tauc_eps = 1.;
PetscErrorCode InvTaucParamExp::toTauc( PetscReal p, PetscReal *value, 
                                           PetscReal *derivative)
{
  if(value != NULL) *value = exp(p);
  if(derivative != NULL) *derivative = exp(p);    
  return 0;
}

PetscErrorCode InvTaucParamExp::fromTauc( PetscReal tauc, PetscReal *OUTPUT)
{
  if(tauc < tauc_eps)
  {
    tauc= tauc_eps;
  }
  *OUTPUT=log(tauc);
  return 0;
}

PetscErrorCode InvTaucParamLinear::toTauc( PetscReal p, PetscReal *value, 
                                           PetscReal *derivative)
{
  *value = m_scale*p;
  *derivative = m_scale;
  return 0;
}

PetscErrorCode InvTaucParamLinear::fromTauc( PetscReal tauc, PetscReal *OUTPUT)
{
  *OUTPUT = tauc/m_scale;
  return 0;
}

