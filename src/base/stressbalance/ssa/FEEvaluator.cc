// Copyright (C) 2012 David Maxwell
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

#include "FEEvaluator.hh"

void FEEvaluator::setElement( PetscInt i, PetscInt j) {
  m_i=i; m_j=j;
  IceGrid &grid = *m_grid;
  m_dofmap.reset(i,j,grid);
  if( (grid.xs <= i) && (i < grid.xs + grid.xm)) {
    if( (grid.ys <= j) && (j < grid.ys + grid.ym)) {
      m_elementOwner = true;
    }
  }
}

void FEEvaluator::testFunctionValues(PetscInt q, PetscInt k, PetscReal *v, PetscReal *dx, PetscReal *dy)
{
  const FEFunctionGerm *values = m_quadrature.testFunctionValues(q,k);
  *v=values->val;
  *dx=values->dx;
  *dy=values->dy;
}

PetscReal FEEvaluator::JxW(PetscInt q) {
  return m_quadrature.quadWeights[q];
}


void FEEvaluator::coordsAt( PetscInt q, PetscReal *x, PetscReal *y )
{
  *x = (m_quadrature.quadPoints[q][0]+1.)/2.*m_grid->dx;
  *y = (m_quadrature.quadPoints[q][1]+1.)/2.*m_grid->dy;
  *x += m_grid->x[m_i];
  *y += m_grid->y[m_j];
}

FEEvaluator2S::~FEEvaluator2S() {
  release();
}

PetscErrorCode FEEvaluator2S::init( IceModelVec2S &vec) {
  release();
  m_vec = &vec;
  PetscErrorCode ierr = m_vec->get_array(m_array ); CHKERRQ(ierr);

  
  m_grid = m_vec->get_grid();
  m_quadrature.init(*m_grid);

  return 0;
}

void FEEvaluator2S::release() {
  if( m_vec != NULL) {
    m_vec->end_access();
    m_vec = NULL;
    m_array = NULL;    
  }
}

void FEEvaluator2S::setElement(PetscInt i,PetscInt j) {
  FEEvaluator::setElement(i,j);
  if( m_elementOwner) {
    m_quadrature.computeTrialFunctionValues( i,j,m_dofmap,m_array,m_values,m_dfdx,m_dfdy);    
  }
}

void FEEvaluator2S::valueAt(int q, PetscReal *v, PetscReal* dfdx, PetscReal *dfdy) {
  PetscReal lv=0;
  PetscReal ldfdx=0; 
  PetscReal ldfdy=0;
  if(m_elementOwner && (q>=0) && (q<4) ) {
    lv = m_values[q];
    ldfdx = m_dfdx[q];
    ldfdy = m_dfdy[q];
  }
  PISMGlobalSum(&lv,v,m_grid->com);
  PISMGlobalSum(&ldfdx,dfdx,m_grid->com);
  PISMGlobalSum(&ldfdy,dfdy,m_grid->com);
}

PetscReal FEEvaluator2S::nodeValue(PetscInt k) {
  PetscInt i,j;
  m_dofmap.localToGlobal(k, &i, &j);
  return m_array[i][j];
}


FEEvaluator2V::~FEEvaluator2V() {
  release();
}

PetscErrorCode FEEvaluator2V::init( IceModelVec2V &vec) {
  release();
  m_vec = &vec;
  PetscErrorCode ierr = m_vec->get_array(m_array ); CHKERRQ(ierr);
  
  m_grid = m_vec->get_grid();
  m_quadrature.init(*m_grid);
  
  return 0;
}

void FEEvaluator2V::release() {
  if( m_vec != NULL) {
    m_vec->end_access();
    m_vec = NULL;
    m_array = NULL;    
  }
}

void FEEvaluator2V::setElement(PetscInt i,PetscInt j) {
  FEEvaluator::setElement(i,j);
  if( m_elementOwner) {
    m_quadrature.computeTrialFunctionValues( i,j,m_dofmap,m_array,m_values,m_Du);    
  }
}

void FEEvaluator2V::valueAt(int q, PetscReal *u, PetscReal *v, PetscReal* Du0, PetscReal *Du1, PetscReal *Du2) {
  PetscReal lu=0;
  PetscReal lv=0;
  PetscReal lDu0=0; 
  PetscReal lDu1=0;
  PetscReal lDu2=0;

  if(m_elementOwner && (q>=0) && (q<4) ) {
    lu = m_values[q].u;
    lv = m_values[q].v;
    lDu0 = m_Du[q][0];
    lDu1 = m_Du[q][1];
    lDu2 = m_Du[q][2];
  }
  PISMGlobalSum(&lu,u,m_grid->com);
  PISMGlobalSum(&lv,v,m_grid->com);

  PISMGlobalSum(&lDu0,Du0,m_grid->com);
  PISMGlobalSum(&lDu1,Du1,m_grid->com);
  PISMGlobalSum(&lDu2,Du2,m_grid->com);
}

PISMVector2 FEEvaluator2V::nodeValue(PetscInt k) {
  PetscInt i,j;
  m_dofmap.localToGlobal(k, &i, &j);
  return m_array[i][j];
}
