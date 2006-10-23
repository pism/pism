// Copyright (C) 2004-2006 Jed Brown and Ed Bueler
//
// This file is part of Pism.
//
// Pism is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// Pism is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with Pism; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

static char help[] =
  "Ice sheet driver for a special purpose: given a stored ice model and \n"
  "a source of mass-balance velocities, compute deformational velocity from SIA\n" 
  "and attribute remaining velocity (i.e. mass-balance minus deformational) \n"
  "to basal drag according to MacAyeal-Morland equations.\n"
  "(Jed Kallen-Brown and Ed Bueler, authors, and Craig Lingle, concepts.).\n";

/* a possible invocation of this executable would be:
 *   $ ./get_drag -if ant_141_101_10yrsmooth.pb -balvel ../ice-input/ant_2006/init.nc \
 *        -gk -d n -verbose -o ant_141_DRAG
 * would produce  ant_141_DRAG.m  containing beta, betax, betay, and balvel
 * where betax,betay have units of  Pa s m^-1 */

#include <cstring>
#include <petscbag.h>
#include "grid.hh"
#include "materials.hh"
#include "iceModel.hh"

class IceDragModel : public IceModel {
public:
  Vec vmagbalvel;
  Vec dragxy;
  Vec vNu[2];
  IceDragModel(IceGrid &g, IceType &i);
#if (WITH_NETCDF)
  PetscErrorCode readBalvelFromFile(const char *fname);
#endif
  PetscErrorCode getEffectiveViscosity();
  PetscErrorCode computeDragFromBalanceVelocity();
  PetscErrorCode writeMatlabDragFile(const char *basename);
  PetscErrorCode destroyDragVecs();
};


IceDragModel::IceDragModel(IceGrid &g, IceType &i)
  : IceModel(g,i) {
  // do nothing (except pass g,i to constructor for IceModel)
}


//FIXME: need to remove C++ style calls to NetCDF library; see iMIOnetcdf.cc
//for now the next thing produces a compile error

#if (WITH_NETCDF)
#include <netcdf.h>

PetscErrorCode IceDragModel::readBalvelFromFile(const char *fname) {
  PetscErrorCode  ierr;
  char filename[PETSC_MAX_PATH_LEN];

  vmagbalvel = vWork2d[2];  // already allocated
  
  if (fname == NULL) { // This can be called from the driver for a default
    SETERRQ(1,"invalid file name");
  } else {
    strcpy(filename, fname);
  }

  NcVar *v_balvel;
  NcFile *f;
  if (grid.rank == 0) {
    f = new NcFile(filename);
    if (! f->is_valid()) {
      ierr = PetscPrintf(grid.com, "Could not open %s for reading.\n", filename); CHKERRQ(ierr);
      return 1;
    }
  
    v_balvel = f->get_var("balvel");   
  }

  // see iMIOnetcdf.cc for original technique
  Vec vzero;
  VecScatter ctx;
  ierr = VecScatterCreateToZero(g2, &ctx, &vzero); CHKERRQ(ierr);  
  ierr = getIndZero(grid.da2, g2, vzero, ctx); CHKERRQ(ierr);

  ierr = ncVarToDAVec(v_balvel, grid.da2, vmagbalvel, g2, vzero);     CHKERRQ(ierr);
  
  ierr = VecDestroy(vzero); CHKERRQ(ierr);
  ierr = VecScatterDestroy(ctx); CHKERRQ(ierr);

  ierr = VecScale(vmagbalvel,1.0/secpera); CHKERRQ(ierr);  // convert to m/s
  
  if (grid.rank == 0) {
    f->close();
    delete f;
  }

  return 0;
}
#endif


PetscErrorCode IceDragModel::getEffectiveViscosity() {
  PetscErrorCode ierr;
  PetscScalar epsilon = 0.0;  // actually get effective viscosity (not bdd below)
  PetscScalar **mask;
  
  // need effective viscosity *everywhere* so mark all grounded points as DRAGGING
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (modMask(mask[i][j]) != MASK_FLOATING)   mask[i][j] = MASK_DRAGGING;
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  // communicate, because computing eff visc requires neighbor vals of mask
  ierr = DALocalToLocalBegin(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);

  // put effective viscosity in working arrays
  vNu[0] = vWork2d[0];
  vNu[1] = vWork2d[1];
  ierr = computeEffectiveViscosity(vNu, epsilon); CHKERRQ(ierr);  

  // view if  "-d n"
  if (nuView[0] != PETSC_NULL && nuView[1] != PETSC_NULL) {
    ierr = DALocalToGlobal(grid.da2, vNu[0], INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecView(g2, nuView[0]); CHKERRQ(ierr);
    ierr = DALocalToGlobal(grid.da2, vNu[1], INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecView(g2, nuView[1]); CHKERRQ(ierr);
    ierr = PetscPrintf(grid.com,"\n final viscosity; pausing for 10 seconds ... \n"); CHKERRQ(ierr);
    ierr = PetscSleep(10); CHKERRQ(ierr);
  }
  return 0;
}

  
// compute drag from velocities and mass-balance-velocity
PetscErrorCode IceDragModel::computeDragFromBalanceVelocity() {
  PetscErrorCode ierr;
  Mat A;
  Vec x, result, rhs;
  PetscScalar **u, **v, **mask, **balvel;

  // allocate Vecs
  const PetscInt M = 2 * grid.p->Mx * grid.p->My;
  ierr = MatCreateMPIAIJ(grid.com, PETSC_DECIDE, PETSC_DECIDE, M, M,
                         13, PETSC_NULL, 13, PETSC_NULL,
                         &A); CHKERRQ(ierr);
  ierr = VecCreateMPI(grid.com, PETSC_DECIDE, M, &x); CHKERRQ(ierr);
  ierr = VecDuplicate(x, &result); CHKERRQ(ierr);
  ierr = VecDuplicate(x, &rhs); CHKERRQ(ierr);
  ierr = VecDuplicate(x, &dragxy); CHKERRQ(ierr);

  // build discrete version of MacAyeal-Morland equations (A x = rhs) at all grounded points
  ierr = assembleMacayealMatrix(vNu, A, rhs); CHKERRQ(ierr);

  // Note rhs contains driving terms  \rho g H \grad h  but  A  contains basal
  // drag term [betax*u betay*v]'.  It must be removed.
  // Also set x = [u v]'  (i.e. interleaved).
  ierr = VecSet(x, 0.0); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vubar, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvbar, &v); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vmagbalvel, &balvel); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscInt    J = 2*j;
      const PetscInt    rowU = i*2*grid.p->My + J;
      const PetscInt    rowV = i*2*grid.p->My + J+1;
      // remove old drag term
      if (intMask(mask[i][j]) == MASK_DRAGGING) {
        ierr = MatSetValue(A, rowU, rowU, - basalDragx(u, v, i, j), ADD_VALUES); CHKERRQ(ierr);
        ierr = MatSetValue(A, rowV, rowV, - basalDragy(u, v, i, j), ADD_VALUES); CHKERRQ(ierr);
      }
      // remove deformational from mass-balance velocities; put in x
      const PetscScalar c = sqrt(PetscSqr(u[i][j]) + PetscSqr(v[i][j]));
      PetscScalar       residualFactor = (balvel[i][j] / c) - 1.0;
      // if (residualFactor < 0.0)   residualFactor = 0.0;  // avoid negative drag coeff here?  probably not
      ierr = VecSetValue(x, rowU, residualFactor * u[i][j], INSERT_VALUES); CHKERRQ(ierr);
      ierr = VecSetValue(x, rowV, residualFactor * v[i][j], INSERT_VALUES); CHKERRQ(ierr);
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vubar, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvbar, &v); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vmagbalvel, &balvel); CHKERRQ(ierr);
  // communicate!
  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(x); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(x); CHKERRQ(ierr);

  // With new A and x, compute  result = Ax - rhs.
  ierr = MatMult(A,x,result); CHKERRQ(ierr);
  ierr = VecAXPY(result,-1.0,rhs); CHKERRQ(ierr);  // note result = 0 if MASK_SHEET
  
  // As result is [betax*u betay*v]', divide by velocities, but only where grounded!
  // where floating, report no drag coeff
  ierr = VecPointwiseDivide(dragxy,result,x); CHKERRQ(ierr);
  ierr = VecScale(dragxy,-1.0); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscInt    J = 2*j;
      const PetscInt    rowU = i*2*grid.p->My + J;
      const PetscInt    rowV = i*2*grid.p->My + J+1;
      if (modMask(mask[i][j]) == MASK_FLOATING) {
        ierr = VecSetValue(dragxy, rowU, 0.0, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValue(dragxy, rowV, 0.0, INSERT_VALUES); CHKERRQ(ierr);
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(dragxy); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(dragxy); CHKERRQ(ierr);
  
  ierr = MatDestroy(A); CHKERRQ(ierr);
  ierr = VecDestroy(x); CHKERRQ(ierr);
  ierr = VecDestroy(result); CHKERRQ(ierr);
  ierr = VecDestroy(rhs); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceDragModel::writeMatlabDragFile(const char *basename) {
  PetscErrorCode  ierr;
  PetscViewer     viewer;
  char b[PETSC_MAX_PATH_LEN];
  char mf[PETSC_MAX_PATH_LEN];  // Matlab format

  strncpy(b, basename, PETSC_MAX_PATH_LEN-4);
  ierr = PetscOptionsGetString(PETSC_NULL, "-o", b, PETSC_MAX_PATH_LEN, PETSC_NULL); CHKERRQ(ierr);
  strcpy(mf, b);
  strcat(mf, ".m");
  ierr = PetscPrintf(grid.com, "writing variables beta, betax, betay, balvel to file %s ...", mf); CHKERRQ(ierr);

  ierr = PetscViewerASCIIOpen(grid.com, mf, &viewer); CHKERRQ(ierr);
  ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB); CHKERRQ(ierr);

  ierr = PetscViewerASCIIPrintf(viewer,"\n\ndisp('iceDragModel output:')\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"echo on\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\nclear year x y xx yy beta betax betay dragxy balvel\n");
  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"echo off\n");  CHKERRQ(ierr);

  ierr = PetscViewerASCIIPrintf(viewer,"year=%10.6f;\n",grid.p->year);  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                                "x=(-%12.3f:%12.3f:%12.3f)/1000.0;\n",
                                grid.p->Lx,grid.p->dx,grid.p->Lx);  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                                "y=(-%12.3f:%12.3f:%12.3f)/1000.0;\n",
                                grid.p->Ly,grid.p->dy,grid.p->Ly);  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"[xx,yy]=meshgrid(x,y);\n\n");  CHKERRQ(ierr);

  ierr=PetscObjectSetName((PetscObject) g2,"balvel"); CHKERRQ(ierr);
  ierr = LVecView(grid.da2, vmagbalvel, g2, viewer); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\nbalvel = reshape(balvel,%d,%d);\n\n",
                                grid.p->Mx,grid.p->My);  CHKERRQ(ierr);

  ierr = PetscObjectSetName((PetscObject) dragxy,"dragxy"); CHKERRQ(ierr);
  ierr = VecView(dragxy, viewer); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"betax = dragxy(1:2:%d);\n",(2*grid.p->Mx*grid.p->My)-1);
        CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\nbetax = reshape(betax,%d,%d);\n\n",
                                grid.p->Mx,grid.p->My);  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"betay = dragxy(2:2:%d);\n",2*grid.p->Mx*grid.p->My);
        CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\nbetay = reshape(betay,%d,%d);\n\n",
                                grid.p->Mx,grid.p->My);  CHKERRQ(ierr);

  /* following block writes Matlab to display */
  ierr = PetscViewerASCIIPrintf(viewer,"echo on\n");  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer, "figure, contourf(x,y,betax',50), colorbar\n"); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"title('drag coeff in x direction (forcex = - betax * u)')\n\n");
        CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer, "figure, contourf(x,y,betay',50), colorbar\n"); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"title('drag coeff in y direction (forcey = - betay * v)')\n\n");
        CHKERRQ(ierr);

  ierr = PetscViewerASCIIPrintf(viewer, "beta = (max(betax,0.0) + max(betay,0.0))/2;\n"); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer, "figure, contour(x,y,balvel',[0 1 5 10 20 50 100 200 500 1000 2000 4000],'k')\n"); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer, "hold on, surf(x,y,log(beta')/log(10)), view(2), colorbar, hold off\n"); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"title('log_{10}(beta) (color) on contours of balance vels (black)')\n\n");
        CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"echo off\n");  CHKERRQ(ierr);

  ierr = PetscViewerPopFormat(viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);
  return 0;
}


PetscErrorCode IceDragModel::destroyDragVecs() {
  PetscErrorCode  ierr;

  ierr = VecDestroy(dragxy); CHKERRQ(ierr);
  // note vmagbalvel, vNu[0], vNu[1] are just IceModel work Vecs
  return 0;
}


int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm    com;
  PetscMPIInt rank, size;

  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);

  com = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  {
    IceGrid    g(com, rank, size);
    IceType*   ice;
    PetscInt   flowlawNumber = 4;  // use hybrid (Goldsby-Kohlstedt) by default
    
    ierr = getFlowLawFromUser(com, ice, flowlawNumber); CHKERRQ(ierr);

    IceDragModel m(g, *ice);

    ierr = m.setFromOptions(); CHKERRQ(ierr);
    ierr = m.initFromOptions(); CHKERRQ(ierr);
    ierr = m.setSoundingFromOptions(); CHKERRQ(ierr);

    ierr = PetscPrintf(com, "computing velocities (with MacAyeal; including eff viscosity iteration) ..."); CHKERRQ(ierr);
    m.setUseMacayealVelocity(PETSC_TRUE);
    ierr = m.velocity(true); CHKERRQ(ierr);
    ierr = PetscPrintf(com, " done \n"); CHKERRQ(ierr);

    ierr = PetscPrintf(com, "computing resulting effective viscosity at all grounded points ..."); CHKERRQ(ierr);
    ierr = m.getEffectiveViscosity(); CHKERRQ(ierr);
    ierr = PetscPrintf(com, " done \n"); CHKERRQ(ierr);

    ierr = PetscPrintf(com, "computing deformational velocities (w/o MacAyeal) ..."); CHKERRQ(ierr);
    m.setUseMacayealVelocity(PETSC_FALSE);
    m.setMuSliding(0.0);  // for deformational, just assume frozen bed
    m.setEnhancementFactor(0.8);  //  reduce amount of deformation to ascribe more
                                  //  of mass-balance velocities to sliding
    ierr = m.velocity(true); CHKERRQ(ierr);
    ierr = PetscPrintf(com, " done \n"); CHKERRQ(ierr);

    ierr = PetscPrintf(com, "seeking file with mass-balance velocities ..."); CHKERRQ(ierr);
    PetscTruth balvelFileSet;
    char balvelFile[PETSC_MAX_PATH_LEN];
    ierr = PetscOptionsGetString(PETSC_NULL, "-balvel", balvelFile,
                                 PETSC_MAX_PATH_LEN, &balvelFileSet); CHKERRQ(ierr);
    if (balvelFileSet == PETSC_TRUE) {
#if (WITH_NETCDF)
      ierr = m.readBalvelFromFile(balvelFile); CHKERRQ(ierr);
#else
      SETERRQ(1,"get_drag requires NetCDF file for balance velocities; NetCDF must be enabled");
#endif
    } else {
      SETERRQ(1,"get_drag requires source for balance velocities; use -balvel foo.nc");
    }
    ierr = PetscPrintf(com, " done \n"); CHKERRQ(ierr);

    ierr = PetscPrintf(com, "computing drag by subtracting deformational from mass balance and using MacAyeal equations ..."); CHKERRQ(ierr);
    ierr = m.computeDragFromBalanceVelocity(); CHKERRQ(ierr);
    ierr = PetscPrintf(com, " done \n"); CHKERRQ(ierr);

    ierr = PetscPrintf(com, "saving Matlab file ..."); CHKERRQ(ierr);
    ierr = m.writeMatlabDragFile("dragfile"); CHKERRQ(ierr);
    ierr = PetscPrintf(com, " done \n"); CHKERRQ(ierr);

    ierr = PetscPrintf(com, "destroying IceDragModel Vecs ..."); CHKERRQ(ierr);
    ierr = m.destroyDragVecs(); CHKERRQ(ierr);
    ierr = PetscPrintf(com, " done\n"); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
