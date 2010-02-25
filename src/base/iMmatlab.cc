// Copyright (C) 2007-2010 Ed Bueler and Constantine Khroulev
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

#include <cstring>
#include <cmath>
#include <petscda.h>
#include "iceModel.hh"

//! Write out the linear system of equations for the last nonlinear iteration (for SSA).
/*!
Writes out the matrix, the right-hand side, and the solution vector in a \em Matlab-readable
format, a <tt>.m</tt> file.
*/
PetscErrorCode IceModel::writeSSAsystemMatlab(IceModelVec2 vNu[2]) {
  PetscErrorCode ierr;
  PetscViewer    viewer;
  char           file_name[PETSC_MAX_PATH_LEN], yearappend[PETSC_MAX_PATH_LEN];

  strcpy(file_name,ssaMatlabFilePrefix);
  snprintf(yearappend, PETSC_MAX_PATH_LEN, "_%.0f.m", grid.year);
  strcat(file_name,yearappend);
  ierr = verbPrintf(2, grid.com, 
             "writing Matlab-readable file for SSA system A xsoln = rhs to file `%s' ...\n",
             file_name); CHKERRQ(ierr);
  ierr = PetscViewerCreate(grid.com, &viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer, PETSC_VIEWER_ASCII);CHKERRQ(ierr);
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, file_name);CHKERRQ(ierr);

  // save linear system; gives system A xsoln = rhs at last (nonlinear) iteration of SSA
  ierr = PetscViewerASCIIPrintf(viewer,
             "%% PISM SSA linear system report.\n"); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
             "%% Writes matrix A (sparse) and vectors xsoln, rhs for the linear system\n"); 
             CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
             "%% (A xsoln = rhs) which was solved at the last step of the nonlinear iteration.\n");
             CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
             "%% Also writes the year, the coordinates x,y and their gridded versions xx,yy.\n"); 
             CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
             "%% Also thickness (H), surface elevation (h), i-offset (staggered grid) \n");
             CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
             "%% vertically-integrated viscosity (nu_0 = nu * H), and j-offset version\n");
             CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
             "%% of same thing (nu_1 = nu * H).\n\n");
             CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
             "%% The thickness can be plotted by\n"); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
             "%%   >> imagesc(x,y,flipud(H')), axis square, colorbar \n");
             CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
             "%% for example. \n");
             CHKERRQ(ierr);

  ierr = PetscViewerASCIIPrintf(viewer,"\n\necho off\n");  CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) SSAStiffnessMatrix,"A"); CHKERRQ(ierr);
  ierr = MatView(SSAStiffnessMatrix, viewer);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) SSARHS,"rhs"); CHKERRQ(ierr);
  ierr = VecView(SSARHS, viewer);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) SSAX,"xsoln"); CHKERRQ(ierr);
  ierr = VecView(SSAX, viewer);CHKERRQ(ierr);

  // save coordinates (for viewing, primarily)
  ierr = PetscViewerASCIIPrintf(viewer,"year=%10.6f;\n",grid.year);  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
            "x=%12.3f + (0:%d)*%12.3f;\n",-grid.Lx,grid.Mx-1,grid.dx);
            CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
            "y=%12.3f + (0:%d)*%12.3f;\n",-grid.Ly,grid.My-1,grid.dy);
            CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"[xx,yy]=meshgrid(x,y);\n\n");  CHKERRQ(ierr);

  // also save thickness and effective viscosity
  ierr = vH.view_matlab(viewer); CHKERRQ(ierr);
  ierr = vh.view_matlab(viewer); CHKERRQ(ierr);

  ierr = vNu[0].set_name("nu_0"); CHKERRQ(ierr);
  ierr = vNu[0].set_attr("long_name", "effective viscosity times thickness (i offset)"); CHKERRQ(ierr);
  ierr = vNu[0].view_matlab(viewer); CHKERRQ(ierr);

  ierr = vNu[1].set_name("nu_1"); CHKERRQ(ierr);
  ierr = vNu[1].set_attr("long_name", "effective viscosity times thickness (j offset)"); CHKERRQ(ierr);
  ierr = vNu[1].view_matlab(viewer); CHKERRQ(ierr);

  ierr = PetscViewerASCIIPrintf(viewer,"echo on\n");  CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
  return 0;
}

