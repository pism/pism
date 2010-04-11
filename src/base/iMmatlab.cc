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

//! Write out the linear system of equations for the SSA linear solve at the last nonlinear iteration in a time step.
/*!
Writes out the matrix, the right-hand side, and the solution vector in a
Matlab/Octave-readable format, a <tt>.m</tt> file.  Also some related
quantities to allow decent plotting.

This method is called by the option \c -ssa_matlab \c foo.  For a \c
-y \c 0 run there will be one system produced in the file \c foo_y0.m.
To use the result in Matlab/Octave type\code
  >> help foo_y0
\endcode

For example, to see the matrix from a low resolution EISMINT-Ross
calculation, the corresponding thickness, and the vector field solution do:\code
$ pross -boot_from ross.nc -Mx 21 -My 21 -Mz 3 -Lz 1e3 -ssaBC ross.nc -riggs riggs.nc -ssa_matlab rosslow
$ octave
...
>> rosslow_y0
>> figure(1), spy(A,'.')
>> figure(2), imagesc(x,y,flipud(thk)), axis square, colorbar
>> figure(3)
>> u = reshape(uv(1:2:end-1),21,21);  v = reshape(uv(2:2:end),21,21);
>> scale = sqrt(max(max(u.^2+v.^2)));
>> quiver(u/scale,v/scale)
\endcode
*/
PetscErrorCode IceModel::writeSSAsystemMatlab(IceModelVec2S vNu[2]) {
  PetscErrorCode ierr;
  PetscViewer    viewer;
  char           file_name[PETSC_MAX_PATH_LEN], yearappend[PETSC_MAX_PATH_LEN];

  strcpy(file_name,ssaMatlabFilePrefix);
  snprintf(yearappend, PETSC_MAX_PATH_LEN, "_y%.0f.m", grid.year);
  strcat(file_name,yearappend);
  ierr = verbPrintf(2, grid.com, 
             "writing Matlab-readable file for SSA system A xsoln = rhs to file `%s' ...\n",
             file_name); CHKERRQ(ierr);
  ierr = PetscViewerCreate(grid.com, &viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer, PETSC_VIEWER_ASCII);CHKERRQ(ierr);
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, file_name);CHKERRQ(ierr);

  // get the command which started the run  [ FIXME: code duplication from
  // IceModel::stampHistoryCommand() ]
  PetscInt argc;
  char **argv;
  ierr = PetscGetArgs(&argc, &argv); CHKERRQ(ierr);
  string cmdstr;  // a string with space-separated command-line arguments:
  for (int j = 0; j < argc; j++)
    cmdstr += string(" ") + argv[j];
  
  // save linear system; gives system A xsoln = rhs at last (nonlinear) iteration of SSA
  ierr = PetscViewerASCIIPrintf(viewer,
    "%% A PISM linear system report for the SSA stress balance from this run:\n"
    "%%   '%s'\n"
    "%% Writes matrix A (sparse), and vectors uv and rhs, for the linear\n"
    "%% system which was solved at the last step of the nonlinear iteration:\n"
    "%%    A * uv = rhs.\n"
    "%% Also writes the year, the coordinates x,y, their gridded versions\n"
    "%% xx,yy, and the thickness (thk) and surface elevation (usurf).\n"
    "%% Also writes i-offsetvalues of vertically-integrated viscosity\n"
    "%% (nuH_0 = nu * H), and j-offset version of same thing (nuH_1 = nu * H);\n"
    "%% these are on the staggered grid.\n",
    cmdstr.c_str());  CHKERRQ(ierr);

  ierr = PetscViewerASCIIPrintf(viewer,"\n\necho off\n");  CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) SSAStiffnessMatrix,"A"); CHKERRQ(ierr);
  ierr = MatView(SSAStiffnessMatrix, viewer);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"clear zzz\n\n");  CHKERRQ(ierr);

  ierr = PetscObjectSetName((PetscObject) SSARHS,"rhs"); CHKERRQ(ierr);
  ierr = VecView(SSARHS, viewer);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) SSAX,"uv"); CHKERRQ(ierr);
  ierr = VecView(SSAX, viewer);CHKERRQ(ierr);

  // save coordinates (for viewing, primarily)
  ierr = PetscViewerASCIIPrintf(viewer,"\nyear=%10.6f;\n",grid.year);  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
            "x=%12.3f + (0:%d)*%12.3f;\n"
            "y=%12.3f + (0:%d)*%12.3f;\n",
            -grid.Lx,grid.Mx-1,grid.dx,-grid.Ly,grid.My-1,grid.dy); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"[xx,yy]=meshgrid(x,y);\n");  CHKERRQ(ierr);

  // also save thickness and effective viscosity
  ierr = vH.view_matlab(viewer); CHKERRQ(ierr);
  ierr = vh.view_matlab(viewer); CHKERRQ(ierr);

  ierr = vNu[0].set_name("nuH_0"); CHKERRQ(ierr);
  ierr = vNu[0].set_attr("long_name", 
    "effective viscosity times thickness (i offset) at current time step"); CHKERRQ(ierr);
  ierr = vNu[0].view_matlab(viewer); CHKERRQ(ierr);

  ierr = vNu[1].set_name("nuH_1"); CHKERRQ(ierr);
  ierr = vNu[1].set_attr("long_name",
    "effective viscosity times thickness (j offset) at current time step"); CHKERRQ(ierr);
  ierr = vNu[1].view_matlab(viewer); CHKERRQ(ierr);

  ierr = PetscViewerASCIIPrintf(viewer,"echo on\n");  CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
  return 0;
}

