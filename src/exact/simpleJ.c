/*
   Copyright (C) 2004-2006 Jed Brown and Ed Bueler
  
   This file is part of Pism.
  
   Pism is free software; you can redistribute it and/or modify it under the
   terms of the GNU General Public License as published by the Free Software
   Foundation; either version 2 of the License, or (at your option) any later
   version.
  
   Pism is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
   FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
   details.
  
   You should have received a copy of the GNU General Public License
   along with Pism; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*  [NEEDS TO BE CONFIRMED AS standard] DIALOGUEs:

user@home:~/pism$ simpleJ
Enter  x, y  separated by space (or newline);
    (units: km, km; e.g. 100 100):
100 100

Results from Test J:
  H   =   500.000 (m)     nu    =     30.000 (MPa a)
  u   =   109.160 (m/a)   v     =     28.921 (m/a)

user@home:~/pism$ simpleJ
Enter  x, y  separated by space (or newline);
    (units: km, km; e.g. 100 100):
0 0

Results from Test J:
  H   =   770.000 (m)     nu    =     19.481 (MPa a)
  u   =     0.000 (m/a)   v     =      0.000 (m/a)

*/

#include <stdio.h>
#include "exactTestIJ.h"

int main() {

  double x, y, u, v, H, nu;
  int    scanret;
  const double secpera=31556926.0;  /* seconds per year; 365.2422 days */
  
  printf("Enter  x, y  separated by space (or newline);\n");
  printf("    (units: km, km; e.g. 100 100):\n");
  scanret = scanf("%lf",&x);
  if (scanret != 1) {
    printf("... input error; exiting\n"); 
    return 1;
  }
  scanret = scanf("%lf",&y);
  if (scanret != 1) {
    printf("... input error; exiting\n"); 
    return 1;
  }

  exactJ(x*1000.0,y*1000.0, &H, &nu, &u, &v);
           
  printf("\nResults from Test J:\n");
  printf("  H   = %9.3f (m)     nu    = %10.3f (MPa a)\n",H,(nu*1.0e-6)/secpera);
  printf("  u   = %9.3f (m/a)   v     = %10.3f (m/a)\n",u*secpera,v*secpera);

  return 0;
}
