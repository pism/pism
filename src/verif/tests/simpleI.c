/*
   Copyright (C) 2004-2006 Jed Brown and Ed Bueler
  
   This file is part of Pism.
  
   Pism is free software; you can redistribute it and/or modify it under the
   terms of the GNU General Public License as published by the Free Software
   Foundation; either version 3 of the License, or (at your option) any later
   version.
  
   Pism is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
   FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
   details.
  
   You should have received a copy of the GNU General Public License
   along with Pism; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*  STANDARD DIALOGUE:

user@home:~/pism$ obj/simpleI
Enter  m, x, y  separated by space (or newline);
    (units: pure, pure, km, km; e.g. 10 100 40):
10 100 40

Results from Test I:
  bed =  -100.000 (m)     tau_c =     17.854 (kPa)
  u   =   252.126 (m/a)   v     =      0.000 (m/a)

*/

#include <stdio.h>
#include "exactTestsIJ.h"

int main() {

  double m, x, y, bed, tauc, u, v;
  int    scanret;
  const double secpera=31556926.0;  /* seconds per year; 365.2422 days */
  
  printf("Enter  m, x, y  separated by space (or newline);\n");
  printf("    (units: pure, km, km; e.g. 10 100 40):\n");
  scanret = scanf("%lf",&m);
  if (scanret != 1) {
    printf("... input error; exiting\n"); 
    return 1;
  }
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

  exactI(m,x*1000.0,y*1000.0, &bed, &tauc, &u, &v);

  printf("\nResults from Test I:\n");
  printf("  bed = %9.3f (m)     tau_c = %10.3f (kPa)\n",bed,tauc/1000.0);
  printf("  u   = %9.3f (m/a)   v     = %10.3f (m/a)\n",u*secpera,v*secpera);

  return 0;
}
