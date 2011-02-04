/*
   Copyright (C) 2007--2008 Ed Bueler
  
   This file is part of PISM.
  
   PISM is free software; you can redistribute it and/or modify it under the
   terms of the GNU General Public License as published by the Free Software
   Foundation; either version 2 of the License, or (at your option) any later
   version.
  
   PISM is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
   FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
   details.
  
   You should have received a copy of the GNU General Public License
   along with PISM; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*  STANDARD DIALOGUEs:

$ ./simpleK 
Enter t in years:
0
Enter z in m (0 < z < 3000 for ice and -1000 < z < 0 for bedrock):
0
Results from Test K:
     T =   260.70501 (K) =   -12.44499 (deg C)   [absolute temperature]

$ ./simpleK 
Enter t in years:
133465
Enter z in m (0 < z < 3000 for ice and -1000 < z < 0 for bedrock):
0
Results from Test K:
     T =   270.55200 (K) =    -2.59800 (deg C)   [absolute temperature]

*/

#include <stdio.h>
#include "exactTestK.h"

int main() {

  double       z, t, TT;
  int          scanret;
  const double secpera=31556926.0;  /* seconds per year; 365.2422 days */
  const double H0 = 3000.0, B0 = -1000.0; 

  /* call to compute alpha_k:  */ 
  /* print_alpha_k(30); */
 
  printf("Enter t in years:\n");
  scanret = scanf("%lf",&t);
  if (scanret != 1) {
    printf("... input error; exiting\n"); 
    return 1;
  }
  printf("Enter z in m (0 < z < %4.0f for ice and %5.0f < z < 0 for bedrock):\n",
         H0,B0);
  scanret = scanf("%lf",&z);
  if (scanret != 1) {
    printf("... input error; exiting\n"); 
    return 1;
  }

  if (t >= 0) {
    exactK(t * secpera, z, &TT, 0);    /* bedrock and ice have different material properties */
  } else {
    exactK(- t * secpera, z, &TT, 1);  /* dumb trick to help test:  use negative t values to choose
                                          bedrock material constants equal to those of ice */
  }

  printf("Results from Test K:\n");
  printf("     T = %11.5f (K) = %11.5f (deg C)   [absolute temperature]\n",TT,TT - 273.15);

  return 0;
}
