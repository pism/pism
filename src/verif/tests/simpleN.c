/*
   Copyright (C) 2010 Ed Bueler
  
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

/*  STANDARD DIALOGUE:  a point near the ELA:

./simpleN 
Enter  x  (in km; 0.0 <= x <= 500.0):   300.0
Results from Test N:
     h = surface elevation    = 1920.00000 (m)
     h_x = surface slope      = -7.20000e-03
     u = ice velocity         =  100.00000 (m a-1)
     M = surface mass balance =   -8.00000 (cm a-1)
     A = ice softness         = 1.29664e-25 (Pa-3 s-1)
*/


#include <stdio.h>
#include "exactTestN.h"

int main() {

  double       H0, L0, xc, x, 
               u, h, hx, M, A;
  const double secpera=31556926.0;  /* seconds per year; 365.2422 days */
  int          scanret, retvalN;
  
  geometry_exactN(&H0, &L0, &xc);

  printf("Enter  x  (in km; 0.0 <= x <= %5.1f):   ", L0 / 1000.0);
  scanret = scanf("%lf",&x);
  if (scanret != 1) {
    printf("... input error; exiting\n"); 
    return 1;
  }

  x = x * 1000.0;

  retvalN = exactN(x,  &h, &hx, &u, &M, &A);

  if (retvalN) {
    printf("SIMPLEN ERROR:  x  out of allowed domain  0.0 <= x <= %5.1f km\n"
           "   ... ending ...\n", L0 / 1000.0);
    return 1;
  }
  if (x > xc) {
    printf("WARNING:  x = %5.1f km is past the calving front at xc = %5.1f km\n",
           x / 1000.0, xc / 1000.0);
  }

  printf("Results from Test N:\n");
  printf(
    "     h = surface elevation    = %10.5f (m)\n"
    "     h_x = surface slope      = %10.5e\n"
    "     u = ice velocity         = %10.5f (m a-1)\n"
    "     M = surface mass balance = %10.5f (cm a-1)\n"
    "     A = ice softness         = %10.5e (Pa-3 s-1)\n",
    h, hx, u * secpera, M * secpera * 100.0, A);
  return 0;
}

