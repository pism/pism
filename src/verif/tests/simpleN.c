/*
   Copyright (C) 2010, 2013 Ed Bueler
  
   This file is part of PISM.
  
   PISM is free software; you can redistribute it and/or modify it under the
   terms of the GNU General Public License as published by the Free Software
   Foundation; either version 3 of the License, or (at your option) any later
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

$ ./simpleN 
Enter  x  (in km; 0.0 <= x <= 500.0):   300.0
Results from Test N:
     H    = ice thickness        =   1920.00000 (m)
     h_x  = surface slope        = -7.20000e-03
     u    = ice velocity         =    300.00000 (m/year)
     M    = surface mass balance =    -24.00000 (cm/year)
     B    = ice hardness         =  1.36989e+08 (Pa s^(1/3))
     beta = ice hardness         =  1.29813e+10 (Pa s m-1)
*/


#include <stdio.h>
#include "exactTestN.h"

int main() {

  double       H0, L0, xc, a, Hela, k, Hc, Tc, x, 
               H, hx, u, M, B, beta;
  const double secpera=31556926.0;  /* seconds per year; 365.2422 days */
  int          scanret, retvalN;
  
  params_exactN(&H0, &L0, &xc, &a, &Hela, &k, &Hc, &Tc);

  printf("Enter  x  (in km; 0.0 <= x <= %5.1f):   ", L0 / 1000.0);
  scanret = scanf("%lf",&x);
  if (scanret != 1) {
    printf("... input error; exiting\n"); 
    return 1;
  }

  x = x * 1000.0;

  retvalN = exactN(x, &H, &hx, &u, &M, &B, &beta);

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
  printf("     H    = ice thickness        = %12.5f (m)\n"
         "     h_x  = surface slope        = %12.5e\n"
         "     u    = ice velocity         = %12.5f (m/year)\n"
         "     M    = surface mass balance = %12.5f (cm/year)\n"
         "     B    = ice hardness         = %12.5e (Pa s^(1/3))\n"
         "     beta = ice hardness         = %12.5e (Pa s m-1)\n",
         H, hx, u * secpera, M * secpera * 100.0, B, beta);
  return 0;
}

