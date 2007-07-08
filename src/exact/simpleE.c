/*
   Copyright (C) 2004-2006 Jed Brown and Ed Bueler
  
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

/*  STANDARD DIALOGUE:

$ ./simpleE
Enter  t  and  r  separated by space (or newline)
    (in yrs and km, resp.; e.g. 15000 500):
15000 500
exact test E also requires cartesian location ...
Enter  x  and  y  separated by space (in km; e.g. 400 200):
400 200
Results for Test E:
    H =  2524.368242 (m)     M =     0.640174 (m/a)    mu = 2.47248e-11 (m/(Pa s))
    u_b =  44.721769 (m/a)   v_b =  22.360884 (m/a)

*/

#include <stdio.h>
#include "exactTestsABCDE.h"

int main() {

  const double SperA=31556926.0;  /* seconds per year; 365.2422 days */

  double year, r, 
         x, y, HE, ME, muE, ubE, vbE;
  int    scanret;
  
  printf("Enter  t  and  r  separated by space (or newline)\n");
  printf("    (in yrs and km, resp.; e.g. 15000 500):\n");
  scanret = scanf("%lf",&year);
  if (scanret != 1) {
    printf("... input error; exiting\n"); 
    return 1;
  }
  scanret = scanf("%lf",&r);
  if (scanret != 1) {
    printf("... input error; exiting\n"); 
    return 1;
  }

  printf("exact test E also requires cartesian location ...\n");
  printf("Enter  x  and  y  separated by space ");
  printf("(in km; e.g. 400 200):\n");
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

  exactE(x*1000.0, y*1000.0, &HE, &ME, &muE, &ubE, &vbE);

  printf("Results for Test E:\n");
  printf("    H = %12.6f (m)     M = %12.6f (m/a)    mu = %11.5e (m/(Pa s))\n",HE,ME*SperA,muE);
  printf("    u_b = %10.6f (m/a)   v_b = %10.6f (m/a)\n",ubE*SperA,vbE*SperA);

  return 0;
}
