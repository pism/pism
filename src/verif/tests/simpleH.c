/*
   Copyright (C) 2004-2006 Jed Brown and Ed Bueler
  
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

/*  STANDARD DIALOGUE:

$ ./simpleH
Enter  f,  t,  and  r  separated by spaces (or newline)
    (in pure, yrs, and km, resp.;  note  f = rho_ice/rho_bedrock;
     e.g. input 0.28 40000 500):
0.28 40000 500
Results for Test H:
    H =  2356.058308 (m)     M =     0.294507 (m/a)    b =  -659.696326 (m)

*/

#include <stdio.h>
#include "exactTestH.h"

int main() {

  const double SperA=31556926.0;  /* seconds per year; 365.2422 days */

  double f, year, r, H, M, b;
  int    scanret;
  
  printf("Enter  f,  t,  and  r  separated by spaces (or newline)\n");
  printf("    (in pure, yrs, and km, resp.;  note  f = rho_ice/rho_bedrock;\n");
  printf("     e.g. input 0.28 40000 500):\n");
  scanret = scanf("%lf",&f);
  if (scanret != 1) {
    printf("... input error; exiting\n"); 
    return 1;
  }
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

  exactH(f, year*SperA, r*1000.0, &H, &M);
  b = -f * H;
  
  printf("Results for Test H:\n");
  printf("    H = %12.6f (m)     M = %12.6f (m/a)    b = %12.6f (m)\n",H,M*SperA,b);

  return 0;
}
