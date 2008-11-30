/*
   Copyright (C) 2008 Ed Bueler
  
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

$ ./simpleM
Enter  r  (in km; e.g. 400.0):   400.0
Results from Test M:
  alpha =  .... (m/a)
  
*/

/* FIXME: add this to dialog?:
(Parameters:
  H0 = 500.0 (m), Rg = 300.0 (km), ug = 100.0 (m/a), barB = 1.900e8 (Pa s^1/3) ) 
*/

#include <stdio.h>
#include "exactTestM.h"

int main() {

  double       r, alpha;
  int          scanret;
  const double secpera=31556926.0;  /* seconds per year; 365.2422 days */
  
  printf("Enter  r  (in km; e.g. 400.0):   ");
  scanret = scanf("%lf",&r);
  if (scanret != 1) {
    printf("... input error; exiting\n"); 
    return 1;
  }

  double EPS_ABS[] = { 1.0e-12, 1.0e-9, 1.0e-7 };
  double EPS_REL[] = { 0.0, 1.0e-14, 1.0e-11 }; 
  int eMresult;
  eMresult = exactM(r*1000.0,&alpha,EPS_ABS[0],EPS_REL[0],1);

  printf("Results from Test M (returned %d):\n", eMresult);
  printf("     alpha = %8.5f (m/a)\n", alpha*secpera);
  return 0;
}
