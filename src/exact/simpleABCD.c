/*
   Copyright (C) 2004-2007 Jed Brown and Ed Bueler
  
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

$ ./simpleABCD
Enter  t  and  r  separated by space (or newline)
    (in yrs and km, resp.; e.g. 15000 500):
15000 500
Results:
    Test A              Test B              Test C              Test D
  H = 2362.8962 (m)   H = 1878.2858 (m)   H = 2386.8463 (m)   H = 2310.3543 (m)
  M =  0.30000 (m/a)  M =  0.00000 (m/a)  M =  0.79562 (m/a)  M =  0.24969 (m/a)

*/

#include <stdio.h>
#include "exactTestsABCDE.h"

int main() {

  const double SperA=31556926.0;  /* seconds per year; 365.2422 days */

  double year, r, 
         HA, MA, HB, MB, HC, MC, HD, MD;
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

  /* evaluate tests */
  exactA(r*1000.0, &HA, &MA);
  exactB(year*SperA, r*1000.0, &HB, &MB);
  exactC(year*SperA, r*1000.0, &HC, &MC);
  exactD(year*SperA, r*1000.0, &HD, &MD);

  printf("Results:\n");
  printf("    Test A              Test B            ");
  printf("  Test C              Test D\n");
  printf("  H = %9.4f (m)   H = %9.4f (m)   H = %9.4f (m)   ",
         HA,HB,HC);
  printf("H = %9.4f (m)\n",HD);
  printf("  M = %8.5f (m/a)  M = %8.5f (m/a)  M = %8.5f (m/a)  ",
         MA*SperA,MB*SperA,MC*SperA);
  printf("M = %8.5f (m/a)\n\n",MD*SperA);

  return 0;
}
