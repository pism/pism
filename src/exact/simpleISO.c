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

/*  STANDARD DIALOGUE:
user@home:~/pism$ obj/simpleISO
Enter  t  and  r  separated by space (or newline)
    (in yrs and km, resp.; e.g. 15000 500):
15000 500

Results:
    Test A              Test B              Test C              Test D
  H = 2362.8962 (m)   H = 1878.2858 (m)   H = 2386.8463 (m)   H = 2310.3543 (m)
  M =  0.30000 (m/a)  M =  0.00000 (m/a)  M =  0.79562 (m/a)  M =  0.24969 (m/a)

exact test E requires cartesian location ...
Enter  x  and  y  separated by space (in km; e.g. 400 200; ctrl-D to exit):
400 200

Results for Test E:
    H =  2524.368242 (m)
    M =     0.640174 (m/a)
    mu = 2.47248e-11 (m/(Pa s))
    u_b =  44.721769 (m/a)
    v_b =  22.360884 (m/a)

*/

#include <stdio.h>
#include "exactTestsABCDE.h"

int main() {

  const double SperA=31556926.0;  /* seconds per year; 365.2422 days */

  double year, r, 
         HA, MA, HB, MB, HC, MC, HD, MD,
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

  /* evaluate tests */
  exactA(r*1000.0, &HA, &MA);
  exactB(year*SperA, r*1000.0, &HB, &MB);
  exactC(year*SperA, r*1000.0, &HC, &MC);
  exactD(year*SperA, r*1000.0, &HD, &MD);

  printf("\nResults:\n");
  printf("    Test A              Test B            ");
  printf("  Test C              Test D\n");
  printf("  H = %9.4f (m)   H = %9.4f (m)   H = %9.4f (m)   ",
         HA,HB,HC);
  printf("H = %9.4f (m)\n",HD);
  printf("  M = %8.5f (m/a)  M = %8.5f (m/a)  M = %8.5f (m/a)  ",
         MA*SperA,MB*SperA,MC*SperA);
  printf("M = %8.5f (m/a)\n\n",MD*SperA);

  printf("exact test E requires cartesian location ...\n");
  printf("Enter  x  and  y  separated by space ");
  printf("(in km; e.g. 400 200; ctrl-D to exit):\n");
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

  printf("\nResults for Test E:\n");
  printf("    H = %12.6f (m)\n",HE);
  printf("    M = %12.6f (m/a)\n",ME*SperA);
  printf("    mu = %11.5e (m/(Pa s))\n",muE);
  printf("    u_b = %10.6f (m/a)\n",ubE*SperA);
  printf("    v_b = %10.6f (m/a)\n\n",vbE*SperA);

  return 0;
}
