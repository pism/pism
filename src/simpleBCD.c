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

/*  STANDARD CASE:
Enter  t and r  separated by newline (in yrs and km, resp.; e.g. 15000 500):
15000
500

Results:
      Test B                   Test C                   Test D
    H =  1878.285815 (m)     H =  2386.846335 (m)     H =  2310.354348 (m)
    M =     0.000000 (m/a)   M =     0.795615 (m/a)   M =     0.249693 (m/a)
*/

#include <stdio.h>
#include "exactTestsBCD.h"

int main() {

  const double SperA=31556926.0;  // seconds per year; 365.2422 days
  double year, r, HB, MB, HC, MC, HD, MD;
  
  printf("Enter  t and r  separated by space (or newline) ");
  printf("(in yrs and km, resp.; e.g. 15000 500):\n");
  scanf("%lf",&year);
  scanf("%lf",&r);

  /* evaluate tests B, C, and D */
  exactB(year*SperA,r*1000.0,&HB,&MB);
  exactC(year*SperA,r*1000.0,&HC,&MC);
  exactD(year*SperA,r*1000.0,&HD,&MD);

  printf("\nResults:\n");
  printf("      Test B                   Test C                   Test D\n");
  printf("    H = %12.6f (m)     H = %12.6f (m)     H = %12.6f (m)\n",HB,HC,HD);
  printf("    M = %12.6f (m/a)   M = %12.6f (m/a)   ",MB*SperA,MC*SperA);
  printf("M = %12.6f (m/a)\n",MD*SperA);

  return 0;
}
