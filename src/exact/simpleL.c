/*
   Copyright (C) 2007 Ed Bueler
  
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

$ ./simpleL
Enter  r  (in km; e.g. 0.0):   0.0
Results from Test L:
  H =  3782.88760 (m)   b = -500.00000 (m)   a =  0.30000 (m/a)

*/

#include <stdio.h>
#include "exactTestL.h"

int main() {

  double       r, H, b, a;
  int          scanret;
  const double secpera=31556926.0;  /* seconds per year; 365.2422 days */
  
  printf("Enter  r  (in km; e.g. 0.0):   ");
  scanret = scanf("%lf",&r);
  if (scanret != 1) {
    printf("... input error; exiting\n"); 
    return 1;
  }

  exactL(r*1000.0,&H,&b,&a);

  printf("Results from Test L:\n");
  printf("  H = %11.5f (m)   b = %10.5f (m)   a = %8.5f (m/a)\n",
         H,b,a*secpera);

/*  produce graphs from this data:
#define N 751
  int i;
  double rr[N];
  double HH[N],bb[N],aa[N];
  for (i = 0; i<N; i++) {
    rr[i] = 750000.0 - 1000.0 * i;
  }
  exactL_list(rr,N,HH,bb,aa);
  for (i = 0; i<N; i++) {
    printf("%10.3f  %11.5f  %10.5f  %8.5f\n",
           rr[i]/1000.0,HH[i],bb[i],aa[i]*secpera);
  }
*/
  return 0;
}
