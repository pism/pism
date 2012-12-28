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

  double EPS_ABS[] = { 1.0e-12, 1.0e-9, 1.0e-7 };
  double EPS_REL[] = { 1.0e-15, 1.0e-14, 1.0e-11 };
  
  printf("Enter  r  (in km; e.g. 0.0):   ");
  scanret = scanf("%lf",&r);
  if (scanret != 1) {
    printf("... input error; exiting\n"); 
    return 1;
  }

  exactL(r*1000.0,&H,&b,&a,EPS_ABS[0],EPS_REL[0],4);

  printf("Results from Test L:\n");
  printf("     H = %11.5f (m)   b = %10.5f (m)   a = %8.5f (m/a)\n",
         H,b,a*secpera);

#define COMMENTARY 0
#if COMMENTARY
  printf("\nAbove were produced with RK Cash-Karp (4,5) method and default tolerances\n");
  printf(" EPS_ABS = %e, EPS_REL = %e.\n",EPS_ABS[0],EPS_REL[0]);
  printf(" Here is a table of values of H using alternative methods and tolerances:\n");
  int method,i,j; 
  for (method=1; method<5; method++) {
    printf("   method = %d  (1=rk8pd,2=rk2,3=rkf45,4=rkck):\n",method);
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
        exactL(r*1000.0,&H,&b,&a,EPS_ABS[i],EPS_REL[j],method);
        printf("     EPS_ABS = %e, EPS_REL = %e:    H = %17.11f\n",EPS_ABS[i],EPS_REL[j],H);
      }
    }
  }
#endif

#define PLOTTABLE 0
#if PLOTTABLE
#define N 751
  printf("\nHere is a plottable table of values r, H, b, a:\n");
  /* run and plot in Matlab:
  $ ./simpleL >& foo.txt   # and then enter r
  $ [EDIT foo.txt TO REMOVE FIRST LINES]
  $ matlab
  >> load foo.txt
  >> plot(foo(:,1),foo(:,3),foo(:,1),foo(:,2)+foo(:,3))
  >> grid on,  xlabel('r  (km)'),  ylabel('elevation  (m)')
  */
  int k;
  const double L = 750.0e3;
  double rr[N], HH[N], bb[N], aa[N];
  for (k= 0; k<N; k++) {
    rr[k] = L * (N - k - 1) / (N-1);
  }
  exactL_list(rr,N,HH,bb,aa);
  for (k = 0; k<N; k++) {
    printf("%10.3f  %11.5f  %10.5f  %8.5f\n",
           rr[k]/1000.0,HH[k],bb[k],aa[k]*secpera);
  }
#endif

  return 0;
}
