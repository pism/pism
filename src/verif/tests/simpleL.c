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

  double EPS_ABS[] = { 1.0e-12, 1.0e-9, 1.0e-7 };
  double EPS_REL[] = { 0.0, 1.0e-14, 1.0e-11 }; 
  exactL(r*1000.0,&H,&b,&a,EPS_ABS[0],EPS_REL[0],1);

  printf("Results from Test L:\n");
  printf("     H = %11.5f (m)   b = %10.5f (m)   a = %8.5f (m/a)\n",
         H,b,a*secpera);

#define COMMENTARY 0
#if COMMENTARY
  printf("(*COMMENTARY*\n Above were produced with RK Cash-Karp method and default (tight) tolerances\n");
  printf(" EPS_ABS = %e, EPS_REL = %e.\n",EPS_ABS[0],EPS_REL[0]);
  printf(" Here is a table of values of H using alternative methods and tolerances:\n");
  int method,i,j; 
  for (method=1; method<5; method++) {
    printf("   method = %d  (1=rkck,2=rk2,3=rk4,4=rk8pd):\n",method);  
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
        if ((method == 2) && (j == 0)) {
          printf("     EPS_ABS = %e, EPS_REL = %e:    <method hangs for this EPS_REL>\n",
                 EPS_ABS[i],EPS_REL[j]);
        } else {
          exactL(r*1000.0,&H,&b,&a,EPS_ABS[i],EPS_REL[j],method);
          printf("     EPS_ABS = %e, EPS_REL = %e:    H = %17.11f\n",EPS_ABS[i],EPS_REL[j],H);
        }
      }
    }
  }
  printf(" *END COMMENTARY*)\n");
#endif


#define GRAPHABLE 0
#if GRAPHABLE
/*
    produce profile and bed graphs from this data.
    e.g. saving output in file "rHba.sce", editing out above stuff, 
    running scilab, do "exec('rHba.sce')" in scilab   
*/
#define N 751
  int k;
  double rr[N];
  double HH[N],bb[N],aa[N];
  for (k= 0; k<N; k++) {
    rr[k] = 750000.0 - 1000.0 * k;
  }
  exactL_list(rr,N,HH,bb,aa);
  printf("\n  rHba = [\n");
  for (k = 0; k<N; k++) {
    printf("%10.3f  %11.5f  %10.5f  %8.5f\n",
           rr[k]/1000.0,HH[k],bb[k],aa[k]*secpera);
  }
  printf("];\n\n plot(rHba(:,1),rHba(:,2)+rHba(:,3),rHba(:,1),rHba(:,3))\n");
#endif

  return 0;
}
