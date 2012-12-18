/*
   Copyright (C) 2012 Ed Bueler

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

$ ./simpleP
FIXME

*/

#include <stdio.h>
#include "exactTestP.h"

int main() {

  double       r, h, magvb, Wcrit, W;
  int          scanret, ierr;
  const double secpera=31556926.0;  /* seconds per year; 365.2422 days */

  double EPS_ABS[] = { 1.0e-12, 1.0e-9, 1.0e-7 };
  double EPS_REL[] = { 0.0, 1.0e-14, 1.0e-11 };

  printf("Enter  r  (in km; e.g. 10.0):   ");
  scanret = scanf("%lf",&r);
  if (scanret != 1) {
    printf("... input error; exiting\n");
    return 1;
  }

  ierr = exactP(r*1000.0,&h,&magvb,&Wcrit,&W,EPS_ABS[0],EPS_REL[0],1);
  if (ierr) {
    printf("simpleP ERROR: exactP() returns %d ... ENDING ...\n", ierr);
    return 1;
  }

  printf("Results from Test P:\n");
  printf("    h = %.4f (m)  Po = %.5f (bar)  |vb| = %.5f (m a-1)  W_c = %.8f (m)  W = %.8f (m)\n",
         h,910.0*9.81*h/1.0e5,magvb*secpera,Wcrit,W);

#define COMMENTARY 1
#if COMMENTARY
  printf("*COMMENTARY*\n Above were produced with RK Cash-Karp method and default tolerances.\n");
  printf(" Here is a table of values of h, |vb|, W using alternative methods and tolerances:\n");
  int method,j;
  for (method=1; method<5; method++) {
    printf("   method = %d  (1=rkck,2=rk2,3=rk4,4=rk8pd):\n",method);
    for (j=0; j<3; j++) {
      exactP(r*1000.0,&h,&magvb,&Wcrit,&W,EPS_ABS[j],EPS_REL[j],method);
      printf("     EPS = [%.1e %.1e]:    h = %11.5f (m)   |vb| = %11.5f (m a-1)   W = %11.8f (m)\n",
             EPS_ABS[j],EPS_REL[j],h,magvb*secpera,W);
    }
  }
  printf("*END COMMENTARY*\n");
#endif

  return 0;
}
