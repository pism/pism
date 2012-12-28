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

  double EPS_ABS[] = { 1.0e-12, 1.0e-9,  1.0e-7  };
  double EPS_REL[] = { 1.0e-15, 1.0e-14, 1.0e-11 };

  printf("Enter  r  (in km; e.g. 10.0):   ");
  scanret = scanf("%lf",&r);
  if (scanret != 1) {
    printf("... input error; exiting\n");
    return 1;
  }

  ierr = exactP(r*1000.0,&h,&magvb,&Wcrit,&W,EPS_ABS[0],EPS_REL[0],1);
  if (ierr) {
    printf("\n\nsimpleP ENDING because of ERROR from exactP():\n");
    error_message_testP(ierr);
    return 1;
  }

  printf("Results from Test P:\n");
  printf("    h = %.4f (m)  Po = %.5f (bar)  |vb| = %.5f (m a-1)\n"
         "    W_c = %.8f (m)  W = %.8f (m)\n",
         h,910.0*9.81*h/1.0e5,magvb*secpera,Wcrit,W);

#define COMMENTARY 0
#if COMMENTARY
  printf("(*COMMENTARY*\n Above were produced with RK Dormand-Prince (8,9) method\n"
         " and default (tight) tolerances EPS_ABS = %.1e, EPS_REL = %.1e.\n",
         EPS_ABS[0],EPS_REL[0]);
  printf("Here is a table of values using alternative methods (1=rk8pd,2=rk2,3=rkf45,4=rkck)\n"
         "and tolerances.\n\n");
  int method,i,j;
  for (method=1; method<5; method++) {
    printf("method = %d:\n",method);
    for (i=0; i<3; i++) {
      printf("    EPS_ABS = %.1e;  EPS_REL = %.1e,  %.1e,  %.1e\n",
             EPS_ABS[i],EPS_REL[0],EPS_REL[1],EPS_REL[2]);
      for (j=0; j<3; j++) {
        ierr = exactP(r*1000.0,&h,&magvb,&Wcrit,&W,EPS_ABS[i],EPS_REL[j],method);
        if (ierr) {
          printf("\n\nsimpleP ENDING because of ERROR from exactP():\n");
          error_message_testP(ierr);
          return 1;
        }
        printf("        W = %.14f (m)\n",W);
      }
    }
  }
  printf(" *END COMMENTARY*)\n");
#endif

  return 0;
}

