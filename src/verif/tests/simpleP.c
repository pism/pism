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
Enter  r  (in km; 0 <= r <= TESTP_L = 22.5):   20.0
Results from Test P:
    h = 180.0000 (m)  Po = 16.0687800 (bar)  |vb| = 46.26644 (m a-1)
    W_c = 0.58184968 (m)  W = 0.67507258 (m)  P = 2.0086731 (bar)

*/

#include <stdio.h>
#include "exactTestP.h"

int main() {

  double       r, h, magvb, Wcrit, W, P;
  int          scanret, ierr;
  const double secpera=31556926.0;  /* seconds per year; 365.2422 days */

  double EPS_ABS[] = { 1.0e-12, 1.0e-9,  1.0e-7  };
  double EPS_REL[] = { 1.0e-15, 1.0e-14, 1.0e-11 };

  printf("Enter  r  (in km; 0 <= r <= TESTP_L = %.1f):   ",TESTP_L/1000.0);
  scanret = scanf("%lf",&r);
  if (scanret != 1) {
    printf("... input error; exiting\n");
    return 1;
  }

  ierr = exactP(r*1000.0,&h,&magvb,&Wcrit,&W,&P,EPS_ABS[0],EPS_REL[0],1);
  if (ierr) {
    printf("\n\nsimpleP ENDING because of ERROR from exactP():\n");
    error_message_testP(ierr);
    return 1;
  }

  printf("Results from Test P:\n");
  printf("    h = %.4f (m)  Po = %.7f (bar)  |vb| = %.5f (m a-1)\n"
         "    W_c = %.8f (m)  W = %.8f (m)  P = %.7f (bar)\n",
         h,910.0*9.81*h/1.0e5,magvb*secpera,Wcrit,W,P/1.0e5);

#define COMMENTARY 0
#if COMMENTARY
  char*  methodnames[4] = { "rk8pd", "rk2", "rkf45", "rkck" };
  printf("\nAbove were produced with RK Dormand-Prince (8,9) method\n"
         "and default (tight) tolerances EPS_ABS = %.1e, EPS_REL = %.1e.\n",
         EPS_ABS[0],EPS_REL[0]);
  printf("Here is a table of values using alternative methods and tolerances.\n\n");
  int method,i,j;
  for (method=1; method<5; method++) {
    printf("method = %d = %s:\n",method,methodnames[method-1]);
    for (i=0; i<3; i++) {
      printf("    EPS_ABS = %.1e\n",EPS_ABS[i]);
      for (j=0; j<3; j++) {
        ierr = exactP(r*1000.0,&h,&magvb,&Wcrit,&W,&P,EPS_ABS[i],EPS_REL[j],method);
        if (ierr) {
          printf("\n\nsimpleP ENDING because of ERROR from exactP():\n");
          error_message_testP(ierr);
          return 1;
        }
        printf("        EPS_REL = %.1e:   W = %.14f\n",EPS_REL[j],W);
      }
    }
  }
#endif

#define PLOTTABLE 0
#if PLOTTABLE
  printf("\nHere is a plottable table of values r, Wcrit, W:\n");
  /* run and plot in Matlab:
  $ ./simpleP >& foo.txt   # and then enter r
  $ [EDIT foo.txt TO REMOVE FIRST LINES]
  $ matlab
  >> load foo.txt
  >> figure(1)      % profile
  >> plot(foo(:,1),910.0*9.81*foo(:,2)/1.0e5,foo(:,1),foo(:,5)), grid on, xlabel('r  (km)')
  >> legend('P_o(r)','P(r)'),  ylabel('pressure  (bar)')
  >> figure(2)      % water amount
  >> plot(foo(:,1),foo(:,3),'-',foo(:,1),foo(:,4),'--'),  grid on,  xlabel('r  (km)')
  >> legend('W_{crit}(r)','W(r)'),  ylabel('water thickness  (m)')
  */
  int k;
  const int N = 300;
  double rr[N], hh[N], magvbvb[N], WWcrit[N], WW[N], PP[N];
  for (k=0; k<N; k++)
    rr[k] = (double)(N - k - 1) * TESTP_L / ((double)N);
  ierr = exactP_list(rr,N,hh,magvbvb,WWcrit,WW,PP,EPS_ABS[0],EPS_REL[0],1);
  if (ierr) {
    printf("\n\nsimpleP ENDING because of ERROR from exactP():\n");
    error_message_testP(ierr);
    return 1;
  }
  printf("      r (km)           h (m)     Wcrit (m)         W (m)       P (bar)\n\n");
  for (k=0; k<N; k++)
    printf("%12.6f %15.9f %13.9f %13.9f %13.9f\n",rr[k]/1000.0,hh[k],WWcrit[k],WW[k],PP[k]/1.0e5);
#endif

  return 0;
}

