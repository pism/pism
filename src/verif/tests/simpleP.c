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
    error_message_testP(ierr);
    printf("simpleP ENDING because of error from exactP() ...\n");
    return 1;
  }

  printf("Results from Test P:\n");
  printf("    h = %.4f (m)  Po = %.5f (bar)  |vb| = %.5f (m a-1)  W_c = %.8f (m)  W = %.8f (m)\n",
         h,910.0*9.81*h/1.0e5,magvb*secpera,Wcrit,W);

  return 0;
}

