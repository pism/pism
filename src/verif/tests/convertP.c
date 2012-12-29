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


#include <stdio.h>
#include <stdlib.h>
#include "exactTestP.h"

int main ( int argc, char *argv[] ) {
  double *rr;
  int    N, k, scanret;
  FILE   *infile;
  /*  
  const double secpera=31556926.0;   seconds per year; 365.2422 days
  double EPS_ABS = 1.0e-12;
  double EPS_REL = 1.0e-15;
  */

  if ( argc != 3 ) {
    printf( "  usage:   convertP infilename outfilename\n");
    return 0;
  } else {
    infile = fopen( argv[1], "r" );
    if ( infile == 0 ) {
      printf( "convertP ERROR:  could not open file %s\n", argv[1]);
      return 1;
    }
  }

  /* read N = number of radius values */
  scanret = fscanf(infile,"%d",&N);
  if (scanret == EOF) {
    printf("can't read N; exiting\n");
    return 1;
  }

  rr = (double *) malloc((size_t)N * sizeof(double));

  for (k=0; k<N; k++) {
    scanret = fscanf(infile,"%lf",&(rr[k]));
    if ((scanret == EOF) && (k<N-1)) {
      printf("reading r fails at k=%d; exiting\n", k);
      return 1;
    }
  }

  printf("N=%d values successfully read from file %s:\n",N,argv[1]);
  for (k=0; k<N; k++) {
    printf("%f ",rr[k]);
  }
  printf("\n");

  /*ierr = exactP(r*1000.0,&h,&magvb,&Wcrit,&W,EPS_ABS[0],EPS_REL[0],1);
  if (ierr) {
    printf("\n\nsimpleP ENDING because of ERROR from exactP():\n");
    error_message_testP(ierr);
    return 1;
  }*/

  free(rr);

  return 0;
}

