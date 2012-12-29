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
  double *rr, *hh, *magvbvb, *WWcrit, *WW;
  int    N, k, scanret;
  FILE   *infile, *outfile;

  const double EPS_ABS = 1.0e-12,
               EPS_REL = 1.0e-15;

  if ( argc != 3 ) {
    printf( "  usage:   convertP infilename outfilename\n");
    return 0;
  } else {
    infile = fopen( argv[1], "r" );
    if ( infile == 0 ) {
      printf( "convertP ERROR:  could not open input file %s\n", argv[1]);
      return 1;
    }
  }

  /* read N = (number of radius values) and allocate */
  scanret = fscanf(infile,"%d",&N);
  if (scanret == EOF) {
    printf("ERROR can't read N; exiting\n");  return 1;
  }
  rr = (double *) malloc((size_t)N * sizeof(double));

  /* read the values */
  for (k=0; k<N; k++) {
    scanret = fscanf(infile,"%lf",&(rr[k]));
    if ((scanret == EOF) && (k<N-1)) {
      printf("ERROR reading r fails at k=%d; exiting\n", k);  return 1;
    }
  }
  fclose(infile);
  printf("N=%d values successfully read from file %s:\n",N,argv[1]);
  /* for (k=0; k<N; k++) { printf("%f\n",rr[k]); } */

  /* check they are decreasing and in  0 <= r < TESTP_L */
  for (k = 0; k<N; k++) {
    if (rr[k] < 0.0) {
      printf("ERROR rr[%d] = %.e is negative; exiting\n",k,rr[k]);  return 1;
    }
    if (rr[k] >= TESTP_L) {
      printf("ERROR rr[%d] = %.e exceeds TESTP_L; exiting\n",k,rr[k]);  return 1;
    }
    if ((k>0) && (rr[k] >= rr[k-1])) {
      printf("ERROR rr[] not decreasing at k=%d; exiting\n",k);  return 1;
    }
  }
  printf("  rr[0] = %.3f > rr[1] > ... > rr[%d] = %.3f\n",rr[0],N-1,rr[N-1]);

  /* now compute h, magvb, Wcrit, W */
  hh = (double *) malloc((size_t)N * sizeof(double));
  magvbvb = (double *) malloc((size_t)N * sizeof(double));
  WWcrit = (double *) malloc((size_t)N * sizeof(double));
  WW = (double *) malloc((size_t)N * sizeof(double));
  k = exactP_list(rr,N,hh,magvbvb,WWcrit,WW,EPS_ABS,EPS_REL,1);  /* use Dormand-Prince (8,9) */
  if (k) {
    printf("ERROR exactP_list() returned an error:\n");
    error_message_testP(k);
    return 1;
  }

  /* write out results */
  outfile = fopen( argv[2], "a" );
  if ( outfile == 0 ) {
    printf( "ERROR could not open output file %s\n", argv[2]);  return 1;
  }
  fprintf(outfile,"%d\n",N);
  for (k=0; k<N; k++) {
    fprintf(outfile,"%17.10f %17.12f %17.10e %17.12f %17.12f\n",
            rr[k],hh[k],magvbvb[k],WWcrit[k],WW[k]);
/*    if (scanret) {
      printf( "ERROR could not write k=%d line of results to %s\n", k, argv[2]);  return 1;
    }*/
  }
  fclose(outfile);
  printf("results successfully written to file %s\n", argv[2]);

  free(rr);  free(hh);  free(magvbvb);  free(WWcrit);  free(WW);

  return 0;
}

