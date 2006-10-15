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
Enter  t and r  separated by newline (in yrs and km, resp.; e.g. 500 500):
500
500
Enter  z  values sep by newline (in m); '-1' to end; e.g. 0 100 500 1500 -1:
0
100
500
1500
-1

Results:
           Test F                         Test G
(functions of r (resp. t and r) only):
      H    =  1925.295290 (m)        H    =  2101.899734 (m)
      M    =    -0.010510 (m/a)      M    =     0.040738 (m/a)
(z=     0.000):
      T    =   265.122620 (K)        T    =   267.835036 (K)
      U    =     0.000000 (m/a)      U    =     0.000000 (m/a)
      w    =     0.000000 (m/a)      w    =     0.000000 (m/a)
      Sig  =     0.264346 (*)        Sig  =     1.215392 (*)
      Sigc =    -0.373726 (*)        Sigc =    -1.323664 (*)
(z=   100.000):
      T    =   263.137595 (K)        T    =   265.849860 (K)
      U    =     0.661716 (m/a)      U    =     2.244496 (m/a)
      w    =     0.000005 (m/a)      w    =    -0.000758 (m/a)
      Sig  =     0.173915 (*)        Sig  =     0.817817 (*)
      Sigc =    -0.306255 (*)        Sigc =    -1.022931 (*)
(z=   500.000):
      T    =   255.486095 (K)        T    =   258.194962 (K)
      U    =     1.785938 (m/a)      U    =     6.217140 (m/a)
      w    =     0.000291 (m/a)      w    =    -0.011984 (m/a)
      Sig  =     0.028439 (*)        Sig  =     0.149934 (*)
      Sigc =    -0.199905 (*)        Sigc =    -0.340039 (*)
(z=  1500.000):
      T    =   238.172200 (K)        T    =   240.856843 (K)
      U    =     2.036372 (m/a)      U    =     7.227603 (m/a)
      w    =     0.002288 (m/a)      w    =    -0.050018 (m/a)
      Sig  =     0.000029 (*)        Sig  =     0.000400 (*)
      Sigc =    -0.193301 (*)        Sigc =     0.365908 (*)
(units: (*) = 10^-3 K/a)
*/

#include <stdio.h>
#include <stdlib.h>
#include "exactTestsFG.h"

int main() {

  const double SperA=31556926.0;  // seconds per year; 365.2422 days
  const double Cp=200.0;     // m;  magnitude of the perturbation in test G
  double year, r, HF, MF, HG, MG;
  double *z, *TF, *UF, *wF, *SigF, *SigcF, *TG, *UG, *wG, *SigG, *SigcG;
  double *mb; /* a block of memory */
  int j, Mz;

  printf("Enter  t and r  separated by newline");
  printf(" (in yrs and km, resp.; e.g. 500 500):\n");
  scanf("%lf",&year);
  scanf("%lf",&r);
  printf("Enter  z  values sep by newline (in m);");
  printf(" '-1' to end; e.g. 0 100 500 1500 -1:\n");

  z = (double *) malloc(501 * sizeof(double));
  if (z == NULL) { 
    fprintf(stderr, "simpleFG: couldn't allocate memory\n"); return -9999; }

  j=0;
  do {
    scanf("%lf",&z[j]);
    j++;
    if (j>490) printf("\n\n\nWARNING: enter -1 to stop soon!!!\n");
  } while (z[j-1]>=0.0);
  Mz=j-1;

  mb = (double *) malloc(10 * Mz * sizeof(double));
  if (mb == NULL) { 
    fprintf(stderr, "simpleFG: couldn't allocate memory\n"); return -9999; }
  TF=mb; UF=mb+Mz*sizeof(double); wF=mb+2*Mz*sizeof(double);
  SigF=mb+3*Mz*sizeof(double); SigcF=mb+4*Mz*sizeof(double);
  TG=mb+5*Mz*sizeof(double); UG=mb+6*Mz*sizeof(double);
  wG=mb+7*Mz*sizeof(double);   
  SigG=mb+8*Mz*sizeof(double); SigcG=mb+9*Mz*sizeof(double);  

  /* evaluate tests F and G */
  bothexact(0.0,r*1000.0,z,Mz,0.0,&HF,&MF,TF,UF,wF,SigF,SigcF);
  bothexact(year*SperA,r*1000.0,z,Mz,Cp,&HG,&MG,TG,UG,wG,SigG,SigcG);

  printf("\nResults:\n           Test F                         Test G\n");
  printf("(functions of r (resp. t and r) only):\n");
  printf("      H    = %12.6f (m)        H    = %12.6f (m)\n",HF,HG);
  printf("      M    = %12.6f (m/a)      M    = %12.6f (m/a)\n",
         MF*SperA,MG*SperA);
  for (j=0; j<Mz; j++) {
    printf("(z=%10.3f):\n",z[j]);
    printf("      T    = %12.6f (K)        T    = %12.6f (K)\n",TF[j],TG[j]);
    printf("      U    = %12.6f (m/a)      U    = %12.6f (m/a)\n",UF[j]*SperA,
           UG[j]*SperA);
    printf("      w    = %12.6f (m/a)      w    = %12.6f (m/a)\n",wF[j]*SperA,
           wG[j]*SperA);
    printf("      Sig  = %12.6f (*)        Sig  = %12.6f (*)\n",
           SigF[j]*SperA*1000.0,SigG[j]*SperA*1000.0);
    printf("      Sigc = %12.6f (*)        Sigc = %12.6f (*)\n",
           SigcF[j]*SperA*1000.0,SigcG[j]*SperA*1000.0);
  }
  printf("(units: (*) = 10^-3 K/a)\n");

  free(mb);
  return 0;
}
