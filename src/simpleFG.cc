// Copyright (C) 2004-2006 Jed Brown and Ed Bueler
//
// This file is part of Pism.
//
// Pism is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// Pism is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with Pism; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include <cstdio>
#include "exactTestsFG.hh"

int main() {

  const double SperA=31556926.0;  // seconds per year; 365.2422 days
  const double Cp=200.0;     // m;  magnitude of the perturbation in test G
  double year, r, HF, MF, HG, MG;
  double *z, *TF, *UF, *wF, *SigF, *SigcF, *TG, *UG, *wG, *SigG, *SigcG;
  int j, Mz;

  printf("Enter  t and r  separated by newline (in yrs and km, resp.):\n");
  scanf("%lf",&year);
  scanf("%lf",&r);
  printf("Enter  z  values separated by newline (in m); enter '-1' to end:\n");
  z=new double[500];
  j=0;
  do {
    scanf("%lf",&z[j]);
    j++;
    if (j>400) printf("\n\n\nERROR: too many z values!!!\n");
  } while (z[j-1]>=0.0);
  Mz=j-1;

  TF=new double[Mz];   UF=new double[Mz];   wF=new double[Mz];
  SigF=new double[Mz];   SigcF=new double[Mz];
  TG=new double[Mz];   UG=new double[Mz];   wG=new double[Mz];
  SigG=new double[Mz];   SigcG=new double[Mz];

  /* evaluate tests F and G */
  bothexact(0.0,r*1000.0,z,Mz,0.0,HF,MF,TF,UF,wF,SigF,SigcF);
  bothexact(year*SperA,r*1000.0,z,Mz,Cp,HG,MG,TG,UG,wG,SigG,SigcG);

  printf("\nResults:\n           Test F                         Test G\n");
  printf("(functions of r (resp. t and r) only):\n");
  printf("      H    = %12.6f (m)        H    = %12.6f (m)\n",HF,HG);
  printf("      M    = %12.6f (m/a)      M    = %12.6f (m/a)\n",MF*SperA,MG*SperA);
  for (j=0; j<Mz; j++) {
    printf("(z=%10.3f):\n",z[j]);
    printf("      T    = %12.6f (K)        T    = %12.6f (K)\n",TF[j],TG[j]);
    printf("      U    = %12.6f (m/a)      U    = %12.6f (m/a)\n",UF[j]*SperA,UG[j]*SperA);
    printf("      w    = %12.6f (m/a)      w    = %12.6f (m/a)\n",wF[j]*SperA,wG[j]*SperA);
    printf("      Sig  = %12.6f (*)        Sig  = %12.6f (*)\n",
           SigF[j]*SperA*1000.0,SigG[j]*SperA*1000.0);
    printf("      Sigc = %12.6f (*)        Sigc = %12.6f (*)\n",
           SigcF[j]*SperA*1000.0,SigcG[j]*SperA*1000.0);
  }
  printf("(units: (*) = 10^-3 K/a)\n");

  delete [] z;
  delete [] TF;  delete [] UF;  delete [] wF;  delete [] SigF;  delete [] SigcF;
  delete [] TG;  delete [] UG;  delete [] wG;  delete [] SigG;  delete [] SigcG;
  return 0;
}
