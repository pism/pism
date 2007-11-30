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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "exactTestL.h"

// heapsort from 
//   http://en.wikibooks.org/wiki/Algorithm_implementation/Sorting/Heapsort
// (found 7/7/07 by ELB;  see also http://en.wikipedia.org/wiki/Heapsort)
//  modified to have two index arrays "follow along" and get rearranged the same way
void heapsort_double_2indfollow(double arr[], int ia[], int ib[], unsigned int N) {
   unsigned int n = N, i = n/2, parent, child;
   double t;
   int tia, tib;

   for (;;) { /* Loops until arr is sorted */
      if (i > 0) { /* First stage - Sorting the heap */
            i--;           /* Save its index to i */
            t = arr[i];    /* Save parent value to t */
            tia = ia[i];
            tib = ib[i];
      } else {     /* Second stage - Extracting elements in-place */
            n--;           /* Make the new heap smaller */
            if (n == 0) return; /* When the heap is empty, we are done */
            t = arr[n];    /* Save last value (it will be overwritten) */
            tia = ia[n];
            tib = ib[n];
            arr[n] = arr[0]; /* Save largest value at the end of arr */
            ia[n] = ia[0];
            ib[n] = ib[0];
      }
      parent = i; /* We will start pushing down t from parent */
      child = i*2 + 1; /* parent's left child */
      while (child < n) { /* Shift operation - pushing the value of t down the heap */
            if (child + 1 < n  &&  arr[child + 1] > arr[child]) {
               child++; /* Choose the largest child */
            }
            if (arr[child] > t) { /* If any child is bigger than the parent */
               arr[parent] = arr[child]; /* Move the largest child up */
               ia[parent] = ia[child];
               ib[parent] = ib[child];
               parent = child; /* Move parent pointer to this child */
               child = parent*2 + 1; /* Find the next child */
            } else { break; /* t's place is found */ }
      }
      arr[parent] = t; /* We save t in the heap */
      ia[parent] = tia;
      ib[parent] = tib;
   }
}


/*  DOESN'T SEEM TO WORK; PRESUMABLY BECAUSE OF ra[1..n] vs ra[0..n-1] issues
//  heapsort adapted from Numerical Recipes in C, 2nd ed.
//  modified to have two index arrays "follow along" and get rearranged the same way
void hpsort2xtra(unsigned long n, double ra[], int inda[], int indb[]) {
  // sorts ra[1],...,ra[n] into ascending order,
  // and reorders inda[1..n],indb[1..n] the same way
  unsigned long i,ir,j,l;
  double rra;
  int dda, ddb;

  if (n < 2) return;
  l = (n >> 1) + 1;
  ir = n;
  for (;;) {
    if (l > 1) {
      rra = ra[--l];
      dda = inda[l];
      ddb = indb[l];
    } else {
      rra = ra[ir];
      dda = inda[ir];
      ddb = indb[ir];
      ra[ir] = ra[1];
      inda[ir] = inda[1];
      indb[ir] = indb[1];
      if (--ir == 1) {
        ra[1] = rra;
        inda[1] = dda;
        indb[1] = ddb;
        break;
      }
    }
    i = l;
    j = l + 1;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1])  j++;
      if (rra < ra[j]) {
        ra[i] = ra[j];
        inda[i] = inda[j];
        indb[i] = indb[j];
        i = j;
        j <<= 1;
      } else  j = ir + 1;
    }
    ra[i] = rra;
    inda[i] = dda;
    indb[i] = ddb;
  }
}
*/

int main() {

// set SMALL 1 to see more info on smaller grid example
#define SMALL 0

#if SMALL
  const int Mx = 11, My = 7;  // good for looking at grid
#else
  const int Mx = 36, My = 36;
#endif
  const double Lx = 700000, Ly = 700000;  // inside the actual margin
  const int imid = (Mx-1)/2, jmid = (My-1)/2;
  const double dx = (2.0*Lx)/(Mx-1), dy = (2.0*Ly)/(My-1);

  // to fill from test L
  double HH[Mx][My], bb[Mx][My], aa[Mx][My];  // slow method
  double HHf[Mx][My], bbf[Mx][My], aaf[Mx][My];  // fast method
  
  // compute radius at each point
  double rr[Mx][My];
  int i,j;
  double x,y;
  printf(" ---- Mx = %d, My = %d ----\n",Mx,My);
  for (i = 0; i < Mx; i++) {
    x = (i - imid) * dx;
    for (j = 0; j < My; j++) {
      y = (j - jmid) * dy;
      rr[i][j] = sqrt(x * x + y * y);
    }
  }

  // slow method: re-solve ODE for each grid point
  printf(" ---- start slow method ----  (takes 10 to 30 seconds for Mx=My=36, perhaps)\n");
  const double EPS_ABS = 1.0e-9, EPS_REL = 0.0; 
  for (i = 0; i < Mx; i++) {
    for (j = 0; j < My; j++) {
      exactL(rr[i][j],&HH[i][j],&bb[i][j],&aa[i][j],EPS_ABS,EPS_REL,1);
    }
  }
  printf(" ---- end slow ----\n");
  
#if SMALL
  // test heapsort
  double myra[8]   = { 1.0, 3.0, 0.0, -3.0, 4.0, 9.0, -2.0, 8.0 };
  int    myinda[8] = { 1, 1, 2, 2, 3, 3, 4, 4 };
  int    myindb[8] = { 1, 2, 1, 2, 1, 2, 1, 2 };
  for (i = 0; i < 8; i++)   { printf(" %6.2f",myra[i]); }   printf("\n");
  for (i = 0; i < 8; i++)   { printf(" %6d",myinda[i]); }   printf("\n");
  for (i = 0; i < 8; i++)   { printf(" %6d",myindb[i]); }   printf("\n");
//  hpsort2xtra(8,myra-1,myinda-1,myindb-1);  // note Num. Recipes arrays are 1,...,N
  heapsort_double_2indfollow(myra,myinda,myindb,8);
  for (i = 0; i < 8; i++)   { printf(" %6.2f",myra[i]); }   printf("\n");
  for (i = 0; i < 8; i++)   { printf(" %6d",myinda[i]); }   printf("\n");
  for (i = 0; i < 8; i++)   { printf(" %6d",myindb[i]); }   printf("\n");
#endif

  // fast method: order r, over all the grid, and then solve ODE once
  const int  MM = Mx * My;
  double  rraa[MM], HHaa[MM], bbaa[MM], aaaa[MM];
  int     iaa[MM], jaa[MM];
  printf(" ---- start fast method ----\n");
  for (i = 0; i < Mx; i++) {
    for (j = 0; j < My; j++) {
      rraa[i*My + j] = - rr[i][j];  // note minus, so ascending --> descending
      iaa[i*My + j] = i;
      jaa[i*My + j] = j;
    }
  }
//  hpsort2xtra(MM,rraa-1,iaa-1,jaa-1);  // sorts into ascending;  O(MM log MM)
  heapsort_double_2indfollow(rraa,iaa,jaa,MM);  // sorts into ascending;  O(MM log MM)
  int k;
  for (k = 0; k < MM; k++)   rraa[k] = - rraa[k];  // now rraa is decreasing
  exactL_list(rraa, MM, HHaa, bbaa, aaaa);  // get soln to test L at these points; solves ODE only once 
  for (k = 0; k < MM; k++) {
    i = iaa[k];   j = jaa[k];
    rr[i][j] = rraa[k];
    HHf[i][j] = HHaa[k];
    bbf[i][j] = bbaa[k];
    aaf[i][j] = aaaa[k];
  }
  printf(" ---- end fast ----\n");

#if SMALL
  // print grid result
  for (i = 0; i < Mx; i++) {
    printf("(i = %d:)\n",i);
    printf("      j: ");
    for (j = 0; j < My; j++) {
      printf("%15d ",j);
    }
    printf("\n      r: ");
    for (j = 0; j < My; j++) {
      printf("%15.3f ",rr[i][j]/1000.0);
    }
    printf("\n  Hslow: ");
    for (j = 0; j < My; j++) {
      printf("%15.10f ",HH[i][j]);
    }
    printf("\n  Hfast: ");
    for (j = 0; j < My; j++) {
      printf("%15.10f ",HHf[i][j]);
    }
    printf("\n");
  }
#endif

  double maxHdiff = 0.0, myHdiff;
  for (i = 0; i < Mx; i++) {
    for (j = 0; j < My; j++) {
      myHdiff = fabs(HH[i][j] - HHf[i][j]);
      maxHdiff = (myHdiff > maxHdiff) ? myHdiff : maxHdiff;
    }
  }
  printf(" (max |Hslow_ij - Hfast_ij|) / (max |Hslow_ij|) = %.16e\n",
         maxHdiff/HH[imid][jmid]);
  return 0;
}
