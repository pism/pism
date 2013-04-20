/*
   Copyright (C) 2011 Ed Bueler
  
   This file is part of PISM.
  
   PISM is free software; you can redistribute it and/or modify it under the
   terms of the GNU General Public License as published by the Free Software
   Foundation; either version 3 of the License, or (at your option) any later
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

$ ./simpleO
Enter z in m (0 <= z < 3000 for ice and -1000 < z < 0 for bedrock):  3000
Results from Test O:
   T(z)  = 223.150 (K)          [absolute temperature]
   T_m   = 271.034 (K)          [pressure-melting (abs.) temperature at base]
   q_i   = 0.033519 (W m-2)     [upward heat flux in ice]
   q_bed = 0.042000 (W m-2)     [        "        in bedrock]
   bmelt = 8.80550e-04 (m a-1)
         = 0.88055 (mm a-1)     [ice-equivalent basal melt rate]

*/

#include <stdio.h>
#include "exactTestO.h"

int main() {

  double       z, TT, Tm, qice, qbed, bmelt;
  int          scanret;
  const double secpera = 31556926.0;  /* seconds per year; 365.2422 days */
  const double H0 = 3000.0, B0 = -1000.0; 
 
  printf("Enter z in m (0 <= z < %4.0f for ice and %5.0f < z < 0 for bedrock):  ",
         H0,B0);
  scanret = scanf("%lf",&z);
  if (scanret != 1) {
    printf("... input error; exiting\n"); 
    return 1;
  }

  exactO(z, &TT, &Tm, &qice, &qbed, &bmelt);

  printf("Results from Test O:\n");
  printf("   T(z)  = %.3f (K)          [absolute temperature]\n",TT);
  printf("   T_m   = %.3f (K)          [pressure-melting (abs.) temperature at base]\n",Tm);
  printf("   q_i   = %f (W m-2)     [upward heat flux in ice]\n",qice);
  printf("   q_bed = %f (W m-2)     [        \"        in bedrock]\n",qbed);
  printf("   bmelt = %.5e (m a-1)\n"
         "         = %.5f (mm a-1)     [ice-equivalent basal melt rate]\n",
         bmelt*secpera,1000.0*bmelt*secpera);

  return 0;
}

