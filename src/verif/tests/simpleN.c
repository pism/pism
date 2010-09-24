/*
   Copyright (C) 2010 Ed Bueler
  
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

$ ./simpleN
Enter  x  (in km; -400.0 <= x <= 400):   200.0
Enter  t  (in years; 0.0 <= t <= ??00):  0.0
Results from Test N:
     H = 
     u =
     tau_d = rho g H h_x = 
     tau_11 = 
  
*/


#include <stdio.h>
#include "exactTestN.h"

int main() {

  double       t, x, u, ux, h, b, hx, taud, taub, tau11, C;
  const double secpera=31556926.0;  /* seconds per year; 365.2422 days */
  int          scanret, retvalN;
  
  printf("Enter  t and x  (in year and km, resp.; t >= 0, -400.0 <= x <= 400;\n"
         "                 e.g. 0.0 300.0):   ");
  scanret = scanf("%lf",&t);
  if (scanret != 1) {
    printf("... input error; exiting\n"); 
    return 1;
  }
  scanret = scanf("%lf",&x);
  if (scanret != 1) {
    printf("... input error; exiting\n"); 
    return 1;
  }

  retvalN = exactN(t*secpera, x*1000.0,
                   &u, &ux, &h, &b, &hx, &taud, &taub, &tau11, &C);

  if (retvalN) {
    printf("SIMPLEN ERROR:  (t,x)  out of allowed domain  t >= 0, -400.0 <= x <= 400\n"
           "   ... ending ...\n");
    return 1;
  } else {
    printf("Results from Test N:\n");
    printf("     u = %10.5f (m/a),     u_x = %10.4e (1/a),\n"
           "     h = %10.5f (m),       b = %10.4f (m),      h_x = %10.5f,\n"
           "     taud = %10.4e (Pa),   taub = %10.4e (Pa),  tau11 = %10.4e (Pa)\n"
           "     C = %10.4e (Pa s m-1)\n",
           u * secpera, ux * secpera, h, b, hx, taud, taub, tau11, C);
    return 0;
  }
}

