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

#include <stdio.h>
#include <math.h>
#include "exactTestI.h"

int exactI(const double m, const double x, const double y, 
           double *bed, double *tauc, double *u, double *v) {
  /* see exact solution for an ice stream sliding over plastic till described
     on pages 237 and 238 of C. Schoof 2006 "A variational approach to ice streams"
     J Fluid Mech 556 pp 227--251 */
  /* const double n = 3.0, p = 1.0 + 1.0/n; // p = 4/3 */
  
  const double L = 40e3;        /* = 40km; y in [-3L,3L]; full width is 6L = 240 km */
  const double aspect = 0.05;
  const double h0 = aspect * L; /* if aspect = 0.05 and L = 40km then h0 = 2000 m */
  const double theta = atan(0.001);   /* a slope of 1/1000, a la Siple streams */
  const double rho = 910, g = 9.81;  /* kg/m^3 and m/s^2, resp. */
  const double f = rho * g * h0 * tan(theta);  /* about 18 kPa given above rho,g,theta,aspect,L */
  const double W = pow(m+1.0,1.0/m) * L;  /* e.g. W = 1.2 * L if m = 10 */
  const double B = 3.7e8;       /* Pa s^{1/3}; hardness given on p. 239 of Schoof; why so big? */  

  const double s = fabs(y/L);

  double C0, C1, C2, C3, C4, z1, z2, z3, z4;
  
  *tauc = f * pow(s,m);

  /* to compute bed, assume bed(x=0)=0 and bed is sloping down (?) for increasing x;
     if tan(theta) = 0.001 and -Lx < x < Lx with Lx = 240km then bed(x=Lx) = -240 m */
  *bed = - x * tan(theta);

  /* formula (4.3) in Schoof; note u is indep of aspect because f/h0 ratio gives C0 */
  if (fabs(y) < W) {
    C0 = 2.0 * pow(f / (B * h0),3.0) * pow(L,4.0);
    /* printf("  C0*secpera = %10.5e\n",C0*31556926.0); */
    C1 = pow(m + 1.0, 4.0/m); C2 = (m+1.0) * C1; C3 = (m+1.0) * C2; C4 = (m+1.0) * C3;
    z1 = ( pow(s,4.0) - C1 ) / 4.0;
    z2 = ( pow(s,m+4.0) - C2 ) / ( (m+1.0) * (m+4.0) );
    z3 = ( pow(s,2.0*m+4.0) - C3 ) / ( (m+1.0)*(m+1.0) * (2.0*m+4.0) );
    z4 = ( pow(s,3.0*m+4.0) - C4 ) / ( pow((m+1.0),3.0) * (3.0*m+4.0) );
    /* printf("  u / C0 = %10.5e\n",- (z1 - 3.0 * z2 + 3.0 * z3 - z4)); */
    *u = - C0 * (z1 - 3.0 * z2 + 3.0 * z3 - z4);  /* comes out positive */
  } else {
    *u = 0.0;
  }
  *v = 0.0;  /* no transverse flow */
  return 0;
}
