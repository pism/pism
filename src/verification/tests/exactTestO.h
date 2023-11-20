/*
   Copyright (C) 2011, 2016 Ed Bueler and Constantine Khroulev

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

#ifndef __exactTestO_h
#define __exactTestO_h 1

#ifdef __cplusplus
extern "C"
{
#endif

/*
Basal-melt rate computation exact solution.  Utterly straightforward arithmetic.

See also src/exact/simpleO.c.

Fills this z-dependent quantity:
     TT    = temperature at z, whether in ice (z >= 0) or in bedrock (z < 0)
Also fills these z-independent quantities:
     Tm    = the pressure-melting temperature (K) at the base of the ice (z=0)
     qice  = upward heat flux (J m-2 s-1) within the ice, 0 <= z < H0
     qbed  = upward heat flux (J m-2 s-1) within the bedrock, -B0 < z < 0
     bmelt = exact solution for melt rate (ice-equivalent m s-1) at base (z=0)
 */

struct TestOParameters {
  double TT;
  double Tm;
  double qice;
  double qbed;
  double bmelt;
};

struct TestOParameters exactO(double z);

#ifdef __cplusplus
}
#endif

#endif  /* __exactTestO_h */
