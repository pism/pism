// Copyright (C) 2007 Ed Bueler
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include "iceModel.hh"

// note that by the ASCII-to-integer conversion, the title is
//        tn[int(singleCharName) - int('0')].title 
// and that the variable name (for NetCDF and Matlab) is 
//        tn[int(singleCharName) - int('0')].name 
const titleNname IceModel::tn[] = {
{"hor. speed at surface (m/a)", "csurf"}, // '0'
{"u at surface (m/a)", "usurf"}, // '1'
{"v at surface (m/a)", "vsurf"}, // '2'
{"w at surface (m/a)", "wsurf"}, // '3'
{}, // '4'
{}, // '5'
{}, // '6'
{}, // '7'
{}, // '8'
{}, // '9'
{}, // ':'
{}, // ';'
{}, // '<'
{}, // '='
{}, // '>'
{}, // '?'
{}, // '@'
{}, // 'A'
{"log(beta) (drag coeff; log_10(Pa s m^-1))", ""}, // 'B'
{"tau_c (till yield stress; bar=10^5Pa)", "tau_c"}, // 'C'
{"D (diffusivity; m^2/s)", "D"}, // 'D'
{"age of ice (years) at kd", "age_kd"}, // 'E'
{"geothermal heat flux (mW/m^2)", "ghf"}, // 'F'
{"grain size (mm) at kd", "G_kd"}, // 'G'
{"H (thickness; m)", "H"}, // 'H'
{}, // 'I'
{}, // 'J'
{}, // 'K'
{"basal melt water thickness (m)", "Hmelt"}, // 'L'
{}, // 'M'
{"SPECIAL; NOT A TITLE", ""}, // 'N' // for both NuView[0] = "(nu*H)_t (I offset)" and NuView[1] = ...; special case!
{}, // 'O'
{}, // 'P'
{}, // 'Q'
{"basal frictional heating (mW/m^2)", ""}, // 'R'
{"Sigma (strain heating; K/a) at kd", "Sigma_kd"}, // 'S'
{"T (temperature; K) at kd", "T_kd"}, // 'T'
{"uvbar[0] (velocity on stag grid; m/a)", ""}, // 'U'
{"uvbar[1] (velocity on stag grid; m/a)", ""}, // 'V'
{}, // 'W'
{"u (velocity; m/a) at kd", "u_kd"}, // 'X'
{"v (velocity; m/a) at kd", "v_kd"}, // 'Y'
{"w (velocity; m/a) at kd", "w_kd"}, // 'Z'
{}, // '['
{}, // '\'
{}, // ']'
{}, // '^'
{}, // '_'
{}, // '`'
{"M (surface accum rate; m/a)", "accum"}, // 'a'
{"b (bed elev; m above sea level)", "bed"}, // 'b'
{"log(speed) (log_10(m/a))", "cbar"}, // 'c'
{}, // 'd'
{"age of ice (years) at id,jd", ""}, // 'e'
{"thickening rate dH/dt (m/a)", "dHdt"}, // 'f'
{"grain size (mm) at id,jd", ""}, // 'g'
{"h (surface elev; m above sea level)", "h"}, // 'h'
{"nu*H (I offset)", ""}, // 'i'
{"nu*H (J offset)", ""}, // 'j'
{}, // 'k'  // currently in use so '-d k' gives same result as -ksp_monitor
{"basal melt rate (m/a)", ""}, // 'l'
{"mask (1=SHEET, 2=DRAG, 3=FLOAT)", "mask"}, // 'm'
{"log_10(nu*H)", ""}, // 'n'
{}, // 'o'
{"bed uplift rate (m/a)", "dbdt"}, // 'p'
{"log(basal sliding speed) (log_10(m/a))", ""}, // 'q'
{"suRface temperature (K)", "Ts"}, // 'r'
{"Sigma (strain heating; K/a) at id,jd", ""}, // 's'
{"T (temperature; K) at id,jd", ""}, // 't'
{"ubar (velocity; m/a)", "ubar"}, // 'u'
{"vbar (velocity; m/a)", "vbar"}, // 'v'
{}, // 'w'
{"u (velocity; m/a) at id,jd", ""}, // 'x'
{"v (velocity; m/a) at id,jd", ""}, // 'y'
{"w (velocity; m/a) at id,jd", ""}, // 'z'
};  // end of static initialization of IceModel::tn


