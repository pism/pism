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


int IceModel::cIndex(const char singleCharName) {
  return int(singleCharName) - int('0');
}


//! Holds short titles for graphical viewers, which are also names used in output Matlab files.
/*! 
This array of strings holds names and titles.  These names and titles are for <i>views</i> of the
internal state of PISM.  That is, they are \em not names and title for <i>parts</i> of the 
internal state of PISM.  Thus, for example, there is \em no one-to-one correspondence between 
the internal variables stored in PETSc Vecs and the entries in this array; it is more complicated 
than that.  

This array is used for both variable names \em and short titles in writeMatlabVars(), among
other places.  It is accessed for short titles \em only in createViewers().

This array is usually referenced by a single character.  Note that by the 
ASCII-to-integer conversion, for instance, the variable name corresponding to the 
single character name 'A' is <tt>tn[int('A') - int('0')].title</tt>.  Such
indexing is facilitated by cIndex().
 */
const titleNname IceModel::tn[] = {
{"hor. speed at surface (m/a)", "csurf"}, // '0'
{"u at ice surface (m/a)", "usurf"}, // '1'
{"v at ice surface (m/a)", "vsurf"}, // '2'
{"w at ice surface (m/a)", "wsurf"}, // '3'
{"u at ice base (m/a)", "ub"}, // '4'
{"v at ice base (m/a)", "vb"}, // '5'
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
{"snow accumulation rate (ice equivalent m/a)", "snow_accum_pera"}, // 'A'
{}, // 'B'
{"tau_c (till yield stress; kPa)", "tau_c"}, // 'C'
{"D (diffusivity; m^2/s)", "D"}, // 'D'
{"age of ice (years) at kd", "age_kd"}, // 'E'
{"geothermal heat flux (mW/m^2)", "ghf"}, // 'F'
{}, // 'G'
{"H (thickness; m)", "H"}, // 'H'
{"till friction angle (deg)", "tillphi"}, // 'I'
{}, // 'J'
{}, // 'K'
{"basal melt water thickness (m)", "Hmelt"}, // 'L'
{}, // 'M'
{"SPECIAL; NOT A TITLE", ""}, // 'N' // for both NuView[0] = "(nu*H)_t (I offset)" and NuView[1] = ...; special case!
{}, // 'O'
{}, // 'P'
{"basal driving stress (kPa)", "f_basal"}, // 'Q'
{"basal frictional heating (mW/m^2)", "basal_heating"}, // 'R'
{"Sigma (strain heating; K/a) at kd", "Sigma_kd"}, // 'S'
{"T (temperature; K) at kd", "T_kd"}, // 'T'
{"uvbar[0] (velocity on stag grid; m/a)", "ubar_staggered"}, // 'U'
{"uvbar[1] (velocity on stag grid; m/a)", "vbar_staggered"}, // 'V'
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
{"ice surface mass balance (m/a)", "accum_pera"}, // 'a'
{"b (bed elev; m above sea level)", "b"}, // 'b'
{"log(speed) (log_10(m/a))", "cbar"}, // 'c'
{}, // 'd'
{"age of ice (years) at id,jd", "age_sounding"}, // 'e'
{"thickening rate dH/dt (m/a)", "dHdt"}, // 'f'
{}, // 'g'
{"h (surface elev; m above sea level)", "h"}, // 'h'
{"nu*H (I offset)", "nu_H_ioffset"}, // 'i'
{"nu*H (J offset)", "nu_H_joffset"}, // 'j'
{}, // 'k'  // currently in use so '-d k' gives same result as -ksp_monitor
{"basal melt rate (m/a)", "basal_melt_rate"}, // 'l'
{"mask (1=SHEET, 2=DRAG, 3=FLOAT)", "mask"}, // 'm'
{"log_10(nu*H)", "log_nu_H"}, // 'n'
{}, // 'o'
{"bed uplift rate (m/a)", "dbdt_pera"}, // 'p'
{"log(basal sliding speed) (log_10(m/a))", "log_basal_sliding_speed"}, // 'q'
{"suRface temperature (K)", "Ts"}, // 'r'
{"Sigma (strain heating; K/a) at id,jd", "strain_heating"}, // 's'
{"T (temperature; K) at id,jd", "T_sounding"}, // 't'
{"ubar (velocity; m/a)", "ubar"}, // 'u'
{"vbar (velocity; m/a)", "vbar"}, // 'v'
{}, // 'w'
{"u (velocity; m/a) at id,jd", "u_sounding"}, // 'x'
{"v (velocity; m/a) at id,jd", "v_sounding"}, // 'y'
{"w (velocity; m/a) at id,jd", "w_sounding"}, // 'z'
};  // end of static initialization of IceModel::tn

