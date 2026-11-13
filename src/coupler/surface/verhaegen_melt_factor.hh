/* Copyright (C) 2026 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

namespace pism {

/*!
 * Implements equation (14) in Verhaegen and Huybrechts (2026),
 * https://doi.org/10.1029/2025JF008748.
 *
 * An empirical fit to the "Ostrem curve" averaged over the 21 debris-covered glacier
 * studies in their Table 1: melt is *enhanced* by a thin, dispersed cover, reaching a
 * maximum of about 1.4 times the clean-ice rate at the effective debris thickness
 * `h_e = 1.5` cm, drops back to the clean-ice rate near the critical thickness
 * `h_c = 4.1` cm, and is suppressed by thicker debris, becoming negligible beyond
 * 1.5 m (Popovnin and Rozova, 2002).
 *
 * @param[in] debris_thickness debris thickness, meters
 *
 * @return the ratio of the sub-debris melt rate to the clean ice melt rate
 */
double verhaegen_melt_factor(double debris_thickness);

} // namespace pism
