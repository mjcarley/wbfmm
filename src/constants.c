/* This file is part of WBFMM, a Wide-Band Fast Multipole Method code
 *
 * Copyright (C) 2019, 2020 Michael Carley
 *
 * WBFMM is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.  WBFMM is distributed in the
 * hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with WBFMM.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <math.h>

#include <glib.h>

#include "wbfmm-private.h"

/*table of \cos m\pi/4 for rotations on upward pass*/
gdouble CmPI_4[] =
  {1, M_SQRT1_2, 0, -M_SQRT1_2, -1, -M_SQRT1_2, 0, M_SQRT1_2, 1} ;
/*table of \cos n\pi/2 for rotations*/
gdouble CnPI_2[] = {1.0, 0.0, -1.0, 0.0} ;

gfloat CmPI_4f[] =
  {1, M_SQRT1_2, 0, -M_SQRT1_2, -1, -M_SQRT1_2, 0, M_SQRT1_2, 1} ;
/*table of \cos n\pi/2 for rotations*/
gfloat CnPI_2f[] = {1.0, 0.0, -1.0, 0.0} ;

gdouble *_wbfmm_SS_coefficients_laplace = NULL,
  *_wbfmm_RR_coefficients_laplace = NULL,
  *_wbfmm_SR_coefficients_laplace = NULL ;
gfloat *_wbfmm_SS_coefficients_laplace_f = NULL,
  *_wbfmm_RR_coefficients_laplace_f = NULL,
  *_wbfmm_SR_coefficients_laplace_f = NULL ;
gint _wbfmm_translation_Nmax = 0 ;

gdouble _wbfmm_shifts_th[49] = {0.0} ;
gdouble _wbfmm_shifts_ph[17] = {0.0} ;
gdouble _wbfmm_shifts_r[WBFMM_SHIFTS_R_NUMBER] = {0.0} ;
gfloat _wbfmm_shifts_th_f[49] = {0.0} ;
gfloat _wbfmm_shifts_ph_f[17] = {0.0} ;
gfloat _wbfmm_shifts_r_f[WBFMM_SHIFTS_R_NUMBER] = {0.0} ;
