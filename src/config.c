/* This file is part of WBFMM, a Wide-Band Fast Multipole Method code
 *
 * Copyright (C) 2019 Michael Carley
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <math.h>
#include <glib.h>
#include <string.h>
#include <stdio.h>

#include <wbfmm.h>

#include "wbfmm-private.h"

gint wbfmm_library_config(wbfmm_library_config_t *c)

{
  c->real_size = sizeof(WBFMM_REAL) ;

  c->switches = WBFMM_COMPILER_FLAGS ;
  
#ifdef HAVE_FMA_INSTRUCTIONS
  c->fma = TRUE ;
#else
  c->fma = FALSE ;
#endif /*HAVE_FMA_INSTRUCTIONS*/

#ifdef HAVE_AVX_INSTRUCTIONS
  c->avx = TRUE ;
#else
  c->avx = FALSE ;
#endif /*HAVE_AVX_INSTRUCTIONS*/

#ifdef HAVE_AVX2_INSTRUCTIONS
  c->avx2 = TRUE ;
#else
  c->avx2 = FALSE ;
#endif /*HAVE_AVX_INSTRUCTIONS*/

#ifdef _OPENMP
  c->openmp = TRUE ;
#else
  c->openmp = FALSE ;
#endif /*_OPENMP*/

  return 0 ;
}

gint wbfmm_library_config_print(wbfmm_library_config_t *c, FILE *f)

{
  fprintf(f,
	  "AVX extensions: %s\n"
	  "FMA extensions: %s\n"
	  "OpenMP:         %s\n"
	  "precision:      %s\n"
	  "compiler flags: %s\n",
	  yes_if_true(c->avx),
	  yes_if_true(c->fma),
	  yes_if_true(c->openmp),
#ifdef WBFMM_SINGLE_PRECISION
	  "single",
#else
	  "double",
#endif /*WBFMM_SINGLE_PRECISION*/
	  c->switches
	  ) ;
  
  return 0 ;
}
