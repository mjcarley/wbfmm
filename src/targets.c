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
#include <stdio.h>
#include <inttypes.h>

#include <glib.h>

#include <wbfmm.h>

#include "wbfmm-private.h"

wbfmm_target_list_t *WBFMM_FUNCTION_NAME(wbfmm_target_list_new)(wbfmm_tree_t *t,
								guint npts,
								gboolean grad)

{
  wbfmm_target_list_t *l ;
  guint nc, nr, nb ;
  
  if ( t->size != sizeof(WBFMM_REAL) )
    g_error("%s: mixed precision not implemented\n"
	    "  (size of tree data type (%lu) not equal to "
	    "size of requested target type (%lu))",
	    __FUNCTION__, t->size, sizeof(WBFMM_REAL)) ;

  g_assert(grad == FALSE) ;
  
  l = (wbfmm_target_list_t *)g_malloc0(sizeof(wbfmm_target_list_t)) ;

  wbfmm_target_tree(l) = t ;
  wbfmm_target_point_number_max(l) = npts ;
  wbfmm_target_point_number(l) = 0 ;
  wbfmm_target_gradient(l) = grad ;

  l->size = sizeof(WBFMM_REAL) ;
  
  l->ip = (guint *)g_malloc0(npts*sizeof(guint)) ;
  l->boxes = (guint32 *)g_malloc0(npts*sizeof(guint32)) ;
  l->points = (gchar *)g_malloc0(npts*3*sizeof(WBFMM_REAL)) ;

  /*size the memory allocation, based on order of leaf-level regular
    expansions*/
  nr = t->order_r[t->depth] ;
  
  switch ( wbfmm_tree_problem(t) ) {
  default:
    g_error("%s: unrecognized problem type %u",
	    __FUNCTION__, wbfmm_tree_problem(t)) ;
    break ;
  case WBFMM_PROBLEM_LAPLACE:
    nc = (nr+1)*(nr+1) ;
    break ;
  case WBFMM_PROBLEM_HELMHOLTZ:
    g_assert_not_reached() ;
    break ;
  }

  if ( grad ) nc *= 4 ;

  l->nc = nc ;

  /*number of leaf boxes*/
  nb = 1 << (3*(t->depth)) ;

  l->cfft = g_malloc0(nc*nb*sizeof(WBFMM_REAL)) ;
  
  return l ;
}
