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
#include <string.h>

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

  wbfmm_target_list_tree(l) = t ;
  wbfmm_target_list_point_number_max(l) = npts ;
  wbfmm_target_list_point_number(l) = 0 ;
  wbfmm_target_list_gradient(l) = grad ;

  l->size = sizeof(WBFMM_REAL) ;
  
  l->ip = (guint *)g_malloc0(npts*sizeof(guint)) ;
  l->boxes = (guint32 *)g_malloc0(npts*sizeof(guint32)) ;

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

static gint compare_morton_indexed(gconstpointer a, gconstpointer b,
				   gpointer data)

{
  guint i, j ;
  wbfmm_tree_t *t = data ;
  guint64 mi, mj ;
  WBFMM_REAL *xi, *xj ;

  i = *((guint *)a) ; j = *((guint *)b) ;

  xi = wbfmm_tree_point_index(t, i) ; 
  xj = wbfmm_tree_point_index(t, j) ;
  /*Morton codes*/
  mi = WBFMM_FUNCTION_NAME(wbfmm_point_index_3d)(xi, wbfmm_tree_origin(t), 
						 wbfmm_tree_width(t)) ;
  mj = WBFMM_FUNCTION_NAME(wbfmm_point_index_3d)(xj, wbfmm_tree_origin(t), 
						 wbfmm_tree_width(t)) ;

  if ( mi < mj ) return -1 ;
  if ( mi > mj ) return  1 ;

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_target_list_add_points)(wbfmm_target_list_t *l,
						       gpointer pts, guint npts,
						       gsize pstr)

{
  gint i, j, k, idx, ns, nsb, nnbr, nb ;
  guint64 neighbours[27] ;
  WBFMM_REAL *x, *xt, D, *xs, *csrc, r ;
  wbfmm_tree_t *t = wbfmm_target_list_tree(l) ;
  guint level = t->depth ;
  wbfmm_box_t box, *boxes ;
  
  if ( l->size != sizeof(WBFMM_REAL) )
    g_error("%s: mixed precision not implemented\n"
	    "  (size of target list data type (%lu) not equal to "
	    "size of requested target type (%lu))",
	    __FUNCTION__, l->size, sizeof(WBFMM_REAL)) ;

  if ( npts > wbfmm_target_list_point_number_max(l) ) 
    g_error("%s: too many points (%u) for target list (%u)",
	    __FUNCTION__, npts, wbfmm_target_list_point_number_max(l)) ;

  wbfmm_target_list_point_number(l) = npts ;
  l->points = (gchar *)pts ; l->pstr = pstr ;

  xt = wbfmm_tree_origin(t) ;
  D = wbfmm_tree_width(t) ;

  for ( i = 0 ; i < npts ; i ++ ) {
    /* x = wbfmm_tree_point_index(t, i) ; */
    x = wbfmm_target_list_point_index(l, i) ;
    if ( x[0] <= xt[0] || x[0] >= xt[0] + D ||
	 x[1] <= xt[1] || x[1] >= xt[1] + D ||
	 x[2] <= xt[2] || x[2] >= xt[2] + D )
      g_error("%s: point (%g,%g,%g) does not lie in box of width %g at "
	      "(%g,%g,%g)",
	      __FUNCTION__, x[0], x[1], x[2], D, xt[0], xt[1], xt[2]) ;
	 
    l->ip[i] = i ;
  }

  /*sort points on the Morton index*/
  g_qsort_with_data(l->ip, npts, sizeof(guint), compare_morton_indexed, 
		    (gpointer)t) ;

  for ( i = 0 ; i < npts ; i ++ ) {
    x = wbfmm_target_list_point_index(l, i) ;
    l->boxes[i] = WBFMM_FUNCTION_NAME(wbfmm_point_box)(t, t->depth, x) ;
  }
  
  nb = 1 << (3*(level)) ;
  nsb = 0 ; boxes = t->boxes[level] ;
  for ( i = 0 ; i < nb ; i ++ ) {
    nnbr = wbfmm_box_neighbours(level, i, neighbours) ;
    for ( j = 0 ; j < nnbr ; j ++ ) {
      box = boxes[neighbours[j]] ;
      nsb += box.n ;
    }
  }

  l->ibox = (gint *)g_malloc0((nb+1)*sizeof(gint)) ;
  l->isrc = (gint *)g_malloc0(nsb*sizeof(gint)) ;

  /*list of indices of source points lying in the near field of each
    box*/
  for ( i = 0 ; i < nb ; i ++ ) {
    l->ibox[i+1] = l->ibox[i] ;
    nnbr = wbfmm_box_neighbours(level, i, neighbours) ;
    for ( j = 0 ; j < nnbr ; j ++ ) {
      box = boxes[neighbours[j]] ;
      for ( k = 0 ; k < box.n ; k ++ ) {
	l->isrc[l->ibox[i+1]] = t->ip[box.i+k] ;
	l->ibox[i+1] ++ ;
      }
    }
  }

  /*number of source coefficients*/
  ns = 0 ;
  for ( i = 0 ; i < npts ; i ++ ) {
    j = l->boxes[i] ;
    ns += l->ibox[j+1] - l->ibox[j] ;
  }

  l->csrc =         g_malloc0(ns*sizeof(WBFMM_REAL)) ;
  l->ics  = (gint *)g_malloc0(npts*sizeof(WBFMM_REAL)) ;

  ns = 0 ;
  for ( i = 0 ; i < npts ; i ++ ) {
    l->ics[i] = ns ;
    csrc = &(((WBFMM_REAL *)(l->csrc))[ns]) ;
    k = l->boxes[i] ;
    x = wbfmm_target_list_point_index(l, i) ;
    for ( j = 0 ; j < l->ibox[k+1] - l->ibox[k] ; j ++ ) {
      idx = l->isrc[l->ibox[k] + j] ;
      xs = wbfmm_tree_point_index(t, idx) ;
      r = (xs[0]-x[0])*(xs[0]-x[0]) + (xs[1]-x[1])*(xs[1]-x[1]) +
	(xs[2]-x[2])*(xs[2]-x[2]) ;
      if ( r > 1e-12 ) {
	r = SQRT(r)*4.0*M_PI ;
	csrc[j] = 1.0/r ;
      } else {
	csrc[j] = 0.0 ;
      }
    }

    ns += l->ibox[k+1] - l->ibox[k] ;
  }
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_target_list_local_coefficients)(wbfmm_target_list_t *l, WBFMM_REAL *work)

{
  gint i, nr, nc ;
  guint level ;
  guint64 b ;
  wbfmm_tree_t *t = wbfmm_target_list_tree(l) ;
  WBFMM_REAL xb[3], wb, xf[3], *x, *cfft ;

  level = t->depth ;
  nr = t->order_r[level] ;
  nc = l->nc ;
  cfft = (WBFMM_REAL *)(l->cfft) ;
  
  for ( i = 0 ; i < wbfmm_target_list_point_number(l) ; i ++ ) {
    b = l->boxes[i] ;
    x = wbfmm_target_list_point_index(l, i) ;
    WBFMM_FUNCTION_NAME(wbfmm_tree_box_centre)(t, level, b, xb, &wb) ;
    xf[0] = x[0] - xb[0] ; xf[1] = x[1] - xb[1] ; xf[2] = x[2] - xb[2] ;
    WBFMM_FUNCTION_NAME(wbfmm_laplace_local_coefficients)(xf, nr, FALSE,
							  &(cfft[i*nc]),
							  work) ;
  }

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_target_list_local_field)(wbfmm_target_list_t *l,
							WBFMM_REAL *src,
							gint sstr,
							WBFMM_REAL *f)

{
  guint level ;
  gint ib, b, j, k, nc, nq, nr, idx ;
  wbfmm_tree_t *t = wbfmm_target_list_tree(l) ;
  wbfmm_box_t *boxes ;
  WBFMM_REAL *cfft, *eval, *csrc ;
  
  g_assert(wbfmm_tree_problem(t) == WBFMM_PROBLEM_LAPLACE) ;
  
  nq = wbfmm_tree_source_size(t) ;

  level = t->depth ;
  nr = t->order_r[level] ;
  nc = l->nc ;
  eval = (WBFMM_REAL *)(l->cfft) ;
  boxes = t->boxes[level] ;

  memset(f, 0, nq*wbfmm_target_list_point_number(l)*sizeof(WBFMM_REAL)) ;
  for ( ib = 0 ; ib < wbfmm_target_list_point_number(l) ; ib ++ ) {
    b = l->boxes[ib] ;
    cfft = boxes[b].mpr ;
    WBFMM_FUNCTION_NAME(wbfmm_laplace_expansion_apply)(cfft,
						       8*nq, nq,
						       &(eval[ib*nc]), nr,
						       &(f[ib*nq])) ;
    /*direct contributions from neighbour boxes*/
    csrc = &(((WBFMM_REAL *)(l->csrc))[l->ics[ib]]) ;
    for ( j = 0 ; j < l->ibox[b+1]-l->ibox[b] ; j ++ ) {
      idx = l->isrc[l->ibox[b]+j] ;
      for ( k = 0 ; k < nq ; k ++ ) {
    	f[ib*nq+k] += src[idx*sstr+k]*csrc[j] ;
      }
    }    
  }
  
  return 0 ;
}
