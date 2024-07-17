/* This file is part of WBFMM, a Wide-Band Fast Multipole Method code
 *
 * Copyright (C) 2020 Michael Carley
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
								guint npts)
  

{
  wbfmm_target_list_t *l ;
  
  if ( t->size != sizeof(WBFMM_REAL) )
    g_error("%s: mixed precision not implemented\n"
	    "  (size of tree data type (%lu) not equal to "
	    "size of requested target type (%lu))",
	    __FUNCTION__, t->size, sizeof(WBFMM_REAL)) ;

  l = (wbfmm_target_list_t *)g_malloc0(sizeof(wbfmm_target_list_t)) ;

  wbfmm_target_list_tree(l) = t ;
  wbfmm_target_list_point_number_max(l) = npts ;
  wbfmm_target_list_point_number(l) = 0 ;

  l->size = sizeof(WBFMM_REAL) ;
  
  l->ip = (guint *)g_malloc0(npts*sizeof(guint)) ;
  l->boxes = (guint32 *)g_malloc0(npts*sizeof(guint32)) ;

  l->cfft = NULL ;
  l->nc = 0 ;

  return l ;
}

static gint compare_morton_indexed(gconstpointer a, gconstpointer b,
				   gpointer data)

{
  guint i, j ;
  wbfmm_target_list_t *l = data ;
  wbfmm_tree_t *t = l->t ;
  guint64 mi, mj ;
  WBFMM_REAL *xi, *xj ;

  i = *((guint *)a) ; j = *((guint *)b) ;

  xi = wbfmm_target_list_point_index(l, i) ;
  xj = wbfmm_target_list_point_index(l, j) ;
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
						       gpointer pts,
						       gsize pstr,
						       guint npts)

{
  gint i, j, k, nsb, nnbr, nb ;
  guint64 neighbours[27] ;
  WBFMM_REAL *x, *xt, D ;
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
  l->points = (char *)pts ; l->pstr = pstr ;
  /* l->normals = (char *)normals ; l->nstr = nstr ; */
  
  xt = wbfmm_tree_origin(t) ;
  D = wbfmm_tree_width(t) ;

  for ( i = 0 ; i < npts ; i ++ ) {
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
		    (gpointer)l) ;

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

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_laplace_target_list_local_coefficients)(wbfmm_target_list_t *l, guint source, WBFMM_REAL *work)

{
  gint i, j, k, npts, nr, nc, idx, tstr ;
  guint level ;
  guint64 b, ns ;
  wbfmm_tree_t *t = wbfmm_target_list_tree(l) ;
  WBFMM_REAL xb[3], wb, xf[3], *x, *cfft ;
  WBFMM_REAL *csrc, r, *xs, *n, nR[3] ;

  g_assert(l->t->problem == WBFMM_PROBLEM_LAPLACE) ;

  g_assert((source == WBFMM_SOURCE_MONOPOLE) ||
	   (source == WBFMM_SOURCE_DIPOLE) ||
	   (source == (WBFMM_SOURCE_MONOPOLE | WBFMM_SOURCE_DIPOLE))) ;

  l->source = source ;
  
  if ( (source == WBFMM_SOURCE_MONOPOLE) ||
       (source == WBFMM_SOURCE_DIPOLE) ) tstr = 1 ;
  if ( source == (WBFMM_SOURCE_MONOPOLE | WBFMM_SOURCE_DIPOLE) ) tstr = 2 ;
  
  level = t->depth ;
  nr = t->order_r[level] ;
  nc = l->nc ;
  cfft = (WBFMM_REAL *)(l->cfft) ;
  l->complex = FALSE ;
  
  for ( i = 0 ; i < wbfmm_target_list_point_number(l) ; i ++ ) {
    b = l->boxes[i] ;
    x = wbfmm_target_list_point_index(l, i) ;
    WBFMM_FUNCTION_NAME(wbfmm_tree_box_centre)(t, level, b, xb, &wb) ;
    xf[0] = x[0] - xb[0] ; xf[1] = x[1] - xb[1] ; xf[2] = x[2] - xb[2] ;
    WBFMM_FUNCTION_NAME(wbfmm_laplace_local_coefficients)(xf, nr, l->field,
							  &(cfft[i*nc]),
							  work) ;
  }

  /*number of source coefficients*/
  npts = wbfmm_target_list_point_number(l) ;
  ns = 0 ;
  for ( i = 0 ; i < npts ; i ++ ) {
    j = l->boxes[i] ;
    /* g_assert(l->ibox[j+1] >= l->ibox[j]) ; */
    ns += l->ibox[j+1] - l->ibox[j] ;
    /* fprintf(stderr, "%d: %lu %d\n", i, ns, l->ibox[j+1] - l->ibox[j]) ; */
  }
  l->ics  = (gint *)g_malloc0(npts*sizeof(WBFMM_REAL)) ;

  /* fprintf(stderr, "ns = %lu; npts = %d\n", ns, npts) ; */
  if ( l->field == WBFMM_FIELD_SCALAR ) {
    l->csrc =         g_malloc0(tstr*ns*sizeof(WBFMM_REAL)) ;

    ns = 0 ;
    if ( source == WBFMM_SOURCE_MONOPOLE ) {
      for ( i = 0 ; i < npts ; i ++ ) {
	l->ics[i] = ns ;
	csrc = &(((WBFMM_REAL *)(l->csrc))[ns*tstr]) ;
	k = l->boxes[i] ;
	x = wbfmm_target_list_point_index(l, i) ;
	for ( j = 0 ; j < l->ibox[k+1] - l->ibox[k] ; j ++ ) {
	  idx = l->isrc[l->ibox[k] + j] ;
	  xs = wbfmm_tree_point_index(t, idx) ;
	  r = (xs[0]-x[0])*(xs[0]-x[0]) + (xs[1]-x[1])*(xs[1]-x[1]) +
	    (xs[2]-x[2])*(xs[2]-x[2]) ;
	  if ( r > WBFMM_LOCAL_CUTOFF_RADIUS*WBFMM_LOCAL_CUTOFF_RADIUS ) {
	    r = SQRT(r)*4.0*M_PI ;
	    csrc[j*tstr+0] = 1.0/r ;
	  } else {
	    csrc[j*tstr+0] = 0.0 ;
	  }
	}
	ns += l->ibox[k+1] - l->ibox[k] ;
      }

      return 0 ;
    }

    if ( source == (WBFMM_SOURCE_MONOPOLE | WBFMM_SOURCE_DIPOLE) ) {
      for ( i = 0 ; i < npts ; i ++ ) {
	/* fprintf(stderr, "%d\n", i) ; */
	l->ics[i] = ns ;
	csrc = &(((WBFMM_REAL *)(l->csrc))[ns*tstr]) ;
	k = l->boxes[i] ;
	x = wbfmm_target_list_point_index(l, i) ;
	for ( j = 0 ; j < l->ibox[k+1] - l->ibox[k] ; j ++ ) {
	  idx = l->isrc[l->ibox[k] + j] ;
	  n = wbfmm_tree_normal_index(t, idx) ;
	  xs = wbfmm_tree_point_index(t, idx) ;
	  r = (xs[0]-x[0])*(xs[0]-x[0]) + (xs[1]-x[1])*(xs[1]-x[1]) +
	    (xs[2]-x[2])*(xs[2]-x[2]) ;
	  if ( r > WBFMM_LOCAL_CUTOFF_RADIUS*WBFMM_LOCAL_CUTOFF_RADIUS ) {
	    r = SQRT(r) ;
	    csrc[j*tstr+0] = 0.25*M_1_PI/r ;
	    csrc[j*tstr+1] =
	      ((x[0]-xs[0])*n[0] +
	       (x[1]-xs[1])*n[1] +
	       (x[2]-xs[2])*n[2])*0.25*M_1_PI/r/r/r ;
	  } else {
	    csrc[j*tstr+0] = csrc[j*tstr+1] = 0.0 ;
	  }
	}
	ns += l->ibox[k+1] - l->ibox[k] ;
      }

  /* fprintf(stderr, "ns = %d; npts = %d\n", ns, npts) ; */
      return 0 ;
    }
    g_assert_not_reached() ;
  }

  if ( l->field == WBFMM_FIELD_GRADIENT ) {
    g_assert_not_reached() ; /*stride information needs to be added*/
    l->csrc =         g_malloc0(3*tstr*ns*sizeof(WBFMM_REAL)) ;

    ns = 0 ;
    for ( i = 0 ; i < npts ; i ++ ) {
      l->ics[i] = ns ;
      csrc = &(((WBFMM_REAL *)(l->csrc))[3*ns]) ;
      k = l->boxes[i] ;
      x = wbfmm_target_list_point_index(l, i) ;
      for ( j = 0 ; j < l->ibox[k+1] - l->ibox[k] ; j ++ ) {
	idx = l->isrc[l->ibox[k] + j] ;
	xs = wbfmm_tree_point_index(t, idx) ;
	r = (xs[0]-x[0])*(xs[0]-x[0]) + (xs[1]-x[1])*(xs[1]-x[1]) +
	  (xs[2]-x[2])*(xs[2]-x[2]) ;
	if ( r > WBFMM_LOCAL_CUTOFF_RADIUS*WBFMM_LOCAL_CUTOFF_RADIUS ) {
	  r = SQRT(r) ;
	  nR[0] = (x[0] - xs[0])/r*0.25*M_1_PI ;
	  nR[1] = (x[1] - xs[1])/r*0.25*M_1_PI ;
	  nR[2] = (x[2] - xs[2])/r*0.25*M_1_PI ;
	  csrc[3*j+0] = -1.0/r/r*nR[0] ;
	  csrc[3*j+1] = -1.0/r/r*nR[1] ;
	  csrc[3*j+2] = -1.0/r/r*nR[2] ;
	} else {
	  csrc[3*j+0] = csrc[3*j+1] = csrc[3*j+2] = 0.0 ;
	}
      }
      ns += l->ibox[k+1] - l->ibox[k] ;
    }
    return 0 ;
  }

  g_assert_not_reached() ;
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_target_list_local_coefficients)(wbfmm_target_list_t *l, WBFMM_REAL k, WBFMM_REAL *work)

{
  gint i, j, kk, npts, nr, nc, ns, idx ;
  guint level ;
  guint64 b ;
  wbfmm_tree_t *t = wbfmm_target_list_tree(l) ;
  WBFMM_REAL xb[3], wb, xf[3], *x, *cfft ;
  WBFMM_REAL *csrc, r, *xs, nR[3], h0[2], h1[2] ;

  g_assert(l->t->problem == WBFMM_PROBLEM_HELMHOLTZ) ;

  level = t->depth ;
  nr = t->order_r[level] ;
  nc = l->nc ;
  cfft = (WBFMM_REAL *)(l->cfft) ;
  l->complex = TRUE ;

  for ( i = 0 ; i < wbfmm_target_list_point_number(l) ; i ++ ) {
    b = l->boxes[i] ;
    x = wbfmm_target_list_point_index(l, i) ;
    WBFMM_FUNCTION_NAME(wbfmm_tree_box_centre)(t, level, b, xb, &wb) ;
    xf[0] = x[0] - xb[0] ; xf[1] = x[1] - xb[1] ; xf[2] = x[2] - xb[2] ;
    WBFMM_FUNCTION_NAME(wbfmm_local_coefficients)(k, xf, nr, l->field,
						  &(cfft[i*nc]),
						  work) ;
  }

  /*number of source coefficients*/
  npts = wbfmm_target_list_point_number(l) ;
  ns = 0 ;
  for ( i = 0 ; i < npts ; i ++ ) {
    j = l->boxes[i] ;
    ns += l->ibox[j+1] - l->ibox[j] ;
  }
  l->ics  = (gint *)g_malloc0(npts*sizeof(WBFMM_REAL)) ;

  if ( l->field == WBFMM_FIELD_SCALAR ) {
    l->csrc =         g_malloc0(2*ns*sizeof(WBFMM_REAL)) ;

    ns = 0 ;
    for ( i = 0 ; i < npts ; i ++ ) {
      l->ics[i] = ns ;
      csrc = &(((WBFMM_REAL *)(l->csrc))[2*ns]) ;
      kk = l->boxes[i] ;
      x = wbfmm_target_list_point_index(l, i) ;
      for ( j = 0 ; j < l->ibox[kk+1] - l->ibox[kk] ; j ++ ) {
	idx = l->isrc[l->ibox[kk] + j] ;
	xs = wbfmm_tree_point_index(t, idx) ;

	r = (xs[0]-x[0])*(xs[0]-x[0]) + (xs[1]-x[1])*(xs[1]-x[1]) +
	  (xs[2]-x[2])*(xs[2]-x[2]) ;
	if ( r > WBFMM_LOCAL_CUTOFF_RADIUS*WBFMM_LOCAL_CUTOFF_RADIUS ) {
	  r = SQRT(r) ;
	  WBFMM_FUNCTION_NAME(wbfmm_bessel_h_init)(k*r, h0, h1) ;
	  h0[0] /= 4.0*M_PI ; h0[1] /= 4.0*M_PI ;
	  csrc[2*j+0] = h0[0] ;
	  csrc[2*j+1] = h0[1] ;
	} else {
	  csrc[2*j+0] = csrc[2*j+1] = 0.0 ;
	}
      }
      ns += l->ibox[kk+1] - l->ibox[kk] ;
    }
    return 0 ;
  }

  if ( l->field == WBFMM_FIELD_GRADIENT ) {
    l->csrc =         g_malloc0(6*ns*sizeof(WBFMM_REAL)) ;

    ns = 0 ;
    for ( i = 0 ; i < npts ; i ++ ) {
      l->ics[i] = ns ;
      csrc = &(((WBFMM_REAL *)(l->csrc))[6*ns]) ;
      kk = l->boxes[i] ;
      x = wbfmm_target_list_point_index(l, i) ;
      for ( j = 0 ; j < l->ibox[kk+1] - l->ibox[kk] ; j ++ ) {
  	idx = l->isrc[l->ibox[kk] + j] ;
  	xs = wbfmm_tree_point_index(t, idx) ;
  	r = (xs[0]-x[0])*(xs[0]-x[0]) + (xs[1]-x[1])*(xs[1]-x[1]) +
  	  (xs[2]-x[2])*(xs[2]-x[2]) ;
  	if ( r > WBFMM_LOCAL_CUTOFF_RADIUS*WBFMM_LOCAL_CUTOFF_RADIUS ) {
  	  r = SQRT(r) ;
  	  nR[0] = (x[0] - xs[0])/r*0.25*M_1_PI ;
  	  nR[1] = (x[1] - xs[1])/r*0.25*M_1_PI ;
  	  nR[2] = (x[2] - xs[2])/r*0.25*M_1_PI ;
	  WBFMM_FUNCTION_NAME(wbfmm_bessel_h_init)(k*r, h0, h1) ;
	  csrc[6*j+0] = -k*h1[0]*nR[0] ;
	  csrc[6*j+1] = -k*h1[1]*nR[0] ;
	  csrc[6*j+2] = -k*h1[0]*nR[1] ;
	  csrc[6*j+3] = -k*h1[1]*nR[1] ;
	  csrc[6*j+4] = -k*h1[0]*nR[2] ;
	  csrc[6*j+5] = -k*h1[1]*nR[2] ;
  	} else {
  	  csrc[6*j+0] = csrc[6*j+1] = csrc[6*j+2] = 
	    csrc[6*j+3] = csrc[6*j+4] = csrc[6*j+5] = 0.0 ;
  	}
      }
      ns += l->ibox[kk+1] - l->ibox[kk] ;
    }
      return 0 ;
  }

  g_assert_not_reached() ;
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_target_list_local_field)(wbfmm_target_list_t *l,
							WBFMM_REAL *src,
							gint sstr,
							WBFMM_REAL *nsrc,
							gint nstr,
							WBFMM_REAL *f,
							gint fstr)

{
  guint level ;
  gint ib, b, j, k, nc, nq, nr, idx, tstr, nf ;
  wbfmm_tree_t *t = wbfmm_target_list_tree(l) ;
  wbfmm_box_t *boxes ;
  WBFMM_REAL *cfft, *eval, *csrc ;

  nq = wbfmm_tree_source_size(t) ;

  level = t->depth ;
  nr = t->order_r[level] ;
  nc = l->nc ;
  eval = (WBFMM_REAL *)(l->cfft) ;
  boxes = t->boxes[level] ;

  if ( l->t->problem == WBFMM_PROBLEM_LAPLACE ) {
    if ( !( l->field  == WBFMM_FIELD_SCALAR ||
	    l->field  == WBFMM_FIELD_GRADIENT ||
	    l->field  == WBFMM_FIELD_CURL ) ) {
      g_error("%s: field type %d not implemented for Laplace problem",
	      __FUNCTION__, l->field) ;
    }
    
    switch ( l->field ) {
    case WBFMM_FIELD_SCALAR:
      nf = 1 ;
      if ( l->source == WBFMM_SOURCE_MONOPOLE ) {
	tstr = 1 ;
	for ( ib = 0 ; ib < wbfmm_target_list_point_number(l) ; ib ++ ) {
	  b = l->boxes[ib] ;
	  /*direct contributions from neighbour boxes*/
	  csrc = &(((WBFMM_REAL *)(l->csrc))[tstr*(l->ics[ib])]) ;
	  for ( j = 0 ; j < l->ibox[b+1]-l->ibox[b] ; j ++ ) {
	    idx = l->isrc[l->ibox[b]+j] ;
	    for ( k = 0 ; k < nq ; k ++ ) {
	      f[ib*fstr+k] += src[idx*sstr+k]*csrc[j*tstr+0] ;
	    }
	  }
	}
      }
      if ( l->source == (WBFMM_SOURCE_MONOPOLE | WBFMM_SOURCE_DIPOLE) ) {
	tstr = 2 ;
	for ( ib = 0 ; ib < wbfmm_target_list_point_number(l) ; ib ++ ) {
	  b = l->boxes[ib] ;
	  /*direct contributions from neighbour boxes*/
	  csrc = &(((WBFMM_REAL *)(l->csrc))[tstr*(l->ics[ib])]) ;
	  for ( j = 0 ; j < l->ibox[b+1]-l->ibox[b] ; j ++ ) {
	    idx = l->isrc[l->ibox[b]+j] ;
	    for ( k = 0 ; k < nq ; k ++ ) {
	      f[ib*fstr+k] +=  src[idx*sstr+k]*csrc[j*tstr+0] ;
	      f[ib*fstr+k] += nsrc[idx*nstr+k]*csrc[j*tstr+1] ;
	    }
	  }
	}
      }
      break ;
    case WBFMM_FIELD_GRADIENT:
      nf = 3 ;
      g_assert_not_reached() ; /*unchecked code*/
      for ( ib = 0 ; ib < wbfmm_target_list_point_number(l) ; ib ++ ) {
      	b = l->boxes[ib] ;
      	/*direct contributions from neighbour boxes*/
      	csrc = &(((WBFMM_REAL *)(l->csrc))[3*(l->ics[ib])]) ;
      	for ( j = 0 ; j < l->ibox[b+1]-l->ibox[b] ; j ++ ) {
      	  idx = l->isrc[l->ibox[b]+j] ;
      	  for ( k = 0 ; k < nq ; k ++ ) {
      	    f[ib*fstr+3*k+0] += src[idx*sstr+k]*csrc[3*j+0] ;
      	    f[ib*fstr+3*k+1] += src[idx*sstr+k]*csrc[3*j+1] ;
      	    f[ib*fstr+3*k+2] += src[idx*sstr+k]*csrc[3*j+2] ;
      	  }
      	}
      }
    case WBFMM_FIELD_CURL:
      nf = 3 ;
      g_assert_not_reached() ; /*unchecked code*/
      for ( ib = 0 ; ib < wbfmm_target_list_point_number(l) ; ib ++ ) {
      	b = l->boxes[ib] ;
      	/*direct contributions from neighbour boxes*/
      	csrc = &(((WBFMM_REAL *)(l->csrc))[3*(l->ics[ib])]) ;
      	for ( j = 0 ; j < l->ibox[b+1]-l->ibox[b] ; j ++ ) {
      	  idx = l->isrc[l->ibox[b]+j] ;
      	  for ( k = 0 ; k < nq ; k ++ ) {
      	    f[ib*fstr+3*k+0] += src[idx*sstr+k]*csrc[3*j+0] ;
      	    f[ib*fstr+3*k+1] += src[idx*sstr+k]*csrc[3*j+1] ;
      	    f[ib*fstr+3*k+2] += src[idx*sstr+k]*csrc[3*j+2] ;
      	  }
      	}
      }      
      break ;
    }

    for ( ib = 0 ; ib < wbfmm_target_list_point_number(l) ; ib ++ ) {
      b = l->boxes[ib] ; cfft = boxes[b].mpr ;
      WBFMM_FUNCTION_NAME(wbfmm_laplace_expansion_apply)(cfft,
							 8*nq, nq,
							 &(eval[ib*nc]), nr,
							 l->field,
							 &(f[ib*fstr]), nf) ;
    }
    
    return 0 ;
  }
  
  if ( l->t->problem == WBFMM_PROBLEM_HELMHOLTZ ) {
    if ( !( l->field  == WBFMM_FIELD_SCALAR ||
	    l->field  == WBFMM_FIELD_GRADIENT ) ) {
      g_error("%s: field type %d not implemented for Helmholtz problem",
	      __FUNCTION__, l->field) ;
    }
    
    switch ( l->field ) {
    case WBFMM_FIELD_SCALAR:
      nf = 2 ;
      for ( ib = 0 ; ib < wbfmm_target_list_point_number(l) ; ib ++ ) {
	b = l->boxes[ib] ;
	/*direct contributions from neighbour boxes*/
	csrc = &(((WBFMM_REAL *)(l->csrc))[2*(l->ics[ib])]) ;
	for ( j = 0 ; j < l->ibox[b+1]-l->ibox[b] ; j ++ ) {
	  idx = l->isrc[l->ibox[b]+j] ;
	  for ( k = 0 ; k < nq ; k ++ ) {
	    WBFMM_REAL sr, si ;
	    sr = src[idx*sstr+2*k+0] ; si = src[idx*sstr+2*k+1] ;
	    f[ib*fstr+2*k+0] += csrc[2*j+0]*sr - csrc[2*j+1]*si ;
	    f[ib*fstr+2*k+1] += csrc[2*j+0]*si + csrc[2*j+1]*sr ;
	  }
	}    
      }
      break ;
    case WBFMM_FIELD_GRADIENT:
      nf = 6 ;
      g_assert_not_reached() ; /*unchecked code*/
      for ( ib = 0 ; ib < wbfmm_target_list_point_number(l) ; ib ++ ) {
      	b = l->boxes[ib] ;
      	/*direct contributions from neighbour boxes*/
      	csrc = &(((WBFMM_REAL *)(l->csrc))[6*(l->ics[ib])]) ;
      	for ( j = 0 ; j < l->ibox[b+1]-l->ibox[b] ; j ++ ) {
      	  idx = l->isrc[l->ibox[b]+j] ;
      	  for ( k = 0 ; k < nq ; k ++ ) {
      	    WBFMM_REAL sr, si ;
      	    sr = src[idx*sstr+2*k+0] ; si = src[idx*sstr+2*k+1] ;
      	    f[ib*fstr+6*k+0] += csrc[6*j+0]*sr - csrc[6*j+1]*si ;
      	    f[ib*fstr+6*k+1] += csrc[6*j+1]*sr + csrc[6*j+0]*si ;

      	    f[ib*fstr+6*k+2] += csrc[6*j+2]*sr - csrc[6*j+3]*si ;
      	    f[ib*fstr+6*k+3] += csrc[6*j+3]*sr + csrc[6*j+2]*si ;

      	    f[ib*fstr+6*k+4] += csrc[6*j+4]*sr - csrc[6*j+5]*si ;
      	    f[ib*fstr+6*k+5] += csrc[6*j+5]*sr + csrc[6*j+4]*si ;
      	  }
      	}
      }
      break ;
    }
    
    for ( ib = 0 ; ib < wbfmm_target_list_point_number(l) ; ib ++ ) {
      b = l->boxes[ib] ;
      cfft = boxes[b].mpr ;
      WBFMM_FUNCTION_NAME(wbfmm_expansion_apply)(cfft,
						 8*nq, nq,
						 &(eval[ib*nc]), nr,
						 l->field,
						 &(f[ib*fstr]), nf) ;
    }

    return 0 ;
  }

  g_error("%s: problem type %d not recognized",
	  __FUNCTION__, l->t->problem) ;
  
  return 0 ;
}
