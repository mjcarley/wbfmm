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

#ifdef _HAVE_CONFIG_H_
#include <config.h>
#endif /*_HAVE_CONFIG_H_*/

#include <math.h>
#include <glib.h>
#include <string.h>
#include <stdio.h>
#include <inttypes.h>

#include <wbfmm.h>

#include "wbfmm-private.h"

#define _DATA_WIDTH       8
#define _DATA_TREE        0
/* #define _DATA_ */

wbfmm_tree_t *FUNCTION_NAME(wbfmm_tree_new)(WBFMM_REAL *x, WBFMM_REAL D,
					     guint maxpoints)

{
  wbfmm_tree_t *t ;
  gint i ;
  WBFMM_REAL *xt ;

  t = (wbfmm_tree_t *)g_malloc(sizeof(wbfmm_tree_t)) ;

  /*maximum number of points in tree*/
  t->maxpoints = maxpoints ;
  /*array for sorted point indices*/
  t->ip = (guint32 *)g_malloc(maxpoints*sizeof(guint32)) ;
  t->npoints = 0 ;

  for ( i = 0 ; i <= WBFMM_TREE_MAX_DEPTH ; i ++ ) t->boxes[i] = NULL ;
  t->boxes[0] = (wbfmm_box_t *)g_malloc(1*sizeof(wbfmm_box_t)) ;

  t->depth = 0 ;

  xt = wbfmm_tree_origin(t) ;
  xt[0] = x[0] ; xt[1] = x[1] ; xt[2] = x[2] ;
  t->D = D ;

  return t ;
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
  mi = FUNCTION_NAME(wbfmm_point_index_3d)(xi, wbfmm_tree_origin(t), 
					   wbfmm_tree_width(t)) ;
  mj = FUNCTION_NAME(wbfmm_point_index_3d)(xj, wbfmm_tree_origin(t), 
					   wbfmm_tree_width(t)) ;

  if ( mi < mj ) return -1 ;
  if ( mi > mj ) return  1 ;

  return 0 ;
}

gint FUNCTION_NAME(wbfmm_tree_add_points)(wbfmm_tree_t *t, 
					  gpointer pts, guint npts, 
					  gsize pstr)

{
  guint i ;
  WBFMM_REAL *x, *xt, D ;

  if ( npts > wbfmm_tree_point_number_max(t) ) 
    g_error("%s: too many points (%u) for tree (%u)",
	    __FUNCTION__, npts, wbfmm_tree_point_number_max(t)) ;

  wbfmm_tree_point_number(t) = npts ;
  t->points = (gchar *)pts ; t->pstr = pstr ;

  xt = wbfmm_tree_origin(t) ;
  D = wbfmm_tree_width(t) ;

  for ( i = 0 ; i < npts ; i ++ ) {
    x = wbfmm_tree_point_index(t, i) ;
    if ( x[0] <= xt[0] || x[0] >= xt[0] + D ||
	 x[1] <= xt[1] || x[1] >= xt[1] + D ||
	 x[2] <= xt[2] || x[2] >= xt[2] + D )
      g_error("%s: point (%g,%g,%g) does not lie in box of width %g at "
	      "(%g,%g,%g)",
	      __FUNCTION__, x[0], x[1], x[2], D, xt[0], xt[1], xt[2]) ;
	 
    t->ip[i] = i ;
  }
  

  /*sort points on the Morton index*/
  g_qsort_with_data(t->ip, npts, sizeof(guint), compare_morton_indexed, 
		    (gpointer)t) ;

  t->boxes[0][0].i = 0 ; 
  t->boxes[0][0].n = npts ; 

  return 0 ;
}


gint FUNCTION_NAME(wbfmm_tree_refine)(wbfmm_tree_t *t)

{
  guint level = wbfmm_tree_depth(t) ;
  wbfmm_box_t *parents, *children ;
  guint np, j ;
  guint64 idx, child, xi, box ;
  WBFMM_REAL *x ;

  wbfmm_tree_add_level(t) ;

  /*number of parent boxes to refine*/
  np = 1 << 3*(level) ;

  parents  = t->boxes[level  ] ;
  children = t->boxes[level+1] ;

  /*this could probably be done with binary searches*/
  for ( idx = 0 ; idx < np ; idx ++ ) {
    /*initialize the first child box*/
    child = wbfmm_box_first_child(idx) ;
    children[child].i = parents[idx].i ;
    children[child].n = 0 ;
    /*start at first parent index*/
    j = parents[idx].i ;
    while ( parents[idx].n != 0 ) {
      /*check if current point is in box*/
      x = wbfmm_tree_point_index(t, t->ip[j]) ;
      xi = FUNCTION_NAME(wbfmm_point_index_3d)(x, wbfmm_tree_origin(t), 
					       wbfmm_tree_width(t)) ;
      box = wbfmm_point_locate_box(xi, level+1) ;
      /* g_assert(box >= child && box < child+8) ; */
      if ( box == child ) {
	parents[idx].n -- ;
	parents[idx].i ++ ;
	children[child].n ++ ;
	j ++ ;
      } else {
	child ++ ;
	children[child].n = 0 ; 
	children[child].i = parents[idx].i ;
      }
    }
    g_assert(parents[idx].n == 0) ;
  }

  return 0 ;
}


static inline guint64 morton_encode(guint32 xsrc, guint32 ysrc, guint32 zsrc)

{
    const guint64  
      mask0 = 
      BITMASK_0000000001000001000001000001000001000001000001000001000001000001,
      mask1 =
      BITMASK_0000001000001000001000001000001000001000001000001000001000001000,
      mask2 =
      BITMASK_0001000000000000000000000000000000000000000000000000000000000000,
      mask3 =
      BITMASK_0000000000000011000000000011000000000011000000000011000000000011,
      mask4 =
      BITMASK_0000000111000000000011000000000011000000000011000000000011000000,
      mask5 =
      BITMASK_0000000000000000000000000000000000001111000000000000000000001111,
      mask6 =
      BITMASK_0000000000000000000000001111000000000000000000001111000000000000,
      mask7 =
      BITMASK_0000000000011111000000000000000000000000000000000000000000000000,
      mask8 =
      BITMASK_0000000000000000000000000000000000000000000000000000000011111111,
      mask9 =
      BITMASK_0000000000000000000000000001111111111111000000000000000000000000;
    guint64  x = xsrc, y = ysrc, z = zsrc ;
    /* 0000000000000000000000000000000000000000000ccccccccccccccccccccc */
    x = (x & mask8) | ((x << 16) & mask9) ;
    y = (y & mask8) | ((y << 16) & mask9) ;
    z = (z & mask8) | ((z << 16) & mask9) ;
    /* 000000000000000000000000000ccccccccccccc0000000000000000cccccccc */
    x = (x & mask5) | ((x << 8) & mask6) | ((x << 16) & mask7) ;
    y = (y & mask5) | ((y << 8) & mask6) | ((y << 16) & mask7) ;
    z = (z & mask5) | ((z << 8) & mask6) | ((z << 16) & mask7) ;
    /* 00000000000ccccc00000000cccc00000000cccc00000000cccc00000000cccc */
    x = (x & mask3) | ((x << 4) & mask4) ;
    y = (y & mask3) | ((y << 4) & mask4) ;
    z = (z & mask3) | ((z << 4) & mask4) ;
    /* 0000000ccc0000cc0000cc0000cc0000cc0000cc0000cc0000cc0000cc0000cc */
    x = (x & mask0) | ((x << 2) & mask1) | ((x << 4) & mask2) ;
    y = (y & mask0) | ((y << 2) & mask1) | ((y << 4) & mask2) ;
    z = (z & mask0) | ((z << 2) & mask1) | ((z << 4) & mask2) ;
    /* 000c00c00c00c00c00c00c00c00c00c00c00c00c00c00c00c00c00c00c00c00c */
    return x | (y << 1) | (z << 2) ;
}

guint64 FUNCTION_NAME(wbfmm_point_index_3d)(WBFMM_REAL *x, WBFMM_REAL *c,
					    WBFMM_REAL D)

{
  guint32 xi, yi, zi ;
  guint64 i ;

  xi = (guint32)((x[0] - c[0])/D*WBFMM_INDEX_SHIFT) ;
  yi = (guint32)((x[1] - c[1])/D*WBFMM_INDEX_SHIFT) ;
  zi = (guint32)((x[2] - c[2])/D*WBFMM_INDEX_SHIFT) ;

  /*interleave*/
  i = morton_encode(xi, yi, zi) ;

  return i ;
}

gint FUNCTION_NAME(wbfmm_tree_leaf_expansions)(wbfmm_tree_t *t, WBFMM_REAL k,
					       WBFMM_REAL *src, gint sstr,
					       WBFMM_REAL *normals, gint nstr,
					       WBFMM_REAL *dipoles, gint dstr,
					       gboolean zero_expansions,
					       WBFMM_REAL *work)

{
  guint32 nb, nc, i, j, ns, d, idx ;
  guint64 im ;
  wbfmm_box_t *boxes ;
  WBFMM_REAL *xs, *q, *fd, *n, xb[3], wb ;

  /*depth of leaves*/
  d = wbfmm_tree_depth(t) ;
  /*order of singular expansions*/
  ns = t->order_s[d] ;
  /*number of boxes*/
  nb = 1 << (3*d) ;
  /*number of coefficients*/
  nc = wbfmm_coefficient_index_nm(ns+1,0) ;

  /*zero the coefficients before accumulating*/
  if ( zero_expansions )
    memset(t->mps[d], 0, nb*nc*2*sizeof(WBFMM_REAL)) ;

  boxes = t->boxes[d] ;

  if ( src == NULL && normals == NULL && dipoles == NULL ) return 0 ;

  if ( normals != NULL && dipoles == NULL ) {
    g_error("%s: normals specified but no dipole strengths (dipoles == NULL)",
	    __FUNCTION__) ;
  }

  if ( normals == NULL && dipoles == NULL ) {
    /* monopoles only */
    for ( i = 0 ; i < nb ; i ++ ) {
      im = (guint64)i ;
      FUNCTION_NAME(wbfmm_box_location_from_index)(im, d, 
						   wbfmm_tree_origin(t), 
						   wbfmm_tree_width(t), xb, 
						   &wb) ;
      xb[0] += 0.5*wb ; xb[1] += 0.5*wb ; xb[2] += 0.5*wb ; 

      for ( j = 0 ; j < boxes[i].n ; j ++ ) {
	idx = t->ip[boxes[i].i+j] ;
	xs = wbfmm_tree_point_index(t,idx) ;
	q = &(src[idx*sstr]) ;
	FUNCTION_NAME(wbfmm_expansion_h_cfft)(k, ns, xb, xs, q,
					      boxes[i].mps, 8, work) ;
      }
    }

    return 0 ;
  }

  if ( src == NULL && normals != NULL ) {
    /*dipoles only, specified as normals and strengths*/
    for ( i = 0 ; i < nb ; i ++ ) {
      im = (guint64)i ;
      FUNCTION_NAME(wbfmm_box_location_from_index)(im, d, 
						   wbfmm_tree_origin(t), 
						   wbfmm_tree_width(t), xb, 
						   &wb) ;
      xb[0] += 0.5*wb ; xb[1] += 0.5*wb ; xb[2] += 0.5*wb ; 

      for ( j = 0 ; j < boxes[i].n ; j ++ ) {
	idx = t->ip[boxes[i].i+j] ;
	xs = wbfmm_tree_point_index(t,idx) ;
	fd = &(dipoles[idx*dstr]) ;
	n  = &(normals[idx*nstr]) ;
	FUNCTION_NAME(wbfmm_expansion_normal_h_cfft)(k, ns, xb, xs, n, fd,
						     boxes[i].mps, 8, work) ;
      }
    }

    return 0 ;
  }

  if ( src != NULL && normals != NULL ) {
    /*mixed sources, dipoles specified as normals and strengths*/
    for ( i = 0 ; i < nb ; i ++ ) {
      im = (guint64)i ;
      FUNCTION_NAME(wbfmm_box_location_from_index)(im, d, 
						   wbfmm_tree_origin(t), 
						   wbfmm_tree_width(t), xb, 
						   &wb) ;
      xb[0] += 0.5*wb ; xb[1] += 0.5*wb ; xb[2] += 0.5*wb ; 

      for ( j = 0 ; j < boxes[i].n ; j ++ ) {
	idx = t->ip[boxes[i].i+j] ;
	xs = wbfmm_tree_point_index(t,idx) ;
	q = &(src[idx*sstr]) ;
	fd = &(dipoles[idx*dstr]) ;
	n  = &(normals[idx*nstr]) ;
	FUNCTION_NAME(wbfmm_expansion_normal_h_cfft)(k, ns, xb, xs, n, fd,
						     boxes[i].mps, 8, work) ;
	FUNCTION_NAME(wbfmm_expansion_h_cfft)(k, ns, xb, xs, q,
					      boxes[i].mps, 8, work) ;
      }
    }

    return 0 ;
  }

  g_assert_not_reached() ;
    
  return 0 ;
}


gint FUNCTION_NAME(wbfmm_tree_box_field)(wbfmm_tree_t *t, guint level,
					 guint b, WBFMM_REAL k,
					 WBFMM_REAL *x, WBFMM_REAL *f,
					 WBFMM_REAL *work)

{
  WBFMM_REAL xb[3], wb, *C ;
  wbfmm_box_t *boxes ;

  boxes = t->boxes[level] ;
  C = boxes[b].mps ;

  FUNCTION_NAME(wbfmm_tree_box_centre)(t, level, b, xb, &wb) ;

  FUNCTION_NAME(wbfmm_expansion_h_evaluate)(k, xb, C, 8,
					    t->order_s[level], x, f, work) ;

  return 0 ;
}

gint FUNCTION_NAME(wbfmm_tree_box_local_field)(wbfmm_tree_t *t, guint level,
					       guint b, WBFMM_REAL k,
					       WBFMM_REAL *x, WBFMM_REAL *f,
					       WBFMM_REAL *src, gint sstr,
					       WBFMM_REAL *normals, gint nstr,
					       WBFMM_REAL *d, gint dstr,
					       gboolean eval_neighbours,
					       WBFMM_REAL *work)

{
  WBFMM_REAL xb[3], wb, *C, *xs, r, h0[2], h1[2], fR[2] ;
  wbfmm_box_t *boxes, box ;
  guint64 neighbours[27] ;
  gint nnbr, i, j, idx ;

  boxes = t->boxes[level] ;
  C = boxes[b].mpr ;

  FUNCTION_NAME(wbfmm_tree_box_centre)(t, level, b, xb, &wb) ;

  FUNCTION_NAME(wbfmm_expansion_j_evaluate)(k, xb, C, 8,
					    t->order_r[level], x, f, work) ;

  if ( !eval_neighbours ) return 0 ;

  if ( src == NULL && normals == NULL && d == NULL ) return 0 ;
  
  if ( normals != NULL && d == NULL ) {
    g_error("%s: normals specified but no dipole strengths (d == NULL)",
	    __FUNCTION__) ;
  }

  /*add the contribution from sources in neighbour boxes*/
  nnbr = wbfmm_box_neighbours(level, b, neighbours) ;
  g_assert(nnbr >= 0 && nnbr < 28) ;

  if ( normals == NULL && d == NULL ) {
    /* monopoles only */
    for ( i = 0 ; i < nnbr ; i ++ ) {
      box = boxes[neighbours[i]] ;
      for ( j = 0 ; j < box.n ; j ++ ) {
	idx = t->ip[box.i+j] ;
	xs = wbfmm_tree_point_index(t, idx) ;
	r = (xs[0]-x[0])*(xs[0]-x[0]) + (xs[1]-x[1])*(xs[1]-x[1]) +
	  (xs[2]-x[2])*(xs[2]-x[2]) ;
	if ( r > 1e-12 ) {
	  r = SQRT(r) ;
	  FUNCTION_NAME(wbfmm_bessel_h_init)(k*r, h0, h1) ;
	  h0[0] /= 4.0*M_PI ; h0[1] /= 4.0*M_PI ; 
	  f[0] += h0[0]*src[idx*sstr+0] - h0[1]*src[idx*sstr+1] ;
	  f[1] += h0[1]*src[idx*sstr+0] + h0[0]*src[idx*sstr+1] ;
	}
      }
    }

    return 0 ;
  }

  /* return 0 ; */
  if ( src == NULL && normals != NULL ) {
    /* dipoles only, specified as normals and strengths */
    for ( i = 0 ; i < nnbr ; i ++ ) {
      box = boxes[neighbours[i]] ;
      for ( j = 0 ; j < box.n ; j ++ ) {
	idx = t->ip[box.i+j] ;
	xs = wbfmm_tree_point_index(t, idx) ;
	r = (xs[0]-x[0])*(xs[0]-x[0]) + (xs[1]-x[1])*(xs[1]-x[1]) +
	  (xs[2]-x[2])*(xs[2]-x[2]) ;
	if ( r > 1e-12 ) {
	  r = SQRT(r) ;
	  FUNCTION_NAME(wbfmm_bessel_h_init)(k*r, h0, h1) ;
	  h1[0] /= 4.0*M_PI ; h1[1] /= 4.0*M_PI ;

	  fR[0] = (normals[idx*nstr+0]*(x[0] - xs[0]) +
		   normals[idx*nstr+1]*(x[1] - xs[1]) +
		   normals[idx*nstr+2]*(x[2] - xs[2]))/r ;
	  fR[1] = d[idx*dstr+1]*fR[0] ; fR[0] *= d[idx*dstr+0] ;

	  f[0] -= k*(h1[0]*fR[0] - h1[1]*fR[1]) ;
	  f[1] -= k*(h1[0]*fR[1] + h1[1]*fR[0]) ;
	}
      }
    }
    
    return 0 ;
  }

  if ( src != NULL && normals != NULL ) {
    /* mixed monopoles and dipoles specified as normals and strengths */
    for ( i = 0 ; i < nnbr ; i ++ ) {
      box = boxes[neighbours[i]] ;
      for ( j = 0 ; j < box.n ; j ++ ) {
	idx = t->ip[box.i+j] ;
	xs = wbfmm_tree_point_index(t, idx) ;
	r = (xs[0]-x[0])*(xs[0]-x[0]) + (xs[1]-x[1])*(xs[1]-x[1]) +
	  (xs[2]-x[2])*(xs[2]-x[2]) ;
	if ( r > 1e-12 ) {
	  r = SQRT(r) ;
	  FUNCTION_NAME(wbfmm_bessel_h_init)(k*r, h0, h1) ;
	  h0[0] /= 4.0*M_PI ; h0[1] /= 4.0*M_PI ;
	  h1[0] /= 4.0*M_PI ; h1[1] /= 4.0*M_PI ;

	  fR[0] = (normals[idx*nstr+0]*(x[0] - xs[0]) +
		   normals[idx*nstr+1]*(x[1] - xs[1]) +
		   normals[idx*nstr+2]*(x[2] - xs[2]))/r ;
	  fR[1] = d[idx*dstr+1]*fR[0] ; fR[0] *= d[idx*dstr+0] ;
	  
	  f[0] += h0[0]*src[idx*sstr+0] - h0[1]*src[idx*sstr+1] ;
	  f[1] += h0[1]*src[idx*sstr+0] + h0[0]*src[idx*sstr+1] ;

	  f[0] -= k*(h1[0]*fR[0] - h1[1]*fR[1]) ;
	  f[1] -= k*(h1[0]*fR[1] + h1[1]*fR[0]) ;
	}
      }
    }
    
    return 0 ;
  }
  
  g_assert_not_reached() ; 
  
  return 0 ;
}

gint FUNCTION_NAME(wbfmm_tree_coefficient_init)(wbfmm_tree_t *t,
						guint l, 
						guint nr, guint ns)

{
  guint nb, nc, i, j ;
  wbfmm_box_t *boxes ;
  WBFMM_REAL *c ;

  g_assert(l <= wbfmm_tree_depth(t)) ;

  /*number of boxes at level l*/
  nb = 1 << (3*l) ;

  t->mps[l] = t->mpr[l] = NULL ;
  t->order_s[l] = ns ; t->order_r[l] = nr ; 

  /*number of coefficients in singular expansions*/
  if ( ns != 0 ) {
    nc = wbfmm_coefficient_index_nm(ns+1,0) ;
    /* nc = wbfmm_coefficient_index_nm(ns+2,0) ; */
    t->mps[l] = g_malloc0(nb*2*nc*sizeof(WBFMM_REAL)) ;
    c = (WBFMM_REAL *)(t->mps[l]) ;
    /*
      set box pointers to start of their coefficients, noting that
      coefficients are packed in groups of eight for shift operations
    */
    boxes = t->boxes[l] ;
    for ( i = 0 ; i < nb ; i += 8 ) {
      for ( j = 0 ; j < 8 ; j ++ ) {
	boxes[i+j].mps = &(c[i*2*nc+2*j]) ;
      }
    }
  }

  if ( nr != 0 ) {
    nc = wbfmm_coefficient_index_nm(nr+1,0) ;
    /* nc = wbfmm_coefficient_index_nm(nr+2,0) ; */
    t->mpr[l] = g_malloc0(nb*2*nc*sizeof(WBFMM_REAL)) ;
    c = (WBFMM_REAL *)(t->mpr[l]) ;
    /*
      set box pointers to start of their coefficients, noting that
      coefficients are packed in groups of eight for shift operations
    */
    boxes = t->boxes[l] ;
    for ( i = 0 ; i < nb ; i += 8 ) {
      for ( j = 0 ; j < 8 ; j ++ ) {
	boxes[i+j].mpr = &(c[i*2*nc+2*j]) ;
      }
    }
  }

  return 0 ;
}

guint64 FUNCTION_NAME(wbfmm_point_box)(wbfmm_tree_t *t, guint level,
				       WBFMM_REAL *x)

{
  guint64 b ;
  guint nb, i, j, k ;
  WBFMM_REAL wb, *x0, dx ;

  nb = 1 << level ;
  wb = t->D/nb ;

  x0 = wbfmm_tree_origin(t) ;
  dx = (x[0] - x0[0])/wb ;
  if ( dx < 0.0 || dx > nb )
    g_error("%s: point (%g,%g,%g) not in octree",
	    __FUNCTION__, x[0], x[1], x[2]) ;
  else
    i = (guint32)floor(dx) ;

  dx = (x[1] - x0[1])/wb ;
  if ( dx < 0.0 || dx > nb )
    g_error("%s: point (%g,%g,%g) not in octree",
	    __FUNCTION__, x[0], x[1], x[2]) ;
  else
    j = (guint32)floor(dx) ;

  dx = (x[2] - x0[2])/wb ;
  if ( dx < 0.0 || dx > nb )
    g_error("%s: point (%g,%g,%g) not in octree",
	    __FUNCTION__, x[0], x[1], x[2]) ;
  else
    k = (guint32)floor(dx) ;

  /* i = (guint32)floor((x[0] - x0[0])/wb) ; */
  /* j = (guint32)floor((x[1] - x0[1])/wb) ; */
  /* k = (guint32)floor((x[2] - x0[2])/wb) ; */

  b = wbfmm_box_index(i, j, k) ;

  return b ;
}
