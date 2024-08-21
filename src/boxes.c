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
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#include <glib.h>

#include <wbfmm.h>

#include "wbfmm-private.h"

#define _DATA_WIDTH       8
#define _DATA_TREE        0
/* #define _DATA_ */

/* #define WBFMM_CHECK_ISNAN */

#ifdef WBFMM_CHECK_ISNAN
#include <stdlib.h>

static gint check_isnan(char *name, WBFMM_REAL *f, gint n)

{
  gint i ;

  for ( i = 0 ; i < n ; i ++ ) {
    if ( isnan(f[i]) ) {
      fprintf(stderr, "%s: NaN at element %d of %d\n", name, i, n) ;
      exit(1) ;
    }
  }
  
  return 0 ;
}
#endif /*WBFMM_CHECK_ISNAN*/

wbfmm_tree_t *WBFMM_FUNCTION_NAME(wbfmm_tree_new)(WBFMM_REAL *x, WBFMM_REAL D,
						  guint maxpoints)

{
  wbfmm_tree_t *t ;
  gint i ;
  WBFMM_REAL *xt ;

  t = (wbfmm_tree_t *)g_malloc0(sizeof(wbfmm_tree_t)) ;

  t->problem = 0 ; t->sorted = FALSE ;
  /*maximum number of points in tree*/
  t->maxpoints = maxpoints ;
  /*number of components in tree sources (set when coefficients are
    initialized)*/
  wbfmm_tree_source_size(t) = 0 ;
  /*array for sorted point indices*/
  t->ip = (guint32 *)g_malloc(maxpoints*sizeof(guint32)) ;
  t->npoints = 0 ;
  t->points = t->normals = NULL ;
  
  t->size = sizeof(WBFMM_REAL) ;
  
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
  mi = WBFMM_FUNCTION_NAME(wbfmm_point_index_3d)(xi, wbfmm_tree_origin(t), 
					   wbfmm_tree_width(t)) ;
  mj = WBFMM_FUNCTION_NAME(wbfmm_point_index_3d)(xj, wbfmm_tree_origin(t), 
					   wbfmm_tree_width(t)) ;

  if ( mi < mj ) return -1 ;
  if ( mi > mj ) return  1 ;

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_tree_add_points)(wbfmm_tree_t *t, 
						gpointer pts, gsize pstr,
						gpointer normals, gsize nstr,
						guint npts, gboolean sorted)

{
  gint i ;
  WBFMM_REAL *x, *xt, D ;

  if ( t->size != sizeof(WBFMM_REAL) )
    g_error("%s: mixed precision not implemented\n"
	    "  (size of tree data type (%lu) not equal to "
	    "size of requested target type (%lu))",
	    __FUNCTION__, t->size, sizeof(WBFMM_REAL)) ;

  if ( npts > wbfmm_tree_point_number_max(t) ) 
    g_error("%s: too many points (%u) for tree (%u)",
	    __FUNCTION__, npts, wbfmm_tree_point_number_max(t)) ;

  wbfmm_tree_point_number(t) = npts ;
  t->points  = (char *)pts ;     t->pstr = pstr ;
  t->normals = (char *)normals ; t->nstr = nstr ;

  t->sorted = sorted ;
  
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

static gint compare_points_morton(gconstpointer a, gconstpointer b,
				  gpointer data)

{
  wbfmm_tree_t *t = data ;
  guint64 mi, mj ;
  WBFMM_REAL *xi, *xj ;
  
  xi = (WBFMM_REAL *)a ; xj = (WBFMM_REAL *)b ;

  /*Morton codes*/
  mi = WBFMM_FUNCTION_NAME(wbfmm_point_index_3d)(xi, wbfmm_tree_origin(t), 
						 wbfmm_tree_width(t)) ;
  mj = WBFMM_FUNCTION_NAME(wbfmm_point_index_3d)(xj, wbfmm_tree_origin(t), 
						 wbfmm_tree_width(t)) ;

  if ( mi < mj ) return -1 ;
  if ( mi > mj ) return  1 ;

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_tree_sort_points)(wbfmm_tree_t *t, 
						 gpointer pts, gsize psize,
						 guint npts)

{
  if ( t->size != sizeof(WBFMM_REAL) )
    g_error("%s: mixed precision not implemented\n"
	    "  (size of tree data type (%lu) not equal to "
	    "size of requested target type (%lu))",
	    __FUNCTION__, t->size, sizeof(WBFMM_REAL)) ;

  /*sort points on the Morton index*/
  g_qsort_with_data((gpointer)pts, npts, psize, compare_points_morton,
		    (gpointer)t) ;

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_tree_refine)(wbfmm_tree_t *t)

{
  guint level = wbfmm_tree_depth(t) ;
  wbfmm_box_t *parents, *children ;
  guint np, nc, j, n ;
  guint64 idx, child, xi, box ;
  WBFMM_REAL *x ;

  /* g_assert(t->problem != 0) ; */
  wbfmm_tree_add_level(t) ;

  /*number of parent boxes to refine*/
  np = 1 << 3*(level) ;

  parents  = t->boxes[level  ] ;
  children = t->boxes[level+1] ;

  /*number of child boxes*/
  nc = 1 << 3*(level + 1) ;
  for ( idx = 0 ; idx < nc ; idx ++ ) {
    children[idx].i = 0 ;
    children[idx].n = 0 ;
  }
  
  /*this could probably be done with binary searches*/
  for ( idx = 0 ; idx < np ; idx ++ ) {
    /*initialize the child boxes*/
    child = wbfmm_box_first_child(idx) ;
    children[child].i = parents[idx].i ;
    children[child].n = 0 ;
    /*start at first parent index*/
    j = parents[idx].i ;
    n = parents[idx].n ;
    /* while ( parents[idx].n != 0 ) { */
    while ( n != 0 ) {
      /*check if current point is in box*/
      x = wbfmm_tree_point_index(t, t->ip[j]) ;
      xi = WBFMM_FUNCTION_NAME(wbfmm_point_index_3d)(x, wbfmm_tree_origin(t), 
					       wbfmm_tree_width(t)) ;
      box = wbfmm_point_locate_box(xi, level+1) ;
      /* g_assert(box >= child && box < child+8) ; */
      /* if ( box < child || box > child + 7 ) { */
      /* 	g_error("%s: box allocation error\n" */
      /* 		"  x            = %lg %lg %lg\n" */
      /* 		"  level        = %u\n" */
      /* 		"  parent index = %lu\n" */
      /* 		"  child index  = %lu\n" */
      /* 		"  box index    = %lu\n", */
      /* 		__FUNCTION__, x[0], x[1], x[2], level, idx, */
      /* 		child, box) ; */
      /* } */
      if ( box == child ) {
	/* parents[idx].n -- ; */
	n -- ;
	parents[idx].i ++ ;
	children[child].n ++ ;
	j ++ ;
      } else {
	child ++ ;
	children[child].n = 0 ; 
	children[child].i = parents[idx].i ;
      }
    }
    /* g_assert(parents[idx].n == 0) ; */
    g_assert(n == 0) ;
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

guint64 WBFMM_FUNCTION_NAME(wbfmm_point_index_3d)(WBFMM_REAL *x, WBFMM_REAL *c,
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

gint WBFMM_FUNCTION_NAME(wbfmm_tree_leaf_expansions)(wbfmm_tree_t *t,
						     WBFMM_REAL k,
						     WBFMM_REAL *src,
						     gint sstr,
						     WBFMM_REAL *dipoles,
						     gint dstr,
						     gboolean zero_expansions,
						     WBFMM_REAL *work)

{
  guint32 nb, nc, i, j, ns, d, idx ;
  guint64 im ;
  wbfmm_box_t *boxes ;
  WBFMM_REAL *xs, *q, *fd, *n, xb[3], wb ;
  gint nq = wbfmm_tree_source_size(t) ;
  
  g_assert(t->problem == WBFMM_PROBLEM_HELMHOLTZ ) ;

  /*depth of leaves*/
  d = wbfmm_tree_depth(t) ;
  /*order of singular expansions*/
  ns = t->order_s[d] ;
  /*number of boxes*/
  nb = 1 << (3*d) ;
  /*number of coefficients*/
  nc = wbfmm_coefficient_number(ns) ;
  
  /*zero the coefficients before accumulating*/
  if ( zero_expansions )
    memset(t->mps[d], 0, nb*nc*2*nq*sizeof(WBFMM_REAL)) ;

  boxes = t->boxes[d] ;

  if ( src == NULL && dipoles == NULL ) return 0 ;

  if ( t->normals == NULL && dipoles != NULL ) {
    g_error("%s: no normals in tree but dipole strengths specified "
	    "(dipoles != NULL)",
	    __FUNCTION__) ;
  }

  if ( dipoles == NULL ) {
    /* monopoles only */
    for ( i = 0 ; i < nb ; i ++ ) {
      im = (guint64)i ;
      WBFMM_FUNCTION_NAME(wbfmm_box_location_from_index)(im, d, 
							 wbfmm_tree_origin(t), 
							 wbfmm_tree_width(t),
							 xb, &wb) ;
      xb[0] += 0.5*wb ; xb[1] += 0.5*wb ; xb[2] += 0.5*wb ; 

      for ( j = 0 ; j < boxes[i].n ; j ++ ) {
	idx = t->ip[boxes[i].i+j] ;
	xs = wbfmm_tree_point_index(t,idx) ;
	q = &(src[idx*sstr]) ;
	WBFMM_FUNCTION_NAME(wbfmm_expansion_h_cfft)(k, ns, xb, xs, q, nq,
					      boxes[i].mps, 8*nq, work) ;
      }
    }

#ifdef WBFMM_CHECK_ISNAN
    check_isnan("leaf coefficients", t->mps[d], nb*nc*2*nq) ;
#endif /*WBFMM_CHECK_ISNAN*/
    
    return 0 ;
  }

  if ( src != NULL && dipoles != NULL ) {
    /*mixed sources, dipoles specified as normals and strengths*/
    for ( i = 0 ; i < nb ; i ++ ) {
      im = (guint64)i ;
      WBFMM_FUNCTION_NAME(wbfmm_box_location_from_index)(im, d, 
						   wbfmm_tree_origin(t), 
						   wbfmm_tree_width(t), xb, 
						   &wb) ;
      xb[0] += 0.5*wb ; xb[1] += 0.5*wb ; xb[2] += 0.5*wb ; 

      for ( j = 0 ; j < boxes[i].n ; j ++ ) {
	idx = t->ip[boxes[i].i+j] ;
	xs = wbfmm_tree_point_index(t,idx) ;
	q = &(src[idx*sstr]) ;
	fd = &(dipoles[idx*dstr]) ;
	n = wbfmm_tree_normal_index(t,idx) ;
	/* n  = &(normals[idx*nstr]) ; */
	WBFMM_FUNCTION_NAME(wbfmm_expansion_normal_h_cfft)(k, ns, xb, xs,
							   n, fd, nq,
							   boxes[i].mps, 8,
							   work) ;
	WBFMM_FUNCTION_NAME(wbfmm_expansion_h_cfft)(k, ns, xb, xs, q, nq,
						    boxes[i].mps, 8, work) ;
      }
    }

    return 0 ;
  }

  /* g_assert_not_reached() ; /\*following code needs modification*\/ */
  
  if ( src == NULL && dipoles != NULL ) {
    /*dipoles only, specified as normals and strengths*/
    for ( i = 0 ; i < nb ; i ++ ) {
      im = (guint64)i ;
      WBFMM_FUNCTION_NAME(wbfmm_box_location_from_index)(im, d, 
						   wbfmm_tree_origin(t), 
						   wbfmm_tree_width(t), xb, 
						   &wb) ;
      xb[0] += 0.5*wb ; xb[1] += 0.5*wb ; xb[2] += 0.5*wb ; 

      for ( j = 0 ; j < boxes[i].n ; j ++ ) {
	idx = t->ip[boxes[i].i+j] ;
	xs = wbfmm_tree_point_index(t,idx) ;
	/* q = &(src[idx*sstr]) ; */
	fd = &(dipoles[idx*dstr]) ;
	/* n  = &(normals[idx*nstr]) ; */
	n = wbfmm_tree_normal_index(t,idx) ;
	WBFMM_FUNCTION_NAME(wbfmm_expansion_normal_h_cfft)(k, ns, xb, xs,
							   n, fd, nq,
							   boxes[i].mps, 8,
							   work) ;
	/* WBFMM_FUNCTION_NAME(wbfmm_expansion_h_cfft)(k, ns, xb, xs, q, nq, */
	/* 					    boxes[i].mps, 8, work) ; */
      }
    }
    /* for ( i = 0 ; i < nb ; i ++ ) { */
    /*   im = (guint64)i ; */
    /*   WBFMM_FUNCTION_NAME(wbfmm_box_location_from_index)(im, d,  */
    /* 						   wbfmm_tree_origin(t),  */
    /* 						   wbfmm_tree_width(t), xb,  */
    /* 						   &wb) ; */
    /*   xb[0] += 0.5*wb ; xb[1] += 0.5*wb ; xb[2] += 0.5*wb ;  */

    /*   for ( j = 0 ; j < boxes[i].n ; j ++ ) { */
    /* 	idx = t->ip[boxes[i].i+j] ; */
    /* 	xs = wbfmm_tree_point_index(t,idx) ; */
    /* 	fd = &(dipoles[idx*dstr]) ; */
    /* 	n  = &(normals[idx*nstr]) ; */
    /* 	WBFMM_FUNCTION_NAME(wbfmm_expansion_normal_h_cfft)(k, ns, xb, xs, */
    /* 							   n, fd, nq, */
    /* 							   boxes[i].mps, */
    /* 							   8, work) ; */
    /*   } */
    /* } */

    return 0 ;
  }

  g_assert_not_reached() ;
    
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_tree_box_field)(wbfmm_tree_t *t, guint level,
					       guint b, WBFMM_REAL k,
					       WBFMM_REAL *x, WBFMM_REAL *f,
					       gint fstr,
					       WBFMM_REAL *work)

{
  WBFMM_REAL xb[3], wb, *C ;
  wbfmm_box_t *boxes ;

  g_assert(t->problem == WBFMM_PROBLEM_HELMHOLTZ ) ;

  boxes = t->boxes[level] ;
  C = boxes[b].mps ;

  WBFMM_FUNCTION_NAME(wbfmm_tree_box_centre)(t, level, b, xb, &wb) ;

  WBFMM_FUNCTION_NAME(wbfmm_expansion_h_evaluate)(k, xb, C, 8,	  
						  t->order_s[level],
						  wbfmm_tree_source_size(t),
						  x, f, fstr,
						  work) ;

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_tree_coefficient_init)(wbfmm_tree_t *t,
						      guint l, 
						      guint nr, guint ns)

{
  guint nb, nc, i, j ;
  gint nq ;
  wbfmm_box_t *boxes ;
  WBFMM_REAL *c ;

  nq = wbfmm_tree_source_size(t) ;
  if ( nq < 1 )
    g_error("%s: tree has invalid number of components in source terms (%d)",
	    __FUNCTION__, nq) ;
  
  g_assert(t->problem == WBFMM_PROBLEM_HELMHOLTZ ) ;
  /* g_assert(t->nq == 1) ; */
  
  g_assert(l <= wbfmm_tree_depth(t)) ;

  /*number of boxes at level l*/
  nb = 1 << (3*l) ;

  t->mps[l] = t->mpr[l] = NULL ;
  t->order_s[l] = ns ; t->order_r[l] = nr ; 

  /*number of coefficients in singular expansions*/
  if ( ns != 0 ) {
    /* nc = wbfmm_coefficient_index_nm(ns+1,0) ; */
    nc = wbfmm_coefficient_number(ns) ;
    /* t->mps[l] = g_malloc0(nb*2*nc*nq*sizeof(WBFMM_REAL)) ; */
    posix_memalign(&(t->mps[l]), 32, nb*2*nc*nq*sizeof(WBFMM_REAL)) ;

    c = (WBFMM_REAL *)(t->mps[l]) ;
    /*
      set box pointers to start of their coefficients, noting that
      coefficients are packed in groups of eight for shift operations
    */
    boxes = t->boxes[l] ;
    for ( i = 0 ; i < nb ; i += 8 ) {
      for ( j = 0 ; j < 8 ; j ++ ) {
	boxes[i+j].mps = &(c[i*2*nc*nq+2*j*nq]) ;
      }
    }
  }

  if ( nr != 0 ) {
    /* nc = wbfmm_coefficient_index_nm(nr+1,0) ; */
    nc = wbfmm_coefficient_number(nr) ;
    /* t->mpr[l] = g_malloc0(nb*2*nc*nq*sizeof(WBFMM_REAL)) ; */
    posix_memalign(&(t->mpr[l]), 32, nb*2*nc*nq*sizeof(WBFMM_REAL)) ;
    c = (WBFMM_REAL *)(t->mpr[l]) ;
    /*
      set box pointers to start of their coefficients, noting that
      coefficients are packed in groups of eight for shift operations
    */
    boxes = t->boxes[l] ;
    for ( i = 0 ; i < nb ; i += 8 ) {
      for ( j = 0 ; j < 8 ; j ++ ) {
	boxes[i+j].mpr = &(c[i*2*nq*nc+2*nq*j]) ;
      }
    }
  }

  return 0 ;
}

guint64 WBFMM_FUNCTION_NAME(wbfmm_point_box)(wbfmm_tree_t *t, guint level,
					     WBFMM_REAL *x)

{
  guint64 b ;
  guint nb, i, j, k ;
  WBFMM_REAL wb, *x0, D ;

  x0 = wbfmm_tree_origin(t) ;
  D = wbfmm_tree_width(t) ;

  if ( x[0] < x0[0] || x[0] > x0[0]+D ||
       x[1] < x0[1] || x[1] > x0[1]+D ||
       x[2] < x0[2] || x[2] > x0[2]+D )
    g_error("%s: point (%g,%g,%g) not in octree",
	    __FUNCTION__, x[0], x[1], x[2]) ;
  
  /* wb = t->D/nb ; */

  nb = 1 << level ;
  wb = D/nb ;
  /* dx = (x[0] - x0[0])/wb ; */
  /* if ( dx < 0.0 || dx > nb ) */
  /*   g_error("%s: point (%g,%g,%g) not in octree", */
  /* 	    __FUNCTION__, x[0], x[1], x[2]) ; */
  /* else */
  /* i = (guint32)floor(dx) ; */

  /* dx = (x[1] - x0[1])/wb ; */
  /* if ( dx < 0.0 || dx > nb ) */
  /*   g_error("%s: point (%g,%g,%g) not in octree", */
  /* 	    __FUNCTION__, x[0], x[1], x[2]) ; */
  /* else */
  /*   j = (guint32)floor(dx) ; */

  /* dx = (x[2] - x0[2])/wb ; */
  /* if ( dx < 0.0 || dx > nb ) */
  /*   g_error("%s: point (%g,%g,%g) not in octree", */
  /* 	    __FUNCTION__, x[0], x[1], x[2]) ; */
  /* else */
    /* k = (guint32)floor(dx) ; */

  i = (guint32)floor((x[0] - x0[0])/wb) ;
  j = (guint32)floor((x[1] - x0[1])/wb) ;
  k = (guint32)floor((x[2] - x0[2])/wb) ;

  b = wbfmm_box_index(i, j, k) ;

  return b ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_tree_coefficient_clear)(wbfmm_tree_t *t,
						       guint l)

{
  guint nb, ncs, ncr, ns, nr, nq ;

  g_assert(l <= wbfmm_tree_depth(t)) ;

  nq = wbfmm_tree_source_size(t) ;
  ns = t->order_s[l] ; nr = t->order_r[l] ; 

  if ( wbfmm_tree_problem(t) == WBFMM_PROBLEM_HELMHOLTZ ) {
    nq *= 2 ; 
    /* ncs = wbfmm_coefficient_index_nm(ns+1,0) ; */
    ncs = wbfmm_coefficient_number(ns) ;
    /* ncr = wbfmm_coefficient_index_nm(nr+1,0) ; */
    ncr = wbfmm_coefficient_number(ns) ;
  } else {
    ncs = (ns+1)*(ns+1) ; 
    ncr = (nr+1)*(nr+1) ; 
  }
  
  /*number of boxes at level l*/
  nb = 1 << (3*l) ;

  /*number of coefficients in singular expansions*/
  if ( ns != 0 ) {
    memset(t->mps[l], 0, nb*nq*ncs*sizeof(WBFMM_REAL)) ;
  }

  /*number of coefficients in regular expansions*/
  if ( nr != 0 ) {
    memset(t->mpr[l], 0, nb*nq*ncr*sizeof(WBFMM_REAL)) ;
  }

  return 0 ;
}

