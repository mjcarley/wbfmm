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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <inttypes.h>

#include <glib.h>

#include <wbfmm.h>

#include "wbfmm-private.h"

static gint _wbfmm_tree_box_local_field_scalar(wbfmm_tree_t *t,
					       guint level,
					       guint b, WBFMM_REAL k,
					       WBFMM_REAL *x,
					       WBFMM_REAL *f,
					       gint fstr,
					       WBFMM_REAL *src, gint sstr,
					       WBFMM_REAL *normals,
					       gint nstr,
					       WBFMM_REAL *d, gint dstr,
					       gboolean eval_neighbours,
					       guint field,
					       WBFMM_REAL *work)

{
  WBFMM_REAL xb[3], wb, *C, *xs, r, h0[2], h1[2], fR[2] ;
  wbfmm_box_t *boxes, box ;
  guint64 neighbours[27] ;
  gint nnbr, i, j, jj, idx, nq ;

  g_assert(field == WBFMM_FIELD_SCALAR ) ;
  
  nq = wbfmm_tree_source_size(t) ;
  
  boxes = t->boxes[level] ;
  C = boxes[b].mpr ;

  WBFMM_FUNCTION_NAME(wbfmm_tree_box_centre)(t, level, b, xb, &wb) ;
  
  WBFMM_FUNCTION_NAME(wbfmm_expansion_j_evaluate)(k, xb, C, 8*nq,  
						  t->order_r[level],
						  wbfmm_tree_source_size(t),
						  x, f, fstr,
						  work) ;


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
	  WBFMM_FUNCTION_NAME(wbfmm_bessel_h_init)(k*r, h0, h1) ;
	  h0[0] /= 4.0*M_PI ; h0[1] /= 4.0*M_PI ;
	  for ( jj = 0 ; jj < nq ; jj ++ ) {
	    f[fstr*jj+0] +=
	      h0[0]*src[idx*sstr+2*jj+0] - h0[1]*src[idx*sstr+2*jj+1] ;
	    f[fstr*jj+1] +=
	      h0[1]*src[idx*sstr+2*jj+0] + h0[0]*src[idx*sstr+2*jj+1] ;
	  }
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
	  WBFMM_FUNCTION_NAME(wbfmm_bessel_h_init)(k*r, h0, h1) ;
	  h1[0] /= 4.0*M_PI ; h1[1] /= 4.0*M_PI ;

	  fR[0] = (normals[idx*nstr+0]*(x[0] - xs[0]) +
		   normals[idx*nstr+1]*(x[1] - xs[1]) +
		   normals[idx*nstr+2]*(x[2] - xs[2]))/r ;
	  fR[1] = d[idx*dstr+1]*fR[0] ; fR[0] *= d[idx*dstr+0] ;

	  f[0] -= k*(h1[0]*fR[0] - h1[1]*fR[1]) ;
	  f[1] -= k*(h1[0]*fR[1] + h1[1]*fR[0]) ;
	  /* f[0] += k*(h1[0]*fR[0] - h1[1]*fR[1]) ; */
	  /* f[1] += k*(h1[0]*fR[1] + h1[1]*fR[0]) ; */
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
	  WBFMM_FUNCTION_NAME(wbfmm_bessel_h_init)(k*r, h0, h1) ;
	  h0[0] /= 4.0*M_PI ; h0[1] /= 4.0*M_PI ;
	  h1[0] /= 4.0*M_PI ; h1[1] /= 4.0*M_PI ;

	  fR[0] = (normals[idx*nstr+0]*(x[0] - xs[0]) +
		   normals[idx*nstr+1]*(x[1] - xs[1]) +
		   normals[idx*nstr+2]*(x[2] - xs[2]))/r ;
	  fR[1] = d[idx*dstr+1]*fR[0] ; fR[0] *= d[idx*dstr+0] ;
	  
	  f[0] += h0[0]*src[idx*sstr+0] - h0[1]*src[idx*sstr+1] ;
	  f[1] += h0[1]*src[idx*sstr+0] + h0[0]*src[idx*sstr+1] ;

	  /* f[0] += k*(h1[0]*fR[0] - h1[1]*fR[1]) ; */
	  /* f[1] += k*(h1[0]*fR[1] + h1[1]*fR[0]) ; */
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

static gint _wbfmm_tree_box_local_field_gradient(wbfmm_tree_t *t,
						 guint level,
						 guint b, WBFMM_REAL k,
						 WBFMM_REAL *x,
						 WBFMM_REAL *f,
						 gint fstr,
						 WBFMM_REAL *src, gint sstr,
						 WBFMM_REAL *normals,
						 gint nstr,
						 WBFMM_REAL *d, gint dstr,
						 gboolean eval_neighbours,
						 guint field,
						 WBFMM_REAL *work)

{
  WBFMM_REAL xb[3], wb, *C, *xs, r, h0[2], h1[2], nR[3], fr, fi ;
  wbfmm_box_t *boxes, box ;
  guint64 neighbours[27] ;
  gint nnbr, i, j, jj, idx, nq ;

  g_assert(field == WBFMM_FIELD_GRADIENT ) ;
  
  nq = wbfmm_tree_source_size(t) ;
  
  boxes = t->boxes[level] ;
  C = boxes[b].mpr ;

  WBFMM_FUNCTION_NAME(wbfmm_tree_box_centre)(t, level, b, xb, &wb) ;
  
  WBFMM_FUNCTION_NAME(wbfmm_expansion_j_grad_evaluate)(k, xb, C, 8*nq,  
						       t->order_r[level],
						       wbfmm_tree_source_size(t),
						       x, f, fstr,
						       work) ;


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
	  nR[0] = (x[0] - xs[0])/r*0.25*M_1_PI ;
	  nR[1] = (x[1] - xs[1])/r*0.25*M_1_PI ;
	  nR[2] = (x[2] - xs[2])/r*0.25*M_1_PI ;
	  WBFMM_FUNCTION_NAME(wbfmm_bessel_h_init)(k*r, h0, h1) ;
	  for ( jj = 0 ; jj < nq ; jj ++ ) {
	    fr = h1[0]*src[idx*sstr+2*jj+0] - h1[1]*src[idx*sstr+2*jj+1] ;
	    fi = h1[1]*src[idx*sstr+2*jj+0] + h1[0]*src[idx*sstr+2*jj+1] ;
	    f[jj*fstr+0] -= k*fr*nR[0] ;
	    f[jj*fstr+1] -= k*fi*nR[0] ;
	    f[jj*fstr+2] -= k*fr*nR[1] ;
	    f[jj*fstr+3] -= k*fi*nR[1] ;
	    f[jj*fstr+4] -= k*fr*nR[2] ;
	    f[jj*fstr+5] -= k*fi*nR[2] ;
	    /* f[fstr*jj+0] += */
	    /*   h0[0]*src[idx*sstr+2*jj+0] - h0[1]*src[idx*sstr+2*jj+1] ; */
	    /* f[fstr*jj+1] += */
	    /*   h0[1]*src[idx*sstr+2*jj+0] + h0[0]*src[idx*sstr+2*jj+1] ; */
	  }
	}
      }
    }

    return 0 ;
  }
  
  g_assert_not_reached() ; 
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_tree_box_local_field)(wbfmm_tree_t *t,
						     guint level,
						     guint b, WBFMM_REAL k,
						     WBFMM_REAL *x,
						     WBFMM_REAL *f,
						     gint fstr,
						     WBFMM_REAL *src, gint sstr,
						     WBFMM_REAL *normals,
						     gint nstr,
						     WBFMM_REAL *d, gint dstr,
						     gboolean eval_neighbours,
						     guint field,
						     WBFMM_REAL *work)

{
  WBFMM_REAL xb[3], *C, *xs, r, h0[2], h1[2], fR[2] ;
  wbfmm_box_t *boxes, box ;
  guint64 neighbours[27] ;
  gint nnbr, i, j, jj, idx, nq ;

  g_assert(t->problem == WBFMM_PROBLEM_HELMHOLTZ ) ;
  /* nq = wbfmm_tree_source_size(t) ; */
  
  /* boxes = t->boxes[level] ; */
  /* C = boxes[b].mpr ; */

  /* WBFMM_FUNCTION_NAME(wbfmm_tree_box_centre)(t, level, b, xb, &wb) ; */

  switch ( field ) {
  default:
    g_error("%s: unrecognized field identifier (%u)",
	    __FUNCTION__, field) ;
    break ;
  case WBFMM_FIELD_SCALAR:
    if ( fstr < 2 )
      g_error("%s: field stride (%d) too small for field calculation",
	      __FUNCTION__, fstr) ;
    return _wbfmm_tree_box_local_field_scalar(t, level, b, k, x, f, fstr,
					      src, sstr, normals, nstr, d,
					      dstr, eval_neighbours, field,
					      work) ;
    break ;
  case WBFMM_FIELD_GRADIENT:
    if ( fstr < 6 )
      g_error("%s: field stride (%d) too small for gradient calculation",
	      __FUNCTION__, fstr) ;
    return _wbfmm_tree_box_local_field_gradient(t, level, b, k, x, f, fstr,
						src, sstr, normals, nstr, d,
						dstr, eval_neighbours, field,
						work) ;
    break ;
  case WBFMM_FIELD_SCALAR | WBFMM_FIELD_GRADIENT :
    g_assert_not_reached() ;
    break ;
  }

  g_assert_not_reached() ;

  return 0 ;
  
  WBFMM_FUNCTION_NAME(wbfmm_expansion_j_evaluate)(k, xb, C, 8*nq,  
						  t->order_r[level],
						  wbfmm_tree_source_size(t),
						  x, f, fstr,
						  work) ;


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
	  WBFMM_FUNCTION_NAME(wbfmm_bessel_h_init)(k*r, h0, h1) ;
	  h0[0] /= 4.0*M_PI ; h0[1] /= 4.0*M_PI ;
	  for ( jj = 0 ; jj < nq ; jj ++ ) {
	    f[fstr*jj+0] +=
	      h0[0]*src[idx*sstr+2*jj+0] - h0[1]*src[idx*sstr+2*jj+1] ;
	    f[fstr*jj+1] +=
	      h0[1]*src[idx*sstr+2*jj+0] + h0[0]*src[idx*sstr+2*jj+1] ;
	  }
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
	  WBFMM_FUNCTION_NAME(wbfmm_bessel_h_init)(k*r, h0, h1) ;
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
	  WBFMM_FUNCTION_NAME(wbfmm_bessel_h_init)(k*r, h0, h1) ;
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
