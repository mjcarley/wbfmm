/* This file is part of WBFMM, a Wide-Band Fast Multipole Method code
 *
 * Copyright (C) 2019, 2021, 2024 Michael Carley
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

#include <glib.h>

#include <wbfmm.h>

#include "wbfmm-private.h"

#ifdef HAVE_AVX_INSTRUCTIONS
#include <immintrin.h>
#endif /*HAVE_AVX_INSTRUCTIONS*/

gint WBFMM_FUNCTION_NAME(wbfmm_laplace_field_grad)(WBFMM_REAL *xs,
						   gint xstride,
						   WBFMM_REAL *src,
						   gint sstride,
						   gint nq,
						   WBFMM_REAL *normals,
						   gint nstr,
						   WBFMM_REAL *dipoles,
						   gint dstr,
						   gint nsrc,
						   WBFMM_REAL *xf,
						   WBFMM_REAL *field,
						   gint fstr)

{
  gint i, j ;
  WBFMM_REAL r, th, ph, nR[3] ;

  if ( src == NULL && normals == NULL && dipoles == NULL ) return 0 ;

  if ( fstr < 3 && nq > 1 )
    g_error("%s: field data stride (%d) must be greater than two",
	    __FUNCTION__, fstr) ;

  g_assert(sstride >= nq) ;
  
  if ( normals != NULL && dipoles == NULL )
    g_error("%s: normals specified but no dipole strengths (dipoles == NULL)",
	    __FUNCTION__) ;

  if ( normals == NULL && dipoles == NULL ) {
    for ( i = 0 ; i < nsrc ; i ++ ) {
      WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(&(xs[i*xstride]), xf, 
							&r, &th, &ph) ;
      nR[0] = (xf[0] - xs[i*xstride+0])/r*0.25*M_1_PI ;
      nR[1] = (xf[1] - xs[i*xstride+1])/r*0.25*M_1_PI ;
      nR[2] = (xf[2] - xs[i*xstride+2])/r*0.25*M_1_PI ;
      for ( j = 0 ; j < nq ; j ++ ) {
	field[fstr*j+0] -= src[i*sstride+j]/r/r*nR[0] ;
	field[fstr*j+1] -= src[i*sstride+j]/r/r*nR[1] ;
	field[fstr*j+2] -= src[i*sstride+j]/r/r*nR[2] ;
      }
    }

    return 0 ;
  }

  g_assert_not_reached() ;
  
  return 0 ;
}

static gint laplace_field_curl(WBFMM_REAL *xs,
			       gint xstride,
			       WBFMM_REAL *src,
			       gint sstride,
			       gint nq,
			       WBFMM_REAL *normals,
			       gint nstr,
			       WBFMM_REAL *dipoles,
			       gint dstr,
			       gint nsrc,
			       WBFMM_REAL *xf,
			       WBFMM_REAL *field,
			       gint fstr)

{
  gint i ;
  WBFMM_REAL r, th, ph, nR[3] ;

  if ( src == NULL && normals == NULL && dipoles == NULL ) return 0 ;

  if ( normals == NULL && dipoles == NULL ) {
    for ( i = 0 ; i < nsrc ; i ++ ) {
      WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(&(xs[i*xstride]), xf, 
							&r, &th, &ph) ;
      if ( r > WBFMM_LOCAL_CUTOFF_RADIUS ) {
	nR[0] = (xf[0] - xs[i*xstride+0])/r/r/r*0.25*M_1_PI ;
	nR[1] = (xf[1] - xs[i*xstride+1])/r/r/r*0.25*M_1_PI ;
	nR[2] = (xf[2] - xs[i*xstride+2])/r/r/r*0.25*M_1_PI ;

	field[0] -= src[i*sstride+2]*nR[1] - src[i*sstride+1]*nR[2] ;
	field[1] -= src[i*sstride+0]*nR[2] - src[i*sstride+2]*nR[0] ;
	field[2] -= src[i*sstride+1]*nR[0] - src[i*sstride+0]*nR[1] ;
      }
    }

    return 0 ;
  }

  g_assert_not_reached() ;
  
  return 0 ;
}

static gint laplace_field_curl_gradient(WBFMM_REAL *xs,
					gint xstride,
					WBFMM_REAL *src,
					gint sstride,
					gint nq,
					WBFMM_REAL *normals,
					gint nstr,
					WBFMM_REAL *dipoles,
					gint dstr,
					gint nsrc,
					WBFMM_REAL *xf,
					WBFMM_REAL *field,
					gint fstr)

{
  gint i ;
  WBFMM_REAL r, th, ph, nR[3] ;

  if ( src == NULL && normals == NULL && dipoles == NULL ) return 0 ;

  if ( normals == NULL && dipoles == NULL ) {
    for ( i = 0 ; i < nsrc ; i ++ ) {
      WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(&(xs[i*xstride]), xf, 
							&r, &th, &ph) ;
      if ( r > WBFMM_LOCAL_CUTOFF_RADIUS ) {
	nR[0] = (xf[0] - xs[i*xstride+0])/r/r/r*0.25*M_1_PI ;
	nR[1] = (xf[1] - xs[i*xstride+1])/r/r/r*0.25*M_1_PI ;
	nR[2] = (xf[2] - xs[i*xstride+2])/r/r/r*0.25*M_1_PI ;

	field[0] -= src[i*sstride+2]*nR[1] - src[i*sstride+1]*nR[2] ;
	field[1] -= src[i*sstride+0]*nR[2] - src[i*sstride+2]*nR[0] ;
	field[2] -= src[i*sstride+1]*nR[0] - src[i*sstride+0]*nR[1] ;
      }
    }

    return 0 ;
  }

  g_assert_not_reached() ;
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_laplace_field)(WBFMM_REAL *xs, gint xstride,
					      WBFMM_REAL *src, gint sstride,
					      gint nq,
					      WBFMM_REAL *normals, gint nstr,
					      WBFMM_REAL *dipoles, gint dstr,
					      gint nsrc,
					      WBFMM_REAL *xf, WBFMM_REAL *field,
					      gint fstr)

{
  gint i, j ;
  WBFMM_REAL r, th, ph, nr ;

  if ( src == NULL && normals == NULL && dipoles == NULL ) return 0 ;

  g_assert(sstride >= nq) ;
  
  if ( normals != NULL && dipoles == NULL )
    g_error("%s: normals specified but no dipole strengths (dipoles == NULL)",
	    __FUNCTION__) ;

  if ( normals == NULL && dipoles == NULL ) {
    for ( i = 0 ; i < nsrc ; i ++ ) {
      WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(&(xs[i*xstride]), xf, 
							&r, &th, &ph) ;
      for ( j = 0 ; j < nq ; j ++ )
	field[j*fstr] += src[i*sstride+j]/r/4.0/M_PI ;
    }

    return 0 ;
  }

  if ( src == NULL && normals != NULL ) {
    /*dipoles only*/
    for ( i = 0 ; i < nsrc ; i ++ ) {
      WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(&(xs[i*xstride]), xf, 
							&r, &th, &ph) ;
      nr =
	(xf[0] - xs[i*xstride+0])*normals[i*nstr+0] +
	(xf[1] - xs[i*xstride+1])*normals[i*nstr+1] + 
	(xf[2] - xs[i*xstride+2])*normals[i*nstr+2] ;
      nr /= r*r*r ;
      for ( j = 0 ; j < nq ; j ++ )
	field[j*fstr] += dipoles[i*dstr+j]*nr/4.0/M_PI ;
    }
    
    /* for ( j = 0 ; j < nq ; j ++ ) field[j*fstr] /= 4.0*M_PI ; */

    return 0 ;
  }

  if ( src != NULL && normals != NULL ) {
    /*monopoles and dipoles*/
    for ( i = 0 ; i < nsrc ; i ++ ) {
      WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(&(xs[i*xstride]), xf, 
							&r, &th, &ph) ;
      nr =
	(xf[0] - xs[i*xstride+0])*normals[i*nstr+0] +
	(xf[1] - xs[i*xstride+1])*normals[i*nstr+1] + 
	(xf[2] - xs[i*xstride+2])*normals[i*nstr+2] ;
      nr /= r*r*r ;
      for ( j = 0 ; j < nq ; j ++ )
	field[j*fstr] +=
	  dipoles[i*dstr+j]*nr + src[i*sstride+j]/r/4.0/M_PI ;
    }
    
    return 0 ;
  }
  
  g_assert_not_reached() ;
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_laplace_field_direct)(WBFMM_REAL *xs,
						     gint xstr,
						     WBFMM_REAL *n,
						     gint nstr,
						     gint nsrc,
						     WBFMM_REAL *src,
						     gint sstr,
						     WBFMM_REAL *d,
						     gint dstr,
						     gint nq,
						     guint field,
						     WBFMM_REAL *xf,
						     WBFMM_REAL *f,
						     gint fstr)

{
  if ( src == NULL && d == NULL ) return 0 ;

  if ( n != NULL && d == NULL )
    g_error("%s: normals specified but no dipole strengths (d == NULL)",
	    __FUNCTION__) ;
  
  switch ( field ) {
  default:
    g_error("%s: unrecognized field type %u\n", __FUNCTION__, field) ;
    break ;
  case WBFMM_FIELD_SCALAR:
    WBFMM_FUNCTION_NAME(wbfmm_laplace_field)(xs, xstr, src, sstr, nq,
					     n, nstr, d, dstr, nsrc,
					     xf, f, fstr) ;
    break ;
  case WBFMM_FIELD_GRADIENT:
  case WBFMM_FIELD_SCALAR | WBFMM_FIELD_GRADIENT:
    WBFMM_FUNCTION_NAME(wbfmm_laplace_field_grad)(xs, xstr, src, sstr, nq,
						  n, nstr, d, dstr, nsrc,
						  xf, f, fstr) ;
    break ;
  case WBFMM_FIELD_CURL:
  case WBFMM_FIELD_SCALAR | WBFMM_FIELD_CURL:
    if ( nq < 3 ) {
      g_error("%s: not enough source components (%d) for curl calculation",
	      __FUNCTION__, nq) ;
    }
    
    laplace_field_curl(xs, xstr, src, sstr, nq,
		       n, nstr, d, dstr, nsrc,
		       xf, f, fstr) ;
    break ;
  case WBFMM_FIELD_CURL | WBFMM_FIELD_GRADIENT:
    if ( nq < 3 ) {
      g_error("%s: not enough source components (%d) for curl calculation",
	      __FUNCTION__, nq) ;
    }
    
    laplace_field_curl_gradient(xs, xstr, src, sstr, nq,
				n, nstr, d, dstr, nsrc,
				xf, f, fstr) ;
    break ;
  }
  
  return 0 ;
}
