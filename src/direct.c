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

gint WBFMM_FUNCTION_NAME(wbfmm_total_field)(WBFMM_REAL k,
					    WBFMM_REAL *xs, gint xstride,
					    WBFMM_REAL *src, gint sstride,
					    WBFMM_REAL *normals, gint nstr,
					    WBFMM_REAL *dipoles, gint dstr,
					    gint nq, gint nsrc,
					    WBFMM_REAL *xf, WBFMM_REAL *field,
					    gint fstr)

{
  gint i, j ;
  WBFMM_REAL r, th, ph, h0[2], h1[2], fR[2], fd[6], norm ;

  if ( fstr < 2 && nq > 1 )
    g_error("%s: field stride (%d) must be greater than 1 for "
	    "multi-component sources (nq=%d)", __FUNCTION__, fstr, nq) ;
  
  if ( src == NULL && normals == NULL && dipoles == NULL ) return 0 ;

  if ( normals != NULL && dipoles == NULL )
    g_error("%s: normals specified but no dipole strengths (dipoles == NULL)",
	    __FUNCTION__) ;

  norm = 0.25*M_1_PI ;
  if ( normals == NULL && dipoles == NULL ) {
    for ( i = 0 ; i < nsrc ; i ++ ) {
      WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(&(xs[i*xstride]), xf, 
						  &r, &th, &ph) ;
      WBFMM_FUNCTION_NAME(wbfmm_bessel_h_init)(k*r, h0, h1) ;
      
      for ( j = 0 ; j < nq ; j ++ ) {
	field[j*fstr + 0] +=
	  (h0[0]*src[i*sstride+2*j+0] - h0[1]*src[i*sstride+2*j+1])*norm ;
	field[j*fstr + 1] +=
	  (h0[1]*src[i*sstride+2*j+0] + h0[0]*src[i*sstride+2*j+1])*norm ;
      }
    }

    return 0 ;
  }

  /* g_assert(nq == 1) ; */
  if ( src != NULL && normals != NULL ) {
    g_assert(sstride >= 2*nq) ; 
    g_assert(dstr >= 2*nq) ; 
    for ( i = 0 ; i < nsrc ; i ++ ) {
      WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(&(xs[i*xstride]), xf, 
						  &r, &th, &ph) ;
      WBFMM_FUNCTION_NAME(wbfmm_bessel_h_init)(k*r, h0, h1) ;

      for ( j = 0 ; j < nq ; j ++ ) {      
	fd[0] = normals[i*nstr+0]*dipoles[i*dstr+2*j+0] ;
	fd[1] = normals[i*nstr+0]*dipoles[i*dstr+2*j+1] ;
	fd[2] = normals[i*nstr+1]*dipoles[i*dstr+2*j+0] ;
	fd[3] = normals[i*nstr+1]*dipoles[i*dstr+2*j+1] ;
	fd[4] = normals[i*nstr+2]*dipoles[i*dstr+2*j+0] ;
	fd[5] = normals[i*nstr+2]*dipoles[i*dstr+2*j+1] ;

	fR[0]  = fd[0]*(xf[0] - xs[i*xstride+0]) ;
	fR[1]  = fd[1]*(xf[0] - xs[i*xstride+0]) ;
	fR[0] += fd[2]*(xf[1] - xs[i*xstride+1]) ;
	fR[1] += fd[3]*(xf[1] - xs[i*xstride+1]) ;
	fR[0] += fd[4]*(xf[2] - xs[i*xstride+2]) ;
	fR[1] += fd[5]*(xf[2] - xs[i*xstride+2]) ;
	
	fR[0] /= r ; fR[1] /= r ;
	
	field[j*fstr+0] += h0[0]*src[i*sstride+2*j+0] -
	  h0[1]*src[i*sstride+2*j+1] ;
	field[j*fstr+1] += h0[1]*src[i*sstride+2*j+0] +
	  h0[0]*src[i*sstride+2*j+1] ;
	
	field[j*fstr+0] -= k*(h1[0]*fR[0] - h1[1]*fR[1]) ;
	field[j*fstr+1] -= k*(h1[0]*fR[1] + h1[1]*fR[0]) ;
      }
	/* fd[0] = normals[i*nstr+0]*dipoles[i*dstr+0] ; */
	/* fd[1] = normals[i*nstr+0]*dipoles[i*dstr+1] ; */
	/* fd[2] = normals[i*nstr+1]*dipoles[i*dstr+0] ; */
	/* fd[3] = normals[i*nstr+1]*dipoles[i*dstr+1] ; */
	/* fd[4] = normals[i*nstr+2]*dipoles[i*dstr+0] ; */
	/* fd[5] = normals[i*nstr+2]*dipoles[i*dstr+1] ; */

	/* fR[0]  = fd[0]*(xf[0] - xs[i*xstride+0]) ; */
	/* fR[1]  = fd[1]*(xf[0] - xs[i*xstride+0]) ; */
	/* fR[0] += fd[2]*(xf[1] - xs[i*xstride+1]) ; */
	/* fR[1] += fd[3]*(xf[1] - xs[i*xstride+1]) ; */
	/* fR[0] += fd[4]*(xf[2] - xs[i*xstride+2]) ; */
	/* fR[1] += fd[5]*(xf[2] - xs[i*xstride+2]) ; */
	
	/* fR[0] /= r ; fR[1] /= r ; */
	
	/* field[2*j+0] += h0[0]*src[i*sstride+0] - h0[1]*src[i*sstride+1] ; */
	/* field[2*j+1] += h0[1]*src[i*sstride+0] + h0[0]*src[i*sstride+1] ; */
	
	/* field[2*j+0] -= k*(h1[0]*fR[0] - h1[1]*fR[1]) ; */
	/* field[2*j+1] -= k*(h1[0]*fR[1] + h1[1]*fR[0]) ; */
    }

    /*G&D normalization of Legendre polynomials*/
    for ( j = 0 ; j < 2*nq ; j ++ ) field[j] /= 4.0*M_PI ;

    return 0 ;
  }

  g_assert_not_reached() ; /*following code needs modification*/
  
  if ( src == NULL && dipoles != NULL ) {
    g_assert(sstride >= 2*nq) ; 
    g_assert(dstr >= 2*nq) ; 
    for ( i = 0 ; i < nsrc ; i ++ ) {
      WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(&(xs[i*xstride]), xf, 
						  &r, &th, &ph) ;
      WBFMM_FUNCTION_NAME(wbfmm_bessel_h_init)(k*r, h0, h1) ;

      fd[0] = normals[i*nstr+0]*dipoles[i*dstr+0] ;
      fd[1] = normals[i*nstr+0]*dipoles[i*dstr+1] ;
      fd[2] = normals[i*nstr+1]*dipoles[i*dstr+0] ;
      fd[3] = normals[i*nstr+1]*dipoles[i*dstr+1] ;
      fd[4] = normals[i*nstr+2]*dipoles[i*dstr+0] ;
      fd[5] = normals[i*nstr+2]*dipoles[i*dstr+1] ;

      fR[0]  = fd[0]*(xf[0] - xs[i*xstride+0]) ;
      fR[1]  = fd[1]*(xf[0] - xs[i*xstride+0]) ;
      fR[0] += fd[2]*(xf[1] - xs[i*xstride+1]) ;
      fR[1] += fd[3]*(xf[1] - xs[i*xstride+1]) ;
      fR[0] += fd[4]*(xf[2] - xs[i*xstride+2]) ;
      fR[1] += fd[5]*(xf[2] - xs[i*xstride+2]) ;
      
      fR[0] /= r ; fR[1] /= r ;
      
      field[0] -= k*(h1[0]*fR[0] - h1[1]*fR[1]) ;
      field[1] -= k*(h1[0]*fR[1] + h1[1]*fR[0]) ;
    }

    /*G&D normalization of Legendre polynomials*/
    field[0] /= 4.0*M_PI ; field[1] /= 4.0*M_PI ;

    return 0 ;
  }

  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_total_field_grad)(WBFMM_REAL k,
						 WBFMM_REAL *xs, gint xstride,
						 WBFMM_REAL *src, gint sstride,
						 WBFMM_REAL *normals, gint nstr,
						 WBFMM_REAL *dipoles, gint dstr,
						 gint nq, gint nsrc,
						 WBFMM_REAL *xf,
						 WBFMM_REAL *field, gint fstr)

{
  gint i, j ;
  WBFMM_REAL r, th, ph, h0[2], h1[2], fR[2], fd[6], nR[3], fr, fi ;

  if ( src == NULL && normals == NULL && dipoles == NULL ) return 0 ;

  if ( normals != NULL && dipoles == NULL )
    g_error("%s: normals specified but no dipole strengths (dipoles == NULL)",
	    __FUNCTION__) ;

  if ( normals == NULL && dipoles == NULL ) {
    for ( i = 0 ; i < nsrc ; i ++ ) {
      WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(&(xs[i*xstride]), xf, 
							&r, &th, &ph) ;
      WBFMM_FUNCTION_NAME(wbfmm_bessel_h_init)(k*r, h0, h1) ;

      nR[0] = (xf[0] - xs[i*xstride+0])/r*0.25*M_1_PI ;
      nR[1] = (xf[1] - xs[i*xstride+1])/r*0.25*M_1_PI ;
      nR[2] = (xf[2] - xs[i*xstride+2])/r*0.25*M_1_PI ;

      for ( j = 0 ; j < nq ; j ++ ) {
	fr = h1[0]*src[i*sstride+2*j+0] - h1[1]*src[i*sstride+2*j+1] ;
	fi = h1[1]*src[i*sstride+2*j+0] + h1[0]*src[i*sstride+2*j+1] ;
	field[j*fstr+0] -= k*fr*nR[0] ;
	field[j*fstr+1] -= k*fi*nR[0] ;
	field[j*fstr+2] -= k*fr*nR[1] ;
	field[j*fstr+3] -= k*fi*nR[1] ;
	field[j*fstr+4] -= k*fr*nR[2] ;
	field[j*fstr+5] -= k*fi*nR[2] ;
      }
    }

    return 0 ;
  }

  g_assert_not_reached() ; /*gradients of dipole fields not
			     implemented yet*/
  
  if ( src == NULL && dipoles != NULL ) {
    for ( i = 0 ; i < nsrc ; i ++ ) {
      WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(&(xs[i*xstride]), xf, 
						  &r, &th, &ph) ;
      WBFMM_FUNCTION_NAME(wbfmm_bessel_h_init)(k*r, h0, h1) ;

      fd[0] = normals[i*nstr+0]*dipoles[i*dstr+0] ;
      fd[1] = normals[i*nstr+0]*dipoles[i*dstr+1] ;
      fd[2] = normals[i*nstr+1]*dipoles[i*dstr+0] ;
      fd[3] = normals[i*nstr+1]*dipoles[i*dstr+1] ;
      fd[4] = normals[i*nstr+2]*dipoles[i*dstr+0] ;
      fd[5] = normals[i*nstr+2]*dipoles[i*dstr+1] ;

      fR[0]  = fd[0]*(xf[0] - xs[i*xstride+0]) ;
      fR[1]  = fd[1]*(xf[0] - xs[i*xstride+0]) ;
      fR[0] += fd[2]*(xf[1] - xs[i*xstride+1]) ;
      fR[1] += fd[3]*(xf[1] - xs[i*xstride+1]) ;
      fR[0] += fd[4]*(xf[2] - xs[i*xstride+2]) ;
      fR[1] += fd[5]*(xf[2] - xs[i*xstride+2]) ;
      
      fR[0] /= r ; fR[1] /= r ;
      
      field[0] -= k*(h1[0]*fR[0] - h1[1]*fR[1]) ;
      field[1] -= k*(h1[0]*fR[1] + h1[1]*fR[0]) ;
    }

    /*G&D normalization of Legendre polynomials*/
    field[0] /= 4.0*M_PI ; field[1] /= 4.0*M_PI ;

    return 0 ;
  }

  if ( src != NULL && normals != NULL ) {
    for ( i = 0 ; i < nsrc ; i ++ ) {
      WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(&(xs[i*xstride]), xf, 
						  &r, &th, &ph) ;
      WBFMM_FUNCTION_NAME(wbfmm_bessel_h_init)(k*r, h0, h1) ;

      fd[0] = normals[i*nstr+0]*dipoles[i*dstr+0] ;
      fd[1] = normals[i*nstr+0]*dipoles[i*dstr+1] ;
      fd[2] = normals[i*nstr+1]*dipoles[i*dstr+0] ;
      fd[3] = normals[i*nstr+1]*dipoles[i*dstr+1] ;
      fd[4] = normals[i*nstr+2]*dipoles[i*dstr+0] ;
      fd[5] = normals[i*nstr+2]*dipoles[i*dstr+1] ;

      fR[0]  = fd[0]*(xf[0] - xs[i*xstride+0]) ;
      fR[1]  = fd[1]*(xf[0] - xs[i*xstride+0]) ;
      fR[0] += fd[2]*(xf[1] - xs[i*xstride+1]) ;
      fR[1] += fd[3]*(xf[1] - xs[i*xstride+1]) ;
      fR[0] += fd[4]*(xf[2] - xs[i*xstride+2]) ;
      fR[1] += fd[5]*(xf[2] - xs[i*xstride+2]) ;
      
      fR[0] /= r ; fR[1] /= r ;
      
      field[0] += h0[0]*src[i*sstride+0] - h0[1]*src[i*sstride+1] ;
      field[1] += h0[1]*src[i*sstride+0] + h0[0]*src[i*sstride+1] ;

      field[0] -= k*(h1[0]*fR[0] - h1[1]*fR[1]) ;
      field[1] -= k*(h1[0]*fR[1] + h1[1]*fR[0]) ;
    }

    /*G&D normalization of Legendre polynomials*/
    field[0] /= 4.0*M_PI ; field[1] /= 4.0*M_PI ;

    return 0 ;
  }
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_total_dipole_field)(WBFMM_REAL k,
						   WBFMM_REAL *xs,
						   gint xstride,
						   WBFMM_REAL *src,
						   gint sstride,
						   gint nsrc,
						   WBFMM_REAL *xf,
						   WBFMM_REAL *field)

{
  gint i ;
  WBFMM_REAL r, th, ph, h0[2], h1[2], fR[2] ;

  /* field[0] = field[1] = 0.0 ; */

  for ( i = 0 ; i < nsrc ; i ++ ) {
    WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(&(xs[i*xstride]), xf, 
						&r, &th, &ph) ;
    WBFMM_FUNCTION_NAME(wbfmm_bessel_h_init)(k*r, h0, h1) ;
    fR[0]  = src[i*sstride+0]*(xf[0] - xs[i*xstride+0]) ;
    fR[1]  = src[i*sstride+1]*(xf[0] - xs[i*xstride+0]) ;
    fR[0] += src[i*sstride+2]*(xf[1] - xs[i*xstride+1]) ;
    fR[1] += src[i*sstride+3]*(xf[1] - xs[i*xstride+1]) ;
    fR[0] += src[i*sstride+4]*(xf[2] - xs[i*xstride+2]) ;
    fR[1] += src[i*sstride+5]*(xf[2] - xs[i*xstride+2]) ;

    fR[0] /= r ; fR[1] /= r ;
    
    field[0] -= k*(h1[0]*fR[0] - h1[1]*fR[1]) ;
    field[1] -= k*(h1[0]*fR[1] + h1[1]*fR[0]) ;
  }

  /*G&D normalization of Legendre polynomials*/
  field[0] /= 4.0*M_PI ; field[1] /= 4.0*M_PI ;

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_total_normal_field)(WBFMM_REAL k,
						   WBFMM_REAL *xs,
						   gint xstride,
						   WBFMM_REAL *ns,
						   gint nstride,
						   WBFMM_REAL *src,
						   gint sstride,
						   gint nsrc,
						   WBFMM_REAL *xf,
						   WBFMM_REAL *field)

{
  gint i ;
  WBFMM_REAL q[6], f[2] = {0.0} ;

  /* field[0] = field[1] = 0.0 ; */
  for ( i = 0 ; i < nsrc ; i ++ ) {
    q[0] = ns[i*nstride+0]*src[i*sstride+0] ;
    q[1] = ns[i*nstride+0]*src[i*sstride+1] ;
    q[2] = ns[i*nstride+1]*src[i*sstride+0] ;
    q[3] = ns[i*nstride+1]*src[i*sstride+1] ;
    q[4] = ns[i*nstride+2]*src[i*sstride+0] ;
    q[5] = ns[i*nstride+2]*src[i*sstride+1] ;

    f[0] = f[1] = 0.0 ;
    WBFMM_FUNCTION_NAME(wbfmm_total_dipole_field)(k, &(xs[i*xstride]), 0,
					    q, 0, 1, xf, f) ;
    field[0] += f[0] ; field[1] += f[1] ;
  }
  
  return 0 ;
}
