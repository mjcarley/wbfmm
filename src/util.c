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

gint WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(WBFMM_REAL *x0,
						       WBFMM_REAL *x,
						       WBFMM_REAL *r,
						       WBFMM_REAL *th,
						       WBFMM_REAL *ph)

{
  *r = 
    (x[0]-x0[0])*(x[0]-x0[0]) +
    (x[1]-x0[1])*(x[1]-x0[1]) +
    (x[2]-x0[2])*(x[2]-x0[2]) ;
  if ( *r == 0.0 ) { *ph = *th = 0.0 ; return 0 ; }

  *r = SQRT((*r)) ;
  *ph = ATAN2(x[1]-x0[1], x[0]-x0[0]) ;

  *th = ACOS((x[2]-x0[2])/(*r)) ;

  return 0 ;
}

/*
  inputs: Pnm1, Pn, normalized Legendre functions for n-1, and n

  \bar{P}_{n}^{m} = (-1)^m 
  \sqrt((n-m)!/(n+m)!\times (2n+1)/4\pi)*P_n^|m|(\cos\theta) 

  C = \cos\theta, S = \sin\theta

  on output: arrays will contain Legendre functions for n, n+1 (note
  this means the arrays are swapped, which is why they are passed as
  pointers to pointers) so they have the same sense, but with n
  incremented

  there is no check on array bounds
*/

gint WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(WBFMM_REAL **Pnm1,
							 WBFMM_REAL **Pn,
							 gint n,
							 WBFMM_REAL C,
							 WBFMM_REAL S)

{
  gint m ;
  WBFMM_REAL *pn = *Pn, *pnm1 = *Pnm1, sq2np1 ;

  /*Cheng, Crutchfield, Gimbutas, Greengard, et al. normalization*/
  /* pnm1[n+1] = -SQRT((2.0*n+3)/(2.0*n+2))*S*pn[n] ; */

  /*Gumerov and Duraiswami normalization*/
  pnm1[n+1] = SQRT((2.0*n+3)/(2.0*n+2))*S*pn[n] ;
  sq2np1 = SQRT(2.0*n+1) ;

  for ( m = 0 ; m <= n ; m ++ ) {
    pnm1[m] = C*sq2np1*pn[m] - SQRT((WBFMM_REAL)(n*n-m*m)/(2*n-1))*pnm1[m] ;
    pnm1[m] /= SQRT((WBFMM_REAL)((n+1)*(n+1)-m*m)/(2*n+3)) ;
  }

  /*swap the pointers*/
  *Pn = pnm1 ; *Pnm1 = pn ;

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_legendre_init)(WBFMM_REAL C, WBFMM_REAL S, 
					      WBFMM_REAL *P0, WBFMM_REAL *P10,
					      WBFMM_REAL *P11)

{
  /*Cheng, Crutchfield, Gimbutas, Greengard, et al. normalization*/
  /* *P0 = 1.0 ; */
  /* *P10 = C*SQRT(3.0) ; */
  /* *P11 = -S*SQRT(1.5) ; */

  /*Gumerov and Duraiswami normalization*/
  *P0 = 1.0 ;
  *P10 = C*SQRT(3.0) ;
  *P11 = S*SQRT(1.5) ;
  *P0 /= sqrt(4*M_PI) ;
  *P10 /= sqrt(4*M_PI) ;
  *P11 /= sqrt(4*M_PI) ;

  return 0 ;
}

static WBFMM_REAL wbfmm_bessel_j_n_series(gint n, WBFMM_REAL x)

{
  WBFMM_REAL jn, x2q = 1.0, x2 = x*x, sc ;
  gint q ;
#ifdef WBFMM_SINGLE_PRECISION
  WBFMM_REAL tol = 1e-7 ;
#else
  WBFMM_REAL tol = 1e-15 ;
#endif /*WBFMM_SINGLE_PRECISION*/

  sc = 1.0/(2*n+1) ;
  for ( q = 0 ; q < n ; q ++ ) sc *= x/(2*q+1) ;

  jn = x2q ;
  for ( q = 0 ; (q < 64) && ( fabs(x2q) > tol) ; q ++ ) {
    x2q *= -x2/((2*q+2)*(2*n+2*(q+1)+1)) ;
    jn += x2q ;
  }

  jn *= sc ;

  return jn ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_bessel_j_init)(WBFMM_REAL x, WBFMM_REAL *j0,
					      WBFMM_REAL *j1)

{
  if ( x > 0.0 ) {
    *j0 = SIN(x)/x ; *j1 = *j0/x - COS(x)/x ;
    return 0 ;
  }

  *j0 = 1.0 ; *j1 = 0.0 ;

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_bessel_j_recursion)(WBFMM_REAL *jnm1,
						   WBFMM_REAL *jn, 
						   WBFMM_REAL x, gint n)

{
  WBFMM_REAL jnp1, cutoff = n+1 ;

  if ( x > cutoff ) {
    jnp1 = (*jn)*(2*n+1)/x - (*jnm1) ;
  
    *jnm1 = *jn ; *jn = jnp1 ;

    return 0 ;
  }

  *jnm1 = *jn ; *jn =  wbfmm_bessel_j_n_series(n+1, x) ;

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_bessel_h_init)(WBFMM_REAL x, WBFMM_REAL *h0,
					      WBFMM_REAL *h1)

{
  h0[0] = SIN(x)/x ; h0[1] = -COS(x)/x ;
  h1[0] = h0[0]/x + h0[1] ;
  h1[1] = h0[1]/x - h0[0] ;

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_bessel_h_recursion)(WBFMM_REAL *hnm1,
						   WBFMM_REAL *hn, 
						   WBFMM_REAL x, gint n)

{
  WBFMM_REAL hnp1 ;

  hnp1 = hn[0]*(2*n+1)/x - hnm1[0] ;
  hnm1[0] = hn[0] ; hn[0] = hnp1 ;

  hnp1 = hn[1]*(2*n+1)/x - hnm1[1] ;
  hnm1[1] = hn[1] ; hn[1] = hnp1 ;

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_total_field)(WBFMM_REAL k,
					    WBFMM_REAL *xs, gint xstride,
					    WBFMM_REAL *src, gint sstride,
					    WBFMM_REAL *normals, gint nstr,
					    WBFMM_REAL *dipoles, gint dstr,
					    gint nsrc,
					    WBFMM_REAL *xf, WBFMM_REAL *field)

{
  gint i ;
  WBFMM_REAL r, th, ph, h0[2], h1[2], fR[2], fd[6] ;

  /* field[0] = field[1] = 0.0 ; */

  if ( src == NULL && normals == NULL && dipoles == NULL ) return 0 ;

  if ( normals != NULL && dipoles == NULL )
    g_error("%s: normals specified but no dipole strengths (dipoles == NULL)",
	    __FUNCTION__) ;

  if ( normals == NULL && dipoles == NULL ) {
    for ( i = 0 ; i < nsrc ; i ++ ) {
      WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(&(xs[i*xstride]), xf, 
						  &r, &th, &ph) ;
      WBFMM_FUNCTION_NAME(wbfmm_bessel_h_init)(k*r, h0, h1) ;
      field[0] += h0[0]*src[i*sstride+0] - h0[1]*src[i*sstride+1] ;
      field[1] += h0[1]*src[i*sstride+0] + h0[0]*src[i*sstride+1] ;
    }
    
    /*G&D normalization of Legendre polynomials*/
    field[0] /= 4.0*M_PI ; field[1] /= 4.0*M_PI ;

    return 0 ;
  }

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

gint WBFMM_FUNCTION_NAME(wbfmm_coordinate_transform)(WBFMM_REAL *x, 
						     WBFMM_REAL *ix,
						     WBFMM_REAL *iy,
						     WBFMM_REAL *iz,
						     WBFMM_REAL *y)

/*
  transform point x to coordinate system (ix,iy,iz) so that
  y = (x.ix, x.iy, x.iz)
*/

{
  y[0] = x[0]*ix[0] + x[1]*ix[1] + x[2]*ix[2] ;
  y[1] = x[0]*iy[0] + x[1]*iy[1] + x[2]*iy[2] ;
  y[2] = x[0]*iz[0] + x[1]*iz[1] + x[2]*iz[2] ;

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_shift_coordinates)(WBFMM_REAL *x, WBFMM_REAL *y,
						  WBFMM_REAL *ix,
						  WBFMM_REAL *iy,
						  WBFMM_REAL *iz,
						  WBFMM_REAL *r)
/*
  coordinate system for shift from x to y, with z axis = (y-x)/|y-x|
*/

{
  WBFMM_REAL l ;

  iz[0] = y[0] - x[0] ; iz[1] = y[1] - x[1] ; iz[2] = y[2] - x[2] ; 

  *r = SQRT(iz[0]*iz[0] + iz[1]*iz[1] + iz[2]*iz[2]) ;

  g_assert((*r) != 0.0) ;

  iz[0] /= (*r) ; iz[1] /= (*r) ; iz[2] /= (*r) ;

  /* if ( iz[0] != 0.0 ) { */
  /*   ix[0] = iz[0] ; ix[1] = -iz[1] ; ix[2] = 0.0 ;  */
  /* } else { */
  /*   if ( iz[1] != 0.0 ) { */
  /*     ix[0] = 1.0 ; ix[1] = iz[1] ; ix[2] = -iz[2] ; */
  /*   } else { */
  /*     ix[0] = 1.0 ; ix[1] = 0.0 ; ix[2] = 0.0 ; */
  /*   } */
  /* } */

  if ( iz[0] != 0.0 ) {
    ix[0] = iz[1] ; ix[1] = -iz[0] ; ix[2] = 0.0 ; 
  } else {
    if ( iz[1] != 0.0 ) {
      ix[0] = 1.0 ; ix[1] = iz[2] ; ix[2] = -iz[1] ;
    } else {
      ix[0] = 1.0 ; ix[1] = 0.0 ; ix[2] = 0.0 ;
    }
  }

  l = SQRT(ix[0]*ix[0] + ix[1]*ix[1] + ix[2]*ix[2]) ;
  ix[0] /= l ; ix[1] /= l ; ix[2] /= l ;

  iy[0] = iz[1]*ix[2] - iz[2]*ix[1] ;
  iy[1] = iz[2]*ix[0] - iz[0]*ix[2] ;
  iy[2] = iz[0]*ix[1] - iz[1]*ix[0] ;

  return 0 ;
}

gint print_bits_uint(FILE *f, guint x)

{
  gint i ;
  guint size = sizeof(guint) ;
  guint maxPow = 1 << (size*8 - 1) ;

  for( i = 0 ; i < 8*size ; ++i ) {
    /* print last bit and shift left. */
    fprintf(f, "%u", (x & maxPow ? 1 : 0)) ;
    x = x << 1 ;
}

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_box_location_from_index)(guint64 idx,
							guint32 level,
							WBFMM_REAL *x0,
							WBFMM_REAL D,
							WBFMM_REAL *x,
							WBFMM_REAL *wb)

{
  guint nb ;
  guint32 i, j, k ;

  nb = 1 << level ;
  if ( !(idx < (1 << 3*level)) ) {
    g_error("%s: box %lu is not at level %u", __FUNCTION__, idx, level) ;
  }

  *wb = D/nb ;

  wbfmm_box_location(idx, &i, &j, &k) ;

  x[0] = x0[0] + (*wb)*i ;
  x[1] = x0[1] + (*wb)*j ;
  x[2] = x0[2] + (*wb)*k ;

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_tree_box_centre)(wbfmm_tree_t *t, guint32 level,
						guint64 b, WBFMM_REAL *xb,
						WBFMM_REAL *wb)

{
  WBFMM_FUNCTION_NAME(wbfmm_box_location_from_index)(b, level, 
					       wbfmm_tree_origin(t), 
					       wbfmm_tree_width(t), xb, 
					       wb) ;
  xb[0] += 0.5*(*wb) ; xb[1] += 0.5*(*wb) ; xb[2] += 0.5*(*wb) ; 

  return 0 ;
}

static gint print_box(FILE *f, wbfmm_tree_t *t, guint64 idx, 
		      guint level, wbfmm_box_t b, gboolean print_empty)

{
  WBFMM_REAL x[3], w ;
  guint i ;

  if ( !print_empty ) {
    if ( b.n == 0 ) return 0 ;
  }

  WBFMM_FUNCTION_NAME(wbfmm_box_location_from_index)(idx, level, 
					       wbfmm_tree_origin(t),
					       wbfmm_tree_width(t),
					       x, &w) ;
    
  fprintf(f, "%lg %lg %lg %lg %u", x[0], x[1], x[2], w, b.n) ;
  for ( i = 0 ; i < b.n ; i ++ ) {
    fprintf(f, " %u", t->ip[b.i+i]) ;
  }
  fprintf(f, "\n") ;
  
  return 0 ;
}

gint wbfmm_tree_print(FILE *f, wbfmm_tree_t *t, guint level, 
		      gboolean print_empty)

{
  guint i, nb, lv ;
  guint64 idx ;
  WBFMM_REAL *x ;
  wbfmm_box_t *boxes ;

  fprintf(f, "points %u\n", wbfmm_tree_point_number(t)) ;

  /*print points in original list order*/
  for ( i = 0 ; i < wbfmm_tree_point_number(t) ; i ++ ) {
    x = wbfmm_tree_point_index(t, i) ;
    fprintf(f, "%u %lg %lg %lg\n", i, x[0], x[1], x[2]) ;
  }

  /*print leaf box data: position, width, point list*/
  lv = wbfmm_tree_depth(t) - 1 ;
  nb = 1 << (3*lv) ;  
  boxes = t->boxes[lv] ;

  for ( idx = 0 ; idx < nb ; idx ++ ) 
    print_box(f, t, idx, lv, boxes[idx], print_empty) ;

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_points_origin_width)(WBFMM_REAL *x,
						    gint str, gint n,
						    WBFMM_REAL *xmin,
						    WBFMM_REAL *xmax,
						    WBFMM_REAL *D,
						    gboolean init_limits)

{
  /* WBFMM_REAL xmin[3] = {G_MAXFLOAT, G_MAXFLOAT, G_MAXFLOAT} ; */
  gint i ;

  if ( init_limits ) {
    xmin[0] =  G_MAXFLOAT ; xmin[1] =  G_MAXFLOAT ; xmin[2] =  G_MAXFLOAT ;
    xmax[0] = -G_MAXFLOAT ; xmax[0] = -G_MAXFLOAT ; xmax[0] = -G_MAXFLOAT ;
  }
  
  for ( i = 0 ; i < n ; i ++ ) {
    xmin[0] = MIN(xmin[0], x[str*i+0]) ;
    xmin[1] = MIN(xmin[1], x[str*i+1]) ;
    xmin[2] = MIN(xmin[2], x[str*i+2]) ;
    xmax[0] = MAX(xmax[0], x[str*i+0]) ;
    xmax[1] = MAX(xmax[1], x[str*i+1]) ;
    xmax[2] = MAX(xmax[2], x[str*i+2]) ;
  }

  *D = xmax[0] - xmin[0] ;
  *D = MAX(*D, xmax[1] - xmin[1]) ;
  *D = MAX(*D, xmax[2] - xmin[2]) ;

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_shift_angles)(WBFMM_REAL *xi, WBFMM_REAL *xj,
					     WBFMM_REAL *th, WBFMM_REAL *ph,
					     WBFMM_REAL *ch, WBFMM_REAL *r)

{
  WBFMM_REAL ix0[3], iy0[3], iz0[3], ix[3], iy[3], iz[3] ;

  ix0[0] = 1.0 ; ix0[1] = 0.0 ; ix0[2] = 0.0 ;
  iy0[0] = 0.0 ; iy0[1] = 1.0 ; iy0[2] = 0.0 ;
  iz0[0] = 0.0 ; iz0[1] = 0.0 ; iz0[2] = 1.0 ;
  
  WBFMM_FUNCTION_NAME(wbfmm_shift_coordinates)(xi, xj, ix, iy, iz, r) ;
  WBFMM_FUNCTION_NAME(wbfmm_rotation_angles)(ix0, iy0, iz0, ix, iy, iz, 
					     th, ph, ch) ;

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_tree_write_sources)(wbfmm_tree_t *t,
						   WBFMM_REAL *q, gint stride,
						   FILE *f)

{
  guint i, idx ;
  WBFMM_REAL *x ;

  if ( q == NULL ) {
    for ( i = 0 ; i < t->npoints ; i ++ ) {
      idx = t->ip[i] ;
      x = wbfmm_tree_point_index(t, idx) ;
      fprintf(f, "%u %u %g %g %g\n", i, idx, x[0], x[1], x[2]) ;
    }

    return 0 ;
  }

  for ( i = 0 ; i < t->npoints ; i ++ ) {
    idx = t->ip[i] ;
    x = wbfmm_tree_point_index(t, idx) ;
    fprintf(f, "%u %u %g %g %g\n", i, idx, x[0], x[1], x[2]) ;
  }

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_rotation_write_coefficients)(WBFMM_REAL *H,
							    gint N, FILE *f)

{
  gint n, m, nu, idx ;

  for ( n = 0 ; n <= N ; n ++ ) {
    for ( nu = 0 ; nu <= n ; nu ++ ) {
      for ( m = 0 ; m <= n ; m ++ ) {
	idx = wbfmm_rotation_index_numn( nu,m,n) ;
#ifdef WBFMM_SINGLE_PRECISION
	fprintf(f, "%d %d %d %g\n", nu, m, n, H[idx]) ;
#else
	fprintf(f, "%d %d %d %lg\n", nu, m, n, H[idx]) ;
#endif /*WBFMM_SINGLE_PRECISION*/
      }
    }
  }
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_truncation_number)(wbfmm_tree_t *t,
						  WBFMM_REAL k, guint level,
						  WBFMM_REAL tol)

/* Gumerov and Duraiswami, J. Acoust. Soc. Am. 2009, p191 */
  
{
  gint p ;
  WBFMM_REAL D, ka, phi, plo, del ;
  guint nb ;

  del = 2.0 ;
  nb = 1 << level ;

  D = wbfmm_tree_width(t) ;
  ka = D/nb ;

  /*radius of circumsphere at level `level'*/
  ka *= 0.5*SQRT(3.0)*k ;

  plo = 1.0 - LOG(tol*(1.0-1.0/del)*SQRT(1.0-1.0/del))/LOG(del) ;

  phi = CBRT(3.0*LOG(1.0/tol)) ;
  phi *= phi ;
  phi = ka + 0.5*phi*CBRT(ka) ;

  p = (gint)ceil(SQRT(SQRT(plo*plo*plo*plo + phi*phi*phi*phi))) ;
  
  return p ;
}

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

  return 0 ;
}

gint wbfmm_library_config_print(wbfmm_library_config_t *c, FILE *f)

{
  fprintf(f,
	  "AVX extensions: %s\n"
	  "FMA extensions: %s\n"
	  "precision:      %s\n"
	  "compiler flags: %s\n",
	  yes_if_true(c->avx),
	  yes_if_true(c->fma),
#ifdef WBFMM_SINGLE_PRECISION
	  "single",
#else
	  "double",
#endif /*WBFMM_SINGLE_PRECISION*/
	  c->switches
	  ) ;
  
  return 0 ;
}
