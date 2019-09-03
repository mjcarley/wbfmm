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

#include <glib.h>

#include <wbfmm.h>

#include "wbfmm-private.h"

#define wbfmm_index_laplace_nm(_n,_m) ((_n)*(_n)+(2*(_m))-1)

/*table of \cos m\pi/4 for rotations on upward pass*/
WBFMM_REAL CmPI_4[] =
  {1, M_SQRT1_2, 0, -M_SQRT1_2, -1, -M_SQRT1_2, 0, M_SQRT1_2, 1} ;
/*table of \cos n\pi/2 for rotations*/
WBFMM_REAL CnPI_2[] = {1.0, 0.0, -1.0, 0.0} ;

#define cos_n_PI_4(_n) (CmPI_4[(_n)%8])
#define sin_n_PI_4(_n) (CmPI_4[((_n)+6)%8])
#define cos_n_PI_2(_n) (CnPI_2[(_n)%4])
#define sin_n_PI_2(_n) (CnPI_2[((_n)+3)%4])

gint WBFMM_FUNCTION_NAME(wbfmm_expansion_laplace_cfft)(gint N,
						       WBFMM_REAL *x0,
						       WBFMM_REAL *xs,
						       WBFMM_REAL *q, gint nq,
						       WBFMM_REAL *cfft,
						       gint cstr,
						       WBFMM_REAL *work)

/*
  generate expansion coefficients for Laplace equation

  inputs

  N:    order of expansion;
  x0:   origin;
  xs:   source location;
  q:    source (real);
  nq:   number of components in q; 
  cfft: expansion coefficients (not zeroed internally to allow accumulation);
  cstr: stride in cfft;
  work: workspace

  packing in coefficient array indexed by
  idx=wbfmm_index_laplace_nm(n,m) for m not equal to zero, and idx=n^2
  for m == 0 (since these are real);

  it is assumed that coefficients are densely packed from cfft[idx*cstr]

  a check is performed to ensure that cstr >= nq
*/
  
{
  WBFMM_REAL r, th, ph, rn ;
  WBFMM_REAL Cth, Sth, *Pn, *Pnm1, Cph, Sph, Cmph[64], Smph[64] ;
  gint n, m, idx, i ;

  g_assert(cstr >= nq) ;
  
  Pnm1 = &(work[0]) ; Pn = &(Pnm1[2*(2*N+1)]) ;

  WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(x0, xs, &r, &th, &ph) ;
  Cth = COS(th) ; Sth = SIN(th) ;
  Cph = COS(ph) ; Sph = SIN(ph) ; 

  WBFMM_FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth,
					   &(Pnm1[0]), &(Pn[0]), &(Pn[1])) ;
  Cmph[0] = 1.0 ; Smph[0] = 0.0 ;
  Cmph[1] = Cph ; Smph[1] = Sph ;
  
  n = 0 ; 
  m = 0 ;
  rn = 1.0 ;
  idx = n*n ;
  for ( i = 0 ; i < nq ; i ++ )
    cfft[cstr*idx+i] += q[i]*rn*Pnm1[0]/(2*n+1) ;
  
  n = 1 ; 
  m = 0 ;
  rn = r ;
  idx = n*n ;
  for ( i = 0 ; i < nq ; i ++ )
    cfft[cstr*idx+i] += q[i]*rn*Pn[0]/(2*n+1) ;
  
  m = 1 ; 
  idx = wbfmm_index_laplace_nm(n,m) ;
  for ( i = 0 ; i < nq ; i ++ ) {
    /* cfft[cstr*idx+2*i+0] += q[i]*rn*Pn[m]*Cmph[m]/(2*n+1) ; */
    /* cfft[cstr*idx+2*i+1] -= q[i]*rn*Pn[m]*Smph[m]/(2*n+1) ; */
    cfft[cstr*(idx+0)+i] += q[i]*rn*Pn[m]*Cmph[m]/(2*n+1) ;
    cfft[cstr*(idx+1)+i] -= q[i]*rn*Pn[m]*Smph[m]/(2*n+1) ;
  }
  
  for ( n = 2 ; n <= N ; n ++ ) {
    rn *= r ;
    WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pnm1, &Pn,
							n-1, Cth, Sth) ;
    Cmph[n] = Cmph[n-1]*Cph - Smph[n-1]*Sph ;
    Smph[n] = Smph[n-1]*Cph + Cmph[n-1]*Sph ;
    
    m = 0 ; 
    idx = n*n ;
    for ( i = 0 ; i < nq ; i ++ )
      cfft[cstr*idx+i] += q[i]*rn*Pn[0]/(2*n+1) ;
    
    for ( m = 1 ; m <= n ; m ++ ) {
      idx = wbfmm_index_laplace_nm(n,m) ;
      for ( i = 0 ; i < nq ; i ++ ) {
	/* cfft[cstr*idx+2*i+0] += q[i]*rn*Pn[m]*Cmph[m]/(2*n+1) ; */
	/* cfft[cstr*idx+2*i+1] -= q[i]*rn*Pn[m]*Smph[m]/(2*n+1) ; */
	cfft[cstr*(idx+0)+i] += q[i]*rn*Pn[m]*Cmph[m]/(2*n+1) ;
	cfft[cstr*(idx+1)+i] -= q[i]*rn*Pn[m]*Smph[m]/(2*n+1) ;
      }
    }
  }

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_expansion_laplace_evaluate)(WBFMM_REAL *x0,
							   WBFMM_REAL *cfft,
							   gint cstr, 
							   gint N, gint nq,
							   WBFMM_REAL *xf,
							   WBFMM_REAL *field,
							   WBFMM_REAL *work)

{
  WBFMM_REAL r, th, ph, rn ;
  WBFMM_REAL Cth, Sth, *Pn, *Pnm1, Cph, Sph, Cmph[64], Smph[64] ;
  gint n, m, idx, i ;

  Pnm1 = &(work[0]) ; Pn = &(Pnm1[2*(2*N+1)]) ;

  WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(x0, xf, &r, &th, &ph) ;
  Cth = COS(th) ; Sth = SIN(th) ; 
  Cph = COS(ph) ; Sph = SIN(ph) ; 

  /*initialize recursions*/
  WBFMM_FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth,
					   &(Pnm1[0]), &(Pn[0]), &(Pn[1])) ;
  Cmph[0] = 1.0 ; Smph[0] = 0.0 ;
  Cmph[1] = Cph ; Smph[1] = Sph ;

  /*first two terms by hand*/
  n = 0 ; 
  m = 0 ;
  rn = 1.0/r ;
  idx = n*n ;
  for ( i = 0 ; i < nq ; i ++ ) field[i] += cfft[cstr*idx+i]*rn*Pnm1[m] ;

  n = 1 ; 
  m = 0 ; 
  rn /= r ;
  idx = n*n ;
  for ( i = 0 ; i < nq ; i ++ ) field[i] += cfft[cstr*idx+i]*rn*Pn[m] ;

  m = 1 ; 
  idx = wbfmm_index_laplace_nm(n,m) ;
  for ( i = 0 ; i < nq ; i ++ ) {
    /* field[i] += 2.0*Pn[m]*rn*(cfft[cstr*idx+2*i+0]*Cmph[m] - */
    /* 			      cfft[cstr*idx+2*i+1]*Smph[m]) ; */
    field[i] += 2.0*Pn[m]*rn*(cfft[cstr*(idx+0)+i]*Cmph[m] -
			      cfft[cstr*(idx+1)+i]*Smph[m]) ;
  }

  for ( n = 2 ; n <= N ; n ++ ) {
    rn /= r ;
    WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pnm1, &Pn,
							n-1, Cth, Sth) ;
    Cmph[n] = Cmph[n-1]*Cph - Smph[n-1]*Sph ;
    Smph[n] = Smph[n-1]*Cph + Cmph[n-1]*Sph ;

    m = 0 ; 
    idx = n*n ;
    for ( i = 0 ; i < nq ; i ++ ) field[i] += cfft[cstr*idx+i]*rn*Pn[0] ;
    
    for ( m = 1 ; m <= n ; m ++ ) {
      idx = wbfmm_index_laplace_nm(n,m) ;
      for ( i = 0 ; i < nq ; i ++ ) {
	/* field[i] += 2.0*Pn[m]*rn*(cfft[cstr*idx+2*i+0]*Cmph[m] - */
	/* 			  cfft[cstr*idx+2*i+1]*Smph[m]) ; */
	field[i] += 2.0*Pn[m]*rn*(cfft[cstr*(idx+0)+i]*Cmph[m] -
				  cfft[cstr*(idx+1)+i]*Smph[m]) ;
      }
    }
  }
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_expansion_laplace_local_evaluate)(WBFMM_REAL *x0,
								 WBFMM_REAL
								 *cfft,
								 gint cstr, 
								 gint N,
								 gint nq,
								 WBFMM_REAL *xf,
								 WBFMM_REAL
								 *field,
								 WBFMM_REAL
								 *work)

{
  WBFMM_REAL r, th, ph, rn ;
  WBFMM_REAL Cth, Sth, *Pn, *Pnm1, Cph, Sph, Cmph[64], Smph[64] ;
  gint n, m, idx, i ;

  Pnm1 = &(work[0]) ; Pn = &(Pnm1[2*(2*N+1)]) ;

  WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(x0, xf, &r, &th, &ph) ;
  Cth = COS(th) ; Sth = SIN(th) ; 
  Cph = COS(ph) ; Sph = SIN(ph) ; 

  /*initialize recursions*/
  WBFMM_FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth,
					   &(Pnm1[0]), &(Pn[0]), &(Pn[1])) ;
  Cmph[0] = 1.0 ; Smph[0] = 0.0 ;
  Cmph[1] = Cph ; Smph[1] = Sph ;

  /*first two terms by hand*/
  n = 0 ; 
  m = 0 ;
  rn = 1.0 ;
  idx = n*n ;
  for ( i = 0 ; i < nq ; i ++ ) field[i] += cfft[cstr*idx+i]*rn*Pnm1[m] ;

  /* fprintf(stderr, "%g\n", field[0]) ; */
  
  n = 1 ; 
  m = 0 ; 
  rn *= r ;
  idx = n*n ;
  for ( i = 0 ; i < nq ; i ++ ) field[i] += cfft[cstr*idx+i]*rn*Pn[m] ;

  /* fprintf(stderr, "%g\n", field[0]) ; */
  m = 1 ; 
  idx = wbfmm_index_laplace_nm(n,m) ;
  for ( i = 0 ; i < nq ; i ++ ) {
    /* field[i] += 2.0*Pn[m]*rn*(cfft[cstr*idx+2*i+0]*Cmph[m] - */
    /* 			      cfft[cstr*idx+2*i+1]*Smph[m]) ; */
    field[i] += 2.0*Pn[m]*rn*(cfft[cstr*(idx+0)+i]*Cmph[m] -
			      cfft[cstr*(idx+1)+i]*Smph[m]) ;
  }

  for ( n = 2 ; n <= N ; n ++ ) {
    rn *= r ;
    WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pnm1, &Pn,
							n-1, Cth, Sth) ;
    Cmph[n] = Cmph[n-1]*Cph - Smph[n-1]*Sph ;
    Smph[n] = Smph[n-1]*Cph + Cmph[n-1]*Sph ;

    m = 0 ; 
    idx = n*n ;
    for ( i = 0 ; i < nq ; i ++ ) field[i] += cfft[cstr*idx+i]*rn*Pn[0] ;
    
    for ( m = 1 ; m <= n ; m ++ ) {
      idx = wbfmm_index_laplace_nm(n,m) ;
      for ( i = 0 ; i < nq ; i ++ ) {
	/* field[i] += 2.0*Pn[m]*rn*(cfft[cstr*idx+2*i+0]*Cmph[m] - */
	/* 			  cfft[cstr*idx+2*i+1]*Smph[m]) ; */
	field[i] += 2.0*Pn[m]*rn*(cfft[cstr*(idx+0)+i]*Cmph[m] -
				  cfft[cstr*(idx+1)+i]*Smph[m]) ;
      }
    }
  /* fprintf(stderr, "%g\n", field[0]) ; */
  }
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_laplace_field)(WBFMM_REAL *xs, gint xstride,
					      WBFMM_REAL *src, gint sstride,
					      gint nq,
					      WBFMM_REAL *normals, gint nstr,
					      WBFMM_REAL *dipoles, gint dstr,
					      gint nsrc,
					      WBFMM_REAL *xf, WBFMM_REAL *field)

{
  gint i, j ;
  WBFMM_REAL r, th, ph, fR[2], fd[6] ;

  if ( src == NULL && normals == NULL && dipoles == NULL ) return 0 ;

  g_assert(sstride >= nq) ;
  
  if ( normals != NULL && dipoles == NULL )
    g_error("%s: normals specified but no dipole strengths (dipoles == NULL)",
	    __FUNCTION__) ;

  if ( normals == NULL && dipoles == NULL ) {
    for ( i = 0 ; i < nsrc ; i ++ ) {
      WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(&(xs[i*xstride]), xf, 
							&r, &th, &ph) ;
      for ( j = 0 ; j < nq ; j ++ ) field[j] += src[i*sstride+j]/r ;
    }
    
    /*G&D normalization of Legendre polynomials*/
    for ( j = 0 ; j < nq ; j ++ ) field[j] /= 4.0*M_PI ;

    return 0 ;
  }

  g_assert_not_reached() ;
  
  return 0 ;
}

static WBFMM_REAL coaxial_translation_SS_cfft(gint n, gint nd, gint m)

{
  WBFMM_REAL c ;

  g_assert(n >= nd) ;
  
  c  = wbfmm_factorial(n  - m)*wbfmm_factorial(n  + m) ;
  c /= wbfmm_factorial(nd - m)*wbfmm_factorial(nd + m) ;
  c *= (2.0*nd+1)/(2.0*n+1) ;
  c  = SQRT(c) ;
  c /= wbfmm_factorial(n-nd) ;

  /*this needs to be multiplied by (-t)^{n-nd} for a translation of
    distance t*/
  
  return c ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_coaxial_translate_SS_laplace)(WBFMM_REAL *Co,
							     gint cstro,
							     gint No,
							     WBFMM_REAL *Ci,
							     gint cstri,
							     gint Ni,
							     gint nq,
							     WBFMM_REAL t)

{
  gint n, m, nd, idxo, idxi, i ;
  WBFMM_REAL c, tn[64] ;

  tn[0] = 1.0 ;
  m = 0 ; 
  for ( idxo = n = 0 ; n <= No ; (n ++), (idxo = n*n) ) {
    for ( idxi = nd = 0 ; nd <= MIN(n,Ni) ; (nd ++), (idxi = nd*nd) ) {
      c = coaxial_translation_SS_cfft(n, nd, m)*tn[n-nd] ;
      for ( i = 0 ; i < nq ; i ++ ) 
	Co[cstro*idxo+i] += c*Ci[cstri*idxi+i] ;
    }
    tn[n+1] = -tn[n]*t ;
  }
  
  for ( m = 1 ; m <= No ; m ++ ) {
    for ( n = m ; n <= No ; n ++ ) {
      idxo = wbfmm_index_laplace_nm(n,m) ;
      for ( nd = m ; nd <= MIN(n, Ni) ; nd ++ ) {
  	idxi = wbfmm_index_laplace_nm(nd,m) ;
	c = coaxial_translation_SS_cfft(n, nd, m)*tn[n-nd] ;
	for ( i = 0 ; i < nq ; i ++ ) {
	  /* Co[cstro*idxo+2*i+0] += c*Ci[cstri*idxi+2*i+0] ; */
	  /* Co[cstro*idxo+2*i+1] += c*Ci[cstri*idxi+2*i+1] ; */
	  Co[cstro*(idxo+0)+i] += c*Ci[cstri*(idxi+0)+i] ;
	  Co[cstro*(idxo+1)+i] += c*Ci[cstri*(idxi+1)+i] ;
	}
      }
    }
  }
  
  return 0 ;
}

static WBFMM_REAL coaxial_translation_RR_cfft(gint n, gint nd, gint m)

{
  WBFMM_REAL c ;

  g_assert(nd >= n) ;
  
  c  = wbfmm_factorial(nd  - m)*wbfmm_factorial(nd  + m) ;
  c /= wbfmm_factorial(n - m)*wbfmm_factorial(n + m) ;
  c *= (2.0*nd+1)/(2.0*n+1) ;
  c  = SQRT(c) ;
  c /= wbfmm_factorial(nd-n) ;

  /*this needs to be multiplied by (t)^{nd-n} for a translation of
    distance t*/
  
  return c ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_coaxial_translate_RR_laplace)(WBFMM_REAL *Co,
							     gint cstro,
							     gint No,
							     WBFMM_REAL *Ci,
							     gint cstri,
							     gint Ni,
							     gint nq,
							     WBFMM_REAL t)

{
  gint n, m, nd, idxo, idxi, i ;
  WBFMM_REAL c, tn[64] ;

  tn[0] = 1.0 ;
  m  = 0 ; 
  for ( idxi = nd = 0 ; nd <= Ni ; (nd ++), (idxi = nd*nd) ) {
    for ( idxo = n = 0 ; n <= nd ; (n ++), (idxo = n*n) ) {
      c = coaxial_translation_RR_cfft(n, nd, m)*tn[nd-n] ;
      for ( i = 0 ; i < nq ; i ++ ) 
	Co[cstro*idxo+i] += c*Ci[cstri*idxi+i] ;
    }
    tn[nd+1] = tn[nd]*t ;
  }
  
  for ( m = 1 ; m <= No ; m ++ ) {
    for ( nd = m ; nd <= No ; nd ++ ) {
      idxi = wbfmm_index_laplace_nm(nd,m) ;
      for ( n = m ; n <= nd ; n ++ ) {
  	idxo = wbfmm_index_laplace_nm(n,m) ;
	c = coaxial_translation_RR_cfft(n, nd, m)*tn[nd-n] ;
	for ( i = 0 ; i < nq ; i ++ ) {
	  /* Co[cstro*idxo+2*i+0] += c*Ci[cstri*idxi+2*i+0] ; */
	  /* Co[cstro*idxo+2*i+1] += c*Ci[cstri*idxi+2*i+1] ; */
  	  Co[cstro*(idxo+0)+i] += c*Ci[cstri*(idxi+0)+i] ;
  	  Co[cstro*(idxo+1)+i] += c*Ci[cstri*(idxi+1)+i] ;
	}
      }
    }
  }
  
  return 0 ;
}

static WBFMM_REAL coaxial_translation_SR_cfft(gint n, gint nd, gint m)

{
  WBFMM_REAL c ;

  c  = wbfmm_factorial(n  - m)*wbfmm_factorial(n  + m) ;
  c *= wbfmm_factorial(nd - m)*wbfmm_factorial(nd + m) ;
  c  = (2.0*nd+1)/(2.0*n+1)/c ;
  c  = SQRT(c) ;
  c *= wbfmm_factorial(n+nd) ;
  c *= pow(-1.0, n+m) ;
  
  /*this needs to be divided by (t)^{n+nd+1} for a translation of
    distance t*/
  
  return c ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_coaxial_translate_SR_laplace)(WBFMM_REAL *Co,
							     gint cstro,
							     gint No,
							     WBFMM_REAL *Ci,
							     gint cstri,
							     gint Ni,
							     gint nq,
							     WBFMM_REAL t)

{
  gint n, m, nd, idxo, idxi, i ;
  WBFMM_REAL c, tn[64] ;

  tn[0] = 1.0 ;
  for ( n = 1 ; n <= No+Ni+2 ; n ++ ) tn[n] = t*tn[n-1] ;

  m = 0 ;
  for ( idxo = n = 0 ; n <= No ; (n ++), (idxo = n*n) ) {
    for ( idxi = nd = 0 ; nd <= Ni ; (nd ++), (idxi = nd*nd) ) {
      c = coaxial_translation_SR_cfft(n, nd, m)/tn[n+nd+1] ;
      for ( i = 0 ; i < nq ; i ++ )
  	Co[cstro*idxo+i] += c*Ci[cstri*idxi+i] ;
    }
  }

  for ( m = 1 ; m <= No ; m ++ ) {
    for ( n = m ; n <= No ; n ++ ) {
      idxo = wbfmm_index_laplace_nm(n,m) ;
      for ( nd = m ; nd <= Ni ; nd ++ ) {
	idxi = wbfmm_index_laplace_nm(nd,m) ;
  	c = coaxial_translation_SR_cfft(n, nd, m)/tn[n+nd+1] ;
  	for ( i = 0 ; i < nq ; i ++ ) {
  	  /* Co[cstro*idxo+2*i+0] += c*Ci[cstri*idxi+2*i+0] ; */
  	  /* Co[cstro*idxo+2*i+1] += c*Ci[cstri*idxi+2*i+1] ; */
  	  Co[cstro*(idxo+0)+i] += c*Ci[cstri*(idxi+0)+i] ;
  	  Co[cstro*(idxo+1)+i] += c*Ci[cstri*(idxi+1)+i] ;
  	}	
      }
    }
  }
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_rotate_H_laplace)(WBFMM_REAL *Co, gint cstro,
						 WBFMM_REAL *Ci, gint cstri,
						 gint N, gint nq,
						 WBFMM_REAL *H,
						 WBFMM_REAL ph, WBFMM_REAL ch)

{
  gint n, m, nu, idxi, idxo, i ;
  WBFMM_REAL Cmch, Smch, Cnph, Snph, Cch, Sch, Cph, Sph, tr[32] ;
  WBFMM_REAL Hp, Hm ;
  gint j ;
  
  g_assert(nq <= cstri) ;
  g_assert(nq <= cstro) ;
  g_assert(nq <= 32) ;
  
  /*initialize recursions*/
  Cph = COS(ph) ; Sph = SIN(ph) ;
  Cch = COS(ch) ; Sch = SIN(ch) ;

  for ( n = 0 ; n <= N ; n ++ ) {
    nu = 0 ; idxo = n*n ;
    m  = 0 ; idxi = n*n ;
    Hp = H[wbfmm_rotation_index_numn(nu,m,n)] ;
    
    for ( i = 0 ; i < nq ; i ++ ) tr[i] = Hp*Ci[cstri*idxi+i] ;

    Cmch = 1.0 ; Smch = 0.0 ;

    for ( m = 1 ; m <= n ; m ++ ) {
      wbfmm_cos_sin_recursion(Cmch,Smch,Cch,Sch) ;

      idxi = wbfmm_index_laplace_nm(n,m) ;
      Hp = H[wbfmm_rotation_index_numn( nu,m,n)] ;
      Hm = H[wbfmm_rotation_index_numn(-nu,m,n)] ;
      for ( i = 0 ; i < nq ; i ++ ) {
	/* tr[i] += (Ci[cstri*idxi+2*i+0]*Cmch - */
	/* 	  Ci[cstri*idxi+2*i+1]*Smch)*(Hp+Hm) ; */
	tr[i] += (Ci[cstri*(idxi+0)+i]*Cmch -
		  Ci[cstri*(idxi+1)+i]*Smch)*(Hp+Hm) ;
      }
    }

    for ( i = 0 ; i < nq ; i ++ ) Co[cstro*idxo+i] += tr[i] ;

    Cnph = 1.0 ; Snph = 0.0 ;
    for ( nu = 1 ; nu <= n ; nu ++ ) {
      WBFMM_REAL ti[32] = {0.0} ;

      wbfmm_cos_sin_recursion(Cnph,Snph,Cph,Sph) ;

      idxo = wbfmm_index_laplace_nm(n,nu) ;
      m = 0 ; idxi = n*n ;
      Hp = H[wbfmm_rotation_index_numn(nu,m,n)] ;

      for ( i = 0 ; i < nq ; i ++ ) tr[i] = Hp*Ci[cstri*idxi+i] ;
      
      Cmch = 1.0 ; Smch = 0.0 ;
      for ( m = 1 ; m <= n ; m ++ ) {
	wbfmm_cos_sin_recursion(Cmch,Smch,Cch,Sch) ;

	idxi = wbfmm_index_laplace_nm(n,m) ;
	Hp = H[wbfmm_rotation_index_numn( nu,m,n)] ;
	Hm = H[wbfmm_rotation_index_numn(-nu,m,n)] ;
	for ( i = 0 ; i < nq ; i ++ ) {
	  /* tr[i] += (Cmch*Ci[cstri*idxi+2*i+0] - */
	  /* 	    Smch*Ci[cstri*idxi+2*i+1])*(Hp+Hm) ; */
	  /* ti[i] += (Smch*Ci[cstri*idxi+2*i+0] + */
	  /* 	    Cmch*Ci[cstri*idxi+2*i+1])*(Hp-Hm) ; */
	  tr[i] += (Ci[cstri*(idxi+0)+i]*Cmch -
		    Ci[cstri*(idxi+1)+i]*Smch)*(Hp+Hm) ;
	  ti[i] += (Ci[cstri*(idxi+0)+i]*Smch +
		    Ci[cstri*(idxi+1)+i]*Cmch)*(Hp-Hm) ;
	}
      }
      for ( i = 0 ; i < nq ; i ++ ) {
	Co[cstro*(idxo+0)+i] += Cnph*tr[i] + Snph*ti[i] ;
	Co[cstro*(idxo+1)+i] -= Snph*tr[i] - Cnph*ti[i] ;
      }
    }
  }
  
  return 0 ;
}

static inline void increment_buf_cp_real(WBFMM_REAL *E0, WBFMM_REAL *E1,
					 WBFMM_REAL *C,
					 WBFMM_REAL H03p, WBFMM_REAL H03m,
					 WBFMM_REAL H47p, WBFMM_REAL H47m,
					 gint nq,
					 WBFMM_REAL *buf)

{
  WBFMM_REAL H03, H47 ;
  gint i ;
  
  H03 = H03p + H03m ; H47 = H47p + H47m ;
  /*terms using \cos (\pm m\pi/4), ...*/
  for ( i = 0 ; i < nq ; i ++ ) {
    buf[0*nq+i] += (C[0*nq+i]*E0[0] - C[8*nq+0*nq+i]*E0[1])*H03 ;
    buf[1*nq+i] += (C[1*nq+i]*E1[0] - C[8*nq+1*nq+i]*E1[1])*H03 ;
    buf[2*nq+i] += (C[2*nq+i]*E0[0] + C[8*nq+2*nq+i]*E0[1])*H03 ;
    buf[3*nq+i] += (C[3*nq+i]*E1[0] + C[8*nq+3*nq+i]*E1[1])*H03 ;

    buf[4*nq+i] += (C[4*nq+i]*E0[0] - C[8*nq+4*nq+i]*E0[1])*H47 ;
    buf[5*nq+i] += (C[5*nq+i]*E1[0] - C[8*nq+5*nq+i]*E1[1])*H47 ;
    buf[6*nq+i] += (C[6*nq+i]*E0[0] + C[8*nq+6*nq+i]*E0[1])*H47 ;
    buf[7*nq+i] += (C[7*nq+i]*E1[0] + C[8*nq+7*nq+i]*E1[1])*H47 ;

  }

  return ;
}

static inline void increment_buf_cp_complex(WBFMM_REAL *E0, WBFMM_REAL *E1,
					    WBFMM_REAL *C,
					    WBFMM_REAL H03p, WBFMM_REAL H03m,
					    WBFMM_REAL H47p, WBFMM_REAL H47m,
					    gint nq,
					    WBFMM_REAL *buf)

{
  WBFMM_REAL sH03, sH47, dH03, dH47 ;
  gint i ;
  
  sH03 = H03p + H03m ; sH47 = H47p + H47m ;
  dH03 = H03p - H03m ; dH47 = H47p - H47m ;
  for ( i = 0 ; i < nq ; i ++ ) {
    buf[0*2*nq+2*i+0] += (C[0*nq+i]*E0[0] - C[8*nq+0*nq+i]*E0[1])*sH03 ;
    buf[0*2*nq+2*i+1] += (C[0*nq+i]*E0[1] + C[8*nq+0*nq+i]*E0[0])*dH03 ;

    buf[1*2*nq+2*i+0] += (C[1*nq+i]*E1[0] - C[8*nq+1*nq+i]*E1[1])*sH03 ;
    buf[1*2*nq+2*i+1] += (C[1*nq+i]*E1[1] + C[8*nq+1*nq+i]*E1[0])*dH03 ;

    buf[2*2*nq+2*i+0] += (C[2*nq+i]*E0[0] + C[8*nq+2*nq+i]*E0[1])*sH03 ;
    buf[2*2*nq+2*i+1] += (C[2*nq+i]*E0[1] - C[8*nq+2*nq+i]*E0[0])*dH03 ;

    buf[3*2*nq+2*i+0] += (C[3*nq+i]*E1[0] + C[8*nq+3*nq+i]*E1[1])*sH03 ;
    buf[3*2*nq+2*i+1] -= (C[3*nq+i]*E1[1] - C[8*nq+3*nq+i]*E1[0])*dH03 ;

    buf[4*2*nq+2*i+0] += (C[4*nq+i]*E0[0] - C[8*nq+4*nq+i]*E0[1])*sH47 ;
    buf[4*2*nq+2*i+1] += (C[4*nq+i]*E0[1] + C[8*nq+4*nq+i]*E0[0])*dH47 ;

    buf[5*2*nq+2*i+0] += (C[5*nq+i]*E1[0] - C[8*nq+5*nq+i]*E1[1])*sH47 ;
    buf[5*2*nq+2*i+1] += (C[5*nq+i]*E1[1] + C[8*nq+5*nq+i]*E1[0])*dH47 ;

    buf[6*2*nq+2*i+0] += (C[6*nq+i]*E0[0] + C[8*nq+6*nq+i]*E0[1])*sH47 ;
    buf[6*2*nq+2*i+1] += (C[6*nq+i]*E0[1] - C[8*nq+6*nq+i]*E0[0])*dH47 ;

    buf[7*2*nq+2*i+0] += (C[7*nq+i]*E1[0] + C[8*nq+7*nq+i]*E1[1])*sH47 ;
    buf[7*2*nq+2*i+1] -= (C[7*nq+i]*E1[1] - C[8*nq+7*nq+i]*E1[0])*dH47 ;
  }

  return ;
}

static inline void increment_cfft_cp_real(WBFMM_REAL *E0, WBFMM_REAL *E1,
					  WBFMM_REAL *Ci,
					  WBFMM_REAL H03p, WBFMM_REAL H03m,
					  WBFMM_REAL H47p, WBFMM_REAL H47m,
					  gint nq,
					  WBFMM_REAL *Cp)

{
  WBFMM_REAL H03, H47 ;
  gint i ;
  
  H03 = H03p + H03m ; H47 = H47p + H47m ;
  for ( i = 0 ; i < nq ; i ++ ) {
    Cp[0*nq+i] += (Ci[0*nq+0*nq+i]*E0[0] - Ci[8*nq+0*nq+i]*E0[1])*H03 ;
    Cp[0*nq+i] += (Ci[0*nq+1*nq+i]*E1[0] - Ci[8*nq+1*nq+i]*E1[1])*H03 ;
    Cp[0*nq+i] += (Ci[0*nq+2*nq+i]*E1[0] + Ci[8*nq+2*nq+i]*E1[1])*H03 ;
    Cp[0*nq+i] += (Ci[0*nq+3*nq+i]*E0[0] - Ci[8*nq+3*nq+i]*E0[1])*H03 ;

    Cp[0*nq+i] += (Ci[0*nq+4*nq+i]*E0[0] - Ci[8*nq+4*nq+i]*E0[1])*H47 ;
    Cp[0*nq+i] += (Ci[0*nq+5*nq+i]*E1[0] - Ci[8*nq+5*nq+i]*E1[1])*H47 ;
    Cp[0*nq+i] += (Ci[0*nq+6*nq+i]*E1[0] + Ci[8*nq+6*nq+i]*E1[1])*H47 ;
    Cp[0*nq+i] += (Ci[0*nq+7*nq+i]*E0[0] - Ci[8*nq+7*nq+i]*E0[1])*H47 ;
  }

  return ;
}

static inline void increment_cfft_cp_complex(WBFMM_REAL *E0, WBFMM_REAL *E1,
					     WBFMM_REAL *E0ph, WBFMM_REAL *E1ph,
					     WBFMM_REAL *Ci,
					     WBFMM_REAL H03p, WBFMM_REAL H03m,
					     WBFMM_REAL H47p, WBFMM_REAL H47m,
					     gint nq,
					     WBFMM_REAL *Cp)

{
  WBFMM_REAL sH03, sH47, dH03, dH47, tr, ti ;
  gint i ;
  
  sH03 = H03p + H03m ; sH47 = H47p + H47m ;
  dH03 = H03p - H03m ; dH47 = H47p - H47m ;
  /*terms using \cos (\pm m\pi/4), ...*/
  for ( i = 0 ; i < nq ; i ++ ) {
    tr = (Ci[0*nq+0*nq+i]*E0[0] - Ci[8*nq+0*nq+i]*E0[1])*sH03 ;
    ti = (Ci[8*nq+0*nq+i]*E0[0] + Ci[0*nq+0*nq+i]*E0[1])*dH03 ;

    Cp[0*nq+i] += E0ph[0]*tr + E0ph[1]*ti ;
    Cp[8*nq+i] += E0ph[0]*ti - E0ph[1]*tr ;

    tr = (Ci[0*nq+1*nq+i]*E1[0] - Ci[8*nq+1*nq+i]*E1[1])*sH03 ;
    ti = (Ci[8*nq+1*nq+i]*E1[0] + Ci[0*nq+1*nq+i]*E1[1])*dH03 ;

    Cp[0*nq+i] += E1ph[0]*tr + E1ph[1]*ti ;
    Cp[8*nq+i] += E1ph[0]*ti - E1ph[1]*tr ;

    tr = (Ci[0*nq+2*nq+i]*E1[0] + Ci[8*nq+2*nq+i]*E1[1])*sH03 ;
    ti = (Ci[8*nq+2*nq+i]*E1[0] - Ci[0*nq+2*nq+i]*E1[1])*dH03 ;

    Cp[0*nq+i] += E0ph[0]*tr - E0ph[1]*ti ;
    Cp[8*nq+i] += E0ph[0]*ti + E0ph[1]*tr ;
    
    tr = (Ci[0*nq+3*nq+i]*E0[0] - Ci[8*nq+3*nq+i]*E0[1])*sH03 ;
    ti = (Ci[8*nq+3*nq+i]*E0[0] + Ci[0*nq+3*nq+i]*E0[1])*dH03 ;

    Cp[0*nq+i] += E1ph[0]*tr - E1ph[1]*ti ;
    Cp[8*nq+i] += E1ph[0]*ti + E1ph[1]*tr ;


    tr = (Ci[0*nq+4*nq+i]*E0[0] - Ci[8*nq+4*nq+i]*E0[1])*sH47 ;
    ti = (Ci[8*nq+4*nq+i]*E0[0] + Ci[0*nq+4*nq+i]*E0[1])*dH47 ;

    Cp[0*nq+i] += E0ph[0]*tr + E0ph[1]*ti ;
    Cp[8*nq+i] += E0ph[0]*ti - E0ph[1]*tr ;

    tr = (Ci[0*nq+5*nq+i]*E1[0] - Ci[8*nq+5*nq+i]*E1[1])*sH47 ;
    ti = (Ci[8*nq+5*nq+i]*E1[0] + Ci[0*nq+5*nq+i]*E1[1])*dH47 ;

    Cp[0*nq+i] += E1ph[0]*tr + E1ph[1]*ti ;
    Cp[8*nq+i] += E1ph[0]*ti - E1ph[1]*tr ;

    tr = (Ci[0*nq+6*nq+i]*E1[0] + Ci[8*nq+6*nq+i]*E1[1])*sH47 ;
    ti = (Ci[8*nq+6*nq+i]*E1[0] - Ci[0*nq+6*nq+i]*E1[1])*dH47 ;

    Cp[0*nq+i] += E0ph[0]*tr - E0ph[1]*ti ;
    Cp[8*nq+i] += E0ph[0]*ti + E0ph[1]*tr ;
    
    tr = (Ci[0*nq+7*nq+i]*E0[0] - Ci[8*nq+7*nq+i]*E0[1])*sH47 ;
    ti = (Ci[8*nq+7*nq+i]*E0[0] + Ci[0*nq+7*nq+i]*E0[1])*dH47 ;

    Cp[0*nq+i] += E1ph[0]*tr - E1ph[1]*ti ;
    Cp[8*nq+i] += E1ph[0]*ti + E1ph[1]*tr ;
  }

  return ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_child_parent_shift_laplace)(WBFMM_REAL *Cp,
							   gint Np,
							   WBFMM_REAL *Cc,
							   gint Nc,
							   gint nq,
							   WBFMM_REAL *H03, 
							   WBFMM_REAL *H47,
							   gint Lh,
							   WBFMM_REAL t,
							   WBFMM_REAL *work)

/*
  child-parent shift for Laplace problem, based on Helmholtz version
  in shift.c, but allowing for vector inputs, and with indexing
  changed to Laplace convention

  03: "lower" boxes (think of Morton index)
  47: "upper" boxes

  assumes all coefficients are densely packed in groups of 8*nq
  coefficients of the same index, in order of child box Morton index
*/

{
  WBFMM_REAL *Cr, buf[128], H, Hp, Hm, tn[64] = {0.0}, c ;
  WBFMM_REAL E0[2], E1[2], E0ph[2], E1ph[2], Cnph, Snph ;
  gint nu, n, m, ic, ip, str, i ;

  g_assert(nq <= 8) ;
  
  /*stride in number of elements per entry*/
  str = 8*nq ;

  /*used to store the rotated child coefficients*/
  Cr = &(work[str*(Np+1)*(Np+1)]) ;
  memset(Cr, 0, str*(Nc+1)*(Nc+1)*sizeof(WBFMM_REAL)) ;

  /*rotate child box coefficients using Cr as temporary storage*/
  for ( n = 0 ; n <= Nc ; n ++ ) {
    memset(buf, 0, 16*nq*sizeof(WBFMM_REAL)) ;
    nu = 0 ; ip = n*n ;
    m  = 0 ; ic = n*n ;

    H = H03[wbfmm_rotation_index_numn(nu,m,n)] ;
    for ( i = 0 ; i < 4*nq ; i ++ ) buf[i] = H*Cc[str*ic+i] ;
    
    H = H47[wbfmm_rotation_index_numn(nu,m,n)] ;
    for ( i = 4*nq ; i < 8*nq ; i ++ ) buf[i] = H*Cc[str*ic+i] ;

    for ( m = 1 ; m <= n ; m ++ ) {
      ic = wbfmm_index_laplace_nm(n,m) ;

      E0[0] = cos_n_PI_4(m)   ; E0[1] = sin_n_PI_4(m) ;
      E1[0] = cos_n_PI_4(3*m) ; E1[1] = sin_n_PI_4(3*m) ;

      increment_buf_cp_real(E0, E1,
      			    &(Cc[str*ic]),
      			    H03[wbfmm_rotation_index_numn( nu,m,n)],
      			    H03[wbfmm_rotation_index_numn(-nu,m,n)],
      			    H47[wbfmm_rotation_index_numn( nu,m,n)],
      			    H47[wbfmm_rotation_index_numn(-nu,m,n)],
      			    nq, buf) ;
    }

    for ( i = 0 ; i < 8*nq ; i ++ ) Cr[str*ip+i] += buf[i] ;

    for ( nu = 1 ; nu <= n ; nu ++ ) {
      memset(buf, 0, 16*nq*sizeof(WBFMM_REAL)) ;

      E0ph[0] = cos_n_PI_2(nu) ; E0ph[1] = -sin_n_PI_2(nu) ;
      E1ph[0] = cos_n_PI_2(nu) ; E1ph[1] =  sin_n_PI_2(nu) ;
      
      ip = wbfmm_index_laplace_nm(n,nu) ;
      m = 0 ; ic = n*n ;
      H = H03[wbfmm_rotation_index_numn(nu,m,n)] ;
      for ( i = 0 ; i < 4*nq ; i ++ ) buf[2*i+0] = H*Cc[str*ic+i] ;
      H = H47[wbfmm_rotation_index_numn(nu,m,n)] ;
      for ( i = 4*nq ; i < 8*nq ; i ++ ) buf[2*i+0] = H*Cc[str*ic+i] ;
      
      for ( m = 1 ; m <= n ; m ++ ) {
	E0[0] = cos_n_PI_4(m)   ; E0[1] = sin_n_PI_4(m) ;
	E1[0] = cos_n_PI_4(3*m) ; E1[1] = sin_n_PI_4(3*m) ;
      
	ic = wbfmm_index_laplace_nm(n,m) ;
	increment_buf_cp_complex(E0, E1, &(Cc[str*ic]),
				 H03[wbfmm_rotation_index_numn( nu,m,n)],
				 H03[wbfmm_rotation_index_numn(-nu,m,n)],
				 H47[wbfmm_rotation_index_numn( nu,m,n)],
				 H47[wbfmm_rotation_index_numn(-nu,m,n)],
				 nq, buf) ;
      }
      for ( i = 0 ; i < nq ; i ++ ) {
	Cr[str*(ip+0)+0*nq+i] +=
	  buf[0*2*nq+2*i+0]*E0ph[0] + buf[0*2*nq+2*i+1]*E0ph[1] ;
	Cr[str*(ip+1)+0*nq+i] +=
	  buf[0*2*nq+2*i+1]*E0ph[0] - buf[0*2*nq+2*i+0]*E0ph[1] ;

	Cr[str*(ip+0)+1*nq+i] +=
	  buf[1*2*nq+2*i+0]*E1ph[0] + buf[1*2*nq+2*i+1]*E1ph[1] ;
	Cr[str*(ip+1)+1*nq+i] +=
	  buf[1*2*nq+2*i+1]*E1ph[0] - buf[1*2*nq+2*i+0]*E1ph[1] ;

	Cr[str*(ip+0)+2*nq+i] +=
	  buf[2*2*nq+2*i+0]*E1ph[0] + buf[2*2*nq+2*i+1]*E1ph[1] ;
	Cr[str*(ip+1)+2*nq+i] -=
	  buf[2*2*nq+2*i+1]*E1ph[0] - buf[2*2*nq+2*i+0]*E1ph[1] ;

	Cr[str*(ip+0)+3*nq+i] +=
	  buf[3*2*nq+2*i+0]*E0ph[0] + buf[3*2*nq+2*i+1]*E0ph[1] ;
	Cr[str*(ip+1)+3*nq+i] +=
	  buf[3*2*nq+2*i+1]*E0ph[0] - buf[3*2*nq+2*i+0]*E0ph[1] ;

	Cr[str*(ip+0)+4*nq+i] +=
	  buf[4*2*nq+2*i+0]*E0ph[0] + buf[4*2*nq+2*i+1]*E0ph[1] ;
	Cr[str*(ip+1)+4*nq+i] +=
	  buf[4*2*nq+2*i+1]*E0ph[0] - buf[4*2*nq+2*i+0]*E0ph[1] ;

	Cr[str*(ip+0)+5*nq+i] +=
	  buf[5*2*nq+2*i+0]*E1ph[0] + buf[5*2*nq+2*i+1]*E1ph[1] ;
	Cr[str*(ip+1)+5*nq+i] +=
	  buf[5*2*nq+2*i+1]*E1ph[0] - buf[5*2*nq+2*i+0]*E1ph[1] ;

	Cr[str*(ip+0)+6*nq+i] +=
	  buf[6*2*nq+2*i+0]*E1ph[0] + buf[6*2*nq+2*i+1]*E1ph[1] ;
	Cr[str*(ip+1)+6*nq+i] -=
	  buf[6*2*nq+2*i+1]*E1ph[0] - buf[6*2*nq+2*i+0]*E1ph[1] ;

	Cr[str*(ip+0)+7*nq+i] +=
	  buf[7*2*nq+2*i+0]*E0ph[0] + buf[7*2*nq+2*i+1]*E0ph[1] ;
	Cr[str*(ip+1)+7*nq+i] +=
	  buf[7*2*nq+2*i+1]*E0ph[0] - buf[7*2*nq+2*i+0]*E0ph[1] ;
      }
    }
  }

  /* { */
  /*   gint off ; */
  /*   off = 5 ; */
  /*   fprintf(stderr, */
  /* 	    "%lg %lg %lg %lg %lg %lg %lg %lg %lg\n", */
  /* 	    Cr[0*str + off], */
  /* 	    Cr[1*str + off], */
  /* 	    Cr[(2+0)*str + off], Cr[(2+1)*str + off], */
  /* 	    Cr[4*str + off], */
  /* 	    Cr[(5+0)*str + off], Cr[(5+1)*str + off], */
  /* 	    Cr[(7+0)*str + off], Cr[(7+1)*str + off]) ; */
  /* } */
  
  /*Cr now contains rotated child coefficients, translate and store in
    work*/
  memset(work, 0, str*(Np+1)*(Np+1)*sizeof(WBFMM_REAL)) ;
  tn[0] = 1.0 ;
  m  = 0 ;
  for ( ip = n = 0 ; n <= Np ; (n ++), (ip = n*n) ) {
    for ( ic = nu = 0 ; nu <= MIN(n, Nc) ; (nu ++), (ic = nu*nu) ) {
      c = coaxial_translation_SS_cfft(n, nu, m)*tn[n-nu] ;
      g_assert(tn[n-nu] != 0.0) ;
      /* for ( i = 0 ; i < str ; i ++ ) */
      for ( i = 0 ; i < 8*nq ; i ++ )
  	work[str*ip+i] += c*Cr[str*ic+i] ;
    }
    tn[n+1] = -tn[n]*t ;
  }
  
  for ( m = 1 ; m <= Np ; m ++ ) {
    for ( n = m ; n <= Np ; n ++ ) {
      ip = wbfmm_index_laplace_nm(n,m) ;
      for ( nu = m ; nu <= n ; nu ++ ) {
  	ic = wbfmm_index_laplace_nm(nu,m) ;
  	c = coaxial_translation_SS_cfft(n, nu, m)*tn[n-nu] ;
	/* g_assert(tn[n-nu] != 0.0) ; */
	for ( i = 0 ; i < 8*nq ; i ++ ) {
	  work[str*(ip+0)+i] += c*Cr[str*(ic+0)+i] ;
	  work[str*(ip+1)+i] += c*Cr[str*(ic+1)+i] ;
	}
      }
    }
  }

  /* { */
  /*   gint off ; */
  /*   off = 7 ; */
  /*   fprintf(stderr, */
  /* 	    "%lg %lg %lg %lg %lg %lg %lg %lg %lg\n", */
  /* 	    work[0*str + off], */
  /* 	    work[1*str + off], */
  /* 	    work[(2+0)*str + off], work[(2+1)*str + off], */
  /* 	    work[4*str + off], */
  /* 	    work[(5+0)*str + off], work[(5+1)*str + off], */
  /* 	    work[(7+0)*str + off], work[(7+1)*str + off]) ; */
  /* } */
  
  /*work now contains rotated and shifted coefficients, perform
    reverse rotation into parent coefficient array*/
  for ( n = 0 ; n <= Np ; n ++ ) {
    nu = 0 ; ip = n*n ;
    m  = 0 ; ic = n*n ;

    H = H03[wbfmm_rotation_index_numn(nu,m,n)] ;
    for ( i = 0 ; i < nq ; i ++ ) {
      Cp[str*ip+i] +=
	H*(work[str*ic+0*nq+i] + work[str*ic+1*nq+i] +
	   work[str*ic+2*nq+i] + work[str*ic+3*nq+i]) ;
    }
    H = H47[wbfmm_rotation_index_numn(nu,m,n)] ;
    for ( i = 0 ; i < nq ; i ++ ) {
      Cp[str*ip+i] +=
	H*(work[str*ic+4*nq+i] + work[str*ic+5*nq+i] +
	   work[str*ic+6*nq+i] + work[str*ic+7*nq+i]) ;
    }

    E0ph[0] = cos_n_PI_4(nu) ; E0ph[1] = sin_n_PI_4(nu) ;
    for ( m = 1 ; m <= n ; m ++ ) {
      E0[0] = cos_n_PI_2(m) ; E0[1] = -sin_n_PI_2(m) ;
      E1[0] = cos_n_PI_2(m) ; E1[1] =  sin_n_PI_2(m) ;

      ic = wbfmm_index_laplace_nm(n,m) ;
      increment_cfft_cp_real(E0, E1, &(work[str*ic]),
			     H03[wbfmm_rotation_index_numn( nu,m,n)],
			     H03[wbfmm_rotation_index_numn(-nu,m,n)],
			     H47[wbfmm_rotation_index_numn( nu,m,n)],
			     H47[wbfmm_rotation_index_numn(-nu,m,n)],
			     nq, &(Cp[str*ip])) ;
    }

    for ( nu = 1 ; nu <= n ; nu ++ ) {
      E0ph[0] = cos_n_PI_4(  nu) ; E0ph[1] = sin_n_PI_4(  nu) ;
      E1ph[0] = cos_n_PI_4(3*nu) ; E1ph[1] = sin_n_PI_4(3*nu) ;

      ip = wbfmm_index_laplace_nm(n,nu) ;
      m = 0 ; ic = n*n ;
      H = H03[wbfmm_rotation_index_numn(nu,m,n)] ;
      for ( i = 0 ; i < nq ; i ++ ) {
      	Cp[str*(ip+0)+i] += work[str*ic+0*nq+i]*H*E0ph[0] ;
      	Cp[str*(ip+1)+i] -= work[str*ic+0*nq+i]*H*E0ph[1] ;
      	Cp[str*(ip+0)+i] += work[str*ic+1*nq+i]*H*E1ph[0] ;
      	Cp[str*(ip+1)+i] -= work[str*ic+1*nq+i]*H*E1ph[1] ;

      	Cp[str*(ip+0)+i] += work[str*ic+2*nq+i]*H*E0ph[0] ;
      	Cp[str*(ip+1)+i] += work[str*ic+2*nq+i]*H*E0ph[1] ;

      	Cp[str*(ip+0)+i] += work[str*ic+3*nq+i]*H*E1ph[0] ;
      	Cp[str*(ip+1)+i] += work[str*ic+3*nq+i]*H*E1ph[1] ;
      }
      H = H47[wbfmm_rotation_index_numn(nu,m,n)] ;
      for ( i = 0 ; i < nq ; i ++ ) {
      	Cp[str*(ip+0)+i] += work[str*ic+4*nq+i]*H*E0ph[0] ;
      	Cp[str*(ip+1)+i] -= work[str*ic+4*nq+i]*H*E0ph[1] ;
      	Cp[str*(ip+0)+i] += work[str*ic+5*nq+i]*H*E1ph[0] ;
      	Cp[str*(ip+1)+i] -= work[str*ic+5*nq+i]*H*E1ph[1] ;

      	Cp[str*(ip+0)+i] += work[str*ic+6*nq+i]*H*E0ph[0] ;
      	Cp[str*(ip+1)+i] += work[str*ic+6*nq+i]*H*E0ph[1] ;

      	Cp[str*(ip+0)+i] += work[str*ic+7*nq+i]*H*E1ph[0] ;
      	Cp[str*(ip+1)+i] += work[str*ic+7*nq+i]*H*E1ph[1] ;
      }

      for ( m = 1 ; m <= n ; m ++ ) {
	E0[0] = cos_n_PI_2(m) ; E0[1] = -sin_n_PI_2(m) ;
	E1[0] = cos_n_PI_2(m) ; E1[1] =  sin_n_PI_2(m) ;

	ic = wbfmm_index_laplace_nm(n,m) ;
	increment_cfft_cp_complex(E0, E1, E0ph, E1ph, &(work[str*ic]),
				  H03[wbfmm_rotation_index_numn( nu,m,n)],
				  H03[wbfmm_rotation_index_numn(-nu,m,n)],
				  H47[wbfmm_rotation_index_numn( nu,m,n)],
				  H47[wbfmm_rotation_index_numn(-nu,m,n)],
				  nq, &(Cp[str*ip])) ;
      }
    }
  }

  /* { */
  /*   gint off ; */
  /*   off = 0 ; */
  /*   fprintf(stderr, */
  /* 	    "%lg %lg %lg %lg %lg %lg %lg %lg %lg\n", */
  /* 	    Cp[0*str + off], */
  /* 	    Cp[1*str + off], */
  /* 	    Cp[(2+0)*str + off], Cp[(2+1)*str + off], */
  /* 	    Cp[4*str + off], */
  /* 	    Cp[(5+0)*str + off], Cp[(5+1)*str + off], */
  /* 	    Cp[(7+0)*str + off], Cp[(7+1)*str + off]) ; */
  /* } */
  
  return 0 ;
}
