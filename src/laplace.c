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

WBFMM_REAL *_wbfmm_SS_coefficients_laplace = NULL,
  *_wbfmm_RR_coefficients_laplace = NULL,
  *_wbfmm_SR_coefficients_laplace = NULL ;
gint _wbfmm_translation_Nmax = 0 ;

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

  n = 1 ; 
  m = 0 ; 
  rn *= r ;
  idx = n*n ;
  for ( i = 0 ; i < nq ; i ++ ) field[i] += cfft[cstr*idx+i]*rn*Pn[m] ;

  m = 1 ; 
  idx = wbfmm_index_laplace_nm(n,m) ;
  for ( i = 0 ; i < nq ; i ++ ) {
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
	field[i] += 2.0*Pn[m]*rn*(cfft[cstr*(idx+0)+i]*Cmph[m] -
				  cfft[cstr*(idx+1)+i]*Smph[m]) ;
      }
    }
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
    
    for ( j = 0 ; j < nq ; j ++ ) field[j] /= 4.0*M_PI ;

    return 0 ;
  }

  g_assert_not_reached() ;
  
  return 0 ;
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
      c = wbfmm_coaxial_translation_SS_cfft(n, nd, m)*tn[n-nd] ;
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
	c = wbfmm_coaxial_translation_SS_cfft(n, nd, m)*tn[n-nd] ;
	for ( i = 0 ; i < nq ; i ++ ) {
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
      /* c = coaxial_translation_RR_cfft(n, nd, m)*tn[nd-n] ; */
      c = wbfmm_coaxial_translation_RR_cfft(n, nd, m)*tn[nd-n] ;
      for ( i = 0 ; i < nq ; i ++ ) Co[cstro*idxo+i] += c*Ci[cstri*idxi+i] ;
    }
    tn[nd+1] = tn[nd]*t ;
  }
  
  for ( m = 1 ; m <= No ; m ++ ) {
    for ( nd = m ; nd <= No ; nd ++ ) {
      idxi = wbfmm_index_laplace_nm(nd,m) ;
      for ( n = m ; n <= nd ; n ++ ) {
  	idxo = wbfmm_index_laplace_nm(n,m) ;
	/* c = coaxial_translation_RR_cfft(n, nd, m)*tn[nd-n] ; */
	c = wbfmm_coaxial_translation_RR_cfft(n, nd, m)*tn[nd-n] ;
	for ( i = 0 ; i < nq ; i ++ ) {
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
      for ( i = 0 ; i < nq ; i ++ ) Co[cstro*idxo+i] += c*Ci[cstri*idxi+i] ;
    }
  }

  for ( m = 1 ; m <= No ; m ++ ) {
    for ( n = m ; n <= No ; n ++ ) {
      idxo = wbfmm_index_laplace_nm(n,m) ;
      for ( nd = m ; nd <= Ni ; nd ++ ) {
	idxi = wbfmm_index_laplace_nm(nd,m) ;
  	c = coaxial_translation_SR_cfft(n, nd, m)/tn[n+nd+1] ;
  	for ( i = 0 ; i < nq ; i ++ ) {
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

gint WBFMM_FUNCTION_NAME(wbfmm_coaxial_translate_laplace_init)(gint N)

{
  gint n, m, nu, i, Nmax, ne ;
  
  /* if ( _wbfmm_SS_coefficients_laplace != NULL ) return 0 ; */

  Nmax = _wbfmm_translation_Nmax = N + 1 ;
  
  ne = _wbfmm_SS_coefficient_index_nmnu(Nmax+1,0,0) ;
  _wbfmm_SS_coefficients_laplace =
    (WBFMM_REAL *)g_malloc0(ne*sizeof(WBFMM_REAL)) ;

  for ( n = 0 ; n <= Nmax ; n ++ ) {
    for ( m = 0 ; m <= n ; m ++ ) {
      for ( nu = m ; nu <= n ; nu ++ ) {
	i = _wbfmm_SS_coefficient_index_nmnu(n, m, nu) ;
	g_assert(i < ne) ;
  	g_assert(_wbfmm_SS_coefficients_laplace[i] == 0.0) ;
	_wbfmm_SS_coefficients_laplace[i] =
	  coaxial_translation_SS_cfft(n, nu, m) ;
      }
    }
  }

  ne = (Nmax+2)*(Nmax+1)*(Nmax+1)/2 ;
  
  _wbfmm_RR_coefficients_laplace =
    (WBFMM_REAL *)g_malloc0(ne*sizeof(WBFMM_REAL)) ;

  for ( m = 0 ; m <= Nmax ; m ++ ) {
    for ( nu = m ; nu <= Nmax ; nu ++ ) {
      for ( n = m ; n <= nu ; n ++ ) {
  	i = _wbfmm_RR_coefficient_index_nmnu(n, m, nu) ;
	g_assert(i < ne) ;
  	g_assert(_wbfmm_RR_coefficients_laplace[i] == 0.0) ;
  	_wbfmm_RR_coefficients_laplace[i] =
  	  coaxial_translation_RR_cfft(n, nu, m) ;
      }
    }
  }
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_tree_laplace_coefficient_init)(wbfmm_tree_t *t,
							      guint l, 
							      guint nr,
							      guint ns)

{
  gint nb, nc, nq, i, j ;
  wbfmm_box_t *boxes ;
  WBFMM_REAL *c ;

  if ( (nq = wbfmm_tree_source_size(t)) < 1 )
    g_error("%s: tree has invalid number of components in source terms (%d)",
	    __FUNCTION__, nq) ;
  
  g_assert(t->problem == WBFMM_PROBLEM_LAPLACE ) ;
  
  g_assert(l <= wbfmm_tree_depth(t)) ;

  /*number of boxes at level l*/
  nb = 1 << (3*l) ;

  t->mps[l] = t->mpr[l] = NULL ;
  t->order_s[l] = ns ; t->order_r[l] = nr ; 

  /*number of coefficients in singular expansions*/
  if ( ns != 0 ) {
    nc = (ns+1)*(ns+1) ;
    t->mps[l] = g_malloc0(nb*nc*nq*sizeof(WBFMM_REAL)) ;
    c = (WBFMM_REAL *)(t->mps[l]) ;
    /*
      set box pointers to start of their coefficients, noting that
      coefficients are packed in groups of eight for shift operations
    */
    boxes = t->boxes[l] ;
    for ( i = 0 ; i < nb ; i += 8 ) {
      for ( j = 0 ; j < 8 ; j ++ ) {
	boxes[i+j].mps = &(c[i*nq*nc+nq*j]) ;
      }
    }
  }

  if ( nr != 0 ) {
    nc = (nr+1)*(nr+1) ;
    t->mpr[l] = g_malloc0(nb*nc*nq*sizeof(WBFMM_REAL)) ;
    /* nc = wbfmm_coefficient_index_nm(nr+1,0) ; */
    /* t->mpr[l] = g_malloc0(nb*2*nc*sizeof(WBFMM_REAL)) ; */
    c = (WBFMM_REAL *)(t->mpr[l]) ;
    /*
      set box pointers to start of their coefficients, noting that
      coefficients are packed in groups of eight for shift operations
    */
    boxes = t->boxes[l] ;
    for ( i = 0 ; i < nb ; i += 8 ) {
      for ( j = 0 ; j < 8 ; j ++ ) {
	boxes[i+j].mpr = &(c[i*nq*nc+nq*j]) ;
	/* boxes[i+j].mpr = &(c[i*2*nc+2*j]) ; */
      }
    }
  }
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_tree_laplace_leaf_expansions)(wbfmm_tree_t *t,
							     WBFMM_REAL *src,
							     gint sstr,
							     WBFMM_REAL *normals,
							     gint nstr,
							     WBFMM_REAL *dipoles,
							     gint dstr,
							     gboolean zero_expansions,
							     WBFMM_REAL *work)

{
  guint32 nb, nc, i, j, ns, d, idx, nq ;
  guint64 im ;
  wbfmm_box_t *boxes ;
  WBFMM_REAL *xs, *q, *fd, *n, xb[3], wb ;

  if ( (nq = wbfmm_tree_source_size(t)) < 1 )
    g_error("%s: tree has invalid number of components in source terms (%d)",
	    __FUNCTION__, nq) ;

  g_assert(t->problem == WBFMM_PROBLEM_LAPLACE ) ;

  /*depth of leaves*/
  d = wbfmm_tree_depth(t) ;
  /*order of singular expansions*/
  ns = t->order_s[d] ;
  /*number of boxes*/
  nb = 1 << (3*d) ;
  /*number of coefficients*/
  nc = (ns+1)*(ns+1) ;
  /* nc = wbfmm_coefficient_index_nm(ns+1,0) ; */

  /*zero the coefficients before accumulating*/
  if ( zero_expansions )
    memset(t->mps[d], 0, nb*nc*nq*sizeof(WBFMM_REAL)) ;

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
      WBFMM_FUNCTION_NAME(wbfmm_box_location_from_index)(im, d, 
							 wbfmm_tree_origin(t), 
							 wbfmm_tree_width(t),
							 xb, &wb) ;
      xb[0] += 0.5*wb ; xb[1] += 0.5*wb ; xb[2] += 0.5*wb ; 

      for ( j = 0 ; j < boxes[i].n ; j ++ ) {
	idx = t->ip[boxes[i].i+j] ;
	xs = wbfmm_tree_point_index(t,idx) ;
	q = &(src[idx*sstr]) ;
	/* WBFMM_FUNCTION_NAME(wbfmm_expansion_h_cfft)(k, ns, xb, xs, q, */
	/* 				      boxes[i].mps, 8, work) ; */
	WBFMM_FUNCTION_NAME(wbfmm_expansion_laplace_cfft)(ns, xb, xs, q, nq,
							  boxes[i].mps, 8*nq,
							  work) ;
      }
    }

    return 0 ;
  }

  g_assert_not_reached() ;
#if 0  
  if ( src == NULL && normals != NULL ) {
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
	fd = &(dipoles[idx*dstr]) ;
	n  = &(normals[idx*nstr]) ;
	WBFMM_FUNCTION_NAME(wbfmm_expansion_normal_h_cfft)(k, ns, xb, xs, n, fd,
						     boxes[i].mps, 8, work) ;
      }
    }

    return 0 ;
  }

  if ( src != NULL && normals != NULL ) {
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
	n  = &(normals[idx*nstr]) ;
	WBFMM_FUNCTION_NAME(wbfmm_expansion_normal_h_cfft)(k, ns, xb, xs, n, fd,
						     boxes[i].mps, 8, work) ;
	WBFMM_FUNCTION_NAME(wbfmm_expansion_h_cfft)(k, ns, xb, xs, q,
					      boxes[i].mps, 8, work) ;
      }
    }

    return 0 ;
  }

  g_assert_not_reached() ;
#endif
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_box_fields_laplace)(wbfmm_tree_t *t,
						   gint level,
						   WBFMM_REAL *xf,
						   WBFMM_REAL *field,
						   WBFMM_REAL *work)

{
  gint nbox, i, nq ;
  WBFMM_REAL xb[3], wb ;
  wbfmm_box_t *boxes ;

  nq = wbfmm_tree_source_size(t) ;
  nbox = 1 << 3*level ;
  
  boxes = t->boxes[level  ] ;
  memset(field, 0, nq*sizeof(WBFMM_REAL)) ;
  
  for ( i = 0 ; i < nbox ; i ++ ) {
    WBFMM_FUNCTION_NAME(wbfmm_tree_box_centre)(t, level, i, xb, &wb) ;
    WBFMM_FUNCTION_NAME(wbfmm_expansion_laplace_evaluate)(xb, boxes[i].mps,
							  8*nq,
							  t->order_s[level],
							  nq, xf, field, work) ;
  }
  
 return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_tree_box_laplace_local_field)(wbfmm_tree_t *t,
							     guint level,
							     guint b,
							     WBFMM_REAL *x,
							     WBFMM_REAL *f,
							     WBFMM_REAL *src,
							     gint sstr,
							     WBFMM_REAL
							     *normals,
							     gint nstr,
							     WBFMM_REAL *d,
							     gint dstr,
							     gboolean
							     eval_neighbours,
							     WBFMM_REAL *work)

{
  WBFMM_REAL xb[3], wb, *C, *xs, r, h0[2], h1[2], fR[2] ;
  wbfmm_box_t *boxes, box ;
  guint64 neighbours[27] ;
  gint nnbr, i, j, idx, nq ;

  g_assert(t->problem == WBFMM_PROBLEM_LAPLACE ) ;

  nq = wbfmm_tree_source_size(t) ;

  boxes = t->boxes[level] ;
  C = boxes[b].mpr ;

  WBFMM_FUNCTION_NAME(wbfmm_tree_box_centre)(t, level, b, xb, &wb) ;

  WBFMM_FUNCTION_NAME(wbfmm_expansion_laplace_local_evaluate)(xb, C, 8*nq,
							t->order_r[level],
							nq, x, f, work) ;

  if ( !eval_neighbours ) return 0 ;

  if ( src == NULL && normals == NULL && d == NULL ) return 0 ;
  
  if ( normals != NULL && d == NULL ) {
    g_error("%s: normals specified but no dipole strengths (d == NULL)",
	    __FUNCTION__) ;
  }

  /* return 0 ; */
  
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
	  /* WBFMM_FUNCTION_NAME(wbfmm_bessel_h_init)(k*r, h0, h1) ; */
	  h0[0] /= 4.0*M_PI ; h0[1] /= 4.0*M_PI ; 
	  f[0] += h0[0]*src[idx*sstr+0] - h0[1]*src[idx*sstr+1] ;
	  f[1] += h0[1]*src[idx*sstr+0] + h0[0]*src[idx*sstr+1] ;
	}
      }
    }

    return 0 ;
  }
  
  g_assert_not_reached() ; 
  
  return 0 ;
}
