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

/**
 * @file   translations.c
 * @author Michael Carley <ensmjc@rpc-ensmjc.bath.ac.uk>
 * @date   Mon Jun 24 11:49:41 2019
 * 
 * @brief  Translation operations
 * 
 * 
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

/* #define CHECK_COEFFICIENTS */

/* #define WBFMM_CHECK_ISNAN */

#ifdef WBFMM_CHECK_ISNAN
#include <stdlib.h>

static gint check_isnan(gchar *name, WBFMM_REAL *f, gint n)

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

#ifdef CHECK_COEFFICIENTS
#include <stdio.h>
#endif

/* #define recursion_anm(_n,_m)						\ */
/*   (((_n) < (ABS(_m))) ? 0 :						\ */
/*    (sqrt((WBFMM_REAL)((_n)+1+ABS(_m))*((_n)+1-ABS(_m))/(2*(_n)+1)/(2*(_n)+3)))) */

/**
 * @defgroup translations Translation coefficients and operations
 *
 * @brief Generation and application of translation operators for 
 * multipole expansions
 *
 * Recursive computation of real and complex translation coefficients
 * using the methods of Gumerov and Duraiswami,
 * http://dx.doi.org/10.1137/S1064827501399705
 *
 */

/* @{ */

#ifdef CHECK_COEFFICIENTS

gint coaxial_index(gint l, gint m, gint n, gint *idx, gint *sgn)

/*
  generate the index into the coefficient matrix for general l,m,n,

  l >= 0, n >= 0, -n <= m <= n

  sgn is the sign multiplier for the coefficient, i.e. use
  sgn*cfft[idx] (see G&D, 4.82)
*/

{
  m = ABS(m) ;

  if ( m > n ) { g_assert_not_reached() ; *idx = *sgn = 0 ; return 0 ; }

  if ( l >= n ) {
    *idx = wbfmm_coaxial_index_lmn(l, m, n) ;
    *sgn = 1 ;
    g_assert(*idx >= 0) ;
    return 0 ;
  }

  *idx = wbfmm_coaxial_index_lmn(n, m, l) ;
  *sgn = minus_one_pow(n+l) ;
  g_assert(*idx >= 0) ;

  return 0 ;
}

static gint coaxial_coefficient_recursion_check(WBFMM_REAL *cfft, gint L)

{
  gint l, m, n, i1, i2, i3, i4, sgn ;
  WBFMM_REAL a1, a2, a3, a4, err, emax ;

  /*fixed m check 4.79*/
  emax = 0.0 ;
  for ( l = 0 ; l < L ; l ++ ) {
    for ( n = 0 ; n <= l ; n ++ ) {
      for ( m = -n ; m <= n ; m ++ ) {
	a1 = recursion_anm(n-1, m) ;
  	a2 = recursion_anm(n  , m) ;
  	a3 = recursion_anm(l  , m) ;
  	a4 = recursion_anm(l-1, m) ;
	coaxial_index(l, m, n-1, &i1, &sgn) ; a1 *= sgn ;
	coaxial_index(l, m, n+1, &i2, &sgn) ; a2 *= sgn ;
	coaxial_index(l+1, m, n, &i3, &sgn) ; a3 *= sgn ;
	coaxial_index(l-1, m, n, &i4, &sgn) ; a4 *= sgn ;
	err =	(a1*cfft[i1] - a2*cfft[i2]) -
	  (a3*cfft[i3] - a4*cfft[i4]) ;
	fprintf(stderr, "%2d %2d %2d %+e (%lg,%lg,%lg,%lg)\n",
		l, m, n, err,
		cfft[i1], cfft[i2],
		cfft[i3], cfft[i4]) ;
	emax = MAX(fabs(err), emax) ;
      }
    }
  }

  fprintf(stderr, "Maximum error in recursion relation (4.79): %e\n", emax) ;

  /*check on 4.80*/
  emax = 0.0 ;
  for ( l = 0 ; l < L ; l ++ ) {
    for ( n = 0 ; n < l ; n ++ ) {
      for ( m = -n+1 ; m <= n-1 ; m ++ ) {
	a1 = recursion_bnm(n  ,  m  ) ;
  	a2 = recursion_bnm(n+1, -m-1) ;
  	a3 = recursion_bnm(l+1,  m  ) ;
  	a4 = recursion_bnm(l  , -m-1) ;
	coaxial_index(l,   m+1, n-1, &i1, &sgn) ; a1 *= sgn ;
	coaxial_index(l,   m+1, n+1, &i2, &sgn) ; a2 *= sgn ;
	coaxial_index(l+1, m  , n  , &i3, &sgn) ; a3 *= sgn ;
	coaxial_index(l-1, m,   n  , &i4, &sgn) ; a4 *= sgn ;
	err = (a1*cfft[i1] - a2*cfft[i2]) -
	  (a3*cfft[i3] - a4*cfft[i4]) ;
	fprintf(stderr, "%2d %2d %2d %+8lg (%lg,%lg,%lg,%lg)\n",
		l, m, n, err,
		cfft[i1], cfft[i2],
		cfft[i3], cfft[i4]) ;
	emax = MAX(fabs(err), emax) ;
      }
    }
  }

  fprintf(stderr, "Maximum error in recursion relation (4.80): %e\n", emax) ;

  return 0 ;
}

#endif /*CHECK_COEFFICIENTS*/

WBFMM_REAL WBFMM_FUNCTION_NAME(recursion_anm)(gint n, gint m)

{
  m = ABS(m) ;

  if ( n < m ) return 0.0 ;

  return SQRT((WBFMM_REAL)(n+1+m)*(n+1-m)/(2*n+1)/(2*n+3)) ;
}

WBFMM_REAL WBFMM_FUNCTION_NAME(recursion_bnm)(gint n, gint m)

{
  /* g_assert(n >= 0) ; */
  if ( n < ABS(m) ) return 0.0 ;

  if ( m >= 0 ) return SQRT((WBFMM_REAL)(n-m-1)*(n-m)/(2*n-1)/(2*n+1)) ;

  g_assert((-n <= m) && (m < 0)) ;

  return -SQRT((WBFMM_REAL)(n-m-1)*(n-m)/(2*n-1)/(2*n+1)) ;
}

static void coefficient_recursions_coaxial(WBFMM_REAL *cfft, gint L, gint nc)

/*
  recursions from G&D for evaluation of coaxial translation
  coefficients; cfft should be initialized as required before function
  call; recursions are real and applied to nc components in cfft
*/

{
  gint l, m, n, i1, i2, i3, i4, i ;
  WBFMM_REAL a1, a2, a3, a4 ;

  /*G&D 4.79*/
  m = 0 ;
  /*fill (E|F)_{ln}^0*/
  for ( n = 0 ; n <= L ; n ++ ) {
    for ( l = n+1 ; l <= 2*L-n ; l ++ ) {
      i1 = wbfmm_coaxial_index_lmn(l,   m, n-1) ;
      i2 = wbfmm_coaxial_index_lmn(l,   m, n+1) ;
      i3 = wbfmm_coaxial_index_lmn(l+1, m, n  ) ;
      i4 = wbfmm_coaxial_index_lmn(l-1, m, n  ) ;
      a1 = WBFMM_FUNCTION_NAME(recursion_anm)(n-1, m) ;
      a2 = WBFMM_FUNCTION_NAME(recursion_anm)(n  , m) ;
      a3 = WBFMM_FUNCTION_NAME(recursion_anm)(l  , m) ;
      a4 = WBFMM_FUNCTION_NAME(recursion_anm)(l-1, m) ;
      for ( i = 0 ; i < nc ; i ++ ) 
	cfft[nc*i2+i] = 
	  (a1*cfft[nc*i1+i] - a3*cfft[nc*i3+i] + a4*cfft[nc*i4+i])/a2 ;
    }
  }

  /*G&D 4.84*/
  /*fill (E|F)_{l,m+1}^{m+1}*/
  for ( m = 0 ; m <= L ; m ++ ) {
    for ( l = m+1 ; l <= 2*L-m ; l ++ ) {
      a1 = WBFMM_FUNCTION_NAME(recursion_bnm)(m+1, -m-1) ;
      a2 = WBFMM_FUNCTION_NAME(recursion_bnm)(l  , -m-1) ;
      a3 = WBFMM_FUNCTION_NAME(recursion_bnm)(l+1,  m  ) ;
      i1 = wbfmm_coaxial_index_lmn(l  , m+1, m+1) ;
      i2 = wbfmm_coaxial_index_lmn(l-1, m  , m  ) ;
      i3 = wbfmm_coaxial_index_lmn(l+1, m  , m  ) ;
      for ( i = 0 ; i < nc ; i ++ ) 
	cfft[nc*i1+i] = (a2*cfft[nc*i2+i] - a3*cfft[nc*i3+i])/a1 ;
    }
  }

  for ( m = 1 ; m <= L ; m ++ ) {
    for ( n = m ; n <= L ; n ++ ) {
      for ( l = n+1 ; l <= 2*L-n ; l ++ ) {
	i1 = wbfmm_coaxial_index_lmn(l  , m, n-1) ;
	i2 = wbfmm_coaxial_index_lmn(l  , m, n+1) ;
	i3 = wbfmm_coaxial_index_lmn(l+1, m, n  ) ;
	i4 = wbfmm_coaxial_index_lmn(l-1, m, n  ) ;
	a1 = WBFMM_FUNCTION_NAME(recursion_anm)(n-1, m) ;
	a2 = WBFMM_FUNCTION_NAME(recursion_anm)(n  , m) ;
	a3 = WBFMM_FUNCTION_NAME(recursion_anm)(l  , m) ;
	a4 = WBFMM_FUNCTION_NAME(recursion_anm)(l-1, m) ;
	for ( i = 0 ; i < nc ; i ++ ) 
	  cfft[nc*i2+i] = 
	    (a1*cfft[nc*i1+i] - a3*cfft[nc*i3+i] + a4*cfft[nc*i4+i])/a2 ;
      }
    }
  }
  
  return ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_coefficients_RR_coaxial)(WBFMM_REAL *cfftRR,
							gint L,
							WBFMM_REAL kr, 
							WBFMM_REAL *work)

/*
  generate coefficients for coaxial shift kr (positive or negative) in
  z direction, G&D section 4.8

  indexing on coaxial coefficients using wbfmm_coaxial_index_lmn(l,m,n) 
  l >= n >= m >= 0 (no check on indexing)

  these are the (R|R) coefficients, which are identical to the (S|S)
  coefficients, G&D (4.62)
*/

{
  gint l, m, n, idx, sgnl, sw = 1 ;
  WBFMM_REAL jlm1, jl, *cfft ;

#ifdef WBFMM_SINGLE_PRECISION
  if ( L > 32 ) 
    g_error("%s: cannot generate single-precision translations for L > 32 (L==%d)",
	    __FUNCTION__, L) ;
#endif /*WBFMM_SINGLE_PRECISION*/
  
  if ( kr < 0 ) { kr = -kr ; sw = -1 ; }

  g_assert(kr > 0.0) ;

  WBFMM_FUNCTION_NAME(wbfmm_bessel_j_init)(kr, &jlm1, &jl) ; 
  cfft = work ;

  /*G&D section 4.8.1, G&D 4.83*/
  /*initialize (R|R)_{l0}^0*/
  /*(-1)^l term reverses sense of translation*/
  l = 0 ; m = 0 ; n = 0 ; sgnl = 1 ;
  idx = wbfmm_coaxial_index_lmn(l, m, n) ;
  cfft[idx] = sgnl*SQRT(2*l+1)*jlm1 ;
  l = 1 ; sgnl *= sw ;
  idx = wbfmm_coaxial_index_lmn(l, m, n) ;
  cfft[idx] = sgnl*SQRT(2*l+1)*jl ;

  /*(-1)^l term reverses sense of translation*/
  for ( l = 2 ; l <= 2*L ; l ++ ) {
    WBFMM_FUNCTION_NAME(wbfmm_bessel_j_recursion)(&jlm1, &jl, kr, l-1) ;
    sgnl *= sw ;
    idx = wbfmm_coaxial_index_lmn(l, m, n) ;
    cfft[idx] = sgnl*SQRT(2*l+1)*jl ;
  }

  /*recursion is the same for either direction of translation*/
  coefficient_recursions_coaxial(cfft, L, 1) ;

#ifdef WBFMM_CHECK_ISNAN
  check_isnan("RR coefficients", cfft, (L+1)*(L+2)*(L+3)/6) ;
#endif /*WBFMM_CHECK_ISNAN*/

  /*copy to output array (probably not ideal, but better for development)*/
  memcpy(cfftRR, cfft, (L+1)*(L+2)*(L+3)/6*sizeof(WBFMM_REAL)) ;

#ifdef CHECK_COEFFICIENTS
  coaxial_coefficient_recursion_check(cfftRR, L) ;
#endif

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_coefficients_SR_coaxial)(WBFMM_REAL *cfftSR,
							gint L,
							WBFMM_REAL kr, 
							WBFMM_REAL *work)

/*
  generate coefficients for coaxial shift kr (positive or negative) in
  z direction, G&D section 4.8

  indexing on coaxial coefficients using wbfmm_coaxial_index_lmn(l,m,n) 
  l >= n >= m >= 0 (no check on indexing)

  these are the (S|R) coefficients, for singular to regular
  translations; note these are complex, with two REALs per coefficient
*/

{
  gint l, m, n, idx, sgnl, sw = 1 ;
  WBFMM_REAL hlm1[2], hl[2], *cfft ;
  
#ifdef WBFMM_SINGLE_PRECISION
  if ( L > 32 ) 
    g_error("%s: cannot generate single-precision translations for L > 32 (L==%d)",
	    __FUNCTION__, L) ;
#endif /*WBFMM_SINGLE_PRECISION*/

  if ( kr < 0 ) { kr = -kr ; sw = -1 ; }

  g_assert(kr > 0.0) ;
  
  WBFMM_FUNCTION_NAME(wbfmm_bessel_h_init)(kr, hlm1, hl) ; 
  cfft = work ;

  /*G&D section 4.8.1, G&D 4.83*/
  /*initialize (S|R)_{l0}^0*/
  l = 0 ; m = 0 ; n = 0 ; sgnl = 1 ;
  idx = wbfmm_coaxial_index_lmn(l, m, n) ;
  cfft[2*idx+0] = sgnl*SQRT(2*l+1)*hlm1[0] ;
  cfft[2*idx+1] = sgnl*SQRT(2*l+1)*hlm1[1] ;
  /*(-1)^l term reverses sense of translation*/
  l = 1 ; sgnl *= sw ;
  idx = wbfmm_coaxial_index_lmn(l, m, n) ;
  cfft[2*idx+0] = sgnl*SQRT(2*l+1)*hl[0] ;
  cfft[2*idx+1] = sgnl*SQRT(2*l+1)*hl[1] ;

  /*(-1)^l term reverses sense of translation*/
  for ( l = 2 ; l <= 2*L ; l ++ ) {
    WBFMM_FUNCTION_NAME(wbfmm_bessel_h_recursion)(hlm1, hl, kr, l-1) ;
    sgnl *= sw ;
    idx = wbfmm_coaxial_index_lmn(l, m, n) ;
    cfft[2*idx+0] = sgnl*SQRT(2*l+1)*hl[0] ;
    cfft[2*idx+1] = sgnl*SQRT(2*l+1)*hl[1] ;
  }
  
  /*recursion is the same for either direction of translation*/
  coefficient_recursions_coaxial(cfft, L, 2) ;

#ifdef WBFMM_CHECK_ISNAN
  check_isnan("SR coefficients", cfft, 2*(L+1)*(L+2)*(L+3)/6) ;
#endif /*WBFMM_CHECK_ISNAN*/
  
  /*copy to output array (probably not ideal, but better for development)*/
  memcpy(cfftSR, cfft, 2*(L+1)*(L+2)*(L+3)/6*sizeof(WBFMM_REAL)) ;

#ifdef CHECK_COEFFICIENTS
  coaxial_coefficient_recursion_check(cfftSR, L) ;
#endif

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_coaxial_translate_ref)(WBFMM_REAL *Co,
						      gint cstro, 
						      gint No,
						      WBFMM_REAL *Ci,
						      gint cstri,
						      gint Ni,
						      gint nq,
						      WBFMM_REAL *cfft,
						      gint L,
						      gboolean complex,
						      WBFMM_REAL sc)

/*
  apply coaxial translation cfft to input coefficients Ci, up to order
  Ni, accumulating output in Co, up to order No, G&D (4.76)

  array bounds and validity of inputs are not checked

  Co is not zeroed (to allow for accumulation in place)

  complex: TRUE for complex shift coefficients (e.g. S|R), FALSE for
  real (e.g. R|R)
*/

{
  gint l, m, n, idxp, idxm, idxc, sgn, sgnn, sgnm, offp, offm, j ;

  g_assert(L >= No) ;
  g_assert(L >= Ni) ;

  if ( cstro < nq )
    g_error("%s: output coefficient stride (%d) less than number of source "
	    " components (%d)", __FUNCTION__, cstro, nq) ;
  if ( cstri < nq )
    g_error("%s: input coefficient stride (%d) less than number of source "
	    " components (%d)", __FUNCTION__, cstro, nq) ;
  
  if ( !complex ) {
    sgnn = 1 ;
    for ( n = 0 ; n <= No ; n ++ ) {
      m = 0 ;
      offp = 2*wbfmm_coefficient_index_nm(n, m)*cstro ;
      for ( j = 0 ; j < 2*nq ; j ++ ) Co[offp+j] *= sc ;
      sgn = sgnn ;
      for ( l = m ; l < n ; (l ++), (sgn = -sgn) ) {
	idxp = wbfmm_coefficient_index_nm(l, m) ;
	idxc = wbfmm_coaxial_index_lmn(n,m,l) ;	
	for ( j = 0 ; j < nq ; j ++ ) {
	  Co[offp+2*j+0] += sgn*cfft[idxc]*Ci[2*idxp*cstri+2*j+0] ;
	  Co[offp+2*j+1] += sgn*cfft[idxc]*Ci[2*idxp*cstri+2*j+1] ;
	}
      }
      sgn = 1 ;
      for ( l = n ; l <= Ni ; l ++ ) {
	idxp = wbfmm_coefficient_index_nm(l, m) ;
	idxc = wbfmm_coaxial_index_lmn(l,m,n) ;
	for ( j = 0 ; j < nq ; j ++ ) {
	  Co[offp+2*j+0] += sgn*cfft[idxc]*Ci[2*idxp*cstri+2*j+0] ;
	  Co[offp+2*j+1] += sgn*cfft[idxc]*Ci[2*idxp*cstri+2*j+1] ;
	}
      }	

      sgnm = -sgnn ;
      for ( m = 1 ; m <= n ; m ++ ) {
	WBFMM_REAL buf[64] = {0.0} ;
	sgn = sgnm ;
	/*loop on input and coefficients*/
	for ( l = m ; l < n ; (l ++), (sgn = -sgn) ) {
	  idxp = wbfmm_coefficient_index_nm(l, m) ;
	  idxm = idxp - 2*m ;
	  idxc = wbfmm_coaxial_index_lmn(n,m,l) ;
	  for ( j = 0 ; j < nq ; j ++ ) {
	    buf[4*j+0] += sgn*cfft[idxc]*Ci[2*idxp*cstri+2*j+0] ;
	    buf[4*j+1] += sgn*cfft[idxc]*Ci[2*idxp*cstri+2*j+1] ;
	    buf[4*j+2] += sgn*cfft[idxc]*Ci[2*idxm*cstri+2*j+0] ;
	    buf[4*j+3] += sgn*cfft[idxc]*Ci[2*idxm*cstri+2*j+1] ;
	  }
	}
	sgn = 1 ;
	for ( l = n ; l <= Ni ; l ++ ) {
	  idxp = wbfmm_coefficient_index_nm(l, m) ;
	  idxm = idxp - 2*m ;
	  idxc = wbfmm_coaxial_index_lmn(l,m,n) ;
	  for ( j = 0 ; j < nq ; j ++ ) {	  
	    buf[4*j+0] += sgn*cfft[idxc]*Ci[2*idxp*cstri+2*j+0] ;
	    buf[4*j+1] += sgn*cfft[idxc]*Ci[2*idxp*cstri+2*j+1] ;
	    buf[4*j+2] += sgn*cfft[idxc]*Ci[2*idxm*cstri+2*j+0] ;
	    buf[4*j+3] += sgn*cfft[idxc]*Ci[2*idxm*cstri+2*j+1] ;
	  }
	}
	sgnm = -sgnm ;
	offp = 2*cstro*wbfmm_coefficient_index_nm(n,  m) ;
	offm = offp - 4*cstro*m ;
	for ( j = 0 ; j < nq ; j ++ ) {
	  Co[offp+2*j+0] = sc*Co[offp+2*j+0] + buf[4*j+0] ;
	  Co[offp+2*j+1] = sc*Co[offp+2*j+1] + buf[4*j+1] ;
	  Co[offm+2*j+0] = sc*Co[offm+2*j+0] + buf[4*j+2] ;
	  Co[offm+2*j+1] = sc*Co[offm+2*j+1] + buf[4*j+3] ;
	}
      }
      sgnn = -sgnn ;
    }
    return 0 ;
  }
  
  sgnn = 1 ;
  for ( n = 0 ; n <= No ; n ++ ) {
    m = 0 ; offp = 2*cstro*wbfmm_coefficient_index_nm(n, m) ;
    sgn = sgnn ;
    for ( j = 0 ; j < 2*nq ; j ++ ) Co[offp+j] *= sc ;
    for ( l = m ; l < n ; (l ++), (sgn = -sgn) ) {
      idxp = wbfmm_coefficient_index_nm(l, m) ;
      idxc = wbfmm_coaxial_index_lmn(n,m,l) ;
      for ( j = 0 ; j < nq ; j ++ ) {
	Co[offp+2*j+0] += sgn*
	  (cfft[2*idxc+0]*Ci[2*idxp*cstri+2*j+0] -
	   cfft[2*idxc+1]*Ci[2*idxp*cstri+2*j+1]) ;
	Co[offp+2*j+1] += sgn*
	  (cfft[2*idxc+1]*Ci[2*idxp*cstri+2*j+0] +
	   cfft[2*idxc+0]*Ci[2*idxp*cstri+2*j+1]) ;
      }
    }
    for ( l = n ; l <= Ni ; l ++ ) {
      idxp = wbfmm_coefficient_index_nm(l, m) ;
      idxc = wbfmm_coaxial_index_lmn(l,m,n) ;
      for ( j = 0 ; j < nq ; j ++ ) {
	Co[offp+2*j+0] += 
	  (cfft[2*idxc+0]*Ci[2*idxp*cstri+2*j+0] -
	   cfft[2*idxc+1]*Ci[2*idxp*cstri+2*j+1]) ;
	Co[offp+2*j+1] += 
	  (cfft[2*idxc+1]*Ci[2*idxp*cstri+2*j+0] +
	   cfft[2*idxc+0]*Ci[2*idxp*cstri+2*j+1]) ;
      }
    }

    sgnm = -sgnn ;
    for ( m = 1 ; m <= n ; m ++ ) {
      WBFMM_REAL buf[64] = {0.0} ;
      sgn = sgnm ;
      /*loop on input and coefficients*/
      for ( l = m ; l < n ; (l ++), (sgn = -sgn) ) {
	idxp = wbfmm_coefficient_index_nm(l, m) ;
	idxm = idxp - 2*m ;
	idxc = wbfmm_coaxial_index_lmn(n,m,l) ;
	for ( j = 0 ; j < nq ; j ++ ) {
	  buf[4*j+0] += sgn*
	    (cfft[2*idxc+0]*Ci[2*idxp*cstri+2*j+0] -
	     cfft[2*idxc+1]*Ci[2*idxp*cstri+2*j+1]) ;
	  buf[4*j+1] += sgn*
	    (cfft[2*idxc+1]*Ci[2*idxp*cstri+2*j+0] +
	     cfft[2*idxc+0]*Ci[2*idxp*cstri+2*j+1]) ;
	  buf[4*j+2] += sgn*
	    (cfft[2*idxc+0]*Ci[2*idxm*cstri+2*j+0] -
	     cfft[2*idxc+1]*Ci[2*idxm*cstri+2*j+1]) ;
	  buf[4*j+3] += sgn*
	    (cfft[2*idxc+1]*Ci[2*idxm*cstri+2*j+0] +
	     cfft[2*idxc+0]*Ci[2*idxm*cstri+2*j+1]) ;
	}
      }
      
      for ( l = n ; l <= Ni ; l ++ ) {
	idxp = wbfmm_coefficient_index_nm(l, m) ;
	idxm = idxp - 2*m ;
	idxc = wbfmm_coaxial_index_lmn(l,m,n) ;
	for ( j = 0 ; j < nq ; j ++ ) {
	  buf[4*j+0] +=
	    (cfft[2*idxc+0]*Ci[2*idxp*cstri+2*j+0] -
	     cfft[2*idxc+1]*Ci[2*idxp*cstri+2*j+1]) ;
	  buf[4*j+1] +=
	    (cfft[2*idxc+1]*Ci[2*idxp*cstri+2*j+0] +
	     cfft[2*idxc+0]*Ci[2*idxp*cstri+2*j+1]) ;
	  buf[4*j+2] +=
	    (cfft[2*idxc+0]*Ci[2*idxm*cstri+2*j+0] -
	     cfft[2*idxc+1]*Ci[2*idxm*cstri+2*j+1]) ;
	  buf[4*j+3] +=
	    (cfft[2*idxc+1]*Ci[2*idxm*cstri+2*j+0] +
	     cfft[2*idxc+0]*Ci[2*idxm*cstri+2*j+1]) ;
	}
      }
      offp = 2*cstro*wbfmm_coefficient_index_nm(n,  m) ;
      offm = offp - 4*cstro*m ;
      for ( j = 0 ; j < nq ; j ++ ) {
	  Co[offp+2*j+0] = sc*Co[offp+2*j+0] + buf[4*j+0] ;
	  Co[offp+2*j+1] = sc*Co[offp+2*j+1] + buf[4*j+1] ;
	  Co[offm+2*j+0] = sc*Co[offm+2*j+0] + buf[4*j+2] ;
	  Co[offm+2*j+1] = sc*Co[offm+2*j+1] + buf[4*j+3] ;
      }
      sgnm = -sgnm ;
    }
    sgnn = -sgnn ;
  }
  
  return 0 ;
}

/* @} */
