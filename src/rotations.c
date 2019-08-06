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

/*rotation operations using Gumerov and Duraiswami,
  http://dx.doi.org/10.1137/S1064827501399705 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <math.h>
#include <string.h>

#include <glib.h>

#include <wbfmm.h>

#include "wbfmm-private.h"

/* #define CHECK_COEFFICIENTS */

#ifndef WBFMM_SINGLE_PRECISION

/*AVX optimization for double precision calculations*/

#ifdef HAVE_AVX_INSTRUCTIONS
#include <immintrin.h>
#define WBFMM_USE_AVX
#endif

#endif /*WBFMM_SINGLE_PRECISION*/

#ifdef CHECK_COEFFICIENTS
#include <stdio.h>
#endif

/**
 * @file   rotations.c
 * @author Michael Carley <ensmjc@rpc-ensmjc.bath.ac.uk>
 * @date   Mon Jun 24 10:03:31 2019
 * 
 * @brief  Rotation coefficients and operations
 * 
 * Recursive computation of rotation coefficients using the methods of
 * Gumerov and Duraiswami, http://dx.doi.org/10.1137/S1064827501399705
 * 
 */

/** 
 * @defgroup rotations Rotation coefficients and operations
 *
 * @brief Computation and application of rotation operators
 *
 * Recursive computation of rotation coefficients using the methods of
 * Gumerov and Duraiswami, http://dx.doi.org/10.1137/S1064827501399705
 */

/* @{ */

#ifdef CHECK_COEFFICIENTS
static gint rotation_coefficients_check_recursion(WBFMM_REAL *H, WBFMM_REAL Cth,
						  WBFMM_REAL Sth, gint N)

{
  gint n, m, nu, idx1, idx2, idx3, idx4 ;
  WBFMM_REAL b1, b2, b3, a4, err, emax ;

  /*recursion (5.55)*/
  emax = 0.0 ;
  for ( n = 2 ; n <= N ; n ++ ) {
    for ( nu = -n+1 ; nu <= n-1 ; nu ++ ) {
      for ( m = 0 ; m <= n-1 ; m ++ ) {
	idx1 = wbfmm_rotation_index_numn(nu  , m+1, n-1) ;
	idx2 = wbfmm_rotation_index_numn(nu+1, m  , n  ) ;
	idx3 = wbfmm_rotation_index_numn(nu-1, m  , n  ) ;
	idx4 = wbfmm_rotation_index_numn(nu  , m  , n  ) ;
	b1 = recursion_bnm(n  ,  m   ) ;
	b2 = recursion_bnm(n  , -nu-1) ;
	b3 = recursion_bnm(n  ,  nu-1) ;
	a4 = recursion_anm(n-1,  nu  ) ;
	err = 2.0*b1*H[idx1] - b2*(1.0-Cth)*H[idx2] + b3*(1.0+Cth)*H[idx3] +
	  2.0*a4*Sth*H[idx4] ;
	emax = MAX(emax, ABS(err)) ;
	fprintf(stderr, "%2d %2d %2d %2d %+e (%lg,%lg,%lg,%lg)\n",
		n, nu, m, idx1, err, H[idx1], H[idx2], H[idx3], H[idx4]) ;
      }
    }
  }

  fprintf(stderr, "Maximum error in recursion relation (5.55): %e\n", emax) ;

  return 0 ;
}

#endif /*CHECK_COEFFICIENTS*/

gint FUNCTION_NAME(wbfmm_rotation_angles)(WBFMM_REAL *ix, WBFMM_REAL *iy,
					  WBFMM_REAL *iz, 
					  WBFMM_REAL *jx, WBFMM_REAL *jy,
					  WBFMM_REAL *jz, 
					  WBFMM_REAL *th, WBFMM_REAL *ph,
					  WBFMM_REAL *ch)

/*
  rotation angles in G&D section 5
  
  (ix,iy,iz) unit vectors on coordinate axes in initial system
  (jx,jy,jz) unit vectors on coordinate axes in rotated system
*/

{
  WBFMM_REAL ex, ey, ez ;

  ex = ix[0]*jz[0] + ix[1]*jz[1] + ix[2]*jz[2] ;
  ey = iy[0]*jz[0] + iy[1]*jz[1] + iy[2]*jz[2] ;
  ez = iz[0]*jz[0] + iz[1]*jz[1] + iz[2]*jz[2] ;

  *th = ACOS(ez) ;
  *ch = ATAN2(ey, ex) ;

  ex = jx[0]*iz[0] + jx[1]*iz[1] + jx[2]*iz[2] ;
  ey = jy[0]*iz[0] + jy[1]*iz[1] + jy[2]*iz[2] ;

  *ph = ATAN2(ey, ex) ;
  
  return 0 ;
}

gint FUNCTION_NAME(wbfmm_coefficients_H_rotation)(WBFMM_REAL *H, gint N, 
						  WBFMM_REAL th, 
						  WBFMM_REAL *work)

{
  WBFMM_REAL Cth, Sth, *Pnm1, *Pn, b1, b2, b, a3, *h ;
  gint idx, n, m, nu, idx1, idx2, idx3 ;

  /*note that with the G&D normalization, the Legendre polynomials in
    (5.48) are related to the normalized Legendre polynomials by a
    scaling (-1)^m*sqrt(2n+1)/sqrt(4\pi)*/
  Pnm1 = &(work[0]) ; Pn = &(Pnm1[2*(2*N+1)]) ;
  h = &(Pn[2*(2*N+1)]) ;

  Cth = COS(th) ; Sth = SIN(th) ;

  /*initialize the m=0 entries of H, G&D (5.48)*/
  FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth, &(Pnm1[0]), &(Pn[0]), &(Pn[1])) ;
  n = 0 ; m = 0 ;
  nu = 0 ;
  idx = wbfmm_rotation_index_numn(nu, m, n) ;
  h[idx] = Pnm1[nu]/sqrt(2*n+1)*sqrt(4.0*M_PI) ;

  n = 1 ;
  for ( nu = 0 ; nu <= n ; nu ++ ) {
    idx  = wbfmm_rotation_index_numn( nu, m, n) ;
    idx1 = wbfmm_rotation_index_numn(-nu, m, n) ;
    h[idx] = h[idx1] = Pn[nu]/sqrt(2*n+1)*sqrt(4.0*M_PI) ;
  }

  for ( n = 2 ; n <= 2*N ; n ++ ) {
    FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pnm1, &Pn, n-1, Cth, Sth) ;
    for ( nu = 0 ; nu <= n ; nu ++ ) {
      idx  = wbfmm_rotation_index_numn( nu, m, n) ;
      idx1 = wbfmm_rotation_index_numn(-nu, m, n) ;
      h[idx] = h[idx1] = Pn[nu]/sqrt(2*n+1)*sqrt(4.0*M_PI) ;
    }
  }

  /*apply recursion 5.55, in this order to advance (nu,n) planes in
    the direction of increasing m*/
  for ( m = 0 ; m <= N ; m ++ ) {
    for ( n = m+2 ; n <= 2*N-m ; n ++ ) {
      for ( nu = -n+1 ; nu <= n-1 ; nu ++ ) {
	b  = FUNCTION_NAME(recursion_bnm)(n, m) ;
	b1 = FUNCTION_NAME(recursion_bnm)(n, -nu-1) ;
	b2 = FUNCTION_NAME(recursion_bnm)(n, nu-1) ;
	a3 = FUNCTION_NAME(recursion_anm)(n-1, nu) ;
	idx  = wbfmm_rotation_index_numn(nu  , m+1, n-1) ;
	idx1 = wbfmm_rotation_index_numn(nu+1, m  , n  ) ;
	idx2 = wbfmm_rotation_index_numn(nu-1, m  , n  ) ;
	idx3 = wbfmm_rotation_index_numn(nu  , m  , n  ) ;
	h[idx] = (0.5*(b1*(1.0 - Cth)*h[idx1] - b2*(1.0 + Cth)*h[idx2]) -
		  a3*Sth*h[idx3])/b ;
      }
    }
  }

  memcpy(H, h, wbfmm_rotation_index_numn(N+1,0,N+1)*sizeof(WBFMM_REAL)) ;

#ifdef CHECK_COEFFICIENTS
  rotation_coefficients_check_recursion(H, Cth, Sth, N) ;

#endif /*CHECK_COEFFICIENTS*/

  return 0 ;
}

#ifdef WBFMM_USE_AVX
gint FUNCTION_NAME(wbfmm_rotate_H)(WBFMM_REAL *Co, gint cstro, 
				   gint N, WBFMM_REAL *Ci, gint cstri,
				   WBFMM_REAL *H,
				   WBFMM_REAL ph, WBFMM_REAL ch)

/*
  apply rotation (matrix H from wbfmm_coefficients_H_rotation) to
  rotate input coefficients Ci into output Co, through angles
  (th,ph,ch), G&D, section 6, and (2.27)
*/

{
  gsize nu, n, m, offp, offm ;
  WBFMM_REAL Cmch, Smch, Cnph, Snph, Cch, Sch, Cph, Sph ;
  WBFMM_REAL tmp, Hp, Hm, CC, SS, CS, SC ;

  /*initialize recursions*/
  Cph = COS(ph) ; Sph = SIN(ph) ;
  Cch = COS(ch) ; Sch = SIN(ch) ;

  /* inside loops, trigonmetric quantities are calculated using
   * recursions and take the following values:
   *
   * Smch = SIN(m*ch) ; Cmch = COS(m*ch) ;
   * Cnph = COS(nu*ph) ; Snph = SIN(nu*ph) ;
   * 
   * Er + j Ei = \exp(j(\pm m\chi - \pm \nu\phi))
   * Er = COS(m*ch-nu*ph) ; Ei = SIN(m*ch-nu*ph)
   *
   * CC = COS(m*ch)*COS(nu*ph) 
   * SC = SIN(m*ch)*COS(nu*ph) 
   * CS = COS(m*ch)*SIN(nu*ph) 
   * SS = SIN(m*ch)*SIN(nu*ph) 
   *
   * offX (X = `p', `m') = offset into array, `p' for `plus' indices,
   * `m' for `minus'
   */

  for ( n = 0 ; n <= N ; n ++ ) {
    {
      __attribute__ ((aligned (32))) WBFMM_REAL tmul[10]={0.0} ;
      __m256d ECp0, ECp1, ECm0, ECm1, op1, En ;
      
      /* ECp0 = _mm256_set1_pd(0.0) ; ECp1 = _mm256_set1_pd(0.0) ; */
      ECm0 = _mm256_set1_pd(0.0) ; ECm1 = _mm256_set1_pd(0.0) ;

      nu = 0 ; Cnph = 1.0 ; Snph = 0.0 ;

      m = 0 ; Cmch = 1.0 ; Smch = 0.0 ;

      offm = offp = 2*cstri*wbfmm_coefficient_index_nm(n,m) ;

      Hp = H[wbfmm_rotation_index_numn(nu,m,n)] ;
      
      CC = 1.0 ; SS = 0.0 ; CS = 0.0 ; SC = 0.0 ;

      En = _mm256_set_pd(0.0, 0.0, 0.0, Hp) ;
#ifdef HAVE_FMA_INSTRUCTIONS
      op1 = _mm256_set1_pd(Ci[offp+0]) ;
      ECp0 = _mm256_mul_pd(op1, En) ;
      /* ECp0 = _mm256_fmadd_pd(op1, En, ECp0) ; */
      
      op1 = _mm256_set1_pd(Ci[offp+1]) ;
      ECp1 = _mm256_mul_pd(op1, En) ;
      /* ECp1 = _mm256_fmadd_pd(op1, En, ECp1) ; */
#else /*HAVE_FMA_INSTRUCTIONS*/
      op1 = _mm256_set1_pd(Ci[offp+0]) ;
      /* op1 = _mm256_mul_pd(op1, En) ; */
      /* ECp0 = _mm256_add_pd(op1, ECp0) ; */
      ECp0 = _mm256_mul_pd(op1, En) ;

      op1 = _mm256_set1_pd(Ci[offp+1]) ;
      ECp1 = _mm256_mul_pd(op1, En) ;
      /* op1 = _mm256_mul_pd(op1, En) ; */
      /* ECp1 = _mm256_add_pd(op1, ECp1) ; */
#endif /*HAVE_FMA_INSTRUCTIONS*/

      for ( m = 1 ; m <= n ; m ++ ) {
	Hp = H[wbfmm_rotation_index_numn( nu,m,n)] ;
	Hm = H[wbfmm_rotation_index_numn(-nu,m,n)] ;

	/* offp = 2*cstri*wbfmm_coefficient_index_nm(n,m) ; */
	/* offm = 2*cstri*wbfmm_coefficient_index_nm(n,-m) ; */
	offp += 2*cstri ; offm -= 2*cstri ;
	
	tmp = Cmch ; 
	Cmch = Cmch*Cch - Smch*Sch ;
	Smch = Smch*Cch + tmp*Sch ;
	
	CC = Cmch*Cnph ; SS = Smch*Snph ;
	CS = Cmch*Snph ; SC = Smch*Cnph ;
	
	En = _mm256_set_pd(Hm*(SC + CS), Hm*(CC - SS),
			   Hp*(SC - CS), Hp*(CC + SS)) ;
#ifdef HAVE_FMA_INSTRUCTIONS
	op1 = _mm256_set1_pd(Ci[offp+0]) ;
	ECp0 = _mm256_fmadd_pd(op1, En, ECp0) ;

	op1 = _mm256_set1_pd(Ci[offp+1]) ;
	ECp1 = _mm256_fmadd_pd(op1, En, ECp1) ;

	op1 = _mm256_set1_pd(Ci[offm+0]) ;
	ECm0 = _mm256_fmadd_pd(op1, En, ECm0) ;

	op1 = _mm256_set1_pd(Ci[offm+1]) ;
	ECm1 = _mm256_fmadd_pd(op1, En, ECm1) ;	
#else /*HAVE_FMA_INSTRUCTIONS*/
	op1 = _mm256_set1_pd(Ci[offp+0]) ;
	op1 = _mm256_mul_pd(op1, En) ;
	ECp0 = _mm256_add_pd(op1, ECp0) ;

	op1 = _mm256_set1_pd(Ci[offp+1]) ;
	op1 = _mm256_mul_pd(op1, En) ;
	ECp1 = _mm256_add_pd(op1, ECp1) ;

	op1 = _mm256_set1_pd(Ci[offm+0]) ;
	op1 = _mm256_mul_pd(op1, En) ;
	ECm0 = _mm256_add_pd(op1, ECm0) ;

	op1 = _mm256_set1_pd(Ci[offm+1]) ;
	op1 = _mm256_mul_pd(op1, En) ;
	ECm1 = _mm256_add_pd(op1, ECm1) ;
#endif /*HAVE_FMA_INSTRUCTIONS*/
      }
      
      offp = 2*cstro*wbfmm_coefficient_index_nm(n, nu) ;
      _mm256_store_pd(&(tmul[0]), ECp0) ;
      _mm256_store_pd(&(tmul[2]), ECp1) ;
      _mm256_store_pd(&(tmul[4]), ECm0) ;
      _mm256_store_pd(&(tmul[6]), ECm1) ;
      Co[offp+0] += tmul[ 0] - tmul[ 3] + tmul[ 4] + tmul[ 7] ;
      Co[offp+1] += tmul[ 2] + tmul[ 1] + tmul[ 6] - tmul[ 5] ;
    }
    
    for ( nu = 1 ; nu <= n ; nu ++ ) {
      __attribute__ ((aligned (32))) WBFMM_REAL tmul[16]={0.0} ;
      __m256d ECp0, ECp1, ECm0, ECm1, op1, En ;
      
      /* ECp0 = _mm256_set1_pd(0.0) ; ECp1 = _mm256_set1_pd(0.0) ; */
      ECm0 = _mm256_set1_pd(0.0) ; ECm1 = _mm256_set1_pd(0.0) ;

      tmp = Cnph ; 
      Cnph = Cnph*Cph - Snph*Sph ;
      Snph = Snph*Cph + tmp*Sph ;

      m = 0 ; Cmch = 1.0 ; Smch = 0.0 ;

      offm = offp = 2*cstri*wbfmm_coefficient_index_nm(n,m) ;

      Hp = H[wbfmm_rotation_index_numn( nu,m,n)] ;
      Hm = H[wbfmm_rotation_index_numn(-nu,m,n)] ;

      /* CC = Cmch*Cnph ; SS = Smch*Snph ; */
      /* CS = Cmch*Snph ; SC = Smch*Cnph ; */
      CC =      Cnph ; SS = 0.0 ;
      CS =      Snph ; SC = 0.0 ;

      En = _mm256_set_pd(Hm*CS, Hm*CC, -Hp*CS, Hp*CC) ;
#ifdef HAVE_FMA_INSTRUCTIONS
      op1 = _mm256_set1_pd(Ci[offp+0]) ;
      ECp0 = _mm256_mul_pd(op1, En) ;
      /* ECp0 = _mm256_fmadd_pd(op1, En, ECp0) ; */
      
      op1 = _mm256_set1_pd(Ci[offp+1]) ;
      ECp1 = _mm256_mul_pd(op1, En) ;
      /* ECp1 = _mm256_fmadd_pd(op1, En, ECp1) ; */
#else /*HAVE_FMA_INSTRUCTIONS*/
      op1 = _mm256_set1_pd(Ci[offp+0]) ;
      ECp0 = _mm256_mul_pd(op1, En) ;
      /* op1 = _mm256_mul_pd(op1, En) ; */
      /* ECp0 = _mm256_add_pd(op1, ECp0) ; */

      op1 = _mm256_set1_pd(Ci[offp+1]) ;
      ECp1 = _mm256_mul_pd(op1, En) ;
      /* op1 = _mm256_mul_pd(op1, En) ; */
      /* ECp1 = _mm256_add_pd(op1, ECp1) ; */
#endif /*HAVE_FMA_INSTRUCTIONS*/
      
      for ( m = 1 ; m <= n ; m ++ ) {
	/*rotation coefficients for \pm\nu*/
	Hp = H[wbfmm_rotation_index_numn( nu,m,n)] ;
	Hm = H[wbfmm_rotation_index_numn(-nu,m,n)] ;

	/* g_assert(wbfmm_rotation_index_numn( nu,m,n) - */
	/* 	 wbfmm_rotation_index_numn(-nu,m,n) == 2*nu*(n+1)) ; */
	
	/* offp = 2*cstri*wbfmm_coefficient_index_nm(n, m) ; */
	/* offm = 2*cstri*wbfmm_coefficient_index_nm(n,-m) ; */
	offp += 2*cstri ; offm -= 2*cstri ;

	tmp = Cmch ; 
	Cmch = Cmch*Cch - Smch*Sch ;
	Smch = Smch*Cch + tmp*Sch ;

	CC = Cmch*Cnph ; SS = Smch*Snph ;
	CS = Cmch*Snph ; SC = Smch*Cnph ;

	En = _mm256_set_pd(Hm*(SC + CS), Hm*(CC - SS),
			   Hp*(SC - CS), Hp*(CC + SS)) ;

#ifdef HAVE_FMA_INSTRUCTIONS
	op1 = _mm256_set1_pd(Ci[offp+0]) ;
	ECp0 = _mm256_fmadd_pd(op1, En, ECp0) ;

	op1 = _mm256_set1_pd(Ci[offp+1]) ;
	ECp1 = _mm256_fmadd_pd(op1, En, ECp1) ;

	op1 = _mm256_set1_pd(Ci[offm+0]) ;
	ECm0 = _mm256_fmadd_pd(op1, En, ECm0) ;

	op1 = _mm256_set1_pd(Ci[offm+1]) ;
	ECm1 = _mm256_fmadd_pd(op1, En, ECm1) ;	
#else /*HAVE_FMA_INSTRUCTIONS*/
	op1 = _mm256_set1_pd(Ci[offp+0]) ;
	op1 = _mm256_mul_pd(op1, En) ;
	ECp0 = _mm256_add_pd(op1, ECp0) ;

	op1 = _mm256_set1_pd(Ci[offp+1]) ;
	op1 = _mm256_mul_pd(op1, En) ;
	ECp1 = _mm256_add_pd(op1, ECp1) ;

	op1 = _mm256_set1_pd(Ci[offm+0]) ;
	op1 = _mm256_mul_pd(op1, En) ;
	ECm0 = _mm256_add_pd(op1, ECm0) ;

	op1 = _mm256_set1_pd(Ci[offm+1]) ;
	op1 = _mm256_mul_pd(op1, En) ;
	ECm1 = _mm256_add_pd(op1, ECm1) ;
#endif /*HAVE_FMA_INSTRUCTIONS*/
      }

      /*put the accumulated results back into tmul*/
      _mm256_store_pd(&(tmul[0]), ECp0) ;
      _mm256_store_pd(&(tmul[4]), ECp1) ;
      _mm256_store_pd(&(tmul[8]), ECm0) ;
      _mm256_store_pd(&(tmul[12]), ECm1) ;
      
      /*output indices for \pm\nu*/
      offp = 2*cstro*wbfmm_coefficient_index_nm(n, nu) ;
      offm = offp - 4*cstro*nu ;
      /* offm = 2*cstro*wbfmm_coefficient_index_nm(n,-nu) ; */
      Co[offp+0] += tmul[ 0] - tmul[ 5] + tmul[10] + tmul[15] ;
      Co[offp+1] += tmul[ 4] + tmul[ 1] + tmul[14] - tmul[11] ;
      Co[offm+0] += tmul[ 8] + tmul[13] + tmul[ 2] - tmul[ 7] ;
      Co[offm+1] += tmul[12] - tmul[ 9] + tmul[ 6] + tmul[ 3] ;
    }
  }

  return 0 ;
}

#else /*WBFMM_USE_AVX*/

gint FUNCTION_NAME(wbfmm_rotate_H)(WBFMM_REAL *Co, gint cstro, 
				   gint N, WBFMM_REAL *Ci, gint cstri,
				   WBFMM_REAL *H,
				   WBFMM_REAL ph, WBFMM_REAL ch)

/*
  apply rotation (matrix H from wbfmm_coefficients_H_rotation) to
  rotate input coefficients Ci into output Co, through angles
  (th,ph,ch), G&D, section 6, and (2.27)
*/

{
  gint nu, n, m, offp, offm ;
  WBFMM_REAL Cmch, Smch, Cnph, Snph, Cch, Sch, Cph, Sph ;
  WBFMM_REAL tmp, Hp, Hm, CC, SS, CS, SC, E[4], tmul[16] ;

  /*initialize recursions*/
  Cph = COS(ph) ; Sph = SIN(ph) ;
  Cch = COS(ch) ; Sch = SIN(ch) ;

  /* inside loops, trigonmetric quantities are calculated using
   * recursions and take the following values:
   *
   * Smch = SIN(m*ch) ; Cmch = COS(m*ch) ;
   * Cnph = COS(nu*ph) ; Snph = SIN(nu*ph) ;
   * 
   * Er + j Ei = \exp(j(\pm m\chi - \pm \nu\phi))
   * Er = COS(m*ch-nu*ph) ; Ei = SIN(m*ch-nu*ph)
   *
   * CC = COS(m*ch)*COS(nu*ph) 
   * SC = SIN(m*ch)*COS(nu*ph) 
   * CS = COS(m*ch)*SIN(nu*ph) 
   * SS = SIN(m*ch)*SIN(nu*ph) 
   *
   * offX (X = `p', `m') = offset into array, `p' for `plus' indices,
   * `m' for `minus'
   */

  for ( n = 0 ; n <= N ; n ++ ) {
    WBFMM_REAL buf[4] ;
    nu = 0 ; Cnph = 1.0 ; Snph = 0.0 ;

    m = 0 ; Cmch = 1.0 ; Smch = 0.0 ;

    offp = 2*cstri*wbfmm_coefficient_index_nm(n,m) ;
    Hp = H[wbfmm_rotation_index_numn(nu,m,n)] ;

    CC = Cmch*Cnph ; SS = Smch*Snph ;
    CS = Cmch*Snph ; SC = Smch*Cnph ;

    E[0] = Hp*CC + Hp*SS ; 
    E[1] = Hp*SC - Hp*CS ; 

    buf[0] = E[0]*Ci[offp+0] - E[1]*Ci[offp+1] ;
    buf[1] = E[0]*Ci[offp+1] + E[1]*Ci[offp+0] ;

    for ( m = 1 ; m <= n ; m ++ ) {
      Hp = H[wbfmm_rotation_index_numn( nu,m,n)] ;
      Hm = H[wbfmm_rotation_index_numn(-nu,m,n)] ;

      offp = 2*cstri*wbfmm_coefficient_index_nm(n,m) ;
      offm = 2*cstri*wbfmm_coefficient_index_nm(n,-m) ;

      tmp = Cmch ; 
      Cmch = Cmch*Cch - Smch*Sch ;
      Smch = Smch*Cch + tmp*Sch ;

      CC = Cmch*Cnph ; SS = Smch*Snph ;
      CS = Cmch*Snph ; SC = Smch*Cnph ;

      E[0] = Hp*CC + Hp*SS ; 
      E[1] = Hp*SC - Hp*CS ; 
      E[2] = Hm*CC - Hm*SS ; 
      E[3] = Hm*SC + Hm*CS ; 

      buf[0] += E[0]*Ci[offp+0] - E[1]*Ci[offp+1] +
	E[2]*Ci[offm+0] + E[3]*Ci[offm+1] ;
      buf[1] += E[0]*Ci[offp+1] + E[1]*Ci[offp+0] +
	E[2]*Ci[offm+1] - E[3]*Ci[offm+0] ;
    }

    offp = 2*cstro*wbfmm_coefficient_index_nm(n, nu) ;
    Co[offp+0] += buf[0] ; Co[offp+1] += buf[1] ;

    for ( nu = 1 ; nu <= n ; nu ++ ) {
      tmp = Cnph ; 
      Cnph = Cnph*Cph - Snph*Sph ;
      Snph = Snph*Cph + tmp*Sph ;

      m = 0 ; Cmch = 1.0 ; Smch = 0.0 ;

      offp = 2*cstri*wbfmm_coefficient_index_nm(n,m) ;

      Hp = H[wbfmm_rotation_index_numn( nu,m,n)] ;
      Hm = H[wbfmm_rotation_index_numn(-nu,m,n)] ;

      CC = Cmch*Cnph ; SS = Smch*Snph ;
      CS = Cmch*Snph ; SC = Smch*Cnph ;

      E[0] = Hp*CC + Hp*SS ; 
      E[1] = Hp*SC - Hp*CS ; 
      E[2] = Hm*CC - Hm*SS ; 
      E[3] = Hm*SC + Hm*CS ; 

      tmul[ 0] = E[0]*Ci[offp+0] ;
      tmul[ 1] = E[1]*Ci[offp+0] ;
      tmul[ 2] = E[2]*Ci[offp+0] ;
      tmul[ 3] = E[3]*Ci[offp+0] ;
      tmul[ 4] = E[0]*Ci[offp+1] ;
      tmul[ 5] = E[1]*Ci[offp+1] ;
      tmul[ 6] = E[2]*Ci[offp+1] ;
      tmul[ 7] = E[3]*Ci[offp+1] ;
      tmul[8] = tmul[9] = tmul[10] = tmul[11] =
      	tmul[12] = tmul[13] = tmul[14] = tmul[15] = 0.0 ;
	
      for ( m = 1 ; m <= n ; m ++ ) {
	/*rotation coefficients for \pm\nu*/
	Hp = H[wbfmm_rotation_index_numn( nu,m,n)] ;
	Hm = H[wbfmm_rotation_index_numn(-nu,m,n)] ;

	offp = 2*cstri*wbfmm_coefficient_index_nm(n, m) ;
	offm = 2*cstri*wbfmm_coefficient_index_nm(n,-m) ;

	tmp = Cmch ; 
	Cmch = Cmch*Cch - Smch*Sch ;
	Smch = Smch*Cch + tmp*Sch ;

	CC = Cmch*Cnph ; SS = Smch*Snph ;
	CS = Cmch*Snph ; SC = Smch*Cnph ;

	E[0] = Hp*CC + Hp*SS ; 
	E[1] = Hp*SC - Hp*CS ; 
	E[2] = Hm*CC - Hm*SS ; 
	E[3] = Hm*SC + Hm*CS ; 

	tmul[ 0] += E[0]*Ci[offp+0] ;
	tmul[ 1] += E[1]*Ci[offp+0] ;
	tmul[ 2] += E[2]*Ci[offp+0] ;
	tmul[ 3] += E[3]*Ci[offp+0] ;
	tmul[ 4] += E[0]*Ci[offp+1] ;
	tmul[ 5] += E[1]*Ci[offp+1] ;
	tmul[ 6] += E[2]*Ci[offp+1] ;
	tmul[ 7] += E[3]*Ci[offp+1] ;
	tmul[ 8] += E[0]*Ci[offm+0] ;
	tmul[ 9] += E[1]*Ci[offm+0] ;
	tmul[10] += E[2]*Ci[offm+0] ;
	tmul[11] += E[3]*Ci[offm+0] ;
	tmul[12] += E[0]*Ci[offm+1] ;
	tmul[13] += E[1]*Ci[offm+1] ;
	tmul[14] += E[2]*Ci[offm+1] ;
	tmul[15] += E[3]*Ci[offm+1] ;
      }
      /*output indices for \pm\nu*/
      offp = 2*cstro*wbfmm_coefficient_index_nm(n, nu) ;
      offm = 2*cstro*wbfmm_coefficient_index_nm(n,-nu) ;
      Co[offp+0] += tmul[ 0] - tmul[ 5] + tmul[10] + tmul[15] ;
      Co[offp+1] += tmul[ 4] + tmul[ 1] + tmul[14] - tmul[11] ;
      Co[offm+0] += tmul[ 8] + tmul[13] + tmul[ 2] - tmul[ 7] ;
      Co[offm+1] += tmul[12] - tmul[ 9] + tmul[ 6] + tmul[ 3] ;
    }
  }

  return 0 ;
}

#endif /*WBFMM_USE_AVX*/

/* @} */
