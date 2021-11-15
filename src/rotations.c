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

#include <blaswrap.h>

#include <wbfmm.h>

#include "wbfmm-private.h"
#include "wbfmm-avx.h"

#ifdef HAVE_AVX_INSTRUCTIONS
#include <immintrin.h>
#endif /*HAVE_AVX_INSTRUCTIONS*/

/* #define CHECK_COEFFICIENTS */


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

gint WBFMM_FUNCTION_NAME(wbfmm_rotation_angles)(WBFMM_REAL *ix, WBFMM_REAL *iy,
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

gint WBFMM_FUNCTION_NAME(wbfmm_coefficients_H_rotation)(WBFMM_REAL *H, gint N, 
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
  WBFMM_FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth,
					   &(Pnm1[0]), &(Pn[0]), &(Pn[1])) ;
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
    WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pnm1, &Pn, n-1,
							Cth, Sth) ;
    for ( nu = 0 ; nu <= n ; nu ++ ) {
      idx  = wbfmm_rotation_index_numn( nu, m, n) ;
      idx1 = wbfmm_rotation_index_numn(-nu, m, n) ;
      h[idx] = h[idx1] = Pn[nu]/sqrt(2*n+1)*sqrt(4.0*M_PI) ;
#ifdef WBFMM_CHECK_ISNAN
      if ( isnan(h[idx]) )
      	g_error("%s: NaN at (nu,m,n)==(%d,%d,%d)", __FUNCTION__, nu, m, n) ;
#endif /*WBFMM_CHECK_ISNAN*/
    }
  }

  /*apply recursion 5.55, in this order to advance (nu,n) planes in
    the direction of increasing m*/
  for ( m = 0 ; m <= N ; m ++ ) {
    for ( n = m+2 ; n <= 2*N-m ; n ++ ) {
      for ( nu = -n+1 ; nu <= n-1 ; nu ++ ) {
	b  = WBFMM_FUNCTION_NAME(recursion_bnm)(n, m) ;
	b1 = WBFMM_FUNCTION_NAME(recursion_bnm)(n, -nu-1) ;
	b2 = WBFMM_FUNCTION_NAME(recursion_bnm)(n, nu-1) ;
	a3 = WBFMM_FUNCTION_NAME(recursion_anm)(n-1, nu) ;
	idx  = wbfmm_rotation_index_numn(nu  , m+1, n-1) ;
	idx1 = wbfmm_rotation_index_numn(nu+1, m  , n  ) ;
	idx2 = wbfmm_rotation_index_numn(nu-1, m  , n  ) ;
	idx3 = wbfmm_rotation_index_numn(nu  , m  , n  ) ;
	h[idx] = (0.5*(b1*(1.0 - Cth)*h[idx1] - b2*(1.0 + Cth)*h[idx2]) -
		  a3*Sth*h[idx3])/b ;
#ifdef WBFMM_CHECK_ISNAN
	if ( isnan(h[idx]) )
	  g_error("%s: NaN at (nu,m,n)==(%d,%d,%d)", __FUNCTION__, nu, m, n) ;
#endif /*WBFMM_CHECK_ISNAN*/
      }
    }
  }

  memcpy(H, h, wbfmm_rotation_index_numn(N+1,0,N+1)*sizeof(WBFMM_REAL)) ;

#ifdef CHECK_COEFFICIENTS
  rotation_coefficients_check_recursion(H, Cth, Sth, N) ;

#endif /*CHECK_COEFFICIENTS*/

  return 0 ;
}

#ifdef HAVE_AVX_INSTRUCTIONS
#ifndef WBFMM_SINGLE_PRECISION
/*only compile AVX for double precision functions, for now anyway*/

gint WBFMM_FUNCTION_NAME(wbfmm_rotate_H_avx)(WBFMM_REAL *Co, gint cstro, 
					     WBFMM_REAL *Ci, gint cstri,
					     gint N,
					     gint nq,
					     WBFMM_REAL *H,
					     WBFMM_REAL ph, WBFMM_REAL ch,
					     WBFMM_REAL sc)

/*
  apply rotation (matrix H from wbfmm_coefficients_H_rotation) to
  rotate input coefficients Ci into output Co, through angles
  (th,ph,ch), G&D, section 6, and (2.27)
*/

{
  gint nu, n, m, offp, offm ;
  WBFMM_REAL Cnph, Snph, Cch, Sch, Cph, Sph ;
  __m256d rCch, rSch, rsgn ;
  
  g_assert(nq == 1) ;
  
  /*initialize recursions*/
  Cph = COS(ph) ; Sph = SIN(ph) ;
  Cch = COS(ch) ; Sch = SIN(ch) ;

  /* inside loops, trigonmetric quantities are calculated using
   * recursions and take the following values:
   *
   * Smch = SIN(m*ch) ; Cmch = COS(m*ch) ;
   * Cnph = COS(nu*ph) ; Snph = SIN(nu*ph) ;
   * 
   * offX (X = `p', `m') = offset into array, `p' for `plus' indices,
   * `m' for `minus'
   *
   * `r' indicates AVX register corresponding to a given `normal'
   * variable
   */
  
  rCch = _mm256_set1_pd(Cch) ; rSch = _mm256_set1_pd(Sch) ;
  rsgn = _mm256_set_pd(-1.0, 1.0, 1.0, -1.0) ;

  for ( n = 0 ; n <= N ; n ++ ) {
    __attribute__ ((aligned (32))) WBFMM_REAL tmul[4] ;
    __m256d rtC, rtS, rCmch, rSmch, rHm, rHp, rtmul ;
#ifndef HAVE_FMA_INSTRUCTIONS
    __m256d op1 ;
#endif /*HAVE_FMA_INSTRUCTIONS*/
    
    nu = 0 ; Cnph = 1.0 ; Snph = 0.0 ;

    m = 0 ; rCmch = _mm256_set1_pd(1.0) ; rSmch = _mm256_setzero_pd() ;

    rHp = _mm256_broadcast_sd(&(H[wbfmm_rotation_index_numn( nu,m,n)])) ;

    offp = 2*cstri*wbfmm_coefficient_index_nm(n,m) ;
    rtC = _mm256_set_pd(0.0, 0.0, Ci[offp+1], Ci[offp+0]) ;

    rtmul = _mm256_mul_pd(rtC, rHp) ;
    
    for ( m = 1 ; m <= n ; m ++ ) {
      rHp = _mm256_broadcast_sd(&(H[wbfmm_rotation_index_numn( nu,m,n)])) ;
      rHm = _mm256_broadcast_sd(&(H[wbfmm_rotation_index_numn(-nu,m,n)])) ;
      wbfmm_cos_sin_recursion_avx(rCmch,rSmch,rCch,rSch) ;
      
      offp = 2*cstri*wbfmm_coefficient_index_nm(n, m) ;
      offm = offp - 4*cstri*m ;

      rtC = _mm256_set_pd(Ci[offm+1], Ci[offm+0], Ci[offp+1], Ci[offp+0]) ;
      rtS = _mm256_permute_pd(rtC, 5) ;
      rtS = _mm256_mul_pd(rtS, rsgn) ;
      rtC = _mm256_mul_pd(rtC, rCmch) ;

#ifdef HAVE_FMA_INSTRUCTIONS
      rtC    = _mm256_fmadd_pd(rtS, rSmch, rtC) ;
      rtmul  = _mm256_fmadd_pd(rtC, rHp, rtmul ) ;      
      rtC    = _mm256_permute2f128_pd(rtC, rtC, 1) ;
      rtmul  = _mm256_fmadd_pd(rtC, rHm, rtmul ) ;      
#else /*HAVE_FMA_INSTRUCTIONS*/
      rtS = _mm256_mul_pd(rtS, rSmch) ;

      op1 = _mm256_mul_pd(rtC, rHp) ;
      rtmul  = _mm256_add_pd(rtmul,  op1) ;
      op1 = _mm256_mul_pd(rtS, rHp) ;
      rtmul  = _mm256_add_pd(rtmul,  op1) ;

      rtC = _mm256_permute2f128_pd(rtC, rtC, 1) ;
      rtS = _mm256_permute2f128_pd(rtS, rtS, 1) ;

      op1 = _mm256_mul_pd(rtC, rHm) ;
      rtmul  = _mm256_add_pd(rtmul,  op1) ;
      op1 = _mm256_mul_pd(rtS, rHm) ;
      rtmul  = _mm256_add_pd(rtmul,  op1) ;
#endif /*HAVE_FMA_INSTRUCTIONS*/
    }

    _mm256_store_pd(tmul, rtmul) ;
    
    offp = 2*cstro*wbfmm_coefficient_index_nm(n, nu) ;
    Co[offp+0] = sc*Co[offp+0] + tmul[0]*Cnph ;
    Co[offp+1] = sc*Co[offp+1] + tmul[1]*Cnph ;
    
    for ( nu = 1 ; nu <= n ; nu ++ ) {
      wbfmm_cos_sin_recursion(Cnph,Snph,Cph,Sph) ;

      m = 0 ; rCmch = _mm256_set1_pd(1.0) ; rSmch = _mm256_setzero_pd() ;

      rHp = _mm256_broadcast_sd(&(H[wbfmm_rotation_index_numn( nu,m,n)])) ;
      rHm = _mm256_broadcast_sd(&(H[wbfmm_rotation_index_numn(-nu,m,n)])) ;

      offp = 2*cstri*wbfmm_coefficient_index_nm(n,m) ;

      rtC = _mm256_set_pd(0.0, 0.0, Ci[offp+1], Ci[offp+0]) ;
      rtC = _mm256_mul_pd(rtC, rCmch) ;

      rtmul = _mm256_mul_pd(rtC, rHp) ;
      rtC = _mm256_permute2f128_pd(rtC, rtC, 1) ;
#ifdef HAVE_FMA_INSTRUCTIONS
      rtmul = _mm256_fmadd_pd(rtC, rHm, rtmul) ;
#else /*HAVE_FMA_INSTRUCTIONS*/
      op1 = _mm256_mul_pd(rtC, rHm) ;
      rtmul = _mm256_add_pd(rtmul, op1) ;
#endif /*HAVE_FMA_INSTRUCTIONS*/
      for ( m = 1 ; m <= n ; m ++ ) {
	/*rotation coefficients for \pm\nu*/
	rHp = _mm256_broadcast_sd(&(H[wbfmm_rotation_index_numn( nu,m,n)])) ;
	rHm = _mm256_broadcast_sd(&(H[wbfmm_rotation_index_numn(-nu,m,n)])) ;

	offp = 2*cstri*wbfmm_coefficient_index_nm(n, m) ;
	offm = offp - 4*cstri*m ;

	wbfmm_cos_sin_recursion_avx(rCmch,rSmch,rCch,rSch) ;

	rtC = _mm256_set_pd(Ci[offm+1], Ci[offm+0], Ci[offp+1], Ci[offp+0]) ;
	rtS = _mm256_permute_pd(rtC, 5) ;
	rtS = _mm256_mul_pd(rtS, rsgn) ;
	rtC = _mm256_mul_pd(rtC, rCmch) ;

#ifdef HAVE_FMA_INSTRUCTIONS
	rtC    = _mm256_fmadd_pd(rtS, rSmch, rtC) ;
	rtmul  = _mm256_fmadd_pd(rtC, rHp, rtmul ) ;      
	rtC    = _mm256_permute2f128_pd(rtC, rtC, 1) ;
	rtmul  = _mm256_fmadd_pd(rtC, rHm, rtmul ) ;      
#else /*HAVE_FMA_INSTRUCTIONS*/
	rtS = _mm256_mul_pd(rtS, rSmch) ;

 	op1 = _mm256_mul_pd(rtC, rHp) ;
	rtmul  = _mm256_add_pd(rtmul,  op1) ;
	op1 = _mm256_mul_pd(rtS, rHp) ;
	rtmul  = _mm256_add_pd(rtmul,  op1) ;

	rtC = _mm256_permute2f128_pd(rtC, rtC, 1) ;
	rtS = _mm256_permute2f128_pd(rtS, rtS, 1) ;

	op1 = _mm256_mul_pd(rtC, rHm) ;
	rtmul  = _mm256_add_pd(rtmul,  op1) ;
	op1 = _mm256_mul_pd(rtS, rHm) ;
	rtmul  = _mm256_add_pd(rtmul,  op1) ;
#endif /*HAVE_FMA_INSTRUCTIONS*/
      }

      /*output indices for \pm\nu*/
      offp = 2*cstro*wbfmm_coefficient_index_nm(n, nu) ;
      offm = offp - 4*cstro*nu ;

      _mm256_store_pd(tmul, rtmul) ;

      Co[offp+0] = sc*Co[offp+0] + tmul[0]*Cnph + tmul[1]*Snph ;
      Co[offp+1] = sc*Co[offp+1] + tmul[1]*Cnph - tmul[0]*Snph ;
      Co[offm+0] = sc*Co[offm+0] + tmul[2]*Cnph - tmul[3]*Snph ;
      Co[offm+1] = sc*Co[offm+1] + tmul[3]*Cnph + tmul[2]*Snph ;
    }
  }
  
  return 0 ;
}

#else /*WBFMM_SINGLE_PRECISION*/
/*make the single-precision AVX function a call to the reference
  function*/

gint WBFMM_FUNCTION_NAME(wbfmm_rotate_H_avx)(WBFMM_REAL *Co, gint cstro, 
					     WBFMM_REAL *Ci, gint cstri,
					     gint N,
					     gint nq,
					     WBFMM_REAL *H,
					     WBFMM_REAL ph, WBFMM_REAL ch,
					     WBFMM_REAL sc)
{
  return WBFMM_FUNCTION_NAME(wbfmm_rotate_H_ref)(Co, cstro, Ci, cstri,
						 N, nq, H, ph, ch, sc) ;
}
#endif /*WBFMM_SINGLE_PRECISION*/
#else /*HAVE_AVX_INSTRUCTIONS*/
gint WBFMM_FUNCTION_NAME(wbfmm_rotate_H_avx)(WBFMM_REAL *Co, gint cstro, 
					     WBFMM_REAL *Ci, gint cstri,
					     gint N,
					     gint nq,
					     WBFMM_REAL *H,
					     WBFMM_REAL ph, WBFMM_REAL ch,
					     WBFMM_REAL sc)
{
  return WBFMM_FUNCTION_NAME(wbfmm_rotate_H_ref)(Co, cstro, Ci, cstri,
						 N, nq, H, ph, ch, sc) ;
}
#endif /*HAVE_AVX_INSTRUCTIONS*/

gint WBFMM_FUNCTION_NAME(wbfmm_rotate_H_ref)(WBFMM_REAL *Co, gint cstro, 
					     WBFMM_REAL *Ci, gint cstri,
					     gint N, gint nq,
					     WBFMM_REAL *H,
					     WBFMM_REAL ph, WBFMM_REAL ch,
					     WBFMM_REAL sc)

/*
  apply rotation (matrix H from wbfmm_coefficients_H_rotation) to
  rotate input coefficients Ci into output Co, through angles
  (th,ph,ch), G&D, section 6, and (2.27)
*/

{
  gint nu, n, m, offp, offm, j ;
  WBFMM_REAL Cmch, Smch, Cnph, Snph, Cch, Sch, Cph, Sph ;
  WBFMM_REAL Hp, Hm ;
  
  if ( cstro < nq )
    g_error("%s: output coefficient stride (%d) less than number of source "
	    " components (%d)", __FUNCTION__, cstro, nq) ;
  if ( cstri < nq )
    g_error("%s: input coefficient stride (%d) less than number of source "
	    " components (%d)", __FUNCTION__, cstri, nq) ;
  
  /*initialize recursions*/
  Cph = COS(ph) ; Sph = SIN(ph) ;
  Cch = COS(ch) ; Sch = SIN(ch) ;

  /* inside loops, trigonmetric quantities are calculated using
   * recursions and take the following values:
   *
   * Smch = SIN( m*ch) ; Cmch = COS( m*ch) ;
   * Cnph = COS(nu*ph) ; Snph = SIN(nu*ph) ;
   * 
   * offX (X = `p', `m') = offset into array, `p' for `plus' indices,
   * `m' for `minus' (+/-m)
   *
   * rotation coefficients for combinations of \nu and m are
   * 
   * Hp = H_{n}^{+\nu, +m} = H_{n}^{-\nu, -m}
   * Hm = H_{n}^{-\nu, +m} = H_{n}^{+\nu, -m}
   * 
   * using the symmetries in G&D
   */

  for ( n = 0 ; n <= N ; n ++ ) {
    WBFMM_REAL tC[4], tS[4] ; /*4*/
    WBFMM_REAL tmul[256] = {0.0} ; /*16*/

    nu = 0 ; Cnph = 1.0 ; Snph = 0.0 ;

    m = 0 ; Cmch = 1.0 ; Smch = 0.0 ;

    offp = 2*cstri*wbfmm_coefficient_index_nm(n,m) ;
    Hp = H[wbfmm_rotation_index_numn(nu,m,n)] ;
    
    for ( j = 0 ; j < nq ; j ++ ) {
      tC[0] = Cmch*Ci[offp+2*j+0] ; tC[1] = Cmch*Ci[offp+2*j+1] ;

      tmul[16*j+ 0] += tC[0]*Hp ; tmul[16*j+ 1] += tC[1]*Hp ;
    }
    
    for ( m = 1 ; m <= n ; m ++ ) {
      Hp = H[wbfmm_rotation_index_numn( nu,m,n)] ;
      Hm = H[wbfmm_rotation_index_numn(-nu,m,n)] ;

      offp = 2*cstri*wbfmm_coefficient_index_nm(n, m) ;
      offm = 2*cstri*wbfmm_coefficient_index_nm(n,-m) ;

      wbfmm_cos_sin_recursion(Cmch,Smch,Cch,Sch) ;

      for ( j = 0 ; j < nq ; j ++ ) {
	tC[0] = Cmch*Ci[offp+2*j+0] ; tC[1] = Cmch*Ci[offp+2*j+1] ;
	tC[2] = Cmch*Ci[offm+2*j+0] ; tC[3] = Cmch*Ci[offm+2*j+1] ;

	tS[0] = Smch*Ci[offp+2*j+0] ; tS[1] = Smch*Ci[offp+2*j+1] ;
	tS[2] = Smch*Ci[offm+2*j+0] ; tS[3] = Smch*Ci[offm+2*j+1] ;
	
	tmul[16*j+ 0] += tC[0]*Hp ; tmul[16*j+ 1] += tC[1]*Hp ;
	
	tmul[16*j+ 4] += tS[0]*Hp ; tmul[16*j+ 5] += tS[1]*Hp ;

	tmul[16*j+10] += tC[2]*Hm ; tmul[16*j+11] += tC[3]*Hm ;

	tmul[16*j+14] += tS[2]*Hm ; tmul[16*j+15] += tS[3]*Hm ;
      }
    }

    offp = 2*cstro*wbfmm_coefficient_index_nm(n, nu) ;
    for ( j = 0 ; j < nq ; j ++ ) {
      Co[offp+2*j+0] = sc*Co[offp+2*j+0] +
	(tmul[16*j+ 0] - tmul[16*j+ 5] + tmul[16*j+10] + tmul[16*j+15])*Cnph ;
      Co[offp+2*j+1] = sc*Co[offp+2*j+1] +
	(tmul[16*j+ 1] + tmul[16*j+ 4] + tmul[16*j+11] - tmul[16*j+14])*Cnph ;
    }
    
    for ( nu = 1 ; nu <= n ; nu ++ ) {
      wbfmm_cos_sin_recursion(Cnph,Snph,Cph,Sph) ;

      m = 0 ; Cmch = 1.0 ; Smch = 0.0 ;

      offp = 2*cstri*wbfmm_coefficient_index_nm(n,m) ;

      Hp = H[wbfmm_rotation_index_numn( nu,m,n)] ;
      Hm = H[wbfmm_rotation_index_numn(-nu,m,n)] ;
      for ( j = 0 ; j < nq ; j ++ ) {
	tC[0] = Cmch*Ci[offp+2*j+0] ; tC[1] = Cmch*Ci[offp+2*j+1] ;
	tC[2] = Smch*Ci[offp+2*j+0] ; tC[3] = Smch*Ci[offp+2*j+1] ;

	tmul[16*j+ 2] = tmul[16*j+ 3] = 0.0 ;
	tmul[16*j+ 6] = tmul[16*j+ 7] = 0.0 ;
	tmul[16*j+10] = tmul[16*j+11] = 0.0 ;
	tmul[16*j+14] = tmul[16*j+15] = 0.0 ;

	tmul[16*j+ 0] = tC[0]*Hp ; tmul[16*j+ 1] = tC[1]*Hp ;
	tmul[16*j+ 4] = tC[2]*Hp ; tmul[16*j+ 5] = tC[3]*Hp ;

	tmul[16*j+ 8] = tC[0]*Hm ; tmul[16*j+ 9] = tC[1]*Hm ;
	tmul[16*j+12] = tC[2]*Hm ; tmul[16*j+13] = tC[3]*Hm ;
      }
      
      for ( m = 1 ; m <= n ; m ++ ) {
	/*rotation coefficients for \pm\nu*/
	Hp = H[wbfmm_rotation_index_numn( nu,m,n)] ;
	Hm = H[wbfmm_rotation_index_numn(-nu,m,n)] ;

	offp = 2*cstri*wbfmm_coefficient_index_nm(n, m) ;
	offm = 2*cstri*wbfmm_coefficient_index_nm(n,-m) ;

	wbfmm_cos_sin_recursion(Cmch,Smch,Cch,Sch) ;

	for ( j = 0 ; j < nq ; j ++ ) {
	  tC[0] = Cmch*Ci[offp+2*j+0] ; tC[1] = Cmch*Ci[offp+2*j+1] ;
	  tC[2] = Cmch*Ci[offm+2*j+0] ; tC[3] = Cmch*Ci[offm+2*j+1] ;

	  tS[0] = Smch*Ci[offp+2*j+0] ; tS[1] = Smch*Ci[offp+2*j+1] ;
	  tS[2] = Smch*Ci[offm+2*j+0] ; tS[3] = Smch*Ci[offm+2*j+1] ;

	  tmul[16*j+ 0] += tC[0]*Hp ; tmul[16*j+ 1] += tC[1]*Hp ;
	  tmul[16*j+ 2] += tC[2]*Hp ; tmul[16*j+ 3] += tC[3]*Hp ;

	  tmul[16*j+ 4] += tS[0]*Hp ; tmul[16*j+ 5] += tS[1]*Hp ;
	  tmul[16*j+ 6] += tS[2]*Hp ; tmul[16*j+ 7] += tS[3]*Hp ;

	  tmul[16*j+ 8] += tC[0]*Hm ; tmul[16*j+ 9] += tC[1]*Hm ;
	  tmul[16*j+10] += tC[2]*Hm ; tmul[16*j+11] += tC[3]*Hm ;

	  tmul[16*j+12] += tS[0]*Hm ; tmul[16*j+13] += tS[1]*Hm ;
	  tmul[16*j+14] += tS[2]*Hm ; tmul[16*j+15] += tS[3]*Hm ;
	}
      }
      /*output indices for \pm\nu*/
      offp = 2*cstro*wbfmm_coefficient_index_nm(n, nu) ;
      offm = offp - 4*cstro*nu ;

      for ( j = 0 ; j < nq ; j ++ ) {
	Co[offp+2*j+0] = sc*Co[offp+2*j+0] +
	  (tmul[16*j+ 0] - tmul[16*j+ 5] + tmul[16*j+10] + tmul[16*j+15])*Cnph +
	  (tmul[16*j+ 1] + tmul[16*j+ 4] + tmul[16*j+11] - tmul[16*j+14])*Snph ;
	Co[offp+2*j+1] = sc*Co[offp+2*j+1] +
	  (tmul[16*j+ 1] + tmul[16*j+ 4] + tmul[16*j+11] - tmul[16*j+14])*Cnph -
	  (tmul[16*j+ 0] - tmul[16*j+ 5] + tmul[16*j+10] + tmul[16*j+15])*Snph ;
	Co[offm+2*j+0] = sc*Co[offm+2*j+0] +
	  (tmul[16*j+ 8] - tmul[16*j+13] + tmul[16*j+ 2] + tmul[16*j+ 7])*Cnph -
	  (tmul[16*j+12] + tmul[16*j+ 9] + tmul[16*j+ 3] - tmul[16*j+ 6])*Snph ;
	Co[offm+2*j+1] = sc*Co[offm+2*j+1] +
	  (tmul[16*j+12] + tmul[16*j+ 9] + tmul[16*j+ 3] - tmul[16*j+ 6])*Cnph +
	  (tmul[16*j+ 8] - tmul[16*j+13] + tmul[16*j+ 2] + tmul[16*j+ 7])*Snph ;
      }
    }
  }

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_coefficients_H_to_T)(WBFMM_REAL *H, gint N,
						    WBFMM_REAL th, 
						    WBFMM_REAL ph, 
						    WBFMM_REAL ch,
						    WBFMM_REAL *T)

/*
 * minimum size of T is 2*wbfmm_T_rotation_matrix_size(N)
 */
  
{
  gint n, m, nu, sizet, offm, offp, offn ;
  WBFMM_REAL Hp, Hm, *Tn ;

  for ( n = 0 ; n <= N ; n ++ ) {
    /*start of T matrix for coefficients of order n is end of
      coefficients for n-1*/
    offn = 2*wbfmm_T_rotation_matrix_size(n-1) ;
    Tn = &(T[offn]) ;
    /*size of matrix*/
    sizet = 2*n + 1 ;

    m = 0 ;
    offp = 2*(n+m)*sizet ;
    for ( nu = -n ; nu <= n ; nu ++ ) {
      /* off = wbfmm_rotation_index_numn(nu, m, n) ; */
      Hp = H[wbfmm_rotation_index_numn(nu, m, n)] ;
      Tn[offp + 2*(nu+n)+0] =
	Hp*(cos(nu*ch)*cos(m*ph) + sin(nu*ch)*sin(m*ph)) ;
      Tn[offp + 2*(nu+n)+1] =
	Hp*(sin(nu*ch)*cos(m*ph) - cos(nu*ch)*sin(m*ph)) ;
    }
    
    for ( m = 1 ; m <= n ; m ++ ) {
      offp = 2*(n+m)*sizet ; offm = 2*(n-m)*sizet ;
      for ( nu = -n ; nu <= n ; nu ++ ) {
	/* off = wbfmm_rotation_index_numn(nu, m, n) ; */
	Hp = H[wbfmm_rotation_index_numn(nu, m, n)] ;
	Tn[offp + 2*(nu+n)+0] =
	  Hp*(cos(nu*ch)*cos(m*ph) + sin(nu*ch)*sin(m*ph)) ;
	Tn[offp + 2*(nu+n)+1] =
	  Hp*(sin(nu*ch)*cos(m*ph) - cos(nu*ch)*sin(m*ph)) ;

	/* off = wbfmm_rotation_index_numn(-nu, m, n) ; */
	Hm = H[wbfmm_rotation_index_numn(-nu, m, n)] ;
	Tn[offm + 2*(nu+n)+0] =
	  Hm*(cos(nu*ch)*cos(m*ph) + sin(nu*ch)*sin(-m*ph)) ;
	Tn[offm + 2*(nu+n)+1] =
	  Hm*(sin(nu*ch)*cos(m*ph) - cos(nu*ch)*sin(-m*ph)) ;
      }
    }
  }
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_rotate_T)(WBFMM_REAL *Co, gint cstro, 
					 WBFMM_REAL *Ci, gint cstri,
					 gint N, gint nq,
					 WBFMM_REAL *T,
					 WBFMM_REAL *sc)


{
  gint n, lda, ldb, ldc, sizet, offn, offi, offo ;
  WBFMM_REAL *Tn, al[] = {1, 0} ;
  
  lda = cstri ; ldc = cstro ;
  for ( n = 0 ; n <= N ; n ++ ) {
    /*select matrix for rotation of order n coefficients*/
    offn = 2*wbfmm_T_rotation_matrix_size(n-1) ;
    Tn = &(T[offn]) ;
    /*size of matrix and leading dimension*/
    ldb = sizet = 2*n + 1 ;
    /*offsets to input and output coefficients*/
    offi = 2*cstri*(wbfmm_coefficient_index_nm(n,-n)) ;
    offo = 2*cstro*(wbfmm_coefficient_index_nm(n,-n)) ;

#ifdef WBFMM_SINGLE_PRECISION
    cgemm_("N", "N", &nq, &sizet, &sizet, al, &(Ci[offi]), &lda,
	   Tn, &ldb, sc, &(Co[offo]), &ldc) ;
#else /*WBFMM_SINGLE_PRECISION*/
    zgemm_("N", "N", &nq, &sizet, &sizet, al, &(Ci[offi]), &lda,
	   Tn, &ldb, sc, &(Co[offo]), &ldc) ;
#endif /*WBFMM_SINGLE_PRECISION*/
  }

  return 0 ;
}
/* @} */
