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
#include "wbfmm-avx.h"

#ifdef HAVE_AVX_INSTRUCTIONS
#include <immintrin.h>
#endif /*HAVE_AVX_INSTRUCTIONS*/

static gint _wbfmm_rotate_H_laplace_ref_nq(WBFMM_REAL *Co, gint cstro,
					   WBFMM_REAL *Ci, gint cstri,
					   gint N, gint nq,
					   WBFMM_REAL *H,
					   WBFMM_REAL ph, WBFMM_REAL ch,
					   WBFMM_REAL sc)

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
      for ( i = 0 ; i < nq ; i ++ ) {
	tr[i] += 2.0*(Ci[cstri*(idxi+0)+i]*Cmch -
		      Ci[cstri*(idxi+1)+i]*Smch)*Hp ;
      }
    }

    for ( i = 0 ; i < nq ; i ++ ) {
      Co[cstro*idxo+i] = tr[i] + sc*Co[cstro*idxo+i] ;
    }

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
	Co[cstro*(idxo+0)+i] =
	  Cnph*tr[i] + Snph*ti[i] + sc*Co[cstro*(idxo+0)+i] ;
	Co[cstro*(idxo+1)+i] =
	  Cnph*ti[i] - Snph*tr[i] + sc*Co[cstro*(idxo+1)+i] ;
      }
    }
  }
  
  return 0 ;
}

static gint _wbfmm_rotate_H_laplace_ref_3(WBFMM_REAL *Co, gint cstro,
					   WBFMM_REAL *Ci, gint cstri,
					   gint N, gint NQ,
					   WBFMM_REAL *H,
					   WBFMM_REAL ph, WBFMM_REAL ch,
					   WBFMM_REAL sc)

{
  gint n, m, nu, idxi, idxo, i, nq ;
  WBFMM_REAL Cmch, Smch, Cnph, Snph, Cch, Sch, Cph, Sph, tr[32] ;
  WBFMM_REAL Hp, Hm ;
  
  nq = 3 ;
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
      for ( i = 0 ; i < nq ; i ++ ) {
	tr[i] += 2.0*(Ci[cstri*(idxi+0)+i]*Cmch -
		      Ci[cstri*(idxi+1)+i]*Smch)*Hp ;
      }
    }

    for ( i = 0 ; i < nq ; i ++ ) {
      Co[cstro*idxo+i] = tr[i] + sc*Co[cstro*idxo+i] ;
    }

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
	Co[cstro*(idxo+0)+i] =
	  Cnph*tr[i] + Snph*ti[i] + sc*Co[cstro*(idxo+0)+i] ;
	Co[cstro*(idxo+1)+i] =
	  Cnph*ti[i] - Snph*tr[i] + sc*Co[cstro*(idxo+1)+i] ;
      }
    }
  }
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_laplace_rotate_H_ref)(WBFMM_REAL *Co, gint cstro,
						     WBFMM_REAL *Ci, gint cstri,
						     gint N, gint nq,
						     WBFMM_REAL *H,
						     WBFMM_REAL ph,
						     WBFMM_REAL ch,
						     WBFMM_REAL sc)

{
  if ( nq == 3 ) 
    return
      _wbfmm_rotate_H_laplace_ref_3(Co, cstro, Ci, cstri, N, nq, H, ph, ch, sc) ;
  
  return
    _wbfmm_rotate_H_laplace_ref_nq(Co, cstro, Ci, cstri, N, nq, H, ph, ch, sc) ;
  
  return 0 ;
}

#ifdef HAVE_AVX_INSTRUCTIONS
#ifndef WBFMM_SINGLE_PRECISION
/*only compile AVX for double precision functions, for now anyway*/

static gint wbfmm_laplace_rotate_H3_avx(WBFMM_REAL *Co, gint cstro,
					WBFMM_REAL *Ci, gint cstri,
					gint N, gint nq,
					WBFMM_REAL *H,
					WBFMM_REAL ph,
					WBFMM_REAL ch,
					WBFMM_REAL sc)

{
  gint n, m, nu, idxi, idxo ;
  WBFMM_REAL Cch, Sch, Cph, Sph ;
  __attribute__ ((aligned (32))) WBFMM_REAL tr[4] ;
  __attribute__ ((aligned (32))) WBFMM_REAL ti[4] ;
  WBFMM_REAL Hp, Hm ;
  __m256d rCch, rSch, rCmch, rSmch, rHm, rHp, rtr, rti, rC, op1 ;
  __m256d rCph, rSph, rCnph, rSnph ; 
  
  g_assert(nq <= cstri) ;
  g_assert(nq <= cstro) ;
  g_assert(nq == 3) ;

  /*initialize recursions*/
  Cph = COS(ph) ; Sph = SIN(ph) ;
  Cch = COS(ch) ; Sch = SIN(ch) ;

  rCch = _mm256_set1_pd(Cch) ; rSch = _mm256_set1_pd(Sch) ;
  rCph = _mm256_set1_pd(Cph) ; rSph = _mm256_set1_pd(Sph) ;

  for ( n = 0 ; n <= N ; n ++ ) {
    nu = 0 ; idxo = n*n ;
    m  = 0 ; idxi = n*n ;

    Hp = H[wbfmm_rotation_index_numn( nu,m,n)] ;
    rHp = _mm256_set1_pd(Hp) ;

    rC    = _mm256_loadu_pd(&(Ci[cstri*idxi])) ;
    rtr   = _mm256_mul_pd(rHp, rC) ;
    rCmch = _mm256_set1_pd(1.0) ; rSmch = _mm256_setzero_pd() ;      

    for ( m = 1 ; m <= n ; m ++ ) {
      idxi = wbfmm_index_laplace_nm(n,m) ;

      wbfmm_cos_sin_recursion_avx(rCmch,rSmch,rCch,rSch) ;
      Hp = H[wbfmm_rotation_index_numn( nu,m,n)] ;
      rHp = _mm256_set1_pd(Hp) ;

      rC  = _mm256_loadu_pd(&(Ci[cstri*(idxi+0)])) ;
      op1 = _mm256_mul_pd(rCmch, rC) ;
      rC  = _mm256_loadu_pd(&(Ci[cstri*(idxi+1)])) ;
      rC  = _mm256_mul_pd(rSmch, rC) ;
      op1 = _mm256_sub_pd(op1, rC) ; /*\cos m\chi C_r - \sin m\chi C_i*/
      op1 = _mm256_mul_pd(op1, rHp) ;
      rtr = _mm256_add_pd(rtr, op1) ;
      rtr = _mm256_add_pd(rtr, op1) ;
    }

    _mm256_store_pd(tr, rtr) ;
    Co[cstro*idxo+0] = tr[0] + sc*Co[cstro*idxo+0] ;
    Co[cstro*idxo+1] = tr[1] + sc*Co[cstro*idxo+1] ;
    Co[cstro*idxo+2] = tr[2] + sc*Co[cstro*idxo+2] ;

    rCnph = _mm256_set1_pd(1.0) ; rSnph = _mm256_setzero_pd() ;      
    for ( nu = 1 ; nu <= n ; nu ++ ) {
      wbfmm_cos_sin_recursion_avx(rCnph,rSnph,rCph,rSph) ;

      idxo = wbfmm_index_laplace_nm(n,nu) ;
      m = 0 ; idxi = n*n ;

      Hp  = H[wbfmm_rotation_index_numn( nu,m,n)] ;
      rHp = _mm256_set1_pd(Hp) ;
      rC  = _mm256_loadu_pd(&(Ci[cstri*idxi])) ;
      rtr = _mm256_mul_pd(rHp, rC) ;

      rti = _mm256_setzero_pd() ;
      
      rCmch = _mm256_set1_pd(1.0) ; rSmch = _mm256_setzero_pd() ;      
      for ( m = 1 ; m <= n ; m ++ ) {
	idxi = wbfmm_index_laplace_nm(n,m) ;
	
	wbfmm_cos_sin_recursion_avx(rCmch,rSmch,rCch,rSch) ;
	Hp  = H[wbfmm_rotation_index_numn( nu,m,n)] ;
	Hm  = H[wbfmm_rotation_index_numn(-nu,m,n)] ;
	rHp = _mm256_set1_pd(Hp+Hm) ;
	rHm = _mm256_set1_pd(Hp-Hm) ;

	rC  = _mm256_loadu_pd(&(Ci[cstri*(idxi+1)])) ;
	op1 = _mm256_mul_pd(rCmch, rC) ;
#ifdef HAVE_FMA_INSTRUCTIONS
	rti = _mm256_fmadd_pd(op1, rHm, rti) ;
#else /*HAVE_FMA_INSTRUCTIONS*/
	op1 = _mm256_mul_pd(op1, rHm) ;
	rti = _mm256_add_pd(rti, op1) ;
#endif /*HAVE_FMA_INSTRUCTIONS*/
	
	op1 = _mm256_mul_pd(rSmch, rC) ;
#ifdef HAVE_FMA_INSTRUCTIONS
	rtr = _mm256_fmsub_pd(op1, rHp, rtr) ;
#else /*HAVE_FMA_INSTRUCTIONS*/	
	op1 = _mm256_mul_pd(op1, rHp) ;
	rtr = _mm256_sub_pd(rtr, op1) ;
#endif /*HAVE_FMA_INSTRUCTIONS*/

	rC  = _mm256_loadu_pd(&(Ci[cstri*(idxi+0)])) ;
	op1 = _mm256_mul_pd(rCmch, rC) ;
#ifdef HAVE_FMA_INSTRUCTIONS
	rtr = _mm256_fmsub_pd(op1, rHp, rtr) ;
#else /*HAVE_FMA_INSTRUCTIONS*/	
	op1 = _mm256_mul_pd(op1, rHp) ;
	rtr = _mm256_add_pd(rtr, op1) ;
#endif /*HAVE_FMA_INSTRUCTIONS*/

	op1 = _mm256_mul_pd(rSmch, rC) ;
#ifdef HAVE_FMA_INSTRUCTIONS
	rti = _mm256_fmadd_pd(op1, rHm, rti) ;
#else /*HAVE_FMA_INSTRUCTIONS*/
	op1 = _mm256_mul_pd(op1, rHm) ;
	rti = _mm256_add_pd(rti, op1) ;
#endif /*HAVE_FMA_INSTRUCTIONS*/
      }

#ifdef HAVE_FMA_INSTRUCTIONS
      op1 = _mm256_mul_pd(rtr,rCnph) ;
      op1 = _mm256_fmadd_pd(rti,rSnph,op1) ;
      _mm256_store_pd(tr, op1) ;
      op1 = _mm256_mul_pd(rtr,rSnph) ;
      op1 = _mm256_fmsub_pd(rti,rCnph,op1) ;
      _mm256_store_pd(ti, op1) ;
#else /*HAVE_FMA_INSTRUCTIONS*/
      _m256d op2 ;
      op1 = _mm256_mul_pd(rtr,rCnph) ; op2 = _mm256_mul_pd(rti,rSnph) ;
      op1 = _mm256_add_pd(op1,op2) ;
      _mm256_store_pd(tr, op1) ;
      op1 = _mm256_mul_pd(rti,rCnph) ; op2 = _mm256_mul_pd(rtr,rSnph) ;
      op1 = _mm256_sub_pd(op1,op2) ;
      _mm256_store_pd(ti, op1) ;
#endif /*HAVE_FMA_INSTRUCTIONS*/
      Co[cstro*(idxo+0)+0] = sc*Co[cstro*(idxo+0)+0] + tr[0] ;
      Co[cstro*(idxo+1)+0] = sc*Co[cstro*(idxo+1)+0] + ti[0] ;
      Co[cstro*(idxo+0)+1] = sc*Co[cstro*(idxo+0)+1] + tr[1] ;
      Co[cstro*(idxo+1)+1] = sc*Co[cstro*(idxo+1)+1] + ti[1] ;
      Co[cstro*(idxo+0)+2] = sc*Co[cstro*(idxo+0)+2] + tr[2] ;
      Co[cstro*(idxo+1)+2] = sc*Co[cstro*(idxo+1)+2] + ti[2] ;
    }
  }
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_laplace_rotate_H_avx)(WBFMM_REAL *Co, gint cstro,
						     WBFMM_REAL *Ci, gint cstri,
						     gint N, gint nq,
						     WBFMM_REAL *H,
						     WBFMM_REAL ph,
						     WBFMM_REAL ch,
						     WBFMM_REAL sc)

{
  gint n, m, nu, idxi, idxo, i ;
  WBFMM_REAL Cch, Sch, Cph, Sph ;
  __attribute__ ((aligned (32))) WBFMM_REAL tr[4] ;
  __attribute__ ((aligned (32))) WBFMM_REAL ti[4] ;
  WBFMM_REAL Hp, Hm ;
  __m256d rCch, rSch, rCmch, rSmch, rHm, rHp, rtr, rti, rC, op1 ;
  __m256d rCph, rSph, rCnph, rSnph ;

  if ( nq == 3 ) return wbfmm_laplace_rotate_H3_avx(Co, cstro, Ci, cstri,
						    N, nq, H, ph, ch, sc) ;
  
  g_assert(nq <= cstri) ;
  g_assert(nq <= cstro) ;
  g_assert(nq <= 4) ;

  /*initialize recursions*/
  Cph = COS(ph) ; Sph = SIN(ph) ;
  Cch = COS(ch) ; Sch = SIN(ch) ;

  rCch = _mm256_set1_pd(Cch) ; rSch = _mm256_set1_pd(Sch) ;
  rCph = _mm256_set1_pd(Cph) ; rSph = _mm256_set1_pd(Sph) ;

  for ( n = 0 ; n <= N ; n ++ ) {
    nu = 0 ; idxo = n*n ;
    m  = 0 ; idxi = n*n ;

    Hp = H[wbfmm_rotation_index_numn( nu,m,n)] ;
    rHp = _mm256_set1_pd(Hp) ;

    rC    = _mm256_loadu_pd(&(Ci[cstri*idxi])) ;
    rtr   = _mm256_mul_pd(rHp, rC) ;
    rCmch = _mm256_set1_pd(1.0) ; rSmch = _mm256_setzero_pd() ;      

    for ( m = 1 ; m <= n ; m ++ ) {
      idxi = wbfmm_index_laplace_nm(n,m) ;

      wbfmm_cos_sin_recursion_avx(rCmch,rSmch,rCch,rSch) ;
      Hp = H[wbfmm_rotation_index_numn( nu,m,n)] ;
      rHp = _mm256_set1_pd(Hp) ;

      rC  = _mm256_loadu_pd(&(Ci[cstri*(idxi+0)])) ;
      op1 = _mm256_mul_pd(rCmch, rC) ;
      rC  = _mm256_loadu_pd(&(Ci[cstri*(idxi+1)])) ;
      rC  = _mm256_mul_pd(rSmch, rC) ;
      op1 = _mm256_sub_pd(op1, rC) ; /*\cos m\chi C_r - \sin m\chi C_i*/
      op1 = _mm256_mul_pd(op1, rHp) ;
      rtr = _mm256_add_pd(rtr, op1) ;
      rtr = _mm256_add_pd(rtr, op1) ;
    }

    _mm256_store_pd(tr, rtr) ;
    for ( i = 0 ; i < nq ; i ++ ) {
      Co[cstro*idxo+i] = tr[i] + sc*Co[cstro*idxo+i] ;
    }

    rCnph = _mm256_set1_pd(1.0) ; rSnph = _mm256_setzero_pd() ;      
    for ( nu = 1 ; nu <= n ; nu ++ ) {
      wbfmm_cos_sin_recursion_avx(rCnph,rSnph,rCph,rSph) ;

      idxo = wbfmm_index_laplace_nm(n,nu) ;
      m = 0 ; idxi = n*n ;

      Hp  = H[wbfmm_rotation_index_numn( nu,m,n)] ;
      rHp = _mm256_set1_pd(Hp) ;
      rC  = _mm256_loadu_pd(&(Ci[cstri*idxi])) ;
      rtr = _mm256_mul_pd(rHp, rC) ;

      rti = _mm256_setzero_pd() ;
      
      rCmch = _mm256_set1_pd(1.0) ; rSmch = _mm256_setzero_pd() ;      
      for ( m = 1 ; m <= n ; m ++ ) {
	idxi = wbfmm_index_laplace_nm(n,m) ;
	
	wbfmm_cos_sin_recursion_avx(rCmch,rSmch,rCch,rSch) ;
	Hp  = H[wbfmm_rotation_index_numn( nu,m,n)] ;
	Hm  = H[wbfmm_rotation_index_numn(-nu,m,n)] ;
	rHp = _mm256_set1_pd(Hp+Hm) ;
	rHm = _mm256_set1_pd(Hp-Hm) ;

	rC  = _mm256_loadu_pd(&(Ci[cstri*(idxi+1)])) ;
	op1 = _mm256_mul_pd(rCmch, rC) ;
#ifdef HAVE_FMA_INSTRUCTIONS
	rti = _mm256_fmadd_pd(op1, rHm, rti) ;
	op1 = _mm256_mul_pd(rSmch, rC) ;
	/*NOTE the fused SUBTRACTION here, which is balanced by a
	  subtraction further down (not the same in the unfused version)*/
	rtr = _mm256_fmsub_pd(op1, rHp, rtr) ;
#else /*HAVE_FMA_INSTRUCTIONS*/
	op1 = _mm256_mul_pd(op1, rHm) ;
	rti = _mm256_add_pd(rti, op1) ;
	op1 = _mm256_mul_pd(rSmch, rC) ;
	op1 = _mm256_mul_pd(op1, rHp) ;
	rtr = _mm256_sub_pd(rtr, op1) ;
#endif /*HAVE_FMA_INSTRUCTIONS*/

	rC  = _mm256_loadu_pd(&(Ci[cstri*(idxi+0)])) ;
	op1 = _mm256_mul_pd(rCmch, rC) ;
#ifdef HAVE_FMA_INSTRUCTIONS
	/*see NOTE above*/
	rtr = _mm256_fmsub_pd(op1, rHp, rtr) ;
	op1 = _mm256_mul_pd(rSmch, rC) ;
	rti = _mm256_fmadd_pd(op1, rHm, rti) ;
#else /*HAVE_FMA_INSTRUCTIONS*/
	op1 = _mm256_mul_pd(op1, rHp) ;
	rtr = _mm256_add_pd(rtr, op1) ;
	op1 = _mm256_mul_pd(rSmch, rC) ;
	op1 = _mm256_mul_pd(op1, rHm) ;
	rti = _mm256_add_pd(rti, op1) ;
#endif /*HAVE_FMA_INSTRUCTIONS*/
      }

#ifdef HAVE_FMA_INSTRUCTIONS
      op1 = _mm256_mul_pd(rtr,rCnph) ;
      op1 = _mm256_fmadd_pd(rti,rSnph,op1) ;
      _mm256_store_pd(tr, op1) ;
      op1 = _mm256_mul_pd(rtr,rSnph) ;
      op1 = _mm256_fmsub_pd(rti,rCnph,op1) ;
      _mm256_store_pd(ti, op1) ;
#else /*HAVE_FMA_INSTRUCTIONS*/
      _m256d op2 ;
      op1 = _mm256_mul_pd(rtr,rCnph) ; op2 = _mm256_mul_pd(rti,rSnph) ;
      op1 = _mm256_add_pd(op1,op2) ;
      _mm256_store_pd(tr, op1) ;
      op1 = _mm256_mul_pd(rti,rCnph) ; op2 = _mm256_mul_pd(rtr,rSnph) ;
      op1 = _mm256_sub_pd(op1,op2) ;
      _mm256_store_pd(ti, op1) ;
#endif /*HAVE_FMA_INSTRUCTIONS*/
      for ( i = 0 ; i < nq ; i ++ ) {
      	Co[cstro*(idxo+0)+i] = sc*Co[cstro*(idxo+0)+i] + tr[i] ;
      	Co[cstro*(idxo+1)+i] = sc*Co[cstro*(idxo+1)+i] + ti[i] ;
      }
    }
  }
  
  return 0 ;
}

#endif
#endif
