/* This file is part of WBFMM, a Wide-Band Fast Multipole Method code
 *
 * Copyright (C) 2019, 2024 Michael Carley
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

#include <blaswrap.h>

#include <wbfmm.h>

#include "wbfmm-private.h"
#include "wbfmm-derivatives.h"

#ifdef HAVE_AVX_INSTRUCTIONS
#include <immintrin.h>
#endif /*HAVE_AVX_INSTRUCTIONS*/

static void box_curl_evaluate(wbfmm_tree_t *t,
			      gint i0, gint i1,
			      WBFMM_REAL *src, gint sstr,
			      WBFMM_REAL *x, WBFMM_REAL *f)

{
  gint idx, j ;
  WBFMM_REAL *xs, r[3], R ;

  for ( j = i0 ; j < i1 ; j ++ ) {
    idx = t->ip[j] ;
    xs = wbfmm_tree_point_index(t, idx) ;
    r[0] = x[0] - xs[0] ; r[1] = x[1] - xs[1] ; r[2] = x[2] - xs[2] ;
    R = r[0]*r[0] + r[1]*r[1] + r[2]*r[2] ;
    if ( R > WBFMM_LOCAL_CUTOFF_RADIUS*WBFMM_LOCAL_CUTOFF_RADIUS ) {
      R *= SQRT(R)*4.0*M_PI ;
      r[0] /= R ; r[1] /= R ; r[2] /= R ; 
      f[0] -= src[idx*sstr+2]*r[1] - src[idx*sstr+1]*r[2] ;
      f[1] -= src[idx*sstr+0]*r[2] - src[idx*sstr+2]*r[0] ;
      f[2] -= src[idx*sstr+1]*r[0] - src[idx*sstr+0]*r[1] ;
    }
  }

  return ;
}

static void gradient_evaluate4(WBFMM_REAL r[12])

/*
 * vectorized evaluation of gradient elements
 * r: on entry contains displacement vectors:
 *
 * |x-x1 x-x2 x-x3 x-x4|
 * |y-y1 y-y2 y-y3 y-y4|
 * |z-z1 z-z2 z-z3 z-z4|
 *
 * and on exit contains
 * |(x-x1)/R_1^3 ...|
 * |(y-y1)/R_1^3 ...|
 * |(z-z1)/R_1^3 ...|
 *
 * equal to \nabla 1/R, R = |r|
 */
  
{
#ifndef WBFMM_SINGLE_PRECISION
#ifdef HAVE_AVX_INSTRUCTIONS
  __m256d rrx, rry, rrz, rR, op1 ;

  rrx = _mm256_loadu_pd(&r[0]) ;
  rry = _mm256_loadu_pd(&r[4]) ;
  rrz = _mm256_loadu_pd(&r[8]) ;

  rR = _mm256_mul_pd(rrx, rrx) ;
  /*this could be done with FMA when I get to a machine that has it*/
  op1 = _mm256_mul_pd(rry, rry) ;
  rR = _mm256_add_pd(rR, op1) ;
  op1 = _mm256_mul_pd(rrz, rrz) ;
  rR = _mm256_add_pd(rR, op1) ;
  /*invsqrt and rsqrt seem not to be available in gcc intrinsics ...*/
  op1 = _mm256_sqrt_pd(rR) ;
  rR = _mm256_mul_pd(rR, op1) ;

  rrx = _mm256_div_pd(rrx, rR) ;
  rry = _mm256_div_pd(rry, rR) ;
  rrz = _mm256_div_pd(rrz, rR) ;
  _mm256_storeu_pd(&(r[0]), rrx) ;  
  _mm256_storeu_pd(&(r[4]), rry) ;  
  _mm256_storeu_pd(&(r[8]), rrz) ;  

#else /*HAVE_AVX_INSTRUCTIONS*/
  WBFMM_REAL R[4] ;
  R[0] = r[4*0+0]*r[4*0+0] + r[4*1+0]*r[4*1+0] + r[4*2+0]*r[4*2+0] ;
  R[1] = r[4*0+1]*r[4*0+1] + r[4*1+1]*r[4*1+1] + r[4*2+1]*r[4*2+1] ;
  R[2] = r[4*0+2]*r[4*0+2] + r[4*1+2]*r[4*1+2] + r[4*2+2]*r[4*2+2] ;
  R[3] = r[4*0+3]*r[4*0+3] + r[4*1+3]*r[4*1+3] + r[4*2+3]*r[4*2+3] ;

  R[0] = 1.0/(R[0]*SQRT(R[0])) ;
  R[1] = 1.0/(R[1]*SQRT(R[1])) ;
  R[2] = 1.0/(R[2]*SQRT(R[2])) ;
  R[3] = 1.0/(R[3]*SQRT(R[3])) ;

  r[4*0+0] *= R[0] ; r[4*1+0] *= R[0] ; r[4*2+0] *= R[0] ; 
  r[4*0+1] *= R[1] ; r[4*1+1] *= R[1] ; r[4*2+1] *= R[1] ; 
  r[4*0+2] *= R[2] ; r[4*1+2] *= R[2] ; r[4*2+2] *= R[2] ; 
  r[4*0+3] *= R[3] ; r[4*1+3] *= R[3] ; r[4*2+3] *= R[3] ; 
#endif /*HAVE_AVX_INSTRUCTIONS*/
#else /*WBFMM_SINGLE_PRECISION*/
  WBFMM_REAL R[4] ;

  R[0] = r[4*0+0]*r[4*0+0] + r[4*1+0]*r[4*1+0] + r[4*2+0]*r[4*2+0] ;
  R[1] = r[4*0+1]*r[4*0+1] + r[4*1+1]*r[4*1+1] + r[4*2+1]*r[4*2+1] ;
  R[2] = r[4*0+2]*r[4*0+2] + r[4*1+2]*r[4*1+2] + r[4*2+2]*r[4*2+2] ;
  R[3] = r[4*0+3]*r[4*0+3] + r[4*1+3]*r[4*1+3] + r[4*2+3]*r[4*2+3] ;

  R[0] = 1.0/(R[0]*SQRT(R[0])) ;
  R[1] = 1.0/(R[1]*SQRT(R[1])) ;
  R[2] = 1.0/(R[2]*SQRT(R[2])) ;
  R[3] = 1.0/(R[3]*SQRT(R[3])) ;

  r[4*0+0] *= R[0] ; r[4*1+0] *= R[0] ; r[4*2+0] *= R[0] ; 
  r[4*0+1] *= R[1] ; r[4*1+1] *= R[1] ; r[4*2+1] *= R[1] ; 
  r[4*0+2] *= R[2] ; r[4*1+2] *= R[2] ; r[4*2+2] *= R[2] ; 
  r[4*0+3] *= R[3] ; r[4*1+3] *= R[3] ; r[4*2+3] *= R[3] ; 
#endif /*WBFMM_SINGLE_PRECISION*/

  return ;
}

static void box_curl_gradient_evaluate(wbfmm_tree_t *t,
				       gint i0, gint i1,
				       WBFMM_REAL *src, gint sstr,
				       WBFMM_REAL *x, WBFMM_REAL *f)

{
  gint idx, j ;
  WBFMM_REAL *xs, dr[3], r, r3, r5, nR[12], *s, df[3] ;

  for ( j = i0 ; j < i1 ; j ++ ) {
    idx = t->ip[j] ;
    xs = wbfmm_tree_point_index(t, idx) ;
    wbfmm_vector_diff(dr,x,xs) ;
    r = wbfmm_vector_length(dr) ;
    if ( r > WBFMM_LOCAL_CUTOFF_RADIUS ) {
      s = &(src[idx*sstr]) ;
      r3 = r*r*r*4.0*M_PI ; r5 = r3*r*r ;
      nR[0] = -dr[0]/r3 ;
      nR[1] = -dr[1]/r3 ;
      nR[2] = -dr[2]/r3 ;
      
      nR[ 3] = 3*dr[0]*dr[0]/r5 - 1.0/r3 ;
      nR[ 4] = 3*dr[0]*dr[1]/r5 ;
      nR[ 5] = 3*dr[0]*dr[2]/r5 ;
      nR[ 6] = 3*dr[1]*dr[0]/r5 ;
      nR[ 7] = 3*dr[1]*dr[1]/r5 - 1.0/r3 ;
      nR[ 8] = 3*dr[1]*dr[2]/r5 ;
      nR[ 9] = 3*dr[2]*dr[0]/r5 ;
      nR[10] = 3*dr[2]*dr[1]/r5 ;
      nR[11] = 3*dr[2]*dr[2]/r5 - 1.0/r3 ;

      /*\nabla(1/R)\times\omega*/
      wbfmm_vector_cross(df,nR,s) ;
      wbfmm_vector_inc(f,df) ;

      /*d/dx field[0,1,2]*/
      wbfmm_vector_cross(df,&(nR[3]),s) ;
      wbfmm_vector_inc(&(f[3]),df) ;

      /*d/dy field[0,1,2]*/
      wbfmm_vector_cross(df,&(nR[6]),s) ;
      wbfmm_vector_inc(&(f[6]),df) ;
      
      /*d/dz field[0,1,2]*/
      wbfmm_vector_cross(df,&(nR[9]),s) ;
      wbfmm_vector_inc(&(f[9]),df) ;
    }
  }

  return ;
}

/* static void box_curl_gradient_evaluate4(wbfmm_tree_t *t, */
/* 					gint i, */
/* 					WBFMM_REAL *src, gint sstr, */
/* 					WBFMM_REAL *x, WBFMM_REAL *f) */

/* { */
/*   /\* gint idx, j ; *\/ */
/*   /\* WBFMM_REAL *xs, dr[3], r, r3, r5, nR[12], *s, df[3] ; *\/ */

/*   gint idx[4], j ; */
/*   WBFMM_REAL *xs[4] ; */
/*   __attribute__ ((aligned (32))) WBFMM_REAL r[12] ; */

/*   for ( j = 0 ; j < 4 ; j ++ ) { */
/*     idx[j] = t->ip[i+j] ; */
/*     xs [j] = wbfmm_tree_point_index(t, idx[j]) ; */
/*     r[4*0+j] = x[0] - xs[j][0] ;  */
/*     r[4*1+j] = x[1] - xs[j][1] ;  */
/*     r[4*2+j] = x[2] - xs[j][2] ;  */
/*   } */

/*   gradient_evaluate4(r) ; */

/*   for ( j = 0 ; j < 4 ; j ++ ) { */
/*     f[0] -= (src[idx[j]*sstr+2]*r[4*1+j] - */
/* 	     src[idx[j]*sstr+1]*r[4*2+j])*0.25*M_1_PI ; */
/*     f[1] -= (src[idx[j]*sstr+0]*r[4*2+j] - */
/* 	     src[idx[j]*sstr+2]*r[4*0+j])*0.25*M_1_PI ; */
/*     f[2] -= (src[idx[j]*sstr+1]*r[4*0+j] - */
/* 	     src[idx[j]*sstr+0]*r[4*1+j])*0.25*M_1_PI ; */
/*   } */
  
/*   /\* for ( j = i0 ; j < i1 ; j ++ ) { *\/ */
/*   /\*   idx = t->ip[j] ; *\/ */
/*   /\*   xs = wbfmm_tree_point_index(t, idx) ; *\/ */
/*   /\*   wbfmm_vector_diff(dr,x,xs) ; *\/ */
/*   /\*   r = wbfmm_vector_length(dr) ; *\/ */
/*   /\*   if ( r > WBFMM_LOCAL_CUTOFF_RADIUS ) { *\/ */
/*   /\*     s = &(src[idx*sstr]) ; *\/ */
/*   /\*     r3 = r*r*r*4.0*M_PI ; r5 = r3*r*r ; *\/ */
/*   /\*     nR[0] = -dr[0]/r3 ; *\/ */
/*   /\*     nR[1] = -dr[1]/r3 ; *\/ */
/*   /\*     nR[2] = -dr[2]/r3 ; *\/ */
      
/*   /\*     nR[ 3] = 3*dr[0]*dr[0]/r5 - 1.0/r3 ; *\/ */
/*   /\*     nR[ 4] = 3*dr[0]*dr[1]/r5 ; *\/ */
/*   /\*     nR[ 5] = 3*dr[0]*dr[2]/r5 ; *\/ */
/*   /\*     nR[ 6] = 3*dr[1]*dr[0]/r5 ; *\/ */
/*   /\*     nR[ 7] = 3*dr[1]*dr[1]/r5 - 1.0/r3 ; *\/ */
/*   /\*     nR[ 8] = 3*dr[1]*dr[2]/r5 ; *\/ */
/*   /\*     nR[ 9] = 3*dr[2]*dr[0]/r5 ; *\/ */
/*   /\*     nR[10] = 3*dr[2]*dr[1]/r5 ; *\/ */
/*   /\*     nR[11] = 3*dr[2]*dr[2]/r5 - 1.0/r3 ; *\/ */

/*   /\*     /\\*\nabla(1/R)\times\omega*\\/ *\/ */
/*   /\*     wbfmm_vector_cross(df,nR,s) ; *\/ */
/*   /\*     wbfmm_vector_inc(f,df) ; *\/ */

/*   /\*     /\\*d/dx field[0,1,2]*\\/ *\/ */
/*   /\*     wbfmm_vector_cross(df,&(nR[3]),s) ; *\/ */
/*   /\*     wbfmm_vector_inc(&(f[3]),df) ; *\/ */

/*   /\*     /\\*d/dy field[0,1,2]*\\/ *\/ */
/*   /\*     wbfmm_vector_cross(df,&(nR[6]),s) ; *\/ */
/*   /\*     wbfmm_vector_inc(&(f[6]),df) ; *\/ */
      
/*   /\*     /\\*d/dz field[0,1,2]*\\/ *\/ */
/*   /\*     wbfmm_vector_cross(df,&(nR[9]),s) ; *\/ */
/*   /\*     wbfmm_vector_inc(&(f[9]),df) ; *\/ */
/*   /\*   } *\/ */
/*   /\* } *\/ */

/*   return ; */
/* } */

static void box_curl_evaluate4(wbfmm_tree_t *t,
			       gint i,
			       WBFMM_REAL *src, gint sstr,
			       WBFMM_REAL *x, WBFMM_REAL *f)

{
  gint idx[4], j ;
  WBFMM_REAL *xs[4] ;
  __attribute__ ((aligned (32))) WBFMM_REAL r[12] ;

  for ( j = 0 ; j < 4 ; j ++ ) {
    idx[j] = t->ip[i+j] ;
    xs [j] = wbfmm_tree_point_index(t, idx[j]) ;
    r[4*0+j] = x[0] - xs[j][0] ; 
    r[4*1+j] = x[1] - xs[j][1] ; 
    r[4*2+j] = x[2] - xs[j][2] ; 
  }

  gradient_evaluate4(r) ;

  for ( j = 0 ; j < 4 ; j ++ ) {
    f[0] -= (src[idx[j]*sstr+2]*r[4*1+j] -
	     src[idx[j]*sstr+1]*r[4*2+j])*0.25*M_1_PI ;
    f[1] -= (src[idx[j]*sstr+0]*r[4*2+j] -
	     src[idx[j]*sstr+2]*r[4*0+j])*0.25*M_1_PI ;
    f[2] -= (src[idx[j]*sstr+1]*r[4*0+j] -
	     src[idx[j]*sstr+0]*r[4*1+j])*0.25*M_1_PI ;
  }

  return ;
}

static void box_curl_evaluate4_sorted(char *y, gsize ysize,
				      WBFMM_REAL *src, gint sstr,
				      WBFMM_REAL *x, WBFMM_REAL *f)

{
  gint j ;
  WBFMM_REAL *xs ;
  __attribute__ ((aligned (32))) WBFMM_REAL r[12] ;

  for ( j = 0 ; j < 4 ; j ++ ) {
    xs = (WBFMM_REAL *)(&(y[j*ysize])) ;
    r[4*0+j] = x[0] - xs[0] ; 
    r[4*1+j] = x[1] - xs[1] ; 
    r[4*2+j] = x[2] - xs[2] ; 
  }

  gradient_evaluate4(r) ;

  for ( j = 0 ; j < 4 ; j ++ ) {
    f[0] -= (src[j*sstr+2]*r[4*1+j] - src[j*sstr+1]*r[4*2+j])*0.25*M_1_PI ;
    f[1] -= (src[j*sstr+0]*r[4*2+j] - src[j*sstr+2]*r[4*0+j])*0.25*M_1_PI ;
    f[2] -= (src[j*sstr+1]*r[4*0+j] - src[j*sstr+0]*r[4*1+j])*0.25*M_1_PI ;
  }

  return ;
}

static gint local_field_evaluate(WBFMM_REAL *x0, WBFMM_REAL *cfft,
				 gint cstr, gint N, gint nq,
				 WBFMM_REAL *xf, WBFMM_REAL *field,
				 WBFMM_REAL *work)

{
  WBFMM_REAL r, th, ph, rn ;
  WBFMM_REAL Cth, Sth, *Pn, *Pnm1 ;
  WBFMM_REAL *Cmph, *Smph ;
  gint n, m, idx, i ;
  
  Pnm1 = &(work[0]) ; Pn = &(Pnm1[N+1]) ;
  Cmph = &(Pn[N+1]) ; Smph = &(Cmph[N+1]) ;

  WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(x0, xf, &r, &th, &ph) ;
  Cth = COS(th) ; Sth = SIN(th) ; 

  /*initialize recursions*/
  WBFMM_FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth,
					   &(Pnm1[0]), &(Pn[0]), &(Pn[1])) ;
  Cmph[0] = 1.0 ; Cmph[1] = COS(ph) ; 
  Smph[0] = 0.0 ; Smph[1] = SIN(ph) ;

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
    /* Cmph[n] = Cmph[n-1]*Cmph[1] - Smph[n-1]*Smph[1] ; */
    /* Smph[n] = Smph[n-1]*Cmph[1] + Cmph[n-1]*Smph[1] ; */
    wbfmm_cos_sin_nph(Cmph,Smph,Cmph[1],Smph[1],n-1) ;
    
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

static gint box_local_field(wbfmm_tree_t *t, guint level, guint b,
			    WBFMM_REAL *x, WBFMM_REAL *f, gint fstr,
			    WBFMM_REAL *src, gint sstr,
			    WBFMM_REAL *d, gint dstr,
			    gboolean eval_neighbours,
			    WBFMM_REAL *work)

{
  WBFMM_REAL xb[3], wb, *C, *xs, r ;
  wbfmm_box_t *boxes, box ;
  guint64 neighbours[27] ;
  gint nnbr, i, j, k, idx, nq ;

  g_assert(t->problem == WBFMM_PROBLEM_LAPLACE ) ;

  nq = wbfmm_tree_source_size(t) ;

  boxes = t->boxes[level] ;
  C = boxes[b].mpr ;

  WBFMM_FUNCTION_NAME(wbfmm_tree_box_centre)(t, level, b, xb, &wb) ;
  
  /* WBFMM_FUNCTION_NAME(wbfmm_laplace_expansion_local_evaluate)(xb, C, 8*nq, */
  /* 							      t->order_r[level], */
  /* 							      nq, x, f, work) ; */
  local_field_evaluate(xb, C, 8*nq, t->order_r[level], nq, x, f, work) ;
  if ( !eval_neighbours ) return 0 ;

  if ( src == NULL && d == NULL ) return 0 ;
  if ( t->normals == NULL && d != NULL ) {
    g_error("%s: no normals in tree but dipole strengths specified "
	    "(d != NULL)",
	    __FUNCTION__) ;
  }

  /*add the contribution from sources in neighbour boxes*/
  nnbr = wbfmm_box_neighbours(level, b, neighbours) ;
  g_assert(nnbr >= 0 && nnbr < 28) ;

  if ( t->normals == NULL && d == NULL ) {
    /* monopoles only */
    for ( i = 0 ; i < nnbr ; i ++ ) {
      box = boxes[neighbours[i]] ;
      for ( j = 0 ; j < box.n ; j ++ ) {
	idx = t->ip[box.i+j] ;
	xs = wbfmm_tree_point_index(t, idx) ;
	r = (xs[0]-x[0])*(xs[0]-x[0]) + (xs[1]-x[1])*(xs[1]-x[1]) +
	  (xs[2]-x[2])*(xs[2]-x[2]) ;
	if ( r > WBFMM_LOCAL_CUTOFF_RADIUS*WBFMM_LOCAL_CUTOFF_RADIUS ) {
	  r = SQRT(r)*4.0*M_PI ;
	  for ( k = 0 ; k < nq ; k ++ ) {
	    f[k*fstr] += src[idx*sstr+k]/r ;
	  }
	}
      }
    }
    
    return 0 ;
  }

  if ( src == NULL && t->normals != NULL ) {
    /*dipoles only*/
    WBFMM_REAL th, ph, nr, *normal ;
    
    for ( i = 0 ; i < nnbr ; i ++ ) {
      box = boxes[neighbours[i]] ;
      for ( j = 0 ; j < box.n ; j ++ ) {
	idx = t->ip[box.i+j] ;
	xs = wbfmm_tree_point_index(t, idx) ;
	normal = wbfmm_tree_normal_index(t,idx) ;

	WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(xs, x, &r, &th, &ph) ;
	if ( r > WBFMM_LOCAL_CUTOFF_RADIUS ) {
	  nr =
	    (x[0] - xs[0])*normal[0] +
	    (x[1] - xs[1])*normal[1] + 
	    (x[2] - xs[2])*normal[2] ;
	  nr /= 4.0*M_PI*r*r*r ;
	  for ( k = 0 ; k < nq ; k ++ ) f[k*fstr] += d[idx*dstr+k]*nr ;
	}
      }
      
    } 

    return 0 ;
  }
  
  if ( src != NULL && t->normals != NULL ) {
    /*sources and dipoles*/
    WBFMM_REAL rr[3], nr, g, *normal ;
    
    for ( i = 0 ; i < nnbr ; i ++ ) {
      box = boxes[neighbours[i]] ;
      for ( j = 0 ; j < box.n ; j ++ ) {
	idx = t->ip[box.i+j] ;
	xs = wbfmm_tree_point_index(t, idx) ;
	normal = wbfmm_tree_normal_index(t,idx) ;

	rr[0] = x[0] - xs[0] ; 
	rr[1] = x[1] - xs[1] ; 
	rr[2] = x[2] - xs[2] ;
	r = rr[0]*rr[0] + rr[1]*rr[1] + rr[2]*rr[2] ;
	if ( r > WBFMM_LOCAL_CUTOFF_RADIUS*WBFMM_LOCAL_CUTOFF_RADIUS ) {
	  r = SQRT(r) ;
	  nr =
	    (x[0] - xs[0])*normal[0] +
	    (x[1] - xs[1])*normal[1] + 
	    (x[2] - xs[2])*normal[2] ;
	  g = 0.25*M_1_PI/r ;
	  for ( k = 0 ; k < nq ; k ++ ) {
	    f[k*fstr] += (d[idx*dstr+k]*nr/r/r + src[idx*sstr+k])*g ;
	  }
	}
      }
      
    } 

    return 0 ;
  }

  g_assert_not_reached() ; 
  
  return 0 ;
}

static gint local_grad_evaluate(WBFMM_REAL *x0, WBFMM_REAL *cfft,
				gint cstr, gint N, gint nq,
				WBFMM_REAL *xf, WBFMM_REAL *field,
				gint fstr, WBFMM_REAL *work)

{
  WBFMM_REAL r, th, ph, rnm1, cr, ci ;
  WBFMM_REAL Cth, Sth, *Pn, *Pnm1 ;
  WBFMM_REAL *Cmph, *Smph ;
  WBFMM_REAL dRnm[6] ;
  gint n, m, idx, i, i1 = 1, i2 = 2, i3 = 3 ;

  if ( fstr < 3 && nq != 1 )
    g_error("%s: field data stride (%d) must be greater than two",
	    __FUNCTION__, fstr) ;

  if ( N == 0 ) return 0 ;

  Pnm1 = &(work[0]) ; Pn = &(Pnm1[N+2]) ;
  Cmph = &(Pn[N+2]) ; Smph = &(Cmph[N+2]) ;

  WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(x0, xf, &r, &th, &ph) ;
  Cth = COS(th) ; Sth = SIN(th) ; 

  /*initialize recursions*/
  WBFMM_FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth,
					   &(Pnm1[0]), &(Pn[0]), &(Pn[1])) ;
  Cmph[0] = 1.0 ; Smph[0] = 0.0 ;
  Cmph[1] = COS(ph) ; Smph[1] = SIN(ph) ;

  /*first two terms by hand; gradient of zero order term is zero*/  
  n = 1 ; 
  m = 0 ; 
  rnm1 = 1.0 ;
  idx = n*n ;
  /* Cmph[n+1] = Cmph[n]*Cmph[1] - Smph[n]*Smph[1] ; */
  /* Smph[n+1] = Smph[n]*Cmph[1] + Cmph[n]*Smph[1] ; */
  /* fprintf(stderr, "Hello\n") ; */
  wbfmm_cos_sin_nph(Cmph,Smph,Cmph[1],Smph[1],n) ;
    
  Rnm_derivatives_1m0(n, m, rnm1, Pnm1, Cmph, Smph, dRnm) ;
  
  for ( i = 0 ; i < nq ; i ++ ) {
    cr = cfft[cstr*idx+i] ;
    
    field[fstr*i+2] += dRnm[WBFMM_DERIVATIVE_Z_R]*cr ;
  }

  m = 1 ; 
  idx = wbfmm_index_laplace_nm(n,m) ;

  Rnm_derivatives_1(n, m, rnm1, Pnm1, Cmph, Smph, dRnm) ;
  for ( i = 0 ; i < nq ; i ++ ) {
    cr = cfft[cstr*(idx+0)+i] ; ci = -cfft[cstr*(idx+1)+i] ;

#ifdef WBFMM_SINGLE_PRECISION
      blaswrap_saxpy(i3, cr, &(dRnm[0]), i2, &(field[i*fstr]), i1) ;
      blaswrap_saxpy(i3, ci, &(dRnm[1]), i2, &(field[i*fstr]), i1) ;
#else /*WBFMM_SINGLE_PRECISION*/
      blaswrap_daxpy(i3, cr, &(dRnm[0]), i2, &(field[i*fstr]), i1) ;
      blaswrap_daxpy(i3, ci, &(dRnm[1]), i2, &(field[i*fstr]), i1) ;
#endif /*WBFMM_SINGLE_PRECISION*/
  }
  
  for ( n = 2 ; n <= N ; n ++ ) {
    rnm1 *= r ;
    WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pnm1, &Pn,
							n-1, Cth, Sth) ;
    /* Cmph[n+1] = Cmph[n]*Cmph[1] - Smph[n]*Smph[1] ; */
    /* Smph[n+1] = Smph[n]*Cmph[1] + Cmph[n]*Smph[1] ; */
    wbfmm_cos_sin_nph(Cmph,Smph,Cmph[1],Smph[1],n) ;

    m = 0 ; 
    idx = n*n ;

    Rnm_derivatives_1m0(n, m, rnm1, Pnm1, Cmph, Smph, dRnm) ;
    for ( i = 0 ; i < nq ; i ++ ) {
      cr = cfft[cstr*idx+i] ;

#ifdef WBFMM_SINGLE_PRECISION
      blaswrap_saxpy(i3, cr, &(dRnm[0]), i2, &(field[i*fstr]), i1) ;
#else /*WBFMM_SINGLE_PRECISION*/
      blaswrap_daxpy(i3, cr, &(dRnm[0]), i2, &(field[i*fstr]), i1) ;
#endif /*WBFMM_SINGLE_PRECISION*/
    }

    for ( m = 1 ; m <= n ; m ++ ) {
      idx = wbfmm_index_laplace_nm(n,m) ;

      Rnm_derivatives_1(n, m, rnm1, Pnm1, Cmph, Smph, dRnm) ;
      for ( i = 0 ; i < nq ; i ++ ) {
	cr = cfft[cstr*(idx+0)+i] ; ci = -cfft[cstr*(idx+1)+i] ;

#ifdef WBFMM_SINGLE_PRECISION
      blaswrap_saxpy(i3, cr, &(dRnm[0]), i2, &(field[i*fstr]), i1) ;
      blaswrap_saxpy(i3, ci, &(dRnm[1]), i2, &(field[i*fstr]), i1) ;
#else /*WBFMM_SINGLE_PRECISION*/
      blaswrap_daxpy(i3, cr, &(dRnm[0]), i2, &(field[i*fstr]), i1) ;
      blaswrap_daxpy(i3, ci, &(dRnm[1]), i2, &(field[i*fstr]), i1) ;
#endif /*WBFMM_SINGLE_PRECISION*/
      }
    }
  }
  
  return 0 ;
}

static gint local_curl_evaluate(WBFMM_REAL *x0, WBFMM_REAL*cfft, gint cstr,
				gint N,	gint nq,
				WBFMM_REAL *xf,	WBFMM_REAL *field, gint fstr,
				WBFMM_REAL *work)

{
  WBFMM_REAL r, th, ph, cr, ci ;
  WBFMM_REAL Cth, Sth, *Pn, *Pnm1 ;
  WBFMM_REAL *Cmph, *Smph ;
  WBFMM_REAL dRnm[6], rnm1 ;
  gint n, m, idx ;

  /*
    fstr is ignored: the curl based on the first three components of
    the source is placed into the first three components of f
   */
  
  if ( nq < 3 )
    g_error("%s: not enough source components (%d) for curl calculation",
	    __FUNCTION__, nq) ;

  if ( N == 0 ) return 0 ;

  Pnm1 = &(work[0]) ; Pn = &(Pnm1[N+2]) ;
  Cmph = &(Pn[N+2]) ; Smph = &(Cmph[N+2]) ;

  WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(x0, xf, &r, &th, &ph) ;
  Cth = COS(th) ; Sth = SIN(th) ; 

  /*initialize recursions*/
  WBFMM_FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth,
					   &(Pnm1[0]), &(Pn[0]), &(Pn[1])) ;
  Cmph[0] = 1.0 ; Smph[0] = 0.0 ;
  Cmph[1] = COS(ph) ; Smph[1] = SIN(ph) ;

  /*first two terms by hand; curl of zero order term is zero*/  
  n = 1 ; 
  m = 0 ; 
  rnm1 = 1.0 ;
  idx = n*n ;
  /* Cmph[n+1] = Cmph[n]*Cmph[1] - Smph[n]*Smph[1] ; */
  /* Smph[n+1] = Smph[n]*Cmph[1] + Cmph[n]*Smph[1] ; */
  /* fprintf(stderr, "Hello\n") ; */
  wbfmm_cos_sin_nph(Cmph,Smph,Cmph[1],Smph[1],n) ;

  Rnm_derivatives_1m0(n,m,rnm1,Pnm1,Cmph,Smph,dRnm) ;
  
  cr = cfft[cstr*idx+1] ;
  field[0] -= dRnm[WBFMM_DERIVATIVE_Z_R]*cr ;
  cr = cfft[cstr*idx+0] ;
  field[1] += dRnm[WBFMM_DERIVATIVE_Z_R]*cr ;
  
  m = 1 ; 
  idx = wbfmm_index_laplace_nm(n,m) ;
  Rnm_derivatives_1(n,m,rnm1,Pnm1,Cmph,Smph,dRnm) ;
  
  cr = cfft[cstr*(idx+0)+2] ; ci = cfft[cstr*(idx+1)+2] ;
  field[0] += dRnm[WBFMM_DERIVATIVE_Y_R]*cr - dRnm[WBFMM_DERIVATIVE_Y_I]*ci ;
  field[1] -= dRnm[WBFMM_DERIVATIVE_X_R]*cr - dRnm[WBFMM_DERIVATIVE_X_I]*ci ;

  cr = cfft[cstr*(idx+0)+1] ; ci = cfft[cstr*(idx+1)+1] ;
  field[0] -= dRnm[WBFMM_DERIVATIVE_Z_R]*cr - dRnm[WBFMM_DERIVATIVE_Z_I]*ci ;
  field[2] += dRnm[WBFMM_DERIVATIVE_X_R]*cr - dRnm[WBFMM_DERIVATIVE_X_I]*ci ;

  cr = cfft[cstr*(idx+0)+0] ; ci = cfft[cstr*(idx+1)+0] ;
  field[1] += dRnm[WBFMM_DERIVATIVE_Z_R]*cr - dRnm[WBFMM_DERIVATIVE_Z_I]*ci ;
  field[2] -= dRnm[WBFMM_DERIVATIVE_Y_R]*cr - dRnm[WBFMM_DERIVATIVE_Y_I]*ci ;
  
  for ( n = 2 ; n <= N ; n ++ ) {
    rnm1 *= r ;
    WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pnm1, &Pn,
							n-1, Cth, Sth) ;
    /* Cmph[n+1] = Cmph[n]*Cmph[1] - Smph[n]*Smph[1] ; */
    /* Smph[n+1] = Smph[n]*Cmph[1] + Cmph[n]*Smph[1] ; */
    wbfmm_cos_sin_nph(Cmph,Smph,Cmph[1],Smph[1],n) ;

    m = 0 ; 
    idx = n*n ;
    Rnm_derivatives_1m0(n,m,rnm1,Pnm1,Cmph,Smph,dRnm) ;

    cr = cfft[cstr*idx+2] ;
    field[0] += dRnm[WBFMM_DERIVATIVE_Y_R]*cr ;
    field[1] -= dRnm[WBFMM_DERIVATIVE_X_R]*cr ;

    cr = cfft[cstr*idx+1] ;
    field[0] -= dRnm[WBFMM_DERIVATIVE_Z_R]*cr ;
    field[2] += dRnm[WBFMM_DERIVATIVE_X_R]*cr ;

    cr = cfft[cstr*idx+0] ;
    field[1] += dRnm[WBFMM_DERIVATIVE_Z_R]*cr ;
    field[2] -= dRnm[WBFMM_DERIVATIVE_Y_R]*cr ;

    for ( m = 1 ; m <= n ; m ++ ) {
      idx = wbfmm_index_laplace_nm(n,m) ;

      Rnm_derivatives_1(n,m,rnm1,Pnm1,Cmph,Smph,dRnm) ;
  
      cr = cfft[cstr*(idx+0)+2] ; ci = cfft[cstr*(idx+1)+2] ;
      field[0] +=
	dRnm[WBFMM_DERIVATIVE_Y_R]*cr - dRnm[WBFMM_DERIVATIVE_Y_I]*ci ;
      field[1] -=
	dRnm[WBFMM_DERIVATIVE_X_R]*cr - dRnm[WBFMM_DERIVATIVE_X_I]*ci ;

      cr = cfft[cstr*(idx+0)+1] ; ci = cfft[cstr*(idx+1)+1] ;
      field[0] -=
	dRnm[WBFMM_DERIVATIVE_Z_R]*cr - dRnm[WBFMM_DERIVATIVE_Z_I]*ci ;
      field[2] +=
	dRnm[WBFMM_DERIVATIVE_X_R]*cr - dRnm[WBFMM_DERIVATIVE_X_I]*ci ;

      cr = cfft[cstr*(idx+0)+0] ; ci = cfft[cstr*(idx+1)+0] ;
      field[1] +=
	dRnm[WBFMM_DERIVATIVE_Z_R]*cr - dRnm[WBFMM_DERIVATIVE_Z_I]*ci ;
      field[2] -=
	dRnm[WBFMM_DERIVATIVE_Y_R]*cr - dRnm[WBFMM_DERIVATIVE_Y_I]*ci ;
    }
  }
  
  return 0 ;
}

static gint local_curl_gradient_evaluate(WBFMM_REAL *x0, WBFMM_REAL*cfft,
					 gint cstr, gint N, gint nq,
					 WBFMM_REAL *xf, WBFMM_REAL *field,
					 gint fstr, WBFMM_REAL *work)

{
  WBFMM_REAL r, th, ph, cr, ci ;
  WBFMM_REAL Cth, Sth, *Pn, *Pnm1, *Pnm2 ;
  WBFMM_REAL *Cmph, *Smph ;
  WBFMM_REAL dRnm[6], d2Rnm[12], rnm1, rnm2 ;
  gint n, m, idx ;

  /*
    fstr is ignored: the curl based on the first three components of
    the source is placed into the first three components of f,
    followed by the nine components of the gradient
   */
  
  if ( nq < 3 )
    g_error("%s: not enough source components (%d) for curl calculation",
	    __FUNCTION__, nq) ;

  if ( N == 0 ) return 0 ;

  Pnm1 = &(work[0]) ; Pn = &(Pnm1[N+2]) ; Pnm2 = &(Pn[N+2]) ;
  memset(Pnm1, 0, (N+2)*sizeof(gdouble)) ;
  memset(Pn  , 0, (N+2)*sizeof(gdouble)) ;
  memset(Pnm2, 0, (N+2)*sizeof(gdouble)) ;
  Cmph = &(Pnm2[N+2]) ; Smph = &(Cmph[N+2]) ;
  memset(Cmph, 0, (N+2)*sizeof(gdouble)) ;
  memset(Smph, 0, (N+2)*sizeof(gdouble)) ;

  WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(x0, xf, &r, &th, &ph) ;
  Cth = COS(th) ; Sth = SIN(th) ; 

  /*initialize recursions*/
  WBFMM_FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth,
					   &(Pnm1[0]), &(Pn[0]), &(Pn[1])) ;
  Cmph[0] = 1.0 ; Smph[0] = 0.0 ;
  Cmph[1] = COS(ph) ; Smph[1] = SIN(ph) ;

  /*first two terms by hand; curl of zero order term is zero*/  
  n = 1 ; 
  m = 0 ; 
  rnm1 = 1.0 ;
  idx = n*n ;
  /* Cmph[n+1] = Cmph[n]*Cmph[1] - Smph[n]*Smph[1] ; */
  /* Smph[n+1] = Smph[n]*Cmph[1] + Cmph[n]*Smph[1] ; */
  wbfmm_cos_sin_nph(Cmph,Smph,Cmph[1],Smph[1],n) ;

  Rnm_derivatives_1m0(n, m, rnm1, Pnm1, Cmph, Smph, dRnm) ;
  
  cr = cfft[cstr*idx+1] ;
  field[0] -= dRnm[WBFMM_DERIVATIVE_Z_R]*cr ;
  cr = cfft[cstr*idx+0] ;
  field[1] += dRnm[WBFMM_DERIVATIVE_Z_R]*cr ;
  
  m = 1 ; 
  idx = wbfmm_index_laplace_nm(n,m) ;
  Rnm_derivatives_1(n,m,rnm1,Pnm1,Cmph,Smph,dRnm) ;
  
  cr = cfft[cstr*(idx+0)+2] ; ci = cfft[cstr*(idx+1)+2] ;
  field[0] += dRnm[WBFMM_DERIVATIVE_Y_R]*cr - dRnm[WBFMM_DERIVATIVE_Y_I]*ci ;
  field[1] -= dRnm[WBFMM_DERIVATIVE_X_R]*cr - dRnm[WBFMM_DERIVATIVE_X_I]*ci ;

  cr = cfft[cstr*(idx+0)+1] ; ci = cfft[cstr*(idx+1)+1] ;
  field[0] -= dRnm[WBFMM_DERIVATIVE_Z_R]*cr - dRnm[WBFMM_DERIVATIVE_Z_I]*ci ;
  field[2] += dRnm[WBFMM_DERIVATIVE_X_R]*cr - dRnm[WBFMM_DERIVATIVE_X_I]*ci ;

  cr = cfft[cstr*(idx+0)+0] ; ci = cfft[cstr*(idx+1)+0] ;
  field[1] += dRnm[WBFMM_DERIVATIVE_Z_R]*cr - dRnm[WBFMM_DERIVATIVE_Z_I]*ci ;
  field[2] -= dRnm[WBFMM_DERIVATIVE_Y_R]*cr - dRnm[WBFMM_DERIVATIVE_Y_I]*ci ;

  rnm2 = 1.0 ;
  for ( n = 2 ; n <= N ; n ++ ) {
    rnm1 *= r ;

    memcpy(Pnm2, Pnm1, (N+2)*sizeof(WBFMM_REAL)) ;
    WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pnm1, &Pn,
							n-1, Cth, Sth) ;
    wbfmm_cos_sin_nph(Cmph,Smph,Cmph[1],Smph[1],n) ;

    m = 0 ; 
    idx = n*n ;
    Rnm_derivatives_1m0(n, m, rnm1, Pnm1, Cmph, Smph, dRnm) ;
    Rnm_derivatives_2m0(n, m, rnm2, Pnm2, Cmph, Smph, d2Rnm) ;

    cr = cfft[cstr*idx+2] ;
    field[ 0] += dRnm[WBFMM_DERIVATIVE_Y_R]*cr ;

    field[ 3] += d2Rnm[WBFMM_DERIVATIVE_XY_R]*cr ;
    field[ 6] += d2Rnm[WBFMM_DERIVATIVE_YY_R]*cr ;
    field[ 9] += d2Rnm[WBFMM_DERIVATIVE_YZ_R]*cr ;
    
    field[ 1] -= dRnm[WBFMM_DERIVATIVE_X_R]*cr ;

    field[ 4] -= d2Rnm[WBFMM_DERIVATIVE_XX_R]*cr ;
    field[ 7] -= d2Rnm[WBFMM_DERIVATIVE_XY_R]*cr ;
    field[10] -= d2Rnm[WBFMM_DERIVATIVE_ZX_R]*cr ;

    cr = cfft[cstr*idx+1] ;
    field[ 0] -= dRnm [WBFMM_DERIVATIVE_Z_R]*cr ;
    field[ 3] -= d2Rnm[WBFMM_DERIVATIVE_ZX_R]*cr ;
    field[ 6] -= d2Rnm[WBFMM_DERIVATIVE_YZ_R]*cr ;
    field[ 9] -= d2Rnm[WBFMM_DERIVATIVE_ZZ_R]*cr ;

    field[ 2] += dRnm[WBFMM_DERIVATIVE_X_R]*cr ;

    field[ 5] += d2Rnm[WBFMM_DERIVATIVE_XX_R]*cr ;
    field[ 8] += d2Rnm[WBFMM_DERIVATIVE_XY_R]*cr ;
    field[11] += d2Rnm[WBFMM_DERIVATIVE_ZX_R]*cr ;

    cr = cfft[cstr*idx+0] ;
    field[ 1] += dRnm[WBFMM_DERIVATIVE_Z_R]*cr ;

    field[ 4] += d2Rnm[WBFMM_DERIVATIVE_ZX_R]*cr ;
    field[ 7] += d2Rnm[WBFMM_DERIVATIVE_YZ_R]*cr ;
    field[10] += d2Rnm[WBFMM_DERIVATIVE_ZZ_R]*cr ;

    field[ 2] -= dRnm[WBFMM_DERIVATIVE_Y_R]*cr ;

    field[ 5] -= d2Rnm[WBFMM_DERIVATIVE_XY_R]*cr ;
    field[ 8] -= d2Rnm[WBFMM_DERIVATIVE_YY_R]*cr ;
    field[11] -= d2Rnm[WBFMM_DERIVATIVE_YZ_R]*cr ;

    m = 1 ;
    idx = wbfmm_index_laplace_nm(n,m) ;

    Rnm_derivatives_1(n, m, rnm1, Pnm1, Cmph, Smph, dRnm) ;
    Rnm_derivatives_2m1(n, m, rnm2, Pnm2, Cmph, Smph, d2Rnm) ;
    
    cr = cfft[cstr*(idx+0)+2] ; ci = cfft[cstr*(idx+1)+2] ;
    field[ 0] +=
      dRnm[WBFMM_DERIVATIVE_Y_R]*cr - dRnm[WBFMM_DERIVATIVE_Y_I]*ci ;
    field[ 3] +=
      d2Rnm[WBFMM_DERIVATIVE_XY_R]*cr - d2Rnm[WBFMM_DERIVATIVE_XY_I]*ci ;
    field[ 6] +=
      d2Rnm[WBFMM_DERIVATIVE_YY_R]*cr - d2Rnm[WBFMM_DERIVATIVE_YY_I]*ci ;
    field[ 9] +=
      d2Rnm[WBFMM_DERIVATIVE_YZ_R]*cr - d2Rnm[WBFMM_DERIVATIVE_YZ_I]*ci ;

    field[ 1] -=
      dRnm[WBFMM_DERIVATIVE_X_R]*cr - dRnm[WBFMM_DERIVATIVE_X_I]*ci ;

    field[ 4] -=
      d2Rnm[WBFMM_DERIVATIVE_XX_R]*cr - d2Rnm[WBFMM_DERIVATIVE_XX_I]*ci ;
    field[ 7] -=
      d2Rnm[WBFMM_DERIVATIVE_XY_R]*cr - d2Rnm[WBFMM_DERIVATIVE_XY_I]*ci ;
    field[10] -=
      d2Rnm[WBFMM_DERIVATIVE_ZX_R]*cr - d2Rnm[WBFMM_DERIVATIVE_ZX_I]*ci ;
    
    cr = cfft[cstr*(idx+0)+1] ; ci = cfft[cstr*(idx+1)+1] ;
    field[ 0] -=
      dRnm[WBFMM_DERIVATIVE_Z_R]*cr - dRnm[WBFMM_DERIVATIVE_Z_I]*ci ;
    field[ 3] -=
      d2Rnm[WBFMM_DERIVATIVE_ZX_R]*cr - d2Rnm[WBFMM_DERIVATIVE_ZX_I]*ci ;
    field[ 6] -=
      d2Rnm[WBFMM_DERIVATIVE_YZ_R]*cr - d2Rnm[WBFMM_DERIVATIVE_YZ_I]*ci ;
    field[ 9] -=
      d2Rnm[WBFMM_DERIVATIVE_ZZ_R]*cr - d2Rnm[WBFMM_DERIVATIVE_ZZ_I]*ci ;

    field[ 2] +=
      dRnm[WBFMM_DERIVATIVE_X_R]*cr - dRnm[WBFMM_DERIVATIVE_X_I]*ci ;
    field[ 5] +=
      d2Rnm[WBFMM_DERIVATIVE_XX_R]*cr - d2Rnm[WBFMM_DERIVATIVE_XX_I]*ci ;
    field[ 8] +=
      d2Rnm[WBFMM_DERIVATIVE_XY_R]*cr - d2Rnm[WBFMM_DERIVATIVE_XY_I]*ci ;
    field[11] +=
      d2Rnm[WBFMM_DERIVATIVE_ZX_R]*cr - d2Rnm[WBFMM_DERIVATIVE_ZX_I]*ci ;
    
    cr = cfft[cstr*(idx+0)+0] ; ci = cfft[cstr*(idx+1)+0] ;
    field[ 1] +=
      dRnm[WBFMM_DERIVATIVE_Z_R]*cr - dRnm[WBFMM_DERIVATIVE_Z_I]*ci ;

    field[ 4] +=
      d2Rnm[WBFMM_DERIVATIVE_ZX_R]*cr - d2Rnm[WBFMM_DERIVATIVE_ZX_I]*ci ;
    field[ 7] +=
      d2Rnm[WBFMM_DERIVATIVE_YZ_R]*cr - d2Rnm[WBFMM_DERIVATIVE_YZ_I]*ci ;
    field[10] +=
      d2Rnm[WBFMM_DERIVATIVE_ZZ_R]*cr - d2Rnm[WBFMM_DERIVATIVE_ZZ_I]*ci ;

    field[ 2] -=
      dRnm[WBFMM_DERIVATIVE_Y_R]*cr - dRnm[WBFMM_DERIVATIVE_Y_I]*ci ;
    field[ 5] -=
      d2Rnm[WBFMM_DERIVATIVE_XY_R]*cr - d2Rnm[WBFMM_DERIVATIVE_XY_I]*ci ;
    field[ 8] -=
      d2Rnm[WBFMM_DERIVATIVE_YY_R]*cr - d2Rnm[WBFMM_DERIVATIVE_YY_I]*ci ;
    field[11] -=
      d2Rnm[WBFMM_DERIVATIVE_YZ_R]*cr - d2Rnm[WBFMM_DERIVATIVE_YZ_I]*ci ;
    
    for ( m = 2 ; m <= n ; m ++ ) {
      idx = wbfmm_index_laplace_nm(n,m) ;

      Rnm_derivatives_1(n, m, rnm1, Pnm1, Cmph, Smph, dRnm) ;
      Rnm_derivatives_2(n, m, rnm2, Pnm2, Cmph, Smph, d2Rnm) ;
  
      cr = cfft[cstr*(idx+0)+2] ; ci = cfft[cstr*(idx+1)+2] ;
      field[ 0] +=
	dRnm[WBFMM_DERIVATIVE_Y_R]*cr - dRnm[WBFMM_DERIVATIVE_Y_I]*ci ;
      field[ 3] +=
	d2Rnm[WBFMM_DERIVATIVE_XY_R]*cr - d2Rnm[WBFMM_DERIVATIVE_XY_I]*ci ;
      field[ 6] +=
	d2Rnm[WBFMM_DERIVATIVE_YY_R]*cr - d2Rnm[WBFMM_DERIVATIVE_YY_I]*ci ;
      field[ 9] +=
	d2Rnm[WBFMM_DERIVATIVE_YZ_R]*cr - d2Rnm[WBFMM_DERIVATIVE_YZ_I]*ci ;

      field[ 1] -=
	dRnm[WBFMM_DERIVATIVE_X_R]*cr - dRnm[WBFMM_DERIVATIVE_X_I]*ci ;

      field[ 4] -=
	d2Rnm[WBFMM_DERIVATIVE_XX_R]*cr - d2Rnm[WBFMM_DERIVATIVE_XX_I]*ci ;
      field[ 7] -=
	d2Rnm[WBFMM_DERIVATIVE_XY_R]*cr - d2Rnm[WBFMM_DERIVATIVE_XY_I]*ci ;
      field[10] -=
	d2Rnm[WBFMM_DERIVATIVE_ZX_R]*cr - d2Rnm[WBFMM_DERIVATIVE_ZX_I]*ci ;

      cr = cfft[cstr*(idx+0)+1] ; ci = cfft[cstr*(idx+1)+1] ;
      field[ 0] -=
	dRnm[WBFMM_DERIVATIVE_Z_R]*cr - dRnm[WBFMM_DERIVATIVE_Z_I]*ci ;
      field[ 3] -=
	d2Rnm[WBFMM_DERIVATIVE_ZX_R]*cr - d2Rnm[WBFMM_DERIVATIVE_ZX_I]*ci ;
      field[ 6] -=
	d2Rnm[WBFMM_DERIVATIVE_YZ_R]*cr - d2Rnm[WBFMM_DERIVATIVE_YZ_I]*ci ;
      field[ 9] -=
	d2Rnm[WBFMM_DERIVATIVE_ZZ_R]*cr - d2Rnm[WBFMM_DERIVATIVE_ZZ_I]*ci ;

      field[ 2] +=
	dRnm[WBFMM_DERIVATIVE_X_R]*cr - dRnm[WBFMM_DERIVATIVE_X_I]*ci ;
      field[ 5] +=
	d2Rnm[WBFMM_DERIVATIVE_XX_R]*cr - d2Rnm[WBFMM_DERIVATIVE_XX_I]*ci ;
      field[ 8] +=
	d2Rnm[WBFMM_DERIVATIVE_XY_R]*cr - d2Rnm[WBFMM_DERIVATIVE_XY_I]*ci ;
      field[11] +=
	d2Rnm[WBFMM_DERIVATIVE_ZX_R]*cr - d2Rnm[WBFMM_DERIVATIVE_ZX_I]*ci ;

      cr = cfft[cstr*(idx+0)+0] ; ci = cfft[cstr*(idx+1)+0] ;
      field[ 1] +=
	dRnm[WBFMM_DERIVATIVE_Z_R]*cr - dRnm[WBFMM_DERIVATIVE_Z_I]*ci ;
      field[ 4] +=
	d2Rnm[WBFMM_DERIVATIVE_ZX_R]*cr - d2Rnm[WBFMM_DERIVATIVE_ZX_I]*ci ;
      field[ 7] +=
	d2Rnm[WBFMM_DERIVATIVE_YZ_R]*cr - d2Rnm[WBFMM_DERIVATIVE_YZ_I]*ci ;
      field[10] +=
	d2Rnm[WBFMM_DERIVATIVE_ZZ_R]*cr - d2Rnm[WBFMM_DERIVATIVE_ZZ_I]*ci ;

      field[ 2] -=
	dRnm[WBFMM_DERIVATIVE_Y_R]*cr - dRnm[WBFMM_DERIVATIVE_Y_I]*ci ;
      field[ 5] -=
	d2Rnm[WBFMM_DERIVATIVE_XY_R]*cr - d2Rnm[WBFMM_DERIVATIVE_XY_I]*ci ;
      field[ 8] -=
	d2Rnm[WBFMM_DERIVATIVE_YY_R]*cr - d2Rnm[WBFMM_DERIVATIVE_YY_I]*ci ;
      field[11] -=
	d2Rnm[WBFMM_DERIVATIVE_YZ_R]*cr - d2Rnm[WBFMM_DERIVATIVE_YZ_I]*ci ;
    }

    rnm2 *= r ;
  }
  
  return 0 ;
}

static gint tree_laplace_box_local_curl(wbfmm_tree_t *t,
					guint level,
					guint b,
					WBFMM_REAL *x,
					WBFMM_REAL *f,
					gint fstr,
					WBFMM_REAL *src,
					gint sstr,
					WBFMM_REAL *d,
					gint dstr,
					gboolean
					eval_neighbours,
					WBFMM_REAL *work)

{
  WBFMM_REAL xb[3], wb, *C ;
  wbfmm_box_t *boxes, *box ;
  guint64 neighbours[27] ;
  gint nnbr, i, j, nq ;

  g_assert(t->problem == WBFMM_PROBLEM_LAPLACE ) ;

  nq = wbfmm_tree_source_size(t) ;

  boxes = t->boxes[level] ;
  C = boxes[b].mpr ;

  WBFMM_FUNCTION_NAME(wbfmm_tree_box_centre)(t, level, b, xb, &wb) ;
  
  local_curl_evaluate(xb, C, 8*nq, t->order_r[level], nq, x, f, fstr, work) ;
  
  if ( !eval_neighbours ) return 0 ;

  if ( src == NULL && d == NULL ) return 0 ;
  
  if ( t->normals == NULL && d != NULL ) {
    g_error("%s: no normals in tree but dipole strengths specified "
	    "(d != NULL)",
	    __FUNCTION__) ;
  }

  /*add the contribution from sources in neighbour boxes*/
  nnbr = wbfmm_box_neighbours(level, b, neighbours) ;
  g_assert(nnbr >= 0 && nnbr < 28) ;

  if ( d == NULL ) {
    /* monopoles only */
    if ( t->sorted ) {
      for ( i = 0 ; i < nnbr ; i ++ ) {
	char *y ;
	gsize ysize = t->pstr ;
	WBFMM_REAL *sy ;
	box = &(boxes[neighbours[i]]) ;
	y = &(t->points[(box->i)*ysize]) ;
	sy = &(src[(box->i)*sstr]) ;
	for ( j = 0 ; j < (gint)(box->n)-4 ; j += 4 )
	  box_curl_evaluate4_sorted(&(y[j*ysize]), ysize,
				    &(sy[j*sstr]), sstr, x, f) ;

	box_curl_evaluate(t, box->i+j, box->i+box->n, src, sstr, x, f) ;
      }
    } else {
      for ( i = 0 ; i < nnbr ; i ++ ) {
	box = &(boxes[neighbours[i]]) ;
	for ( j = 0 ; j < (gint)(box->n)-4 ; j += 4 ) 
	  box_curl_evaluate4(t, box->i+j, src, sstr, x, f) ;
	box_curl_evaluate(t, box->i+j, box->i+box->n, src, sstr, x, f) ;      
      }
    }

    return 0 ;
  }

  g_assert_not_reached() ;
  
  return 0 ;
}

static gint tree_laplace_box_local_grad(wbfmm_tree_t *t,
					guint level,
					guint b,
					WBFMM_REAL *x,
					WBFMM_REAL *f,
					gint fstr,
					WBFMM_REAL *src,
					gint sstr,
					WBFMM_REAL *d,
					gint dstr,
					gboolean
					eval_neighbours,
					WBFMM_REAL *work)

{
  WBFMM_REAL xb[3], wb, *C, *xs, r, nR[3] ;
  wbfmm_box_t *boxes, box ;
  guint64 neighbours[27] ;
  gint nnbr, i, j, k, idx, nq ;

  g_assert(t->problem == WBFMM_PROBLEM_LAPLACE ) ;

  nq = wbfmm_tree_source_size(t) ;

  boxes = t->boxes[level] ;
  C = boxes[b].mpr ;

  WBFMM_FUNCTION_NAME(wbfmm_tree_box_centre)(t, level, b, xb, &wb) ;

  local_grad_evaluate(xb, C, 8*nq, t->order_r[level], nq, x, f, fstr, work) ;
  
  if ( !eval_neighbours ) return 0 ;

  if ( src == NULL && d == NULL ) return 0 ;

  if ( t->normals == NULL && d != NULL ) {
    g_error("%s: no normals in tree but dipole strengths specified "
	    "(d != NULL)",
	    __FUNCTION__) ;
  }

  /*add the contribution from sources in neighbour boxes*/
  nnbr = wbfmm_box_neighbours(level, b, neighbours) ;
  g_assert(nnbr >= 0 && nnbr < 28) ;

  if ( d == NULL ) {
    /* monopoles only */
    for ( i = 0 ; i < nnbr ; i ++ ) {
      box = boxes[neighbours[i]] ;
      for ( j = 0 ; j < box.n ; j ++ ) {
	idx = t->ip[box.i+j] ;
	xs = wbfmm_tree_point_index(t, idx) ;
	r = (xs[0]-x[0])*(xs[0]-x[0]) + (xs[1]-x[1])*(xs[1]-x[1]) +
	  (xs[2]-x[2])*(xs[2]-x[2]) ;
	if ( r > WBFMM_LOCAL_CUTOFF_RADIUS*WBFMM_LOCAL_CUTOFF_RADIUS ) {
	  nR[0] = (x[0] - xs[0])/r ;
	  nR[1] = (x[1] - xs[1])/r ;
	  nR[2] = (x[2] - xs[2])/r ;
	  r = SQRT(r)*4.0*M_PI ;
	  nR[0] /= r ; nR[1] /= r ; nR[2] /= r ; 
	  for ( k = 0 ; k < nq ; k ++ ) {
	    f[k*fstr+0] -= src[idx*sstr+k]*nR[0] ;
	    f[k*fstr+1] -= src[idx*sstr+k]*nR[1] ;
	    f[k*fstr+2] -= src[idx*sstr+k]*nR[2] ;
	  }
	}
      }
    }
    
    return 0 ;
  }

  g_assert_not_reached() ;
  
  if ( src == NULL && d != NULL ) {
    /*dipoles only*/
    /* g_assert_not_reached() ; */
    WBFMM_REAL th, ph, nr ;
    WBFMM_REAL *normal ;
    
    for ( i = 0 ; i < nnbr ; i ++ ) {
      box = boxes[neighbours[i]] ;
      for ( j = 0 ; j < box.n ; j ++ ) {
	idx = t->ip[box.i+j] ;
	xs = wbfmm_tree_point_index(t, idx) ;
	normal = wbfmm_tree_normal_index(t,idx) ;

	WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(xs, x, &r, &th, &ph) ;
	if ( r > WBFMM_LOCAL_CUTOFF_RADIUS ) {
	  nr =
	    (x[0] - xs[0])*normal[0] +
	    (x[1] - xs[1])*normal[1] + 
	    (x[2] - xs[2])*normal[2] ;
	  nr /= 4.0*M_PI*r*r*r ;
	  for ( k = 0 ; k < nq ; k ++ ) f[k] += d[idx*dstr+k]*nr ;
	}
      }
      
    } 

    return 0 ;
  }
  
  if ( src != NULL && d != NULL ) {
    /*sources and dipoles*/
    WBFMM_REAL th, ph, nr, g ;
    WBFMM_REAL *normal ;
    
    for ( i = 0 ; i < nnbr ; i ++ ) {
      box = boxes[neighbours[i]] ;
      for ( j = 0 ; j < box.n ; j ++ ) {
	idx = t->ip[box.i+j] ;
	xs = wbfmm_tree_point_index(t, idx) ;
	normal = wbfmm_tree_normal_index(t,idx) ;

	WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(xs, x, &r, &th, &ph) ;
	if ( r > WBFMM_LOCAL_CUTOFF_RADIUS ) {
	  nr =
	    (x[0] - xs[0])*normal[0] +
	    (x[1] - xs[1])*normal[1] + 
	    (x[2] - xs[2])*normal[2] ;
	  g = 0.25*M_1_PI/r ;
	  for ( k = 0 ; k < nq ; k ++ ) {
	    f[k] += (d[idx*dstr+k]*nr/r/r + src[idx*sstr+k])*g ;
	  }
	}
      }
      
    } 

    return 0 ;
  }

  g_assert_not_reached() ; 
  
  return 0 ;
}

static gint tree_laplace_box_local_curl_gradient(wbfmm_tree_t *t,
						 guint level,
						 guint b,
						 WBFMM_REAL *x,
						 WBFMM_REAL *f,
						 gint fstr,
						 WBFMM_REAL *src,
						 gint sstr,
						 WBFMM_REAL *d,
						 gint dstr,
						 gboolean
						 eval_neighbours,
						 WBFMM_REAL *work)

{
  WBFMM_REAL xb[3], wb, *C ;
  wbfmm_box_t *boxes, *box ;
  guint64 neighbours[27] ;
  gint nnbr, i, j, nq ;

  g_assert(t->problem == WBFMM_PROBLEM_LAPLACE ) ;

  nq = wbfmm_tree_source_size(t) ;

  boxes = t->boxes[level] ;
  C = boxes[b].mpr ;

  WBFMM_FUNCTION_NAME(wbfmm_tree_box_centre)(t, level, b, xb, &wb) ;
  
  local_curl_gradient_evaluate(xb, C, 8*nq, t->order_r[level], nq,
			       x, f, fstr, work) ;
  
  if ( !eval_neighbours ) return 0 ;

  if ( src == NULL && d == NULL ) return 0 ;
  
  /* if ( t->normals == NULL && d != NULL ) { */
  /*   g_error("%s: no normals in tree but dipole strengths specified " */
  /* 	    "(d != NULL)", */
  /* 	    __FUNCTION__) ; */
  /* } */

  /*add the contribution from sources in neighbour boxes*/
  nnbr = wbfmm_box_neighbours(level, b, neighbours) ;
  /* g_assert(nnbr >= 0 && nnbr < 28) ; */

  /* if ( d == NULL ) { */
    /* monopoles only */
    for ( i = 0 ; i < nnbr ; i ++ ) {
      box = &(boxes[neighbours[i]]) ;
      /* for ( j = 0 ; j < (gint)(box->n)-4 ; j += 4 )  */
      /* 	box_curl_gradient_evaluate4(t, box->i+j, src, sstr, x, f) ;       */
      box_curl_gradient_evaluate(t, box->i, box->i+box->n,
				 src, sstr, x, f) ;      
    }

    return 0 ;
  /* } */

  /* g_assert_not_reached() ; */
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_laplace_box_field)(wbfmm_tree_t *t,
						  guint level, guint b,
						  WBFMM_REAL *src, gint sstr,
						  WBFMM_REAL *d, gint dstr,
						  guint field,
						  gboolean eval_neighbours,
						  WBFMM_REAL *x,
						  WBFMM_REAL *f, gint fstr,
						  WBFMM_REAL *work)

{
  g_assert(t->problem == WBFMM_PROBLEM_LAPLACE ) ;

  switch ( field ) {
  default:
    g_error("%s: unrecognized field type %u\n", __FUNCTION__, field) ;
    break ;
  case WBFMM_FIELD_POTENTIAL:
    box_local_field(t, level, b, x, f, fstr,
		    src, sstr, d, dstr,
		    eval_neighbours, work) ;
    break ;
  case WBFMM_FIELD_GRADIENT | WBFMM_FIELD_POTENTIAL:
  case WBFMM_FIELD_GRADIENT:
    tree_laplace_box_local_grad(t, level, b, x, f, fstr,
				src, sstr, d, dstr,
				eval_neighbours, work) ;
    break ;
  case WBFMM_FIELD_CURL | WBFMM_FIELD_POTENTIAL:
  case WBFMM_FIELD_CURL:
    tree_laplace_box_local_curl(t, level, b, x, f, fstr,
				src, sstr, d, dstr,
				eval_neighbours, work) ;    
    break ;
  case WBFMM_FIELD_CURL | WBFMM_FIELD_GRADIENT | WBFMM_FIELD_POTENTIAL:
  case WBFMM_FIELD_CURL | WBFMM_FIELD_GRADIENT:
    tree_laplace_box_local_curl_gradient(t, level, b, x, f, fstr,
					 src, sstr, d, dstr,
					 eval_neighbours, work) ;    
    break ;    
  }
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_laplace_expansion_local_eval)(WBFMM_REAL *x0,
							     WBFMM_REAL
							     *cfft,
							     gint cstr, 
							     gint N,
							     gint nq,
							     guint field,
							     WBFMM_REAL *xf,
							     WBFMM_REAL *f,
							     gint fstr,
							     WBFMM_REAL
							     *work)

{  
  switch ( field ) {
  default:
    g_error("%s: unrecognized field type %u\n", __FUNCTION__, field) ;
    break ;
  case WBFMM_FIELD_POTENTIAL:
    local_field_evaluate(x0, cfft, cstr, N, nq,	xf, f, work) ;
    break ;
  case WBFMM_FIELD_GRADIENT | WBFMM_FIELD_POTENTIAL:
  case WBFMM_FIELD_GRADIENT:
    local_grad_evaluate(x0, cfft, cstr, N, nq,	xf, f, fstr, work) ;
    break ;
  case WBFMM_FIELD_CURL | WBFMM_FIELD_POTENTIAL:
  case WBFMM_FIELD_CURL:
    local_curl_evaluate(x0, cfft, cstr, N, nq, xf, f, fstr, work) ;
    break ;
  case WBFMM_FIELD_CURL | WBFMM_FIELD_GRADIENT | WBFMM_FIELD_POTENTIAL:
  case WBFMM_FIELD_CURL | WBFMM_FIELD_GRADIENT:
    local_curl_gradient_evaluate(x0, cfft, cstr, N, nq, xf, f, fstr, work) ;
    break ;
  }
  
  return 0 ;
}

