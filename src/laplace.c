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

#include <wbfmm.h>

#include "wbfmm-private.h"

#ifdef HAVE_AVX_INSTRUCTIONS
#include <immintrin.h>
#endif /*HAVE_AVX_INSTRUCTIONS*/

gint WBFMM_FUNCTION_NAME(wbfmm_laplace_expansion_cfft)(gint N,
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
  WBFMM_REAL Cth, Sth, *Pn, *Pnm1, Cmph[64], Smph[64] ;
  gint n, m, idx, i ;

  g_assert(cstr >= nq) ;
  
  if ( N == 0 ) { return 0 ; }

  Pnm1 = &(work[0]) ; Pn = &(Pnm1[2*(2*N+1)]) ;

  WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(x0, xs, &r, &th, &ph) ;
  Cth = COS(th) ; Sth = SIN(th) ;

  WBFMM_FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth,
					   &(Pnm1[0]), &(Pn[0]), &(Pn[1])) ;
  Cmph[0] = 1.0 ; Cmph[1] = COS(ph) ; 
  Smph[0] = 0.0 ; Smph[1] = SIN(ph) ;
  
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
    /* Cmph[n] = Cmph[n-1]*Cmph[1] - Smph[n-1]*Smph[1] ; */
    /* Smph[n] = Smph[n-1]*Cmph[1] + Cmph[n-1]*Smph[1] ; */
    wbfmm_cos_sin_nph(Cmph,Smph,Cmph[1],Smph[1],n-1) ;
    
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

gint WBFMM_FUNCTION_NAME(wbfmm_laplace_expansion_dipole_cfft)(gint N,
							      WBFMM_REAL *x0,
							      WBFMM_REAL *xs,
							      WBFMM_REAL *fx,
							      WBFMM_REAL *fy,
							      WBFMM_REAL *fz,
							      gint nq,
							      WBFMM_REAL *cfft,
							      gint cstr,
							      WBFMM_REAL *work)

/*
  generate expansion coefficients for Laplace equation for dipole
  source specified as three components of vector

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
  /* WBFMM_REAL r, th, ph, rn ; */
  /* WBFMM_REAL Cth, Sth, *Pn, *Pnm1, Cmph[64], Smph[64] ; */
  /* WBFMM_REAL anm, b1, b2, Rnmm1, Rnm, Rnmp1 ; */
  /* gint n, m, idx, i ; */

  /* g_assert(nq == 1) ; */
  /* g_assert(cstr >= nq) ; */
  
  /* Pnm1 = &(work[0]) ; Pn = &(Pnm1[N+1]) ; */

  /* WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(x0, xs, &r, &th, &ph) ; */
  /* Cth = COS(th) ; Sth = SIN(th) ; */

  /* WBFMM_FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth, */
  /* 					   &(Pnm1[0]), &(Pn[0]), &(Pn[1])) ; */
  /* Cmph[0] = 1.0 ; Cmph[1] = COS(ph) ;  */
  /* Smph[0] = 0.0 ; Smph[1] = SIN(ph) ; */

  return 0 ;
}

static gint expansion_dipole_increment_cfft(gint N, gint n, gint m,
					    WBFMM_REAL rn,
					    WBFMM_REAL *cfft, gint cstr,
					    WBFMM_REAL *Pn,
					    WBFMM_REAL *Cmph,
					    WBFMM_REAL *Smph,
					    WBFMM_REAL *fx, 
					    WBFMM_REAL *fy, 
					    WBFMM_REAL *fz,
					    gint nq)

{
  WBFMM_REAL Rnm, a0, a1, a2 ;
  gint idx, idxp, idxm, i ;

  /* g_assert(n < N) ; */
  /* g_assert(m != 0) ; */
  /* g_assert(m != 1) ; */
  
  Rnm = rn*Pn[m]/(2*n+1) ;

  a0 =  SQRT((WBFMM_REAL)(2*n+1)/(2*n+3)*(n+m+1)*(n+m+2)) ;
  a1 = -SQRT((WBFMM_REAL)(2*n+1)/(2*n+3)*(n-m+1)*(n-m+2)) ;
  a2 =  SQRT(((WBFMM_REAL)(n+1)*(n+1)-m*m)*(2*n+1)/(2*n+3)) ;

  idxp = wbfmm_index_laplace_nm(n+1,m+1) ;
  idxm = wbfmm_index_laplace_nm(n+1,m-1) ;
  idx = wbfmm_index_laplace_nm(n+1,m) ;
  for ( i = 0 ; i < nq ; i ++ ) {
    /*d/dx and d/dy*/
    cfft[cstr*(idxp+0)+i] += fx[i]*0.5*a0*Rnm*Cmph[m] ;
    cfft[cstr*(idxp+1)+i] -= fx[i]*0.5*a0*Rnm*Smph[m] ;
    cfft[cstr*(idxp+0)+i] -= fy[i]*0.5*a0*Rnm*Smph[m] ;
    cfft[cstr*(idxp+1)+i] -= fy[i]*0.5*a0*Rnm*Cmph[m] ;
    cfft[cstr*(idxm+0)+i] += fx[i]*0.5*a1*Rnm*Cmph[m] ;
    cfft[cstr*(idxm+1)+i] -= fx[i]*0.5*a1*Rnm*Smph[m] ;
    cfft[cstr*(idxm+0)+i] += fy[i]*0.5*a1*Rnm*Smph[m] ;
    cfft[cstr*(idxm+1)+i] += fy[i]*0.5*a1*Rnm*Cmph[m] ;
    /*d/dz*/
    cfft[cstr*(idx+0)+i] += fz[i]*a2*Rnm*Cmph[m] ;
    cfft[cstr*(idx+1)+i] -= fz[i]*a2*Rnm*Smph[m] ;
  }
  
  return 0 ;
}

static gint expansion_dipole_increment_cfft_m0(gint N, gint n,
					       WBFMM_REAL rn,
					       WBFMM_REAL *cfft, gint cstr,
					       WBFMM_REAL *Pn,
					       WBFMM_REAL *Cmph,
					       WBFMM_REAL *Smph,
					       WBFMM_REAL *fx, 
					       WBFMM_REAL *fy, 
					       WBFMM_REAL *fz,
					       gint nq)

{
  WBFMM_REAL Rnm, a0, a2 ;
  gint idx, idxp, i, m = 0 ;
  
  /* g_assert(n < N) ; */
  /* g_assert(m == 0) ; */
  
  Rnm = rn*Pn[m]/(2*n+1) ;

  a0 =  SQRT((WBFMM_REAL)(2*n+1)/(2*n+3)*(n+m+1)*(n+m+2)) ;
  a2 =  SQRT(((WBFMM_REAL)(n+1)*(n+1)-m*m)*(2*n+1)/(2*n+3)) ;
  /*d/dx and d/dy*/
  /*only fill in the m==1 case, leave m==-1 alone (picked up by the
    conjugate)*/
  idxp = wbfmm_index_laplace_nm(n+1,m+1) ;
  idx  = (n+1)*(n+1) ;
  for ( i = 0 ; i < nq ; i ++ ) {
    cfft[cstr*(idxp+0)+i] += fx[i]*0.5*a0*Rnm*Cmph[m] ;
    cfft[cstr*(idxp+1)+i] -= fx[i]*0.5*a0*Rnm*Smph[m] ;
    cfft[cstr*(idxp+0)+i] -= fy[i]*0.5*a0*Rnm*Smph[m] ;
    cfft[cstr*(idxp+1)+i] -= fy[i]*0.5*a0*Rnm*Cmph[m] ;
    cfft[cstr*idx     +i] += fz[i]*a2*Rnm ;
  }

  return 0 ;
}

static gint expansion_dipole_increment_cfft_m1(gint N, gint n,
					       WBFMM_REAL rn,
					       WBFMM_REAL *cfft, gint cstr,
					       WBFMM_REAL *Pn,
					       WBFMM_REAL *Cmph,
					       WBFMM_REAL *Smph,
					       WBFMM_REAL *fx, 
					       WBFMM_REAL *fy, 
					       WBFMM_REAL *fz,
					       gint nq)

{
  WBFMM_REAL Rnm, a0, a1, a2 ;
  gint idx, idx0, idxp, i, m = 1 ;
  
  /* if ( n >= N ) return 0 ; */
  /* g_assert(m == 1) ; */
  
  Rnm = rn*Pn[m]/(2*n+1) ;

  a0 =  SQRT((WBFMM_REAL)(2*n+1)/(2*n+3)*(n+m+1)*(n+m+2)) ;
  a1 = -SQRT((WBFMM_REAL)(2*n+1)/(2*n+3)*(n-m+1)*(n-m+2)) ;
  a2 =  SQRT(((WBFMM_REAL)(n+1)*(n+1)-m*m)*(2*n+1)/(2*n+3)) ;
  idxp = wbfmm_index_laplace_nm(n+1,m+1) ;
  idx0 = (n+1)*(n+1) ;
  idx = wbfmm_index_laplace_nm(n+1,m) ;
  for ( i = 0 ; i < nq ; i ++ ) {
    /*d/dx and d/dy*/
    cfft[cstr*(idxp+0)+i] += fx[i]*0.5*a0*Rnm*Cmph[m] ;
    cfft[cstr*(idxp+1)+i] -= fx[i]*0.5*a0*Rnm*Smph[m] ;
    cfft[cstr*(idxp+0)+i] -= fy[i]*0.5*a0*Rnm*Smph[m] ;
    cfft[cstr*(idxp+1)+i] -= fy[i]*0.5*a0*Rnm*Cmph[m] ;
    /*m=0 coefficient is real*/
    cfft[cstr*(idx0+0)+i] += fx[i]*a1*Rnm*Cmph[m] ;
    cfft[cstr*(idx0+0)+i] += fy[i]*a1*Rnm*Smph[m] ;
    /*d/dz*/
    cfft[cstr*(idx+0)+i] += fz[i]*a2*Rnm*Cmph[m] ;
    cfft[cstr*(idx+1)+i] -= fz[i]*a2*Rnm*Smph[m] ;
  }
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_laplace_expansion_normal_cfft)(gint N,
							      WBFMM_REAL *x0,
							      WBFMM_REAL *xs,
							      WBFMM_REAL
							      *normal,
							      WBFMM_REAL *q,
							      gint nq,
							      WBFMM_REAL *cfft,
							      gint cstr,
							      WBFMM_REAL *work)

/*
  generate expansion coefficients for Laplace equation for dipole
  source specified as normal and strength

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
  WBFMM_REAL Cth, Sth, *Pn, *Pnm1, Cmph[64], Smph[64] ;
  WBFMM_REAL fx[8], fy[8], fz[8] ;
  gint n, m, i ;

  g_assert(nq < 9) ;
  g_assert(cstr >= nq) ;
  
  if ( N == 0 ) { return 0 ; }
  
  Pnm1 = &(work[0]) ; Pn = &(Pnm1[N+1]) ;

  WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(x0, xs, &r, &th, &ph) ;
  Cth = COS(th) ; Sth = SIN(th) ;

  WBFMM_FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth,
					   &(Pnm1[0]), &(Pn[0]), &(Pn[1])) ;
  Cmph[0] = 1.0 ; Cmph[1] = COS(ph) ; 
  Smph[0] = 0.0 ; Smph[1] = SIN(ph) ;

  for ( i = 0 ; i < nq ; i ++ ) {
    fx[i] = normal[0]*q[i] ; fy[i] = normal[1]*q[i] ; fz[i] = normal[2]*q[i] ;
  }
  
  n = 0 ; 
  m = 0 ;
  rn = 1.0 ;
  expansion_dipole_increment_cfft_m0(N, n, rn, cfft, cstr, Pnm1, Cmph, Smph,
				     fx, fy, fz, nq) ;
  
  if ( N == 1 ) return 0 ;
  /* fprintf(stderr, "Hello\n") ; */
  
  n = 1 ; 
  m = 0 ;
  rn *= r ;
  /* Cmph[n+1] = Cmph[n]*Cmph[1] - Smph[n]*Smph[1] ; */
  /* Smph[n+1] = Smph[n]*Cmph[1] + Cmph[n]*Smph[1] ; */
  wbfmm_cos_sin_nph(Cmph,Smph,Cmph[1],Smph[1],n) ;
    
  expansion_dipole_increment_cfft_m0(N, n, rn, cfft, cstr, Pn, Cmph, Smph,
				     fx, fy, fz, nq) ;
  
  m = 1 ; 

  expansion_dipole_increment_cfft_m1(N, n, rn, cfft, cstr, Pn, Cmph, Smph,
				     fx, fy, fz, nq) ;

  if ( N == 2 ) return 0 ;

  for ( n = 2 ; n < N ; n ++ ) {
    rn *= r ;
    WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pnm1, &Pn,
  							n-1, Cth, Sth) ;
    /* Cmph[n+1] = Cmph[n]*Cmph[1] - Smph[n]*Smph[1] ; */
    /* Smph[n+1] = Smph[n]*Cmph[1] + Cmph[n]*Smph[1] ; */
    wbfmm_cos_sin_nph(Cmph,Smph,Cmph[1],Smph[1],n) ;
    
    m = 0 ;
    expansion_dipole_increment_cfft_m0(N, n, rn, cfft, cstr, Pn, Cmph, Smph,
				       fx, fy, fz, nq) ;
    m = 1  ;
    expansion_dipole_increment_cfft_m1(N, n, rn, cfft, cstr, Pn, Cmph, Smph,
				       fx, fy, fz, nq) ;
    
    for ( m = 2 ; m <= n ; m ++ ) {
      expansion_dipole_increment_cfft(N, n, m, rn, cfft, cstr, Pn, Cmph, Smph,
				      fx, fy, fz, nq) ;
    }
  }

  return 0 ;
}


gint WBFMM_FUNCTION_NAME(wbfmm_laplace_expansion_evaluate)(WBFMM_REAL *x0,
							   WBFMM_REAL *cfft,
							   gint cstr, 
							   gint N, gint nq,
							   WBFMM_REAL *xf,
							   WBFMM_REAL *field,
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

gint WBFMM_FUNCTION_NAME(wbfmm_laplace_coaxial_translate_SS)(WBFMM_REAL *Co,
							     gint cstro,
							     gint No,
							     WBFMM_REAL *Ci,
							     gint cstri,
							     gint Ni,
							     gint nq,
							     WBFMM_REAL t,
							     WBFMM_REAL sc)

{
  gint n, m, nd, idxo, idxi, i ;
  WBFMM_REAL c, tn[64] ;

  tn[0] = 1.0 ;
  m = 0 ; 
  for ( idxo = n = 0 ; n <= No ; (n ++), (idxo = n*n) ) {
    for ( i = 0 ; i < nq ; i ++ ) Co[cstro*idxo+i] *= sc ;
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
      for ( i = 0 ; i < nq ; i ++ ) {
	Co[cstro*(idxo+0)+i] *= sc ; Co[cstro*(idxo+1)+i] *= sc ;
      }
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

  c =
    lgamma(nd - m + 1) + lgamma(nd + m + 1) -
    lgamma(n  - m + 1) - lgamma(n  + m + 1) -
    2.0*lgamma(nd - n + 1) ;
  c += log((2.0*nd+1)/(2.0*n+1)) ;

  c = exp(0.5*c) ;
  
  return c ;
  
  c  = wbfmm_factorial(nd  - m)*wbfmm_factorial(nd  + m) ;
  c /= wbfmm_factorial(n - m)*wbfmm_factorial(n + m) ;
  c *= (2.0*nd+1)/(2.0*n+1) ;
  c  = SQRT(c) ;
  c /= wbfmm_factorial(nd-n) ;

  /*this needs to be multiplied by (t)^{nd-n} for a translation of
    distance t*/
  
  return c ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_laplace_coaxial_translate_RR)(WBFMM_REAL *Co,
							     gint cstro,
							     gint No,
							     WBFMM_REAL *Ci,
							     gint cstri,
							     gint Ni,
							     gint nq,
							     WBFMM_REAL t,
							     WBFMM_REAL sc)
  
{
  gint n, m, nd, idxo, idxi, i ;
  WBFMM_REAL c, tn[64] ;

  for ( idxi = nd = 0 ; nd <= Ni ; (nd ++), (idxi = nd*nd) ) {
    for ( idxo = n = 0 ; n <= MIN(nd,Ni) ; (n ++), (idxo = n*n) ) {
      for ( i = 0 ; i < nq ; i ++ ) Co[cstro*idxo+i] *= sc ;
    }
  }
  for ( m = 1 ; m <= No ; m ++ ) {
    for ( nd = m ; nd <= No ; nd ++ ) {
      for ( n = m ; n <= MIN(nd,Ni) ; n ++ ) {
  	idxo = wbfmm_index_laplace_nm(n,m) ;
	for ( i = 0 ; i < nq ; i ++ ) {
	  Co[cstro*(idxo+0)+i] *= sc ; Co[cstro*(idxo+1)+i] *= sc ;
	}
      }
    }
  }

  tn[0] = 1.0 ;
  m  = 0 ; 
  for ( idxi = nd = 0 ; nd <= Ni ; (nd ++), (idxi = nd*nd) ) {
    for ( idxo = n = 0 ; n <= MIN(nd,Ni) ; (n ++), (idxo = n*n) ) {
      /* for ( i = 0 ; i < nq ; i ++ ) Co[cstro*idxo+i] *= sc ; */
      c = wbfmm_coaxial_translation_RR_cfft(n, nd, m)*tn[nd-n] ;
      for ( i = 0 ; i < nq ; i ++ ) Co[cstro*idxo+i] += c*Ci[cstri*idxi+i] ;
    }
    tn[nd+1] = tn[nd]*t ;
  }
  
  for ( m = 1 ; m <= No ; m ++ ) {
    for ( nd = m ; nd <= No ; nd ++ ) {
      idxi = wbfmm_index_laplace_nm(nd,m) ;
      for ( n = m ; n <= MIN(nd,Ni) ; n ++ ) {
  	idxo = wbfmm_index_laplace_nm(n,m) ;
	/* for ( i = 0 ; i < nq ; i ++ ) Co[cstro*idxo+i] *= sc ; */
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

  c =
    2.0*lgamma(n + nd +1) -
    lgamma(n  - m + 1) - lgamma(n  + m + 1) -
    lgamma(nd - m + 1) - lgamma(nd + m + 1) ;
  c += log((2.0*nd+1)/(2.0*n+1)) ;
  c = exp(0.5*c) ;
  c *= pow(-1.0, n+m) ;

  return c ;

  /*this needs to be divided by (t)^{n+nd+1} for a translation of
    distance t*/
  
  return c ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_laplace_coaxial_translate_SR)(WBFMM_REAL *Co,
							     gint cstro,
							     gint No,
							     WBFMM_REAL *Ci,
							     gint cstri,
							     gint Ni,
							     gint nq,
							     WBFMM_REAL t,
							     WBFMM_REAL sc)
  
{
  gint n, m, nd, idxo, idxi, i ;
  WBFMM_REAL c, tn[128] ;

  g_assert(No + Ni < 126) ;
  tn[0] = ( t < 0 ? -1.0 : 1.0) ;
  for ( n = 1 ; n <= No+Ni+2 ; n ++ ) tn[n] = t*tn[n-1] ;

  m = 0 ;
  for ( idxo = n = 0 ; n <= No ; (n ++), (idxo = n*n) ) {
    for ( i = 0 ; i < nq ; i ++ ) Co[cstro*idxo+i] *= sc ;
    for ( idxi = nd = 0 ; nd <= Ni ; (nd ++), (idxi = nd*nd) ) {
      c = wbfmm_coaxial_translation_SR_cfft(n, nd, m)/tn[n+nd+1] ;
      for ( i = 0 ; i < nq ; i ++ ) Co[cstro*idxo+i] += c*Ci[cstri*idxi+i] ;
    }
  }

  for ( m = 1 ; m <= No ; m ++ ) {
    for ( n = m ; n <= No ; n ++ ) {
      idxo = wbfmm_index_laplace_nm(n,m) ;
      for ( i = 0 ; i < nq ; i ++ ) {
	Co[cstro*(idxo+0)+i] *= sc ; Co[cstro*(idxo+1)+i] *= sc ;
      }
      for ( nd = m ; nd <= Ni ; nd ++ ) {
	idxi = wbfmm_index_laplace_nm(nd,m) ;
  	c = wbfmm_coaxial_translation_SR_cfft(n, nd, m)/tn[n+nd+1] ;
  	for ( i = 0 ; i < nq ; i ++ ) {
  	  Co[cstro*(idxo+0)+i] += c*Ci[cstri*(idxi+0)+i] ;
  	  Co[cstro*(idxo+1)+i] += c*Ci[cstri*(idxi+1)+i] ;
  	}	
      }
    }
  }
  
  return 0 ;
}

static WBFMM_REAL coaxial_translation_SS_cfft(gint n, gint nd, gint m)

{
  WBFMM_REAL c ;

  g_assert(n >= nd) ;

  c =
    lgamma(n  - m + 1) + lgamma(n  + m + 1) -
    lgamma(nd - m + 1) - lgamma(nd + m + 1) -
    2.0*lgamma(n - nd + 1) ;

  c += log((2.0*nd+1)/(2.0*n+1)) ;
  c = exp(0.5*c) ;
  
  return c ;
  
  c  = wbfmm_factorial(n  - m)*wbfmm_factorial(n  + m) ;
  c /= wbfmm_factorial(nd - m)*wbfmm_factorial(nd + m) ;
  c *= (2.0*nd+1)/(2.0*n+1) ;
  c  = SQRT(c) ;
  c /= wbfmm_factorial(n-nd) ;

  /*this needs to be multiplied by (-t)^{n-nd} for a translation of
    distance t*/
  
  return c ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_laplace_coaxial_translate_init)(gint N)

{
  gint n, m, nu, i, Nmax, ne ;

  /* fprintf(stderr, "N = %d\n", N) ; */
  
  Nmax = _wbfmm_translation_Nmax = N + 1 ;

  ne = _wbfmm_SS_coefficient_index_nmnu(Nmax+1,0,0) ;

  WBFMM_FUNCTION_NAME(_wbfmm_SS_coefficients_laplace) =
    (WBFMM_REAL *)g_malloc0(ne*sizeof(WBFMM_REAL)) ;

  for ( n = 0 ; n <= Nmax ; n ++ ) {
    for ( m = 0 ; m <= n ; m ++ ) {
      for ( nu = m ; nu <= n ; nu ++ ) {
	i = _wbfmm_SS_coefficient_index_nmnu(n, m, nu) ;
	g_assert(i < ne) ;
  	g_assert(WBFMM_FUNCTION_NAME(_wbfmm_SS_coefficients_laplace)[i]
		 == 0.0) ;
	WBFMM_FUNCTION_NAME(_wbfmm_SS_coefficients_laplace)[i] =
	  coaxial_translation_SS_cfft(n, nu, m) ;
      }
    }
  }

  ne = (Nmax+2)*(Nmax+1)*(Nmax+1)/2 ;
  
  WBFMM_FUNCTION_NAME(_wbfmm_RR_coefficients_laplace) =
    (WBFMM_REAL *)g_malloc0(ne*sizeof(WBFMM_REAL)) ;

  for ( m = 0 ; m <= Nmax ; m ++ ) {
    for ( nu = m ; nu <= Nmax ; nu ++ ) {
      for ( n = m ; n <= nu ; n ++ ) {
  	i = _wbfmm_RR_coefficient_index_nmnu(n, m, nu) ;
	g_assert(i < ne) ;
  	g_assert(WBFMM_FUNCTION_NAME(_wbfmm_RR_coefficients_laplace)[i]
		 == 0.0) ;
  	WBFMM_FUNCTION_NAME(_wbfmm_RR_coefficients_laplace)[i] =
  	  coaxial_translation_RR_cfft(n, nu, m) ;
      }
    }
  }

  ne = Nmax*(6*(Nmax+1)*(Nmax+2) + 1 - 3*(Nmax+1)*(2*Nmax+3) +
  	     2*(Nmax+1)*(Nmax+1)) ;
  
  WBFMM_FUNCTION_NAME(_wbfmm_SR_coefficients_laplace) =
    (WBFMM_REAL *)g_malloc0(ne*sizeof(WBFMM_REAL)) ;
  
  for ( m = 0 ; m <= Nmax ; m ++ ) {
    for ( n = m ; n <= Nmax ; n ++ ) {
      for ( nu = m ; nu <= Nmax ; nu ++ ) {
  	i = _wbfmm_SR_coefficient_index_nmnu(n, m, nu) ;
	g_assert(i < ne) ;
  	g_assert(WBFMM_FUNCTION_NAME(_wbfmm_SR_coefficients_laplace)[i]
		 == 0.0) ;
  	WBFMM_FUNCTION_NAME(_wbfmm_SR_coefficients_laplace)[i] =
  	  coaxial_translation_SR_cfft(n, nu, m) ;
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
    c = (WBFMM_REAL *)(t->mpr[l]) ;
    /*
      set box pointers to start of their coefficients, noting that
      coefficients are packed in groups of eight for shift operations
    */
    boxes = t->boxes[l] ;
    for ( i = 0 ; i < nb ; i += 8 ) {
      for ( j = 0 ; j < 8 ; j ++ ) {
	boxes[i+j].mpr = &(c[i*nq*nc+nq*j]) ;
      }
    }
  }
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_tree_laplace_leaf_expansions)(wbfmm_tree_t *t,
							     WBFMM_REAL *src,
							     gint sstr,
							     WBFMM_REAL
							     *dipoles,
							     gint dstr,
							     gboolean
							     zero_expansions,
							     WBFMM_REAL *work)

{
  guint32 nb, nc, i, j, ns, d, idx, nq ;
  guint64 im ;
  wbfmm_box_t *boxes ;
  WBFMM_REAL *xs, *q, xb[3], wb ;

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

  if ( src == NULL && dipoles == NULL ) return 0 ;
  if ( t->normals == NULL && dipoles != NULL ) {
    g_error("%s: no normals in tree but dipole strengths specified "
	    "(dipoles != NULL)",
	    __FUNCTION__) ;
  }

  /* if ( normals == NULL && dipoles == NULL ) { */
  if ( dipoles == NULL ) {
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
	WBFMM_FUNCTION_NAME(wbfmm_laplace_expansion_cfft)(ns, xb, xs, q, nq,
							  boxes[i].mps, 8*nq,
							  work) ;
      }
    }

    return 0 ;
  }

  if ( src == NULL && dipoles != NULL ) {
    /* dipoles only */
    WBFMM_REAL *normal ;
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
	/* normal = &(normals[nstr*idx]) ; */
	normal = wbfmm_tree_normal_index(t,idx) ;
	q = &(dipoles[idx*dstr]) ;
	WBFMM_FUNCTION_NAME(wbfmm_laplace_expansion_normal_cfft)(ns, xb, xs,
								 normal, q,
								 nq,
								 boxes[i].mps,
								 8*nq,
								 work) ;
      }
    }

    return 0 ;
  }  

  if ( src != NULL && dipoles != NULL ) {
    /* sources and dipoles */
    WBFMM_REAL *normal ;
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
	/* normal = &(normals[nstr*idx]) ; */
	normal = wbfmm_tree_normal_index(t,idx) ;
	q = &(dipoles[idx*dstr]) ;
	WBFMM_FUNCTION_NAME(wbfmm_laplace_expansion_normal_cfft)(ns, xb, xs,
								 normal, q,
								 nq,
								 boxes[i].mps,
								 8*nq,
								 work) ;
	/*this needs a combined expansion generation*/
	q = &(src[idx*sstr]) ;
	WBFMM_FUNCTION_NAME(wbfmm_laplace_expansion_cfft)(ns, xb, xs, q, nq,
							  boxes[i].mps, 8*nq,
							  work) ;
      }
    }

    return 0 ;
  }  
  
  g_assert_not_reached() ;
  
  return 0 ;
}

static gint _wbfmm_laplace_local_coefficients_scalar(WBFMM_REAL *cfft,
						     gint N,
						     WBFMM_REAL r,
						     WBFMM_REAL Cth,
						     WBFMM_REAL Sth,
						     WBFMM_REAL Cph,
						     WBFMM_REAL Sph,
						     WBFMM_REAL *work)

{
  gint n, m, idx ;
  WBFMM_REAL *Pn, *Pnm1 ;
  WBFMM_REAL *Cmph, *Smph ;
  WBFMM_REAL rn ;

  Pnm1 = &(work[0]) ; Pn = &(Pnm1[N+1]) ;
  Cmph = &(Pn[N+1]) ; Smph = &(Cmph[N+1]) ;

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
  cfft[idx] = rn*Pnm1[m] ;

  n = 1 ; 
  m = 0 ; 
  rn *= r ;
  idx = n*n ;
  cfft[idx] = rn*Pn[m] ;

  m = 1 ; 
  idx = wbfmm_index_laplace_nm(n,m) ;
  cfft[idx+0] = 2.0*Pn[m]*rn*Cmph[m] ;
  cfft[idx+1] = 2.0*Pn[m]*rn*Smph[m] ;

  for ( n = 2 ; n <= N ; n ++ ) {
    rn *= r ;
    WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pnm1, &Pn,
							n-1, Cth, Sth) ;
    /* Cmph[n] = Cmph[n-1]*Cph - Smph[n-1]*Sph ; */
    /* Smph[n] = Smph[n-1]*Cph + Cmph[n-1]*Sph ; */
    wbfmm_cos_sin_nph(Cmph,Smph,Cmph[1],Smph[1],n-1) ;

    m = 0 ; 
    idx = n*n ;
    cfft[idx] = rn*Pn[0] ;

    for ( m = 1 ; m <= n ; m ++ ) {
      idx = wbfmm_index_laplace_nm(n,m) ;
      cfft[idx+0] = 2.0*Pn[m]*rn*Cmph[m] ;
      cfft[idx+1] = 2.0*Pn[m]*rn*Smph[m] ;
    }
  }
  
  return 0 ;
}

static gint _wbfmm_laplace_local_coefficients_gradient(WBFMM_REAL *cfft,
						       gint N,
						       WBFMM_REAL r,
						       WBFMM_REAL Cth,
						       WBFMM_REAL Sth,
						       WBFMM_REAL Cph,
						       WBFMM_REAL Sph,
						       WBFMM_REAL *work)

{
  gint n, m, idx ;
  WBFMM_REAL *Pn, *Pnm1, *Cmph, *Smph, rn, anm, b1, b2, Rnmm1, Rnm, Rnmp1 ;

  Pnm1 = &(work[0]) ; Pn = &(Pnm1[N+1]) ;
  Cmph = &(Pn[N+1]) ; Smph = &(Cmph[N+1]) ;

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
  cfft[3*idx+0] = cfft[3*idx+1] = cfft[3*idx+2] = 0.0 ;

  n = 1 ; 
  m = 0 ;
  rn = 1.0 ;
  idx = n*n ;
  /* Cmph[n+1] = Cmph[n]*Cmph[1] - Smph[n]*Smph[1] ; */
  /* Smph[n+1] = Smph[n]*Cmph[1] + Cmph[n]*Smph[1] ; */
  wbfmm_cos_sin_nph(Cmph,Smph,Cmph[1],Smph[1],n) ;

  anm = SQRT((WBFMM_REAL)(2*n+1)/(2*n-1)*(n-m)*(n+m)) ;
  Rnm = rn*Pnm1[m]*anm ;
  cfft[3*idx+0] = cfft[3*idx+1] = 0.0 ;
  cfft[3*idx+2] = Rnm ;

  m = 1 ; 
  idx = wbfmm_index_laplace_nm(n,m) ;
  anm = SQRT((WBFMM_REAL)(2*n+1)/(2*n-1)*(n-m)*(n+m)) ;
  b1  = SQRT((WBFMM_REAL)(2*n+1)/(2*n-1)*(n-m)*(n-m-1)) ;
  b2  = SQRT((WBFMM_REAL)(2*n+1)/(2*n-1)*(n+m)*(n+m-1)) ;
  Rnmm1 = rn*Pnm1[m-1]*b2 ;
  Rnm   = rn*Pnm1[m+0]*anm*2.0 ;
  Rnmp1 = rn*Pnm1[m+1]*b1 ;
  cfft[3*idx+0] =  Rnmm1*Cmph[m-1] - Rnmp1*Cmph[m+1] ;
  cfft[3*idx+1] =  Rnmp1*Smph[m+1] - Rnmm1*Smph[m-1] ;
  cfft[3*idx+2] = -Rnmp1*Smph[m+1] - Rnmm1*Smph[m-1] ;
  cfft[3*idx+3] = -Rnmm1*Cmph[m-1] - Rnmp1*Cmph[m+1] ;
  cfft[3*idx+4] =  Rnm  *Cmph[m+0] ;
  cfft[3*idx+5] = -Rnm  *Smph[m+0] ;
  
  for ( n = 2 ; n <= N ; n ++ ) {
    rn *= r ;
    WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pnm1, &Pn,
							n-1, Cth, Sth) ;
    /* Cmph[n+1] = Cmph[n]*Cmph[1] - Smph[n]*Smph[1] ; */
    /* Smph[n+1] = Smph[n]*Cmph[1] + Cmph[n]*Smph[1] ; */
    wbfmm_cos_sin_nph(Cmph,Smph,Cmph[1],Smph[1],n) ;

    m = 0 ; 
    idx = n*n ;
    anm = SQRT((WBFMM_REAL)(2*n+1)/(2*n-1)*(n-m)*(n+m)) ;
    b1  = SQRT((WBFMM_REAL)(2*n+1)/(2*n-1)*(n-m)*(n-m-1)) ;
    b2  = SQRT((WBFMM_REAL)(2*n+1)/(2*n-1)*(n+m)*(n+m-1)) ;
    Rnmm1 = rn*Pnm1[m+1]*b2 ;
    Rnm   = rn*Pnm1[m]*anm ;
    Rnmp1 = rn*Pnm1[m+1]*b1 ;

    cfft[3*idx+0] = -Rnmp1*Cmph[m+1] ;
    cfft[3*idx+1] = -Rnmp1*Smph[m+1] ;
    cfft[3*idx+2] =  Rnm ;

    for ( m = 1 ; m <= n ; m ++ ) {
      idx = wbfmm_index_laplace_nm(n,m) ;
      anm = SQRT((WBFMM_REAL)(2*n+1)/(2*n-1)*(n-m)*(n+m)) ;
      b1  = SQRT((WBFMM_REAL)(2*n+1)/(2*n-1)*(n-m)*(n-m-1)) ;
      b2  = SQRT((WBFMM_REAL)(2*n+1)/(2*n-1)*(n+m)*(n+m-1)) ;
      Rnmm1 = rn*Pnm1[m-1]*b2 ;
      Rnm   = rn*Pnm1[m+0]*anm*2.0 ;
      Rnmp1 = rn*Pnm1[m+1]*b1 ;
      cfft[3*idx+0] =  Rnmm1*Cmph[m-1] - Rnmp1*Cmph[m+1] ;
      cfft[3*idx+1] =  Rnmp1*Smph[m+1] - Rnmm1*Smph[m-1] ;
      cfft[3*idx+2] = -Rnmp1*Smph[m+1] - Rnmm1*Smph[m-1] ;
      cfft[3*idx+3] = -Rnmm1*Cmph[m-1] - Rnmp1*Cmph[m+1] ;
      cfft[3*idx+4] =  Rnm  *Cmph[m+0] ;
      cfft[3*idx+5] = -Rnm  *Smph[m+0] ;
    }
  }
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_laplace_local_coefficients)(WBFMM_REAL *x,
							   gint N,
							   guint field,
							   WBFMM_REAL *cfft,
							   WBFMM_REAL *work)

{
  WBFMM_REAL r, th, ph, x0[3] = {0.0} ;
  WBFMM_REAL Cth, Sth, Cph, Sph ;
  
  WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(x0, x, &r, &th, &ph) ;
  Cth = COS(th) ; Sth = SIN(th) ; 
  Cph = COS(ph) ; Sph = SIN(ph) ; 

  switch ( field ) {
  default:
    g_error("%s: unrecognized field definition (%u)", __FUNCTION__, field) ;
    break ;
  case WBFMM_FIELD_POTENTIAL:
    return _wbfmm_laplace_local_coefficients_scalar(cfft, N, r, Cth, Sth,
						    Cph, Sph, work) ;
    break ;
  case WBFMM_FIELD_GRADIENT:
    return _wbfmm_laplace_local_coefficients_gradient(cfft, N, r, Cth, Sth,
						      Cph, Sph, work) ;
    break ;
  }

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_laplace_field_coefficients)(WBFMM_REAL *x,
							   gint N,
							   guint field,
							   WBFMM_REAL *cfft,
							   WBFMM_REAL *work)

{
  WBFMM_REAL r, th, ph, rn, x0[3] = {0.0} ;
  WBFMM_REAL Cth, Sth, *Pn, *Pnm1, Cph, Sph ;
  WBFMM_REAL *Cmph, *Smph ;
  gint n, m, idx ;
  
  Pnm1 = &(work[0]) ; Pn = &(Pnm1[N+1]) ;
  Cmph = &(Pn[N+1]) ; Smph = &(Cmph[N+1]) ;

  WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(x0, x, &r, &th, &ph) ;
  Cth = COS(th) ; Sth = SIN(th) ; 
  Cph = COS(ph) ; Sph = SIN(ph) ; 

  /*initialize recursions*/
  WBFMM_FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth,
					   &(Pnm1[0]), &(Pn[0]), &(Pn[1])) ;
  Cmph[0] = 1.0 ; Cmph[1] = Cph ; 
  Smph[0] = 0.0 ; Smph[1] = Sph ;

  /*first two terms by hand*/
  n = 0 ; 
  m = 0 ;
  rn = 1.0/r ;
  idx = n*n ;
  cfft[idx] = rn*Pnm1[m] ;

  n = 1 ; 
  m = 0 ; 
  rn /= r ;
  idx = n*n ;
  cfft[idx] = rn*Pn[m] ;

  m = 1 ; 
  idx = wbfmm_index_laplace_nm(n,m) ;
  cfft[idx+0] = 2.0*Pn[m]*rn*Cmph[m] ;
  cfft[idx+1] = 2.0*Pn[m]*rn*Smph[m] ;

  for ( n = 2 ; n <= N ; n ++ ) {
    rn /= r ;
    WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pnm1, &Pn,
							n-1, Cth, Sth) ;
    /* Cmph[n] = Cmph[n-1]*Cph - Smph[n-1]*Sph ; */
    /* Smph[n] = Smph[n-1]*Cph + Cmph[n-1]*Sph ; */
    wbfmm_cos_sin_nph(Cmph,Smph,Cmph[1],Smph[1],n-1) ;

    m = 0 ; 
    idx = n*n ;
    cfft[idx] = rn*Pn[0] ;

    for ( m = 1 ; m <= n ; m ++ ) {
      idx = wbfmm_index_laplace_nm(n,m) ;
      cfft[idx+0] = 2.0*Pn[m]*rn*Cmph[m] ;
      cfft[idx+1] = 2.0*Pn[m]*rn*Smph[m] ;
    }
  }
  
  return 0 ;
}

static gint _wbfmm_laplace_expansion_apply_scalar(WBFMM_REAL *C,
						  gint cstr,
						  gint nq,
						  WBFMM_REAL *ec,
						  gint N,
						  WBFMM_REAL *f, gint fstr)
{
  gint n, m, idx, i ;
  
  /*first two terms by hand*/
  n = 0 ; 
  m = 0 ;
  idx = n*n ;
  for ( i = 0 ; i < nq ; i ++ ) f[i] += C[cstr*idx+i]*ec[idx] ;

  n = 1 ; 
  m = 0 ; 
  idx = n*n ;
  for ( i = 0 ; i < nq ; i ++ ) f[i] += C[cstr*idx+i]*ec[idx] ;

  m = 1 ; 
  idx = wbfmm_index_laplace_nm(n,m) ;
  for ( i = 0 ; i < nq ; i ++ ) {
    f[i] += C[cstr*(idx+0)+i]*ec[idx+0] - C[cstr*(idx+1)+i]*ec[idx+1] ;
  }

  for ( n = 2 ; n <= N ; n ++ ) {
    m = 0 ; 
    idx = n*n ;
    for ( i = 0 ; i < nq ; i ++ ) f[i] += C[cstr*idx+i]*ec[idx] ;
    
    for ( m = 1 ; m <= n ; m ++ ) {
      idx = wbfmm_index_laplace_nm(n,m) ;
      for ( i = 0 ; i < nq ; i ++ ) {
	f[i] += C[cstr*(idx+0)+i]*ec[idx+0] - C[cstr*(idx+1)+i]*ec[idx+1] ;
      }
    }
  }
  
  return 0 ;
}

static gint _wbfmm_laplace_expansion_apply_gradient(WBFMM_REAL *C,
						    gint cstr,
						    gint nq,
						    WBFMM_REAL *ec,
						    gint N,
						    WBFMM_REAL *f,
						    gint fstr)
{
  gint n, m, idx, i ;
  WBFMM_REAL cr, ci ;
  
  /*first two terms by hand*/
  n = 0 ; 
  m = 0 ;
  idx = n*n ;
  for ( i = 0 ; i < nq ; i ++ ) {
    f[fstr*i+0] += C[cstr*idx+i]*ec[3*idx+0] ;
    f[fstr*i+1] += C[cstr*idx+i]*ec[3*idx+1] ;
    f[fstr*i+2] += C[cstr*idx+i]*ec[3*idx+2] ;
  }
  
  n = 1 ; 
  m = 0 ; 
  idx = n*n ;
  for ( i = 0 ; i < nq ; i ++ ) {
    f[fstr*i+0] += C[cstr*idx+i]*ec[3*idx+0] ;
    f[fstr*i+1] += C[cstr*idx+i]*ec[3*idx+1] ;
    f[fstr*i+2] += C[cstr*idx+i]*ec[3*idx+2] ;
  }

  m = 1 ; 
  idx = wbfmm_index_laplace_nm(n,m) ;
  for ( i = 0 ; i < nq ; i ++ ) {
    cr = C[cstr*(idx+0)+i] ; ci = C[cstr*(idx+1)+i] ;
    f[fstr*i+0] += cr*ec[3*idx+0] + ci*ec[3*idx+1] ;
    f[fstr*i+1] += cr*ec[3*idx+2] + ci*ec[3*idx+3] ;
    f[fstr*i+2] += cr*ec[3*idx+4] + ci*ec[3*idx+5] ;
  }
  
  for ( n = 2 ; n <= N ; n ++ ) {
    m = 0 ; 
    idx = n*n ;
    for ( i = 0 ; i < nq ; i ++ ) {
      f[fstr*i+0] += C[cstr*idx+i]*ec[3*idx+0] ;
      f[fstr*i+1] += C[cstr*idx+i]*ec[3*idx+1] ;
      f[fstr*i+2] += C[cstr*idx+i]*ec[3*idx+2] ;
    }

    for ( m = 1 ; m <= n ; m ++ ) {
      idx = wbfmm_index_laplace_nm(n,m) ;
      for ( i = 0 ; i < nq ; i ++ ) {
	cr = C[cstr*(idx+0)+i] ; ci = C[cstr*(idx+1)+i] ;
	f[fstr*i+0] += cr*ec[3*idx+0] + ci*ec[3*idx+1] ;
	f[fstr*i+1] += cr*ec[3*idx+2] + ci*ec[3*idx+3] ;
	f[fstr*i+2] += cr*ec[3*idx+4] + ci*ec[3*idx+5] ;
      }  
    }
  }
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_laplace_expansion_apply)(WBFMM_REAL *C,
							gint cstr,
							gint nq,
							WBFMM_REAL *ec,
							gint N,
							guint field,
							WBFMM_REAL *f,
							gint fstr)
  
{
  switch ( field ) {
  default:
    g_error("%s: unrecognized field definition (%u)", __FUNCTION__, field) ;
    break ;
  case WBFMM_FIELD_POTENTIAL:
    return _wbfmm_laplace_expansion_apply_scalar(C, cstr, nq, ec, N, f, fstr) ;
    break ;
  case WBFMM_FIELD_GRADIENT:
    return _wbfmm_laplace_expansion_apply_gradient(C, cstr, nq, ec, N, f,
						   fstr) ;
    break ;
  }

  return 0 ;
}
