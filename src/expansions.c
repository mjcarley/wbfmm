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

#include <glib.h>

#include <wbfmm.h>

#include "wbfmm-private.h"

static gint expansion_h_increment_cfft(gint n, gint m, gint sgn,
				       WBFMM_REAL jn,
				       WBFMM_REAL *cfft, gint cstr,
				       WBFMM_REAL *Pn,
				       WBFMM_REAL *Cmph, WBFMM_REAL *Smph,
				       WBFMM_REAL *q, gint nq)
{
  gint idx ;
  gint j ;
  
  idx = wbfmm_coefficient_index_nm(n,sgn*m) ;  

  for ( j = 0 ; j < nq ; j ++ ) {
    cfft[2*idx*cstr+2*j+0] +=
      jn*Pn[m]*(q[2*j+0]*Cmph[m] + sgn*q[2*j+1]*Smph[m]) ;
    cfft[2*idx*cstr+2*j+1] +=
      jn*Pn[m]*(q[2*j+1]*Cmph[m] - sgn*q[2*j+0]*Smph[m]) ;
  }
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_expansion_h_cfft)(WBFMM_REAL k, gint N, 
						 WBFMM_REAL *x0, WBFMM_REAL *xs,
						 WBFMM_REAL *q, gint nq,
						 WBFMM_REAL *cfft, gint cstr,
						 WBFMM_REAL *work)

/*
  workspace size: 4*(2*N+1)
*/

{
  WBFMM_REAL jn, jnm1, r, th, ph, kr ;
  WBFMM_REAL Cth, Sth, *Pn, *Pnm1, Cmph[64], Smph[64] ;
  gint n, m ;

  if ( cstr < nq )
    g_error("%s: coefficient stride (cstr=%d) less than number of source"
	    "elements (nq=%d)", __FUNCTION__, cstr, nq) ;
  
  if ( N == 0 ) { return 0 ; }

  /*G&D (2.17)*/

  Pnm1 = &(work[0]) ; Pn = &(Pnm1[2*(2*N+1)]) ;

  WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(x0, xs, &r, &th, &ph) ;
  Cth = COS(th) ; Sth = SIN(th) ; kr = k*r ;

  /*initialize recursions*/
  WBFMM_FUNCTION_NAME(wbfmm_bessel_j_init)(kr, &jnm1, &jn) ;
  WBFMM_FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth,
					   &(Pnm1[0]), &(Pn[0]), &(Pn[1])) ;
  Cmph[0] = 1.0     ; Smph[0] = 0.0 ;
  Cmph[1] = COS(ph) ; Smph[1] = SIN(ph) ;
  
  /*first entries are done by hand*/
  n = 0 ; 
  m = 0 ; 
  expansion_h_increment_cfft(n, m,  1, jnm1, cfft, cstr, Pnm1, Cmph, Smph,
			     q, nq) ;
  
  if ( N == 1 ) { return 0 ; }

  n = 1 ; 
  m = 0 ; 
  expansion_h_increment_cfft(n, m,  1, jn, cfft, cstr, Pn, Cmph, Smph, q, nq) ;
  
  m = 1 ; 
  expansion_h_increment_cfft(n, m,  1, jn, cfft, cstr, Pn, Cmph, Smph, q, nq) ;
  expansion_h_increment_cfft(n, m, -1, jn, cfft, cstr, Pn, Cmph, Smph, q, nq) ;
  
  for ( n = 2 ; n <= N ; n ++ ) {
    WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pnm1, &Pn,
							n-1, Cth, Sth) ;
    WBFMM_FUNCTION_NAME(wbfmm_bessel_j_recursion)(&jnm1, &jn, kr, n-1) ;
    Cmph[n] = Cmph[n-1]*Cmph[1] - Smph[n-1]*Smph[1] ;
    Smph[n] = Smph[n-1]*Cmph[1] + Cmph[n-1]*Smph[1] ;
    
    m = 0 ; 
    expansion_h_increment_cfft(n, m,  1, jn, cfft, cstr, Pn, Cmph, Smph,
			       q, nq) ;
    
    for ( m = 1 ; m <= n ; m ++ ) {
      expansion_h_increment_cfft(n, m,  1, jn, cfft, cstr, Pn, Cmph, Smph,
				 q, nq) ;
      expansion_h_increment_cfft(n, m, -1, jn, cfft, cstr, Pn, Cmph, Smph,
				 q, nq) ;
    }
  }

  return 0 ;
}

static gint expansion_h_increment(gint n, gint m, gint sgn,
				  WBFMM_REAL *hn,
				  WBFMM_REAL *cfft, gint cstr, 
				  WBFMM_REAL *Pn,
				  WBFMM_REAL Cmph, WBFMM_REAL Smph,
				  WBFMM_REAL *field, gint nq, gint fstr)

{
  gint idx, j ;
  WBFMM_REAL ar, ai, tr, ti ;

  idx = wbfmm_coefficient_index_nm(n,sgn*m) ;
  for ( j = 0 ; j < nq ; j ++ ) {
    ar = cfft[2*idx*cstr+2*j+0] ; ai = cfft[2*idx*cstr+2*j+1] ; 
    tr = ar*hn[0] - ai*hn[1] ; ti = ai*hn[0] + ar*hn[1] ;
    field[j*fstr+0] += (Cmph*tr - sgn*Smph*ti)*Pn[m] ;
    field[j*fstr+1] += (Cmph*ti + sgn*Smph*tr)*Pn[m] ;
  }
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_expansion_h_evaluate)(WBFMM_REAL k,
						     WBFMM_REAL *x0,
						     WBFMM_REAL *cfft,
						     gint cstr,
						     gint N, 
						     gint nq,
						     WBFMM_REAL *xf, 
						     WBFMM_REAL *field,
						     gint fstr,
						     WBFMM_REAL *work)

/*
  cstr stride by element (multiply by two to get to complex entry)
*/
		       
{
  WBFMM_REAL hn[2], hnm1[2], r, th, ph, kr ;
  WBFMM_REAL Cth, Sth, *Pn, *Pnm1, Cmph[64], Smph[64] ;
  gint n, m ;

  if ( cstr < nq )
    g_error("%s: coefficient stride (cstr=%d) less than number of source"
	    "elements (nq=%d)", __FUNCTION__, cstr, nq) ;
  if ( fstr < 2 && nq > 1 )
    g_error("%s: field stride (%d) must be greater than 1 for "
	    "multi-component sources (nq=%d)", __FUNCTION__, fstr, nq) ;

  if ( N == 0 ) { return 0 ; }

  Pnm1 = &(work[0]) ; Pn = &(Pnm1[2*(2*N+1)]) ;

  WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(x0, xf, &r, &th, &ph) ;
  Cth = COS(th) ; Sth = SIN(th) ; 
  kr = k*r ;

  /*initialize recursions*/
  WBFMM_FUNCTION_NAME(wbfmm_bessel_h_init)(kr, hnm1, hn) ;
  WBFMM_FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth,
					   &(Pnm1[0]), &(Pn[0]), &(Pn[1])) ;
  Cmph[0] = 1.0     ; Smph[0] = 0.0 ;
  Cmph[1] = COS(ph) ; Smph[1] = SIN(ph) ;

  /*first two terms by hand*/
  n = 0 ; 
  m = 0 ; 
  expansion_h_increment(n, m,  1, hnm1, cfft, cstr, 
			Pnm1, Cmph[m], Smph[m], field, nq, fstr) ;

  if ( N == 1 ) { return 0 ; }

  n = 1 ; 
  m = 0 ; 
  expansion_h_increment(n, m,  1, hn, cfft, cstr, Pn, Cmph[m], Smph[m],
			field, nq, fstr) ;

  m = 1 ; 
  expansion_h_increment(n, m,  1, hn, cfft, cstr, Pn, Cmph[m], Smph[m],
			field, nq, fstr) ;
  expansion_h_increment(n, m, -1, hn, cfft, cstr, Pn, Cmph[m], Smph[m],
			field, nq, fstr) ;

  for ( n = 2 ; n <= N ; n ++ ) {
    WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pnm1, &Pn,
							n-1, Cth, Sth) ;
    WBFMM_FUNCTION_NAME(wbfmm_bessel_h_recursion)(hnm1, hn, kr, n-1) ;
    Cmph[n] = Cmph[n-1]*Cmph[1] - Smph[n-1]*Smph[1] ;
    Smph[n] = Smph[n-1]*Cmph[1] + Cmph[n-1]*Smph[1] ;

    m = 0 ; 
    expansion_h_increment(n, m,  1, hn, cfft, cstr, 
			  Pn, Cmph[m], Smph[m], field, nq, fstr) ;

    for ( m = 1 ; m <= n ; m ++ ) {
      expansion_h_increment(n, m,  1, hn, cfft, cstr, 
			    Pn, Cmph[m], Smph[m], field, nq, fstr) ;
      expansion_h_increment(n, m, -1, hn, cfft, cstr, 
			    Pn, Cmph[m], Smph[m], field, nq, fstr) ;
    }
  }

  return 0 ;
}

static gint expansion_j_increment0(gint n, gint m,
				   WBFMM_REAL jn,
				   WBFMM_REAL *cfft, gint cstr, 
				   WBFMM_REAL *Pn,
				   WBFMM_REAL Cmph, WBFMM_REAL Smph,
				   WBFMM_REAL *field, gint nq, gint fstr)

{
  gint idx, j ;
  WBFMM_REAL ar, ai, jPnC ;

  idx = wbfmm_coefficient_index_nm(n,m) ;
  jPnC = Cmph*jn*Pn[m] ;
  for ( j = 0 ; j < nq ; j ++ ) {
    ar = cfft[2*idx*cstr+2*j+0] ; ai = cfft[2*idx*cstr+2*j+1] ; 
    field[j*fstr+0] += ar*jPnC ; field[j*fstr+1] += ai*jPnC ;
  }

  return 0 ;
}

static gint expansion_j_increment(gint n, gint m,
				  WBFMM_REAL jn,
				  WBFMM_REAL *cfft, gint cstr, 
				  WBFMM_REAL *Pn,
				  WBFMM_REAL Cmph, WBFMM_REAL Smph,
				  WBFMM_REAL *field, gint nq, gint fstr)

{
  gint idx, j, jdx ;
  WBFMM_REAL ar, ai, br, bi, jPnC, jPnS ;

  idx = wbfmm_coefficient_index_nm(n, m) ;
  jdx = wbfmm_coefficient_index_nm(n,-m) ;
  jPnC = jn*Pn[m]*Cmph ; jPnS = jn*Pn[m]*Smph ;
  for ( j = 0 ; j < nq ; j ++ ) {
    ar = cfft[2*idx*cstr+2*j+0] ; ai = cfft[2*idx*cstr+2*j+1] ; 
    br = cfft[2*jdx*cstr+2*j+0] ; bi = cfft[2*jdx*cstr+2*j+1] ;
    field[j*fstr+0] += jPnC*(ar+br) - jPnS*(ai-bi) ;
    field[j*fstr+1] += jPnC*(ai+bi) + jPnS*(ar-br) ;
  }

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_expansion_j_evaluate)(WBFMM_REAL k,
						     WBFMM_REAL *x0,
						     WBFMM_REAL *cfft,
						     gint cstr,
						     gint N, gint nq,
						     WBFMM_REAL *xf, 
						     WBFMM_REAL *field,
						     gint fstr, 
						     WBFMM_REAL *work)

/*
  cstr stride by element (multiply by two to get to complex entry)
*/
		       
{
  WBFMM_REAL jn, jnm1, r, th, ph, kr ;
  WBFMM_REAL Cth, Sth, *Pn, *Pnm1, Cmph[64], Smph[64] ;
  gint n, m ;

  if ( cstr < nq )
    g_error("%s: coefficient stride (cstr=%d) less than number of source"
	    "elements (nq=%d)", __FUNCTION__, cstr, nq) ;

  Pnm1 = &(work[0]) ; Pn = &(Pnm1[2*(2*N+1)]) ;

  WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(x0, xf, &r, &th, &ph) ;
  Cth = COS(th) ; Sth = SIN(th) ; 
  kr = k*r ;

  /*initialize recursions*/
  WBFMM_FUNCTION_NAME(wbfmm_bessel_j_init)(kr, &jnm1, &jn) ;
  WBFMM_FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth,
					   &(Pnm1[0]), &(Pn[0]), &(Pn[1])) ;
  Cmph[0] = 1.0     ; Smph[0] = 0.0 ;
  Cmph[1] = COS(ph) ; Smph[1] = SIN(ph) ;

  /*first two terms by hand*/
  n = 0 ; 
  m = 0 ; 
  expansion_j_increment0(n, m, jnm1, cfft, cstr, 
			Pnm1, Cmph[m], Smph[m], field, nq, fstr) ;

  n = 1 ; 
  m = 0 ; 
  expansion_j_increment0(n, m, jn, cfft, cstr, Pn, Cmph[m], Smph[m],
			 field, nq, fstr) ;

  m = 1 ; 
  expansion_j_increment(n, m, jn, cfft, cstr, Pn, Cmph[m], Smph[m],
			field, nq, fstr) ;

  for ( n = 2 ; n <= N ; n ++ ) {
    WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pnm1, &Pn,
							n-1, Cth, Sth) ;
    WBFMM_FUNCTION_NAME(wbfmm_bessel_j_recursion)(&jnm1, &jn, kr, n-1) ;
    Cmph[n] = Cmph[n-1]*Cmph[1] - Smph[n-1]*Smph[1] ;
    Smph[n] = Smph[n-1]*Cmph[1] + Cmph[n-1]*Smph[1] ;

    m = 0 ; 
    expansion_j_increment0(n, m, jn, cfft, cstr, 
			   Pn, Cmph[m], Smph[m], field, nq, fstr) ;

    for ( m = 1 ; m <= n ; m ++ ) {
      expansion_j_increment(n, m, jn, cfft, cstr, 
			    Pn, Cmph[m], Smph[m], field, nq, fstr) ;
    }
  }

  return 0 ;
}

static gint expansion_dipole_h_increment_cfft(gint N, gint n, gint m, gint sgn,
					      WBFMM_REAL jn,
					      WBFMM_REAL *cfft, gint cstr,
					      WBFMM_REAL *Pn,
					      WBFMM_REAL *Cmph,
					      WBFMM_REAL *Smph,
					      WBFMM_REAL *fm,
					      WBFMM_REAL *fp,
					      WBFMM_REAL *fz,
					      gint nq)
  
{
  gint idx ;
  WBFMM_REAL Rnm, ab ;
  gint j ;
  
  if ( n >= N ) return 0 ;
  
  Rnm = jn*Pn[m] ;

  /*S_{n-1}^{m+1}*/
  idx = wbfmm_coefficient_index_nm(n-1,sgn*m+1) ;
  ab = WBFMM_FUNCTION_NAME(recursion_bnm)(n, sgn*m) ;
  for ( j = 0 ; j < nq ; j ++ ) {
    cfft[2*idx*cstr+2*j+0] -= ab*Rnm*(fm[2*j+0]*Cmph[m] +
				      sgn*fm[2*j+1]*Smph[m]) ;
    cfft[2*idx*cstr+2*j+1] -= ab*Rnm*(fm[2*j+1]*Cmph[m] -
				      sgn*fm[2*j+0]*Smph[m]) ;
  }
  
  /*S_{n-1}^{m-1}*/
  idx = wbfmm_coefficient_index_nm(n-1,sgn*m-1) ;
  ab = WBFMM_FUNCTION_NAME(recursion_bnm)(n, -sgn*m) ;
  for ( j = 0 ; j < nq ; j ++ ) {
    cfft[2*idx*cstr+2*j+0] -= ab*Rnm*(fp[2*j+0]*Cmph[m] +
				      sgn*fp[2*j+1]*Smph[m]) ;
    cfft[2*idx*cstr+2*j+1] -= ab*Rnm*(fp[2*j+1]*Cmph[m] -
				      sgn*fp[2*j+0]*Smph[m]) ;
  }
  
  /*S_{n-1}^{m}*/
  idx = wbfmm_coefficient_index_nm(n-1,sgn*m) ;
  ab = WBFMM_FUNCTION_NAME(recursion_anm)(n-1, sgn*m) ;
  for ( j = 0 ; j < nq ; j ++ ) {
    cfft[2*idx*cstr+2*j+0] += ab*Rnm*(fz[2*j+0]*Cmph[m] +
				      sgn*fz[2*j+1]*Smph[m]) ;
    cfft[2*idx*cstr+2*j+1] += ab*Rnm*(fz[2*j+1]*Cmph[m] -
				      sgn*fz[2*j+0]*Smph[m]) ;
  }
  
  /*S_{n+1}^{m+1}*/
  idx = wbfmm_coefficient_index_nm(n+1,sgn*m+1) ;
  ab = WBFMM_FUNCTION_NAME(recursion_bnm)(n+1, -sgn*m-1) ;
  for ( j = 0 ; j < nq ; j ++ ) {
    cfft[2*idx*cstr+2*j+0] += ab*Rnm*(fm[2*j+0]*Cmph[m] +
				      sgn*fm[2*j+1]*Smph[m]) ;
    cfft[2*idx*cstr+2*j+1] += ab*Rnm*(fm[2*j+1]*Cmph[m] -
				      sgn*fm[2*j+0]*Smph[m]) ;
  }
  
  /*S_{n+1}^{m-1}*/
  idx = wbfmm_coefficient_index_nm(n+1,sgn*m-1) ;
  ab = WBFMM_FUNCTION_NAME(recursion_bnm)(n+1, sgn*m-1) ;
  for ( j = 0 ; j < nq ; j ++ ) {
    cfft[2*idx*cstr+2*j+0] += ab*Rnm*(fp[2*j+0]*Cmph[m] +
				      sgn*fp[2*j+1]*Smph[m]) ;
    cfft[2*idx*cstr+2*j+1] += ab*Rnm*(fp[2*j+1]*Cmph[m] -
				      sgn*fp[2*j+0]*Smph[m]) ;
  }
  
  /*S_{n+1}^{m}*/
  idx = wbfmm_coefficient_index_nm(n+1,sgn*m) ;
  ab = WBFMM_FUNCTION_NAME(recursion_anm)(n, sgn*m) ;
  for ( j = 0 ; j < nq ; j ++ ) {
    cfft[2*idx*cstr+2*j+0] -= ab*Rnm*(fz[2*j+0]*Cmph[m] +
				      sgn*fz[2*j+1]*Smph[m]) ;
    cfft[2*idx*cstr+2*j+1] -= ab*Rnm*(fz[2*j+1]*Cmph[m] -
				      sgn*fz[2*j+0]*Smph[m]) ;
  }
  
  /* /\*S_{n-1}^{m+1}*\/ */
  /* idx = wbfmm_coefficient_index_nm(n-1,sgn*m+1) ; */
  /* ab = WBFMM_FUNCTION_NAME(recursion_bnm)(n, sgn*m) ; */
  /* cfft[2*idx*cstr+0] -= ab*Rnm*(fm[0]*Cmph[m] + sgn*fm[1]*Smph[m]) ; */
  /* cfft[2*idx*cstr+1] -= ab*Rnm*(fm[1]*Cmph[m] - sgn*fm[0]*Smph[m]) ; */

  /* /\*S_{n-1}^{m-1}*\/ */
  /* idx = wbfmm_coefficient_index_nm(n-1,sgn*m-1) ; */
  /* ab = WBFMM_FUNCTION_NAME(recursion_bnm)(n, -sgn*m) ; */
  /* cfft[2*idx*cstr+0] -= ab*Rnm*(fp[0]*Cmph[m] + sgn*fp[1]*Smph[m]) ; */
  /* cfft[2*idx*cstr+1] -= ab*Rnm*(fp[1]*Cmph[m] - sgn*fp[0]*Smph[m]) ; */

  /* /\*S_{n-1}^{m}*\/ */
  /* idx = wbfmm_coefficient_index_nm(n-1,sgn*m) ; */
  /* ab = WBFMM_FUNCTION_NAME(recursion_anm)(n-1, sgn*m) ; */
  /* cfft[2*idx*cstr+0] += ab*Rnm*(fz[0]*Cmph[m] + sgn*fz[1]*Smph[m]) ; */
  /* cfft[2*idx*cstr+1] += ab*Rnm*(fz[1]*Cmph[m] - sgn*fz[0]*Smph[m]) ; */

  /* /\*S_{n+1}^{m+1}*\/ */
  /* idx = wbfmm_coefficient_index_nm(n+1,sgn*m+1) ; */
  /* ab = WBFMM_FUNCTION_NAME(recursion_bnm)(n+1, -sgn*m-1) ; */
  /* cfft[2*idx*cstr+0] += ab*Rnm*(fm[0]*Cmph[m] + sgn*fm[1]*Smph[m]) ; */
  /* cfft[2*idx*cstr+1] += ab*Rnm*(fm[1]*Cmph[m] - sgn*fm[0]*Smph[m]) ; */

  /* /\*S_{n+1}^{m-1}*\/ */
  /* idx = wbfmm_coefficient_index_nm(n+1,sgn*m-1) ; */
  /* ab = WBFMM_FUNCTION_NAME(recursion_bnm)(n+1, sgn*m-1) ; */
  /* cfft[2*idx*cstr+0] += ab*Rnm*(fp[0]*Cmph[m] + sgn*fp[1]*Smph[m]) ; */
  /* cfft[2*idx*cstr+1] += ab*Rnm*(fp[1]*Cmph[m] - sgn*fp[0]*Smph[m]) ; */

  /* /\*S_{n+1}^{m}*\/ */
  /* idx = wbfmm_coefficient_index_nm(n+1,sgn*m) ; */
  /* ab = WBFMM_FUNCTION_NAME(recursion_anm)(n, sgn*m) ; */
  /* cfft[2*idx*cstr+0] -= ab*Rnm*(fz[0]*Cmph[m] + sgn*fz[1]*Smph[m]) ; */
  /* cfft[2*idx*cstr+1] -= ab*Rnm*(fz[1]*Cmph[m] - sgn*fz[0]*Smph[m]) ; */
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_expansion_dipole_h_cfft)(WBFMM_REAL k, gint N, 
							WBFMM_REAL *x0,
							WBFMM_REAL *xs,
							WBFMM_REAL *fxi,
							WBFMM_REAL *fyi,
							WBFMM_REAL *fzi,
							gint nq,
							WBFMM_REAL *cfft,
							gint cstr,
							WBFMM_REAL *work)

/*
  workspace size: 4*(2*N+1)
*/

{
  WBFMM_REAL jn, jnm1, r, th, ph, kr, fm[2], fp[2], fz[2] ;
  WBFMM_REAL Cth, Sth, *Pn, *Pnm1, Cmph[64], Smph[64] ;
  gint n, m ;

  g_assert(nq == 1) ;
  
  /*G&D (2.17) combined with derivatives (3.2)--(3.7)*/

  /*f_{\pm} = (k/2) \mathbf{f}.(i_{x} \pm j i_{y})*/
  fm[0] = 0.5*k*(fxi[0] + fyi[1]) ; fm[1] = 0.5*k*(fxi[1] - fyi[0]) ;
  fp[0] = 0.5*k*(fxi[0] - fyi[1]) ; fp[1] = 0.5*k*(fxi[1] + fyi[0]) ;
  fz[0] = k*fzi[0] ; fz[1] = k*fzi[1] ;
  
  Pnm1 = &(work[0]) ; Pn = &(Pnm1[2*(2*N+1)]) ;

  WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(x0, xs, &r, &th, &ph) ;
  Cth = COS(th) ; Sth = SIN(th) ; kr = k*r ;

  /*initialize recursions*/
  WBFMM_FUNCTION_NAME(wbfmm_bessel_j_init)(kr, &jnm1, &jn) ;
  WBFMM_FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth,
					   &(Pnm1[0]), &(Pn[0]), &(Pn[1])) ;
  Cmph[0] = 1.0     ; Smph[0] = 0.0 ;
  Cmph[1] = COS(ph) ; Smph[1] = SIN(ph) ;
  
  /*first entries are done by hand*/
  n = 0 ; 
  m = 0 ; 
  expansion_dipole_h_increment_cfft(N, n, m,  1, jnm1, cfft, cstr, Pnm1,
				    Cmph, Smph, fm, fp, fz, nq) ;
  
  n = 1 ; 
  m = 0 ; 
  expansion_dipole_h_increment_cfft(N, n, m,  1, jn, cfft, cstr, Pn,
				    Cmph, Smph, fm, fp, fz, nq) ;
  
  m = 1 ; 
  expansion_dipole_h_increment_cfft(N, n, m,  1, jn, cfft, cstr, Pn,
				    Cmph, Smph, fm, fp, fz, nq) ;
  expansion_dipole_h_increment_cfft(N, n, m, -1, jn, cfft, cstr, Pn,
				    Cmph, Smph, fm, fp, fz, nq) ;
  
  for ( n = 2 ; n <= N ; n ++ ) {
    WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pnm1, &Pn,
							n-1, Cth, Sth) ;
    WBFMM_FUNCTION_NAME(wbfmm_bessel_j_recursion)(&jnm1, &jn, kr, n-1) ;
    Cmph[n] = Cmph[n-1]*Cmph[1] - Smph[n-1]*Smph[1] ;
    Smph[n] = Smph[n-1]*Cmph[1] + Cmph[n-1]*Smph[1] ;
    
    m = 0 ; 
    expansion_dipole_h_increment_cfft(N, n, m,  1, jn, cfft, cstr, Pn,
				      Cmph, Smph, fm, fp, fz, nq) ;
    
    for ( m = 1 ; m <= n ; m ++ ) {
      expansion_dipole_h_increment_cfft(N, n, m,  1, jn, cfft, cstr, Pn,
					Cmph, Smph, fm, fp, fz, nq) ;
      expansion_dipole_h_increment_cfft(N, n, m, -1, jn, cfft, cstr, Pn,
					Cmph, Smph, fm, fp, fz, nq) ;
    }
  }

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_expansion_normal_h_cfft)(WBFMM_REAL k, gint N, 
							WBFMM_REAL *x0,
							WBFMM_REAL *xs,
							WBFMM_REAL *normal,
							WBFMM_REAL *q,
							gint nq,
							WBFMM_REAL *cfft,
							gint cstr, 
							WBFMM_REAL *work)

/*
  workspace size: 4*(2*N+1)
*/

{
  WBFMM_REAL jn, jnm1, r, th, ph, kr, fm[32], fp[32], fz[32] ;
  WBFMM_REAL Cth, Sth, *Pn, *Pnm1, Cmph[64], Smph[64] ;
  gint n, m, j ;

  /* g_assert(nq == 1) ; */
  /*G&D (2.17) combined with derivatives (3.2)--(3.7)*/
  g_assert(nq < 16) ;
  
  /*f_{\pm} = (k/2) \mathbf{f}.(i_{x} \pm j i_{y})*/
  /* fm[0] = 0.5*k*(q[0]*normal[0] + q[1]*normal[1]) ; */
  /* fm[1] = 0.5*k*(q[1]*normal[0] - q[0]*normal[1]) ; */
  /* fp[0] = 0.5*k*(q[0]*normal[0] - q[1]*normal[1]) ; */
  /* fp[1] = 0.5*k*(q[1]*normal[0] + q[0]*normal[1]) ; */
  /* fz[0] = k*q[0]*normal[2] ; fz[1] = k*q[1]*normal[2] ; */
  for ( j = 0 ; j < nq ; j ++ ) {
    fm[2*j+0] = 0.5*k*(q[2*j+0]*normal[0] + q[2*j+1]*normal[1]) ;
    fm[2*j+1] = 0.5*k*(q[2*j+1]*normal[0] - q[2*j+0]*normal[1]) ;
    fp[2*j+0] = 0.5*k*(q[2*j+0]*normal[0] - q[2*j+1]*normal[1]) ;
    fp[2*j+1] = 0.5*k*(q[2*j+1]*normal[0] + q[2*j+0]*normal[1]) ;
    fz[2*j+0] = k*q[2*j+0]*normal[2] ; fz[2*j+1] = k*q[2*j+1]*normal[2] ;
  }
  
  Pnm1 = &(work[0]) ; Pn = &(Pnm1[2*(2*N+1)]) ;

  WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(x0, xs, &r, &th, &ph) ;
  Cth = COS(th) ; Sth = SIN(th) ; kr = k*r ;

  /*initialize recursions*/
  WBFMM_FUNCTION_NAME(wbfmm_bessel_j_init)(kr, &jnm1, &jn) ;
  WBFMM_FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth,
					   &(Pnm1[0]), &(Pn[0]), &(Pn[1])) ;
  Cmph[0] = 1.0     ; Smph[0] = 0.0 ;
  Cmph[1] = COS(ph) ; Smph[1] = SIN(ph) ;
  
  /*first entries are done by hand*/
  n = 0 ; 
  m = 0 ; 
  expansion_dipole_h_increment_cfft(N, n, m,  1, jnm1, cfft, cstr, Pnm1,
				    Cmph, Smph, fm, fp, fz, nq) ;
  
  n = 1 ; 
  m = 0 ; 
  expansion_dipole_h_increment_cfft(N, n, m,  1, jn, cfft, cstr, Pn,
				    Cmph, Smph, fm, fp, fz, nq) ;
  
  m = 1 ; 
  expansion_dipole_h_increment_cfft(N, n, m,  1, jn, cfft, cstr, Pn,
				    Cmph, Smph, fm, fp, fz, nq) ;
  expansion_dipole_h_increment_cfft(N, n, m, -1, jn, cfft, cstr, Pn,
				    Cmph, Smph, fm, fp, fz, nq) ;
  
  for ( n = 2 ; n <= N ; n ++ ) {
    WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pnm1, &Pn, n-1,
							Cth, Sth) ;
    WBFMM_FUNCTION_NAME(wbfmm_bessel_j_recursion)(&jnm1, &jn, kr, n-1) ;
    Cmph[n] = Cmph[n-1]*Cmph[1] - Smph[n-1]*Smph[1] ;
    Smph[n] = Smph[n-1]*Cmph[1] + Cmph[n-1]*Smph[1] ;
    
    m = 0 ; 
    expansion_dipole_h_increment_cfft(N, n, m,  1, jn, cfft, cstr, Pn,
				      Cmph, Smph, fm, fp, fz, nq) ;
    
    for ( m = 1 ; m <= n ; m ++ ) {
      expansion_dipole_h_increment_cfft(N, n, m,  1, jn, cfft, cstr, Pn,
					Cmph, Smph, fm, fp, fz, nq) ;
      expansion_dipole_h_increment_cfft(N, n, m, -1, jn, cfft, cstr, Pn,
					Cmph, Smph, fm, fp, fz, nq) ;
    }
  }

  return 0 ;
}

static gint coefficients_j_increment(gint n, gint m,
				     WBFMM_REAL jn,
				     WBFMM_REAL *Pn,
				     WBFMM_REAL Cmph, WBFMM_REAL Smph,
				     WBFMM_REAL *cfft)

{
  gint idx ;

  idx = wbfmm_conjugate_index_nm(n,m) ;
  cfft[2*idx+0] = Cmph*jn*Pn[m] ;
  cfft[2*idx+1] = Smph*jn*Pn[m] ;
  
  return 0 ;
}

static gint _wbfmm_local_coefficients_scalar(WBFMM_REAL kr,
					     WBFMM_REAL *cfft,
					     gint N,
					     WBFMM_REAL Cth,
					     WBFMM_REAL Sth,
					     WBFMM_REAL Cph,
					     WBFMM_REAL Sph,
					     WBFMM_REAL *work)

{
  WBFMM_REAL jn, jnm1 ;
  WBFMM_REAL *Pn, *Pnm1, Cmph[64], Smph[64] ;
  gint n, m ;

  Pnm1 = &(work[0]) ; Pn = &(Pnm1[2*(2*N+1)]) ;

  /*initialize recursions*/
  WBFMM_FUNCTION_NAME(wbfmm_bessel_j_init)(kr, &jnm1, &jn) ;
  WBFMM_FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth,
					   &(Pnm1[0]), &(Pn[0]), &(Pn[1])) ;
  Cmph[0] = 1.0 ; Smph[0] = 0.0 ;
  Cmph[1] = Cph ; Smph[1] = Sph ;

  /*first two terms by hand*/
  n = 0 ; 
  m = 0 ;
  coefficients_j_increment(n, m, jnm1, Pnm1, Cmph[m], Smph[m], cfft) ;

  n = 1 ;  
  for ( m = 0 ; m <= n ; m ++ ) 
    coefficients_j_increment(n, m, jn, Pn, Cmph[m], Smph[m], cfft) ;

  for ( n = 2 ; n <= N ; n ++ ) {
    WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pnm1, &Pn,
							n-1, Cth, Sth) ;
    WBFMM_FUNCTION_NAME(wbfmm_bessel_j_recursion)(&jnm1, &jn, kr, n-1) ;
    Cmph[n] = Cmph[n-1]*Cph - Smph[n-1]*Sph ;
    Smph[n] = Smph[n-1]*Cph + Cmph[n-1]*Sph ;

    for ( m = 0 ; m <= n ; m ++ ) 
      coefficients_j_increment(n, m, jn, Pn, Cmph[m], Smph[m], cfft) ;
  }
  
  return 0 ;
}

static gint coefficients_j_grad_increment(gint n, gint m, gint sgn,
					  WBFMM_REAL k,
					  WBFMM_REAL jnm1, WBFMM_REAL jnp1,
					  WBFMM_REAL *Pnm1, WBFMM_REAL *Pnp1,
					  WBFMM_REAL Cmphm1, WBFMM_REAL Smphm1,
					  WBFMM_REAL Cmph, WBFMM_REAL Smph,
					  WBFMM_REAL Cmphp1, WBFMM_REAL Smphp1,
					  WBFMM_REAL *cfft)

{
  gint idx, mm1, mp1 ;
  WBFMM_REAL tm, tp, a1, a2, b1, b2, *c ;
  
  /*application of G&D (2004) equation 3.7*/  
  idx = wbfmm_coefficient_index_nm(n,sgn*m) ;
  c = &(cfft[6*idx]) ;
  
  /*coefficient times j_{n-1}(kr), j_{n+1}(kr)*/
  tm = k*jnm1 ; tp = k*jnp1 ;

  mm1 = ABS(sgn*m-1) ; mp1 = ABS(sgn*m+1) ;
  
  /*x and y derivatives*/
  b1  = WBFMM_FUNCTION_NAME(recursion_bnm)(n+1, -sgn*m-1)/2.0 ;
  b2  = WBFMM_FUNCTION_NAME(recursion_bnm)(n  ,  sgn*m  )/2.0 ;
  c[0] = c[4] = Cmphp1*(b1*tp*Pnp1[mp1] - b2*tm*Pnm1[mp1]) ;
  c[1] = c[5] = Smphp1*(b1*tp*Pnp1[mp1] - b2*tm*Pnm1[mp1]) ;
  b1  = WBFMM_FUNCTION_NAME(recursion_bnm)(n+1,  sgn*m-1)/2.0 ;
  b2  = WBFMM_FUNCTION_NAME(recursion_bnm)(n  , -sgn*m  )/2.0 ;
  c[2] = -Smphm1*(b1*tp*Pnp1[mm1] - b2*tm*Pnm1[mm1]) ;
  c[3] =  Cmphm1*(b1*tp*Pnp1[mm1] - b2*tm*Pnm1[mm1]) ;

  c[0] += c[3] ; c[1] -= c[2] ;
  c[2] += c[5] ; c[3] -= c[4] ;

  /*z derivative*/
  a1 = WBFMM_FUNCTION_NAME(recursion_anm)(n-1, m) ;
  a2 = WBFMM_FUNCTION_NAME(recursion_anm)(n  , m) ;

  c[4] = Cmph*(a1*tm*Pnm1[m] - a2*tp*Pnp1[m]) ;
  c[5] = Smph*(a1*tm*Pnm1[m] - a2*tp*Pnp1[m]) ;

  return 0 ;
}

static gint _wbfmm_local_coefficients_gradient(WBFMM_REAL k, WBFMM_REAL r,
					       WBFMM_REAL *cfft,
					       gint N,
					       WBFMM_REAL Cth,
					       WBFMM_REAL Sth,
					       WBFMM_REAL Cph,
					       WBFMM_REAL Sph,
					       WBFMM_REAL *work)

{
  WBFMM_REAL jn, jnm1, jnp1, kr ;
  WBFMM_REAL *Pn, *Pnm1, *Pnp1, Cmph[64], Smph[64] ;
  gint n, m ;

  Pnm1 = &(work[0]) ;
  Pn   = &(Pnm1[N+1]) ;
  Pnp1 = &(Pn[N+3]) ;

  kr = k*r ;

  /*initialize recursions*/
  WBFMM_FUNCTION_NAME(wbfmm_bessel_j_init)(kr, &jn, &jnp1) ;
  WBFMM_FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth,
					   &(Pn[0]), &(Pnp1[0]), &(Pnp1[1])) ;

  Cmph[0] = 1.0 ; Smph[0] = 0.0 ;
  Cmph[1] = Cph ; Smph[1] = Sph ;
  jnm1 = 0.0 ;
  /*first two terms by hand*/
  n = 0 ; 
  m = 0 ; 
  coefficients_j_grad_increment(n, m,  1, k, jnm1, jnp1,
				Pnm1, Pnp1,
				Cmph[m+1], -Smph[m+1],  
				Cmph[m  ],  Smph[m  ],
				Cmph[m+1],  Smph[m+1],
				cfft) ;
  
  jnm1 = jn ;
  WBFMM_FUNCTION_NAME(wbfmm_bessel_j_recursion)(&jn, &jnp1, kr, 1) ;
  n = 1 ; 
  m = 0 ;
  memcpy(Pnm1, Pn, (n+1)*sizeof(WBFMM_REAL)) ;
  WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pn, &Pnp1,
						      1, Cth, Sth) ;
  Cmph[n+1] = Cmph[n]*Cph - Smph[n]*Sph ;
  Smph[n+1] = Smph[n]*Cph + Cmph[n]*Sph ;

  coefficients_j_grad_increment(n, m,  1, k, jnm1, jnp1,
				Pnm1, Pnp1,
				Cmph[m+1], -Smph[m+1],
				Cmph[m  ],  Smph[m  ],
				Cmph[m+1],  Smph[m+1],
				cfft) ;

  m = 1 ;
  coefficients_j_grad_increment(n, m,  1, k, jnm1, jnp1,
				Pnm1, Pnp1,
				Cmph[m-1],  Smph[m-1],  
				Cmph[m  ],  Smph[m  ],
				Cmph[m+1],  Smph[m+1],
				cfft) ;
  coefficients_j_grad_increment(n, m, -1, k, jnm1, jnp1,
				Pnm1, Pnp1,
				Cmph[m+1], -Smph[m+1],  
				Cmph[m  ], -Smph[m  ],
				Cmph[m-1], -Smph[m-1],
				cfft) ;

  for ( n = 2 ; n <= N ; n ++ ) {
    memcpy(Pnm1, Pn, (n+1)*sizeof(WBFMM_REAL)) ;
    WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pn, &Pnp1,
							n, Cth, Sth) ;
    jnm1 = jn ;
    WBFMM_FUNCTION_NAME(wbfmm_bessel_j_recursion)(&jn, &jnp1, kr, n) ;

    Cmph[n+1] = Cmph[n]*Cph - Smph[n]*Sph ;
    Smph[n+1] = Smph[n]*Cph + Cmph[n]*Sph ;

    m = 0 ; 
    coefficients_j_grad_increment(n, m,  1, k, jnm1, jnp1,
				  Pnm1, Pnp1,
				  Cmph[m+1], -Smph[m+1],
				  Cmph[m  ],  Smph[m  ],
				  Cmph[m+1],  Smph[m+1],
				  cfft) ;

    for ( m = 1 ; m <= n ; m ++ ) {
      coefficients_j_grad_increment(n, m,  1, k, jnm1, jnp1,
				    Pnm1, Pnp1,
				    Cmph[m-1],  Smph[m-1],  
				    Cmph[m  ],  Smph[m  ],
				    Cmph[m+1],  Smph[m+1],
				    cfft) ;
      coefficients_j_grad_increment(n, m, -1, k, jnm1, jnp1,
				    Pnm1, Pnp1,
				    Cmph[m+1], -Smph[m+1],
				    Cmph[m  ], -Smph[m  ],
				    Cmph[m-1], -Smph[m-1],
				    cfft) ;
    }
  }
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_local_coefficients)(WBFMM_REAL k,
						   WBFMM_REAL *x,
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
    return _wbfmm_local_coefficients_scalar(k*r, cfft, N, Cth, Sth,
					    Cph, Sph, work) ;
    break ;
  case WBFMM_FIELD_GRADIENT:
    return _wbfmm_local_coefficients_gradient(k, r, cfft, N, Cth, Sth,
					      Cph, Sph, work) ;
    break ;
  }

  return 0 ;
}

static gint _wbfmm_expansion_apply_scalar(WBFMM_REAL *C,
					  gint cstr,
					  gint nq,
					  WBFMM_REAL *ec,
					  gint N,
					  WBFMM_REAL *f,
					  gint fstr)
{
  gint n, m, idx, j, jdx ;
  WBFMM_REAL er, ei ;
  
  for ( n = 0 ; n <= N ; n ++ ) {
    m = 0 ;
    idx = 2*wbfmm_conjugate_index_nm(n,m) ;      
    er = ec[idx+0] ; ei = ec[idx+1] ;
    idx = 2*wbfmm_coefficient_index_nm(n,m) ;  

    for ( j = 0 ; j < nq ; j ++ ) {
      f[j*fstr+0] += C[idx*cstr+2*j+0]*er - C[idx*cstr+2*j+1]*ei ;
      f[j*fstr+1] += C[idx*cstr+2*j+0]*ei + C[idx*cstr+2*j+1]*er ;
    }
    for ( m = 1 ; m <= n ; m ++ ) {
      idx = 2*wbfmm_conjugate_index_nm(n,m) ;      
      er = ec[idx+0] ; ei = ec[idx+1] ;
      idx = 2*wbfmm_coefficient_index_nm(n,  m) ;  
      jdx = 2*wbfmm_coefficient_index_nm(n, -m) ;  
      for ( j = 0 ; j < nq ; j ++ ) {
	f[j*fstr+0] +=
	  er*(C[idx*cstr+2*j+0] + C[jdx*cstr+2*j+0]) +
	  ei*(C[jdx*cstr+2*j+1] - C[idx*cstr+2*j+1]) ;
	f[j*fstr+1] +=
	  ei*(C[idx*cstr+2*j+0] - C[jdx*cstr+2*j+0]) +
	  er*(C[idx*cstr+2*j+1] + C[jdx*cstr+2*j+1]) ;
      }
    }
  }
  
  return 0 ;
}

static gint _wbfmm_expansion_apply_gradient(WBFMM_REAL *C,
					    gint cstr,
					    gint nq,
					    WBFMM_REAL *ec,
					    gint N,
					    WBFMM_REAL *f,
					    gint fstr)
{
  gint n, m, idx, j, jdx ;
  WBFMM_REAL *cp, *cm ;
  
  for ( n = 0 ; n <= N ; n ++ ) {
    m = 0 ;
    idx = 6*wbfmm_coefficient_index_nm(n,m) ;
    cp = &(ec[idx]) ;
    idx = 2*wbfmm_coefficient_index_nm(n,m) ;  

    for ( j = 0 ; j < nq ; j ++ ) {
      f[j*fstr+0] += C[idx*cstr+2*j+0]*cp[0] - C[idx*cstr+2*j+1]*cp[1] ;
      f[j*fstr+1] += C[idx*cstr+2*j+0]*cp[1] + C[idx*cstr+2*j+1]*cp[0] ;
      f[j*fstr+2] += C[idx*cstr+2*j+0]*cp[2] - C[idx*cstr+2*j+1]*cp[3] ;
      f[j*fstr+3] += C[idx*cstr+2*j+0]*cp[3] + C[idx*cstr+2*j+1]*cp[2] ;
      f[j*fstr+4] += C[idx*cstr+2*j+0]*cp[4] - C[idx*cstr+2*j+1]*cp[5] ;
      f[j*fstr+5] += C[idx*cstr+2*j+0]*cp[5] + C[idx*cstr+2*j+1]*cp[4] ;
    }
    for ( m = 1 ; m <= n ; m ++ ) {
      idx = 6*wbfmm_coefficient_index_nm(n, m) ;
      cp = &(ec[idx]) ;
      idx = 6*wbfmm_coefficient_index_nm(n,-m) ;
      cm = &(ec[idx]) ;
      idx = 2*wbfmm_coefficient_index_nm(n,  m) ;  
      jdx = 2*wbfmm_coefficient_index_nm(n, -m) ;  
      for ( j = 0 ; j < nq ; j ++ ) {
	f[j*fstr+0] +=
	  C[idx*cstr+2*j+0]*cp[0] - C[idx*cstr+2*j+1]*cp[1] +
	  C[jdx*cstr+2*j+0]*cm[0] - C[jdx*cstr+2*j+1]*cm[1] ;
	f[j*fstr+1] +=
	  C[idx*cstr+2*j+0]*cp[1] + C[idx*cstr+2*j+1]*cp[0] +
	  C[jdx*cstr+2*j+0]*cm[1] + C[jdx*cstr+2*j+1]*cm[0] ;

	f[j*fstr+2] +=
	  C[idx*cstr+2*j+0]*cp[2] - C[idx*cstr+2*j+1]*cp[3] +
	  C[jdx*cstr+2*j+0]*cm[2] - C[jdx*cstr+2*j+1]*cm[3] ;
	f[j*fstr+3] +=
	  C[idx*cstr+2*j+0]*cp[3] + C[idx*cstr+2*j+1]*cp[2] +
	  C[jdx*cstr+2*j+0]*cm[3] + C[jdx*cstr+2*j+1]*cm[2] ;

	f[j*fstr+4] +=
	  C[idx*cstr+2*j+0]*cp[4] - C[idx*cstr+2*j+1]*cp[5] +
	  C[jdx*cstr+2*j+0]*cm[4] - C[jdx*cstr+2*j+1]*cm[5] ;
	f[j*fstr+5] +=
	  C[idx*cstr+2*j+0]*cp[5] + C[idx*cstr+2*j+1]*cp[4] +
	  C[jdx*cstr+2*j+0]*cm[5] + C[jdx*cstr+2*j+1]*cm[4] ;
      }
    }
  }
  
  return 0 ;
}
 
gint WBFMM_FUNCTION_NAME(wbfmm_expansion_apply)(WBFMM_REAL *C,
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
    return _wbfmm_expansion_apply_scalar(C, cstr, nq, ec, N, f, fstr) ;
    break ;
  case WBFMM_FIELD_GRADIENT:
    return _wbfmm_expansion_apply_gradient(C, cstr, nq, ec, N, f, fstr) ;
    break ;
  }

  return 0 ;
}

