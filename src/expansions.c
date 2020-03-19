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
				       WBFMM_REAL *q)
{
  gint idx ;

  idx = wbfmm_coefficient_index_nm(n,sgn*m) ;
  cfft[2*idx*cstr+0] += jn*Pn[m]*(q[0]*Cmph[m] + sgn*q[1]*Smph[m]) ;
  cfft[2*idx*cstr+1] += jn*Pn[m]*(q[1]*Cmph[m] - sgn*q[0]*Smph[m]) ;

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_expansion_h_cfft)(WBFMM_REAL k, gint N, 
						 WBFMM_REAL *x0, WBFMM_REAL *xs,
						 WBFMM_REAL *q, 
						 WBFMM_REAL *cfft, gint cstr,
						 WBFMM_REAL *work)

/*
  workspace size: 4*(2*N+1)
*/

{
  WBFMM_REAL jn, jnm1, r, th, ph, kr ;
  WBFMM_REAL Cth, Sth, *Pn, *Pnm1, Cph, Sph, Cmph[64], Smph[64] ;
  gint n, m ;

  /*G&D (2.17)*/

  Pnm1 = &(work[0]) ; Pn = &(Pnm1[2*(2*N+1)]) ;

  WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(x0, xs, &r, &th, &ph) ;
  Cth = COS(th) ; Sth = SIN(th) ; kr = k*r ;
  Cph = COS(ph) ; Sph = SIN(ph) ; 

  /*initialize recursions*/
  WBFMM_FUNCTION_NAME(wbfmm_bessel_j_init)(kr, &jnm1, &jn) ;
  WBFMM_FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth,
					   &(Pnm1[0]), &(Pn[0]), &(Pn[1])) ;
  Cmph[0] = 1.0 ; Smph[0] = 0.0 ;
  Cmph[1] = Cph ; Smph[1] = Sph ;
  
  /*first entries are done by hand*/
  n = 0 ; 
  m = 0 ; 
  expansion_h_increment_cfft(n, m,  1, jnm1, cfft, cstr, Pnm1, Cmph, Smph, q) ;
  
  n = 1 ; 
  m = 0 ; 
  expansion_h_increment_cfft(n, m,  1, jn, cfft, cstr, Pn, Cmph, Smph, q) ;
  
  m = 1 ; 
  expansion_h_increment_cfft(n, m,  1, jn, cfft, cstr, Pn, Cmph, Smph, q) ;
  expansion_h_increment_cfft(n, m, -1, jn, cfft, cstr, Pn, Cmph, Smph, q) ;
  
  for ( n = 2 ; n <= N ; n ++ ) {
    WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pnm1, &Pn,
							n-1, Cth, Sth) ;
    WBFMM_FUNCTION_NAME(wbfmm_bessel_j_recursion)(&jnm1, &jn, kr, n-1) ;
    Cmph[n] = Cmph[n-1]*Cph - Smph[n-1]*Sph ;
    Smph[n] = Smph[n-1]*Cph + Cmph[n-1]*Sph ;
    
    m = 0 ; 
    expansion_h_increment_cfft(n, m,  1, jn, cfft, cstr, Pn, Cmph, Smph, q) ;
    
    for ( m = 1 ; m <= n ; m ++ ) {
      expansion_h_increment_cfft(n, m,  1, jn, cfft, cstr, Pn, Cmph, Smph, q) ;
      expansion_h_increment_cfft(n, m, -1, jn, cfft, cstr, Pn, Cmph, Smph, q) ;
    }
  }

  return 0 ;
}

static gint expansion_h_increment(gint n, gint m, gint sgn,
				  WBFMM_REAL *hn,
				  WBFMM_REAL *cfft, gint cstr, 
				  WBFMM_REAL *Pn,
				  WBFMM_REAL Cmph, WBFMM_REAL Smph,
				  WBFMM_REAL *field)

{
  gint idx ;
  WBFMM_REAL ar, ai, tr, ti ;

  idx = wbfmm_coefficient_index_nm(n,sgn*m) ;
  ar = cfft[2*idx*cstr+0] ; ai = cfft[2*idx*cstr+1] ; 
  tr = ar*hn[0] - ai*hn[1] ; ti = ai*hn[0] + ar*hn[1] ;
  field[0] += (Cmph*tr - sgn*Smph*ti)*Pn[m] ;
  field[1] += (Cmph*ti + sgn*Smph*tr)*Pn[m] ;
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_expansion_h_evaluate)(WBFMM_REAL k,
						     WBFMM_REAL *x0,
						     WBFMM_REAL *cfft,
						     gint cstr,
						     gint N, 
						     WBFMM_REAL *xf, 
						     WBFMM_REAL *field,
						     WBFMM_REAL *work)

/*
  cstr stride by element (multiply by two to get to complex entry)
*/
		       
{
  WBFMM_REAL hn[2], hnm1[2], r, th, ph, kr ;
  WBFMM_REAL Cth, Sth, *Pn, *Pnm1, Cph, Sph, Cmph[64], Smph[64] ;
  gint n, m ;

  Pnm1 = &(work[0]) ; Pn = &(Pnm1[2*(2*N+1)]) ;

  WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(x0, xf, &r, &th, &ph) ;
  Cth = COS(th) ; Sth = SIN(th) ; 
  Cph = COS(ph) ; Sph = SIN(ph) ; 
  kr = k*r ;

  /*initialize recursions*/
  WBFMM_FUNCTION_NAME(wbfmm_bessel_h_init)(kr, hnm1, hn) ;
  WBFMM_FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth,
					   &(Pnm1[0]), &(Pn[0]), &(Pn[1])) ;
  Cmph[0] = 1.0 ; Smph[0] = 0.0 ;
  Cmph[1] = Cph ; Smph[1] = Sph ;

  /*first two terms by hand*/
  n = 0 ; 
  m = 0 ; 
  expansion_h_increment(n, m,  1, hnm1, cfft, cstr, 
			Pnm1, Cmph[m], Smph[m], field) ;

  n = 1 ; 
  m = 0 ; 
  expansion_h_increment(n, m,  1, hn, cfft, cstr, Pn, Cmph[m], Smph[m], field) ;

  m = 1 ; 
  expansion_h_increment(n, m,  1, hn, cfft, cstr, Pn, Cmph[m], Smph[m], field) ;
  expansion_h_increment(n, m, -1, hn, cfft, cstr, Pn, Cmph[m], Smph[m], field) ;

  for ( n = 2 ; n <= N ; n ++ ) {
    WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pnm1, &Pn,
							n-1, Cth, Sth) ;
    WBFMM_FUNCTION_NAME(wbfmm_bessel_h_recursion)(hnm1, hn, kr, n-1) ;
    Cmph[n] = Cmph[n-1]*Cph - Smph[n-1]*Sph ;
    Smph[n] = Smph[n-1]*Cph + Cmph[n-1]*Sph ;

    m = 0 ; 
    expansion_h_increment(n, m,  1, hn, cfft, cstr, 
			  Pn, Cmph[m], Smph[m], field) ;

    for ( m = 1 ; m <= n ; m ++ ) {
      expansion_h_increment(n, m,  1, hn, cfft, cstr, 
			    Pn, Cmph[m], Smph[m], field) ;
      expansion_h_increment(n, m, -1, hn, cfft, cstr, 
			    Pn, Cmph[m], Smph[m], field) ;
    }
  }

  return 0 ;
}

static gint expansion_j_increment(gint n, gint m, gint sgn,
				  WBFMM_REAL jn,
				  WBFMM_REAL *cfft, gint cstr, 
				  WBFMM_REAL *Pn,
				  WBFMM_REAL Cmph, WBFMM_REAL Smph,
				  WBFMM_REAL *field)

{
  gint idx ;
  WBFMM_REAL ar, ai, tr, ti ;

  idx = wbfmm_coefficient_index_nm(n,sgn*m) ;
  ar = cfft[2*idx*cstr+0] ; ai = cfft[2*idx*cstr+1] ; 
  tr = ar*jn ; ti = ai*jn ;
  field[0] += (Cmph*tr - sgn*Smph*ti)*Pn[m] ;
  field[1] += (Cmph*ti + sgn*Smph*tr)*Pn[m] ;
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_expansion_j_evaluate)(WBFMM_REAL k,
						     WBFMM_REAL *x0,
						     WBFMM_REAL *cfft,
						     gint cstr,
						     gint N, 
						     WBFMM_REAL *xf, 
						     WBFMM_REAL *field,
						     WBFMM_REAL *work)

/*
  cstr stride by element (multiply by two to get to complex entry)
*/
		       
{
  WBFMM_REAL jn, jnm1, r, th, ph, kr ;
  WBFMM_REAL Cth, Sth, *Pn, *Pnm1, Cph, Sph, Cmph[64], Smph[64] ;
  gint n, m ;

  Pnm1 = &(work[0]) ; Pn = &(Pnm1[2*(2*N+1)]) ;

  WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(x0, xf, &r, &th, &ph) ;
  Cth = COS(th) ; Sth = SIN(th) ; 
  Cph = COS(ph) ; Sph = SIN(ph) ; 
  kr = k*r ;

  /*initialize recursions*/
  WBFMM_FUNCTION_NAME(wbfmm_bessel_j_init)(kr, &jnm1, &jn) ;
  WBFMM_FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth,
					   &(Pnm1[0]), &(Pn[0]), &(Pn[1])) ;
  Cmph[0] = 1.0 ; Smph[0] = 0.0 ;
  Cmph[1] = Cph ; Smph[1] = Sph ;

  /*first two terms by hand*/
  n = 0 ; 
  m = 0 ; 
  expansion_j_increment(n, m,  1, jnm1, cfft, cstr, 
			Pnm1, Cmph[m], Smph[m], field) ;

  n = 1 ; 
  m = 0 ; 
  expansion_j_increment(n, m,  1, jn, cfft, cstr, Pn, Cmph[m], Smph[m], field) ;

  m = 1 ; 
  expansion_j_increment(n, m,  1, jn, cfft, cstr, Pn, Cmph[m], Smph[m], field) ;
  expansion_j_increment(n, m, -1, jn, cfft, cstr, Pn, Cmph[m], Smph[m], field) ;

  for ( n = 2 ; n <= N ; n ++ ) {
    WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pnm1, &Pn,
							n-1, Cth, Sth) ;
    WBFMM_FUNCTION_NAME(wbfmm_bessel_j_recursion)(&jnm1, &jn, kr, n-1) ;
    Cmph[n] = Cmph[n-1]*Cph - Smph[n-1]*Sph ;
    Smph[n] = Smph[n-1]*Cph + Cmph[n-1]*Sph ;

    m = 0 ; 
    expansion_j_increment(n, m,  1, jn, cfft, cstr, 
			  Pn, Cmph[m], Smph[m], field) ;

    for ( m = 1 ; m <= n ; m ++ ) {
      expansion_j_increment(n, m,  1, jn, cfft, cstr, 
			    Pn, Cmph[m], Smph[m], field) ;
      expansion_j_increment(n, m, -1, jn, cfft, cstr, 
			    Pn, Cmph[m], Smph[m], field) ;
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
					      WBFMM_REAL *fz)
  
{
  gint idx ;
  WBFMM_REAL Rnm, ab ;

  if ( n >= N ) return 0 ;
  
  Rnm = jn*Pn[m] ;

  /*S_{n-1}^{m+1}*/
  idx = wbfmm_coefficient_index_nm(n-1,sgn*m+1) ;
  ab = WBFMM_FUNCTION_NAME(recursion_bnm)(n, sgn*m) ;
  cfft[2*idx*cstr+0] -= ab*Rnm*(fm[0]*Cmph[m] + sgn*fm[1]*Smph[m]) ;
  cfft[2*idx*cstr+1] -= ab*Rnm*(fm[1]*Cmph[m] - sgn*fm[0]*Smph[m]) ;

  /*S_{n-1}^{m-1}*/
  idx = wbfmm_coefficient_index_nm(n-1,sgn*m-1) ;
  ab = WBFMM_FUNCTION_NAME(recursion_bnm)(n, -sgn*m) ;
  cfft[2*idx*cstr+0] -= ab*Rnm*(fp[0]*Cmph[m] + sgn*fp[1]*Smph[m]) ;
  cfft[2*idx*cstr+1] -= ab*Rnm*(fp[1]*Cmph[m] - sgn*fp[0]*Smph[m]) ;

  /*S_{n-1}^{m}*/
  idx = wbfmm_coefficient_index_nm(n-1,sgn*m) ;
  ab = WBFMM_FUNCTION_NAME(recursion_anm)(n-1, sgn*m) ;
  cfft[2*idx*cstr+0] += ab*Rnm*(fz[0]*Cmph[m] + sgn*fz[1]*Smph[m]) ;
  cfft[2*idx*cstr+1] += ab*Rnm*(fz[1]*Cmph[m] - sgn*fz[0]*Smph[m]) ;

  /*S_{n+1}^{m+1}*/
  idx = wbfmm_coefficient_index_nm(n+1,sgn*m+1) ;
  ab = WBFMM_FUNCTION_NAME(recursion_bnm)(n+1, -sgn*m-1) ;
  cfft[2*idx*cstr+0] += ab*Rnm*(fm[0]*Cmph[m] + sgn*fm[1]*Smph[m]) ;
  cfft[2*idx*cstr+1] += ab*Rnm*(fm[1]*Cmph[m] - sgn*fm[0]*Smph[m]) ;

  /*S_{n+1}^{m-1}*/
  idx = wbfmm_coefficient_index_nm(n+1,sgn*m-1) ;
  ab = WBFMM_FUNCTION_NAME(recursion_bnm)(n+1, sgn*m-1) ;
  cfft[2*idx*cstr+0] += ab*Rnm*(fp[0]*Cmph[m] + sgn*fp[1]*Smph[m]) ;
  cfft[2*idx*cstr+1] += ab*Rnm*(fp[1]*Cmph[m] - sgn*fp[0]*Smph[m]) ;

  /*S_{n+1}^{m}*/
  idx = wbfmm_coefficient_index_nm(n+1,sgn*m) ;
  ab = WBFMM_FUNCTION_NAME(recursion_anm)(n, sgn*m) ;
  cfft[2*idx*cstr+0] -= ab*Rnm*(fz[0]*Cmph[m] + sgn*fz[1]*Smph[m]) ;
  cfft[2*idx*cstr+1] -= ab*Rnm*(fz[1]*Cmph[m] - sgn*fz[0]*Smph[m]) ;
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_expansion_dipole_h_cfft)(WBFMM_REAL k, gint N, 
							WBFMM_REAL *x0,
							WBFMM_REAL *xs,
							WBFMM_REAL *fxi,
							WBFMM_REAL *fyi,
							WBFMM_REAL *fzi,
							WBFMM_REAL *cfft,
							gint cstr,
							WBFMM_REAL *work)

/*
  workspace size: 4*(2*N+1)
*/

{
  WBFMM_REAL jn, jnm1, r, th, ph, kr, fm[2], fp[2], fz[2] ;
  WBFMM_REAL Cth, Sth, *Pn, *Pnm1, Cph, Sph, Cmph[64], Smph[64] ;
  gint n, m ;

  /*G&D (2.17) combined with derivatives (3.2)--(3.7)*/

  /*f_{\pm} = (k/2) \mathbf{f}.(i_{x} \pm j i_{y})*/
  fm[0] = 0.5*k*(fxi[0] + fyi[1]) ; fm[1] = 0.5*k*(fxi[1] - fyi[0]) ;
  fp[0] = 0.5*k*(fxi[0] - fyi[1]) ; fp[1] = 0.5*k*(fxi[1] + fyi[0]) ;
  fz[0] = k*fzi[0] ; fz[1] = k*fzi[1] ;
  
  Pnm1 = &(work[0]) ; Pn = &(Pnm1[2*(2*N+1)]) ;

  WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(x0, xs, &r, &th, &ph) ;
  Cth = COS(th) ; Sth = SIN(th) ; kr = k*r ;
  Cph = COS(ph) ; Sph = SIN(ph) ; 

  /*initialize recursions*/
  WBFMM_FUNCTION_NAME(wbfmm_bessel_j_init)(kr, &jnm1, &jn) ;
  WBFMM_FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth,
					   &(Pnm1[0]), &(Pn[0]), &(Pn[1])) ;
  Cmph[0] = 1.0 ; Smph[0] = 0.0 ;
  Cmph[1] = Cph ; Smph[1] = Sph ;
  
  /*first entries are done by hand*/
  n = 0 ; 
  m = 0 ; 
  expansion_dipole_h_increment_cfft(N, n, m,  1, jnm1, cfft, cstr, Pnm1,
				    Cmph, Smph, fm, fp, fz) ;
  
  n = 1 ; 
  m = 0 ; 
  expansion_dipole_h_increment_cfft(N, n, m,  1, jn, cfft, cstr, Pn,
				    Cmph, Smph, fm, fp, fz) ;
  
  m = 1 ; 
  expansion_dipole_h_increment_cfft(N, n, m,  1, jn, cfft, cstr, Pn,
				    Cmph, Smph, fm, fp, fz) ;
  expansion_dipole_h_increment_cfft(N, n, m, -1, jn, cfft, cstr, Pn,
				    Cmph, Smph, fm, fp, fz) ;
  
  for ( n = 2 ; n <= N ; n ++ ) {
    WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pnm1, &Pn,
							n-1, Cth, Sth) ;
    WBFMM_FUNCTION_NAME(wbfmm_bessel_j_recursion)(&jnm1, &jn, kr, n-1) ;
    Cmph[n] = Cmph[n-1]*Cph - Smph[n-1]*Sph ;
    Smph[n] = Smph[n-1]*Cph + Cmph[n-1]*Sph ;
    
    m = 0 ; 
    expansion_dipole_h_increment_cfft(N, n, m,  1, jn, cfft, cstr, Pn,
				      Cmph, Smph, fm, fp, fz) ;
    
    for ( m = 1 ; m <= n ; m ++ ) {
      expansion_dipole_h_increment_cfft(N, n, m,  1, jn, cfft, cstr, Pn,
					Cmph, Smph, fm, fp, fz) ;
      expansion_dipole_h_increment_cfft(N, n, m, -1, jn, cfft, cstr, Pn,
					Cmph, Smph, fm, fp, fz) ;
    }
  }

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_expansion_normal_h_cfft)(WBFMM_REAL k, gint N, 
							WBFMM_REAL *x0,
							WBFMM_REAL *xs,
							WBFMM_REAL *normal,
							WBFMM_REAL *q,
							WBFMM_REAL *cfft,
							gint cstr,
							WBFMM_REAL *work)

/*
  workspace size: 4*(2*N+1)
*/

{
  WBFMM_REAL jn, jnm1, r, th, ph, kr, fm[2], fp[2], fz[2] ;
  WBFMM_REAL Cth, Sth, *Pn, *Pnm1, Cph, Sph, Cmph[64], Smph[64] ;
  gint n, m ;

  /*G&D (2.17) combined with derivatives (3.2)--(3.7)*/

  /*f_{\pm} = (k/2) \mathbf{f}.(i_{x} \pm j i_{y})*/
  fm[0] = 0.5*k*(q[0]*normal[0] + q[1]*normal[1]) ;
  fm[1] = 0.5*k*(q[1]*normal[0] - q[0]*normal[1]) ;
  fp[0] = 0.5*k*(q[0]*normal[0] - q[1]*normal[1]) ;
  fp[1] = 0.5*k*(q[1]*normal[0] + q[0]*normal[1]) ;
  fz[0] = k*q[0]*normal[2] ; fz[1] = k*q[1]*normal[2] ;
    
  Pnm1 = &(work[0]) ; Pn = &(Pnm1[2*(2*N+1)]) ;

  WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(x0, xs, &r, &th, &ph) ;
  Cth = COS(th) ; Sth = SIN(th) ; kr = k*r ;
  Cph = COS(ph) ; Sph = SIN(ph) ; 

  /*initialize recursions*/
  WBFMM_FUNCTION_NAME(wbfmm_bessel_j_init)(kr, &jnm1, &jn) ;
  WBFMM_FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth,
					   &(Pnm1[0]), &(Pn[0]), &(Pn[1])) ;
  Cmph[0] = 1.0 ; Smph[0] = 0.0 ;
  Cmph[1] = Cph ; Smph[1] = Sph ;
  
  /*first entries are done by hand*/
  n = 0 ; 
  m = 0 ; 
  expansion_dipole_h_increment_cfft(N, n, m,  1, jnm1, cfft, cstr, Pnm1,
				    Cmph, Smph, fm, fp, fz) ;
  
  n = 1 ; 
  m = 0 ; 
  expansion_dipole_h_increment_cfft(N, n, m,  1, jn, cfft, cstr, Pn,
				    Cmph, Smph, fm, fp, fz) ;
  
  m = 1 ; 
  expansion_dipole_h_increment_cfft(N, n, m,  1, jn, cfft, cstr, Pn,
				    Cmph, Smph, fm, fp, fz) ;
  expansion_dipole_h_increment_cfft(N, n, m, -1, jn, cfft, cstr, Pn,
				    Cmph, Smph, fm, fp, fz) ;
  
  for ( n = 2 ; n <= N ; n ++ ) {
    WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pnm1, &Pn, n-1,
							Cth, Sth) ;
    WBFMM_FUNCTION_NAME(wbfmm_bessel_j_recursion)(&jnm1, &jn, kr, n-1) ;
    Cmph[n] = Cmph[n-1]*Cph - Smph[n-1]*Sph ;
    Smph[n] = Smph[n-1]*Cph + Cmph[n-1]*Sph ;
    
    m = 0 ; 
    expansion_dipole_h_increment_cfft(N, n, m,  1, jn, cfft, cstr, Pn,
				      Cmph, Smph, fm, fp, fz) ;
    
    for ( m = 1 ; m <= n ; m ++ ) {
      expansion_dipole_h_increment_cfft(N, n, m,  1, jn, cfft, cstr, Pn,
					Cmph, Smph, fm, fp, fz) ;
      expansion_dipole_h_increment_cfft(N, n, m, -1, jn, cfft, cstr, Pn,
					Cmph, Smph, fm, fp, fz) ;
    }
  }

  return 0 ;
}

static gint expansion_h_grad_increment(gint n, gint m, gint sgn,
				       WBFMM_REAL k,
				       WBFMM_REAL *hnm1, WBFMM_REAL *hnp1,
				       WBFMM_REAL *cfft, gint cstr, 
				       WBFMM_REAL *Pnm1, WBFMM_REAL *Pnp1,
				       WBFMM_REAL Cmphm1, WBFMM_REAL Smphm1,
				       WBFMM_REAL Cmph, WBFMM_REAL Smph,
				       WBFMM_REAL Cmphp1, WBFMM_REAL Smphp1,
				       WBFMM_REAL *field)

{
  gint idx, mm1, mp1 ;
  WBFMM_REAL ar, ai, tm1r, tm1i, tp1r, tp1i, a1, a2, b1, b2 ;
  WBFMM_REAL d1r, d1i, d2r, d2i ;
  
  /*application of G&D (2004) equation 3.7*/
  
  idx = wbfmm_coefficient_index_nm(n,sgn*m) ;
  ar = k*cfft[2*idx*cstr+0] ; ai = k*cfft[2*idx*cstr+1] ;

  /*coefficient times h_{n-1}(kr), h_{n+1}(kr)*/
  tm1r = ar*hnm1[0] - ai*hnm1[1] ; tm1i = ai*hnm1[0] + ar*hnm1[1] ;
  tp1r = ar*hnp1[0] - ai*hnp1[1] ; tp1i = ai*hnp1[0] + ar*hnp1[1] ;

  /*z derivative*/
  a1 = WBFMM_FUNCTION_NAME(recursion_anm)(n-1, m) ;
  a2 = WBFMM_FUNCTION_NAME(recursion_anm)(n  , m) ;
  d1r = a1*(Cmph*tm1r - Smph*tm1i)*Pnm1[m] -
    a2*(Cmph*tp1r - Smph*tp1i)*Pnp1[m] ;
  d1i = a1*(Cmph*tm1i + Smph*tm1r)*Pnm1[m] -
    a2*(Cmph*tp1i + Smph*tp1r)*Pnp1[m] ;
  field[4] += d1r ; field[5] += d1i ;

  mm1 = ABS(sgn*m-1) ; mp1 = ABS(sgn*m+1) ;
  
  /*x and y derivatives*/
  b1 = WBFMM_FUNCTION_NAME(recursion_bnm)(n+1, -sgn*m-1)/2.0 ;
  b2 = WBFMM_FUNCTION_NAME(recursion_bnm)(n  ,  sgn*m  )/2.0 ;
  d1r =
    b1*(Cmphp1*tp1r - Smphp1*tp1i)*Pnp1[mp1] -
    b2*(Cmphp1*tm1r - Smphp1*tm1i)*Pnm1[mp1] ;
  d1i =
    b1*(Cmphp1*tp1i + Smphp1*tp1r)*Pnp1[mp1] -
    b2*(Cmphp1*tm1i + Smphp1*tm1r)*Pnm1[mp1] ;
  b1 = WBFMM_FUNCTION_NAME(recursion_bnm)(n+1,  sgn*m-1)/2.0 ;
  b2 = WBFMM_FUNCTION_NAME(recursion_bnm)(n  , -sgn*m  )/2.0 ;
  d2r =
    b1*(Cmphm1*tp1r - Smphm1*tp1i)*Pnp1[mm1] -
    b2*(Cmphm1*tm1r - Smphm1*tm1i)*Pnm1[mm1] ;
  d2i =
    b1*(Cmphm1*tp1i + Smphm1*tp1r)*Pnp1[mm1] -
    b2*(Cmphm1*tm1i + Smphm1*tm1r)*Pnm1[mm1] ;

  field[0] += d2r + d1r ;
  field[1] += d2i + d1i ;
  field[2] -= d2i - d1i ;
  field[3] += d2r - d1r ;
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_expansion_h_grad_evaluate)(WBFMM_REAL k,
							  WBFMM_REAL *x0,
							  WBFMM_REAL *cfft,
							  gint cstr,
							  gint N, 
							  WBFMM_REAL *xf, 
							  WBFMM_REAL *field,
							  WBFMM_REAL *work)

/*
  cstr stride by element (multiply by two to get to complex entry)
*/
		       
{
  WBFMM_REAL hn[2], hnm1[2], hnp1[2], r, th, ph, kr ;
  WBFMM_REAL Cth, Sth, *Pn, *Pnm1, *Pnp1, Cph, Sph, Cmph[64], Smph[64] ;
  gint n, m ;

  Pnm1 = &(work[0]) ;
  Pn   = &(Pnm1[2*(2*N+1)]) ;
  Pnp1 = &(Pn[2*(2*N+3)]) ;

  memset(Pnm1, 0, 2*(2*N+1)*sizeof(WBFMM_REAL)) ;
  memset(Pn  , 0, 2*(2*N+3)*sizeof(WBFMM_REAL)) ;
  memset(Pnp1, 0, 2*(2*N+5)*sizeof(WBFMM_REAL)) ;
  
  WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(x0, xf, &r, &th, &ph) ;
  Cth = COS(th) ; Sth = SIN(th) ; 
  Cph = COS(ph) ; Sph = SIN(ph) ; 
  kr = k*r ;

  /*initialize recursions*/
  WBFMM_FUNCTION_NAME(wbfmm_bessel_h_init)(kr, hn, hnp1) ;
  WBFMM_FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth,
					   &(Pn[0]), &(Pnp1[0]), &(Pnp1[1])) ;
  Cmph[0] = 1.0 ; Smph[0] = 0.0 ;
  Cmph[1] = Cph ; Smph[1] = Sph ;
  
  /*first two terms by hand*/
  n = 0 ; 
  m = 0 ; 
  expansion_h_grad_increment(n, m,  1, k, hnm1, hnp1, cfft, cstr, 
			     Pnm1, Pnp1,
			     Cmph[m+1], -Smph[m+1],  
			     Cmph[m  ],  Smph[m  ],
			     Cmph[m+1],  Smph[m+1],
			     field) ;
  
  hnm1[0] = hn[0] ; hnm1[1] = hn[1] ; 
  WBFMM_FUNCTION_NAME(wbfmm_bessel_h_recursion)(hn, hnp1, kr, 1) ;

  n = 1 ; 
  m = 0 ; 
  memcpy(Pnm1, Pn, (n+1)*sizeof(WBFMM_REAL)) ;
  WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pn, &Pnp1,
						      1, Cth, Sth) ;
  Cmph[n+1] = Cmph[n]*Cph - Smph[n]*Sph ;
  Smph[n+1] = Smph[n]*Cph + Cmph[n]*Sph ;

  expansion_h_grad_increment(n, m,  1, k, hnm1, hnp1, cfft, cstr, 
			     Pnm1, Pnp1,
			     Cmph[m+1], -Smph[m+1],
			     Cmph[m  ],  Smph[m  ],
			     Cmph[m+1],  Smph[m+1],
			     field) ;

  m = 1 ;
  expansion_h_grad_increment(n, m,  1, k, hnm1, hnp1, cfft, cstr, 
			     Pnm1, Pnp1,
			     Cmph[m-1],  Smph[m-1],  
			     Cmph[m  ],  Smph[m  ],
			     Cmph[m+1],  Smph[m+1],
			     field) ;
  expansion_h_grad_increment(n, m, -1, k, hnm1, hnp1, cfft, cstr, 
			     Pnm1, Pnp1,
			     Cmph[m+1], -Smph[m+1],  
			     Cmph[m  ], -Smph[m  ],
			     Cmph[m-1], -Smph[m-1],
			     field) ;

  for ( n = 2 ; n <= N ; n ++ ) {
    memcpy(Pnm1, Pn, (n+1)*sizeof(WBFMM_REAL)) ;
    WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pn, &Pnp1,
							n, Cth, Sth) ;
    hnm1[0] = hn[0] ; hnm1[1] = hn[1] ; 
    WBFMM_FUNCTION_NAME(wbfmm_bessel_h_recursion)(hn, hnp1, kr, n) ;

    Cmph[n+1] = Cmph[n]*Cph - Smph[n]*Sph ;
    Smph[n+1] = Smph[n]*Cph + Cmph[n]*Sph ;

    m = 0 ; 
    expansion_h_grad_increment(n, m,  1, k, hnm1, hnp1, cfft, cstr, 
			       Pnm1, Pnp1,
			       Cmph[m+1], -Smph[m+1],
			       Cmph[m  ],  Smph[m  ],
			       Cmph[m+1],  Smph[m+1],
			       field) ;

    for ( m = 1 ; m <= n ; m ++ ) {
      expansion_h_grad_increment(n, m,  1, k, hnm1, hnp1, cfft, cstr, 
				 Pnm1, Pnp1,
				 Cmph[m-1],  Smph[m-1],  
				 Cmph[m  ],  Smph[m  ],
				 Cmph[m+1],  Smph[m+1],
				 field) ;
      expansion_h_grad_increment(n, m, -1, k, hnm1, hnp1, cfft, cstr, 
				 Pnm1, Pnp1,
				 Cmph[m+1], -Smph[m+1],
				 Cmph[m  ], -Smph[m  ],
				 Cmph[m-1], -Smph[m-1],
				 field) ;
    }
  }
  
  return 0 ;
}
