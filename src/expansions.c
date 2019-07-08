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

#ifdef _HAVE_CONFIG_H_
#include <config.h>
#endif /*_HAVE_CONFIG_H_*/

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

gint FUNCTION_NAME(wbfmm_expansion_h_cfft)(WBFMM_REAL k, gint N, 
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

  FUNCTION_NAME(wbfmm_cartesian_to_spherical)(x0, xs, &r, &th, &ph) ;
  Cth = COS(th) ; Sth = SIN(th) ; kr = k*r ;
  Cph = COS(ph) ; Sph = SIN(ph) ; 

  /*initialize recursions*/
  FUNCTION_NAME(wbfmm_bessel_j_init)(kr, &jnm1, &jn) ;
  FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth, &(Pnm1[0]), &(Pn[0]), &(Pn[1])) ;
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
    FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pnm1, &Pn, n-1, Cth, Sth) ;
    FUNCTION_NAME(wbfmm_bessel_j_recursion)(&jnm1, &jn, kr, n-1) ;
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

gint FUNCTION_NAME(wbfmm_expansion_h_evaluate)(WBFMM_REAL k, WBFMM_REAL *x0,
						WBFMM_REAL *cfft, gint cstr,
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

  FUNCTION_NAME(wbfmm_cartesian_to_spherical)(x0, xf, &r, &th, &ph) ;
  Cth = COS(th) ; Sth = SIN(th) ; 
  Cph = COS(ph) ; Sph = SIN(ph) ; 
  kr = k*r ;

  /*initialize recursions*/
  FUNCTION_NAME(wbfmm_bessel_h_init)(kr, hnm1, hn) ;
  FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth, &(Pnm1[0]), &(Pn[0]), &(Pn[1])) ;
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
    FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pnm1, &Pn, n-1, Cth, Sth) ;
    FUNCTION_NAME(wbfmm_bessel_h_recursion)(hnm1, hn, kr, n-1) ;
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
  /* field[0] += jn*Pn[m]*cfft[2*idx*cstr+0]*Cmph ;  */
  /* field[1] += jn*Pn[m]*cfft[2*idx*cstr+1]*Smph*sgn ;  */
  ar = cfft[2*idx*cstr+0] ; ai = cfft[2*idx*cstr+1] ; 
  tr = ar*jn ; ti = ai*jn ;
  field[0] += (Cmph*tr - sgn*Smph*ti)*Pn[m] ;
  field[1] += (Cmph*ti + sgn*Smph*tr)*Pn[m] ;
  
  return 0 ;
}

gint FUNCTION_NAME(wbfmm_expansion_j_evaluate)(WBFMM_REAL k, WBFMM_REAL *x0,
						WBFMM_REAL *cfft, gint cstr,
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

  FUNCTION_NAME(wbfmm_cartesian_to_spherical)(x0, xf, &r, &th, &ph) ;
  Cth = COS(th) ; Sth = SIN(th) ; 
  Cph = COS(ph) ; Sph = SIN(ph) ; 
  kr = k*r ;

  /*initialize recursions*/
  FUNCTION_NAME(wbfmm_bessel_j_init)(kr, &jnm1, &jn) ;
  FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth, &(Pnm1[0]), &(Pn[0]), &(Pn[1])) ;
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
    FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pnm1, &Pn, n-1, Cth, Sth) ;
    FUNCTION_NAME(wbfmm_bessel_j_recursion)(&jnm1, &jn, kr, n-1) ;
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

static gint expansion_dipole_h_increment_cfft(gint n, gint m, gint sgn,
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

  Rnm = jn*Pn[m] ;

  /*S_{n+1}^{m+1}*/
  idx = wbfmm_coefficient_index_nm(n+1,sgn*m+1) ;
  ab = FUNCTION_NAME(recursion_bnm)(n+1, -sgn*m-1) ;
  cfft[2*idx*cstr+0] += ab*Rnm*(fm[0]*Cmph[m] + sgn*fm[1]*Smph[m]) ;
  cfft[2*idx*cstr+1] += ab*Rnm*(fm[1]*Cmph[m] - sgn*fm[0]*Smph[m]) ;

  /*S_{n-1}^{m+1}*/
  idx = wbfmm_coefficient_index_nm(n-1,sgn*m+1) ;
  ab = FUNCTION_NAME(recursion_bnm)(n, sgn*m) ;
  cfft[2*idx*cstr+0] -= ab*Rnm*(fm[0]*Cmph[m] + sgn*fm[1]*Smph[m]) ;
  cfft[2*idx*cstr+1] -= ab*Rnm*(fm[1]*Cmph[m] - sgn*fm[0]*Smph[m]) ;

  /*S_{n+1}^{m-1}*/
  idx = wbfmm_coefficient_index_nm(n+1,sgn*m-1) ;
  ab = FUNCTION_NAME(recursion_bnm)(n+1, sgn*m-1) ;
  cfft[2*idx*cstr+0] += ab*Rnm*(fp[0]*Cmph[m] + sgn*fp[1]*Smph[m]) ;
  cfft[2*idx*cstr+1] += ab*Rnm*(fp[1]*Cmph[m] - sgn*fp[0]*Smph[m]) ;

  /*S_{n-1}^{m-1}*/
  idx = wbfmm_coefficient_index_nm(n-1,sgn*m-1) ;
  ab = FUNCTION_NAME(recursion_bnm)(n, -sgn*m) ;
  cfft[2*idx*cstr+0] -= ab*Rnm*(fp[0]*Cmph[m] + sgn*fp[1]*Smph[m]) ;
  cfft[2*idx*cstr+1] -= ab*Rnm*(fp[1]*Cmph[m] - sgn*fp[0]*Smph[m]) ;

  /*S_{n-1}^{m}*/
  idx = wbfmm_coefficient_index_nm(n-1,sgn*m) ;
  ab = FUNCTION_NAME(recursion_anm)(n-1, sgn*m) ;
  cfft[2*idx*cstr+0] += ab*Rnm*(fz[0]*Cmph[m] + sgn*fz[1]*Smph[m]) ;
  cfft[2*idx*cstr+1] += ab*Rnm*(fz[1]*Cmph[m] - sgn*fz[0]*Smph[m]) ;

  /*S_{n+1}^{m}*/
  idx = wbfmm_coefficient_index_nm(n+1,sgn*m) ;
  ab = FUNCTION_NAME(recursion_anm)(n, sgn*m) ;
  cfft[2*idx*cstr+0] -= ab*Rnm*(fz[0]*Cmph[m] + sgn*fz[1]*Smph[m]) ;
  cfft[2*idx*cstr+1] -= ab*Rnm*(fz[1]*Cmph[m] - sgn*fz[0]*Smph[m]) ;
  
  return 0 ;
}

gint FUNCTION_NAME(wbfmm_expansion_dipole_h_cfft)(WBFMM_REAL k, gint N, 
						  WBFMM_REAL *x0,
						  WBFMM_REAL *xs,
						  WBFMM_REAL *fx,
						  WBFMM_REAL *fy,
						  WBFMM_REAL *fz,
						  WBFMM_REAL *cfft, gint cstr,
						  WBFMM_REAL *work)

/*
  workspace size: 4*(2*N+1)
*/

{
  WBFMM_REAL jn, jnm1, r, th, ph, kr, fm[2], fp[2] ;
  WBFMM_REAL Cth, Sth, *Pn, *Pnm1, Cph, Sph, Cmph[64], Smph[64] ;
  gint n, m ;

  /*G&D (2.17) combined with derivatives (3.2)--(3.7)*/

  /*f_{\pm} = (k/2) \mathbf{f}.(i_{x} \pm j i_{y})*/
  fm[0] = 0.5*k*(fx[0] + fy[1]) ; fm[1] = 0.5*k*(fx[1] - fy[0]) ;
  fp[0] = 0.5*k*(fx[0] - fy[1]) ; fp[1] = 0.5*k*(fx[1] + fy[0]) ;
  
  Pnm1 = &(work[0]) ; Pn = &(Pnm1[2*(2*N+1)]) ;

  FUNCTION_NAME(wbfmm_cartesian_to_spherical)(x0, xs, &r, &th, &ph) ;
  Cth = COS(th) ; Sth = SIN(th) ; kr = k*r ;
  Cph = COS(ph) ; Sph = SIN(ph) ; 

  /*initialize recursions*/
  FUNCTION_NAME(wbfmm_bessel_j_init)(kr, &jnm1, &jn) ;
  FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth, &(Pnm1[0]), &(Pn[0]), &(Pn[1])) ;
  Cmph[0] = 1.0 ; Smph[0] = 0.0 ;
  Cmph[1] = Cph ; Smph[1] = Sph ;
  
  /*first entries are done by hand*/
  n = 0 ; 
  m = 0 ; 
  expansion_dipole_h_increment_cfft(n, m,  1, jnm1, cfft, cstr, Pnm1,
				    Cmph, Smph, fm, fp, fz) ;
  
  n = 1 ; 
  m = 0 ; 
  expansion_dipole_h_increment_cfft(n, m,  1, jn, cfft, cstr, Pn,
				    Cmph, Smph, fm, fp, fz) ;
  
  m = 1 ; 
  expansion_dipole_h_increment_cfft(n, m,  1, jn, cfft, cstr, Pn,
				    Cmph, Smph, fm, fp, fz) ;
  expansion_dipole_h_increment_cfft(n, m, -1, jn, cfft, cstr, Pn,
				    Cmph, Smph, fm, fp, fz) ;
  
  for ( n = 2 ; n <= N ; n ++ ) {
    FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pnm1, &Pn, n-1, Cth, Sth) ;
    FUNCTION_NAME(wbfmm_bessel_j_recursion)(&jnm1, &jn, kr, n-1) ;
    Cmph[n] = Cmph[n-1]*Cph - Smph[n-1]*Sph ;
    Smph[n] = Smph[n-1]*Cph + Cmph[n-1]*Sph ;
    
    m = 0 ; 
    expansion_dipole_h_increment_cfft(n, m,  1, jn, cfft, cstr, Pn,
				      Cmph, Smph, fm, fp, fz) ;
    
    for ( m = 1 ; m <= n ; m ++ ) {
      expansion_dipole_h_increment_cfft(n, m,  1, jn, cfft, cstr, Pn,
					Cmph, Smph, fm, fp, fz) ;
      expansion_dipole_h_increment_cfft(n, m, -1, jn, cfft, cstr, Pn,
					Cmph, Smph, fm, fp, fz) ;
    }
  }

  return 0 ;
}
