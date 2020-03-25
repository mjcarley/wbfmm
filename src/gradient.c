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
							  gint fstr,
							  WBFMM_REAL *work)

/*
  cstr stride by element (multiply by two to get to complex entry)
*/
		       
{
  WBFMM_REAL hn[2], hnm1[2], hnp1[2], r, th, ph, kr ;
  WBFMM_REAL Cth, Sth, *Pn, *Pnm1, *Pnp1, Cph, Sph, Cmph[64], Smph[64] ;
  gint n, m ;

  if ( fstr < 3 )
    g_error("%s: field data stride (%d) must be greater than two",
	    __FUNCTION__, fstr) ;

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

/* static gint expansion_j_grad_increment(gint n, gint m, gint sgn, */
/* 				       WBFMM_REAL jn, */
/* 				       WBFMM_REAL *cfft, gint cstr,  */
/* 				       WBFMM_REAL *Pn, */
/* 				       WBFMM_REAL Cmph, WBFMM_REAL Smph, */
/* 				       WBFMM_REAL *field) */

/* { */
/*   gint idx ; */
/*   WBFMM_REAL ar, ai, tr, ti ; */

/*   idx = wbfmm_coefficient_index_nm(n,sgn*m) ; */
/*   ar = cfft[2*idx*cstr+0] ; ai = cfft[2*idx*cstr+1] ;  */
/*   tr = ar*jn ; ti = ai*jn ; */
/*   field[0] += (Cmph*tr - sgn*Smph*ti)*Pn[m] ; */
/*   field[1] += (Cmph*ti + sgn*Smph*tr)*Pn[m] ; */
  
/*   return 0 ; */
/* } */

static gint expansion_j_grad_increment(gint n, gint m, gint sgn,
				       WBFMM_REAL k,
				       WBFMM_REAL jnm1, WBFMM_REAL jnp1,
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
  tm1r = ar*jnm1 ; tm1i = ai*jnm1 ;
  tp1r = ar*jnp1 ; tp1i = ai*jnp1 ;

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

gint WBFMM_FUNCTION_NAME(wbfmm_expansion_j_grad_evaluate)(WBFMM_REAL k,
							  WBFMM_REAL *x0,
							  WBFMM_REAL *cfft,
							  gint cstr,
							  gint N, 
							  WBFMM_REAL *xf, 
							  WBFMM_REAL *field,
							  gint fstr,
							  WBFMM_REAL *work)

/*
  cstr stride by element (multiply by two to get to complex entry)
*/
		       
{
  WBFMM_REAL jn, jnm1, jnp1, r, th, ph, kr ;
  WBFMM_REAL Cth, Sth, *Pn, *Pnm1, *Pnp1, Cph, Sph, Cmph[64], Smph[64] ;
  gint n, m ;

  if ( fstr < 3 )
    g_error("%s: field data stride (%d) must be greater than two",
	    __FUNCTION__, fstr) ;

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
  WBFMM_FUNCTION_NAME(wbfmm_bessel_j_init)(kr, &jn, &jnp1) ;
  WBFMM_FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth,
					   &(Pn[0]), &(Pnp1[0]), &(Pnp1[1])) ;

  Cmph[0] = 1.0 ; Smph[0] = 0.0 ;
  Cmph[1] = Cph ; Smph[1] = Sph ;

  /*first two terms by hand*/
  n = 0 ; 
  m = 0 ; 
  expansion_j_grad_increment(n, m,  1, k, jnm1, jnp1, cfft, cstr, 
			     Pnm1, Pnp1,
			     Cmph[m+1], -Smph[m+1],  
			     Cmph[m  ],  Smph[m  ],
			     Cmph[m+1],  Smph[m+1],
			     field) ;
  
  jnm1 = jn ;
  WBFMM_FUNCTION_NAME(wbfmm_bessel_j_recursion)(&jn, &jnp1, kr, 1) ;
  n = 1 ; 
  m = 0 ;
  memcpy(Pnm1, Pn, (n+1)*sizeof(WBFMM_REAL)) ;
  WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pn, &Pnp1,
						      1, Cth, Sth) ;
  Cmph[n+1] = Cmph[n]*Cph - Smph[n]*Sph ;
  Smph[n+1] = Smph[n]*Cph + Cmph[n]*Sph ;

  expansion_j_grad_increment(n, m,  1, k, jnm1, jnp1, cfft, cstr, 
			     Pnm1, Pnp1,
			     Cmph[m+1], -Smph[m+1],
			     Cmph[m  ],  Smph[m  ],
			     Cmph[m+1],  Smph[m+1],
			     field) ;

  m = 1 ;
  expansion_j_grad_increment(n, m,  1, k, jnm1, jnp1, cfft, cstr, 
			     Pnm1, Pnp1,
			     Cmph[m-1],  Smph[m-1],  
			     Cmph[m  ],  Smph[m  ],
			     Cmph[m+1],  Smph[m+1],
			     field) ;
  expansion_j_grad_increment(n, m, -1, k, jnm1, jnp1, cfft, cstr, 
			     Pnm1, Pnp1,
			     Cmph[m+1], -Smph[m+1],  
			     Cmph[m  ], -Smph[m  ],
			     Cmph[m-1], -Smph[m-1],
			     field) ;

  for ( n = 2 ; n <= N ; n ++ ) {
    memcpy(Pnm1, Pn, (n+1)*sizeof(WBFMM_REAL)) ;
    WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pn, &Pnp1,
							n, Cth, Sth) ;
    jnm1 = jn ;
    WBFMM_FUNCTION_NAME(wbfmm_bessel_j_recursion)(&jn, &jnp1, kr, n) ;

    Cmph[n+1] = Cmph[n]*Cph - Smph[n]*Sph ;
    Smph[n+1] = Smph[n]*Cph + Cmph[n]*Sph ;

    m = 0 ; 
    expansion_j_grad_increment(n, m,  1, k, jnm1, jnp1, cfft, cstr, 
			       Pnm1, Pnp1,
			       Cmph[m+1], -Smph[m+1],
			       Cmph[m  ],  Smph[m  ],
			       Cmph[m+1],  Smph[m+1],
			       field) ;

    for ( m = 1 ; m <= n ; m ++ ) {
      expansion_j_grad_increment(n, m,  1, k, jnm1, jnp1, cfft, cstr, 
				 Pnm1, Pnp1,
				 Cmph[m-1],  Smph[m-1],  
				 Cmph[m  ],  Smph[m  ],
				 Cmph[m+1],  Smph[m+1],
				 field) ;
      expansion_j_grad_increment(n, m, -1, k, jnm1, jnp1, cfft, cstr, 
				 Pnm1, Pnp1,
				 Cmph[m+1], -Smph[m+1],
				 Cmph[m  ], -Smph[m  ],
				 Cmph[m-1], -Smph[m-1],
				 field) ;
    }
  }
#if 0
#endif
  
  return 0 ;
}
