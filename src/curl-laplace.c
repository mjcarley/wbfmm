/* This file is part of WBFMM, a Wide-Band Fast Multipole Method code
 *
 * Copyright (C) 2019, 2021 Michael Carley
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

gint
WBFMM_FUNCTION_NAME(wbfmm_laplace_expansion_local_curl_evaluate)(WBFMM_REAL *x0,
								 WBFMM_REAL
								 *cfft,
								 gint cstr, 
								 gint N,
								 gint nq,
								 WBFMM_REAL *xf,
								 WBFMM_REAL
								 *field,
								 gint fstr,
								 WBFMM_REAL
								 *work)

{
  WBFMM_REAL r, th, ph, rn, anm, b1, b2, cr, ci ;
  WBFMM_REAL Cth, Sth, *Pn, *Pnm1 ;
  WBFMM_REAL *Cmph, *Smph, Rnmm1, Rnm, Rnmp1 ;
  WBFMM_REAL ddxr, ddxi, ddyr, ddyi, ddzr, ddzi ;
  gint n, m, idx ;

  /*
    fstr is ignored: the curl based on the first three components of
    the source is placed into the first three components of f
   */

  
  if ( nq < 3 )
    g_error("%s: not enough source components (%d) for curl calculation",
	    __FUNCTION__, nq) ;
  
  /* if ( fstr < 3 ) */
  /*   g_error("%s: field data stride (%d) must be greater than two", */
  /* 	    __FUNCTION__, fstr) ; */

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
  rn = 1.0 ;
  idx = n*n ;
  Cmph[n+1] = Cmph[n]*Cmph[1] - Smph[n]*Smph[1] ;
  Smph[n+1] = Smph[n]*Cmph[1] + Cmph[n]*Smph[1] ;

  anm = SQRT((WBFMM_REAL)(2*n+1)/(2*n-1)*(n-m)*(n+m)) ;
  Rnm = rn*Pnm1[m]*anm ;
  
  ddzr = Rnm ;
  
  cr = cfft[cstr*idx+1] ;
  field[0] -= ddzr*cr ;
  cr = cfft[cstr*idx+0] ;
  field[1] += ddzr*cr ;
  
  m = 1 ; 
  idx = wbfmm_index_laplace_nm(n,m) ;
  anm = SQRT((WBFMM_REAL)(2*n+1)/(2*n-1)*(n-m)*(n+m)) ;
  b1  = SQRT((WBFMM_REAL)(2*n+1)/(2*n-1)*(n-m)*(n-m-1)) ;
  b2  = SQRT((WBFMM_REAL)(2*n+1)/(2*n-1)*(n+m)*(n+m-1)) ;
  Rnmm1 = rn*Pnm1[m-1]*b2 ;
  Rnm   = rn*Pnm1[m+0]*anm*2.0 ;
  Rnmp1 = rn*Pnm1[m+1]*b1 ;

  ddxr = -Rnmp1*Cmph[m+1] + Rnmm1*Cmph[m-1] ;
  ddxi = +Rnmp1*Smph[m+1] - Rnmm1*Smph[m-1] ;

  ddyr = -Rnmp1*Smph[m+1] - Rnmm1*Smph[m-1] ;
  ddyi = -Rnmp1*Cmph[m+1] - Rnmm1*Cmph[m-1] ;

  ddzr = +Rnm*Cmph[m+0] ;
  ddzi = -Rnm*Smph[m+0] ;
  
  cr = cfft[cstr*(idx+0)+2] ; ci = cfft[cstr*(idx+1)+2] ;
  field[0] += ddyr*cr + ddyi*ci ;
  cr = cfft[cstr*(idx+0)+1] ; ci = cfft[cstr*(idx+1)+1] ;
  field[0] -= ddzr*cr + ddzi*ci ;

  cr = cfft[cstr*(idx+0)+0] ; ci = cfft[cstr*(idx+1)+0] ;
  field[1] += ddzr*cr + ddzi*ci ;
  cr = cfft[cstr*(idx+0)+2] ; ci = cfft[cstr*(idx+1)+2] ;
  field[1] -= ddxr*cr + ddxi*ci ;
  
  cr = cfft[cstr*(idx+0)+1] ; ci = cfft[cstr*(idx+1)+1] ;
  field[2] += ddxr*cr + ddxi*ci ;
  cr = cfft[cstr*(idx+0)+0] ; ci = cfft[cstr*(idx+1)+0] ;
  field[2] -= ddyr*cr + ddyi*ci ;
  
  for ( n = 2 ; n <= N ; n ++ ) {
    rn *= r ;
    WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pnm1, &Pn,
							n-1, Cth, Sth) ;
    Cmph[n+1] = Cmph[n]*Cmph[1] - Smph[n]*Smph[1] ;
    Smph[n+1] = Smph[n]*Cmph[1] + Cmph[n]*Smph[1] ;

    m = 0 ; 
    idx = n*n ;
    anm = SQRT((WBFMM_REAL)(2*n+1)/(2*n-1)*(n-m)*(n+m)) ;
    b1  = SQRT((WBFMM_REAL)(2*n+1)/(2*n-1)*(n-m)*(n-m-1)) ;
    b2  = SQRT((WBFMM_REAL)(2*n+1)/(2*n-1)*(n+m)*(n+m-1)) ;
    Rnmm1 = rn*Pnm1[m+1]*b2 ;
    Rnm = rn*Pnm1[m]*anm ;
    Rnmp1 = rn*Pnm1[m+1]*b1 ;

    ddxr = -Rnmp1*Cmph[m+1] ;
    ddyr = -Rnmp1*Smph[m+1] ;
    ddzr = Rnm ;

    cr = cfft[cstr*idx+2] ;
    field[0] += ddyr*cr ;
    cr = cfft[cstr*idx+1] ;
    field[0] -= ddzr*cr ;

    cr = cfft[cstr*idx+0] ;
    field[1] += ddzr*cr ;
    cr = cfft[cstr*idx+2] ;
    field[1] -= ddxr*cr ;
    
    cr = cfft[cstr*idx+1] ;
    field[2] += ddxr*cr ;
    cr = cfft[cstr*idx+0] ;
    field[2] -= ddyr*cr ;

    for ( m = 1 ; m <= n ; m ++ ) {
      idx = wbfmm_index_laplace_nm(n,m) ;
      anm = SQRT((WBFMM_REAL)(2*n+1)/(2*n-1)*(n-m)*(n+m)) ;
      b1  = SQRT((WBFMM_REAL)(2*n+1)/(2*n-1)*(n-m)*(n-m-1)) ;
      b2  = SQRT((WBFMM_REAL)(2*n+1)/(2*n-1)*(n+m)*(n+m-1)) ;
      Rnmm1 = rn*Pnm1[m-1]*b2 ;
      Rnm   = rn*Pnm1[m+0]*anm*2.0 ;
      Rnmp1 = rn*Pnm1[m+1]*b1 ;
      ddxr = -Rnmp1*Cmph[m+1] + Rnmm1*Cmph[m-1] ;
      ddxi = +Rnmp1*Smph[m+1] - Rnmm1*Smph[m-1] ;
      
      ddyr = -Rnmp1*Smph[m+1] - Rnmm1*Smph[m-1] ;
      ddyi = -Rnmp1*Cmph[m+1] - Rnmm1*Cmph[m-1] ;
      
      ddzr = +Rnm*Cmph[m+0] ;
      ddzi = -Rnm*Smph[m+0] ;

      cr = cfft[cstr*(idx+0)+2] ; ci = cfft[cstr*(idx+1)+2] ;
      field[0] += ddyr*cr + ddyi*ci ;
      cr = cfft[cstr*(idx+0)+1] ; ci = cfft[cstr*(idx+1)+1] ;
      field[0] -= ddzr*cr + ddzi*ci ;

      cr = cfft[cstr*(idx+0)+0] ; ci = cfft[cstr*(idx+1)+0] ;
      field[1] += ddzr*cr + ddzi*ci ;
      cr = cfft[cstr*(idx+0)+2] ; ci = cfft[cstr*(idx+1)+2] ;
      field[1] -= ddxr*cr + ddxi*ci ;
  
      cr = cfft[cstr*(idx+0)+1] ; ci = cfft[cstr*(idx+1)+1] ;
      field[2] += ddxr*cr + ddxi*ci ;
      cr = cfft[cstr*(idx+0)+0] ; ci = cfft[cstr*(idx+1)+0] ;
      field[2] -= ddyr*cr + ddyi*ci ;
    }
  }
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_laplace_expansion_curl_evaluate)(WBFMM_REAL *x0,
								WBFMM_REAL
								*cfft,
								gint cstr, 
								gint N, gint nq,
								WBFMM_REAL *xf,
								WBFMM_REAL
								*field,
								gint fstr,
								WBFMM_REAL
								*work)

/*
  z derivative of singular functions:
  \partial S_{n}^{m}/\partial z = -a_{n}^{m}S_{n+1}^{m}
  
  a_{n}^{m} = 
  \left
  (2n+1)(n+m+1)(n-m+1)/(2n+3)
  \right]^{1/2}
 */
  
{
  WBFMM_REAL r, th, ph, rn, anm, b1, b2, Snmp1, Snmm1, Snm ;
  WBFMM_REAL Cth, Sth, *Pn, *Pnp1 ;
  WBFMM_REAL *Cmph, *Smph, cr, ci ;
  gint n, m, idx, i ;

  g_assert_not_reached() ; /*not yet implemented, the code below is
			     for gradients*/
  
  if ( fstr < 3 )
    g_error("%s: field data stride (%d) must be greater than two",
	    __FUNCTION__, fstr) ;
  
  Pn = &(work[0]) ; Pnp1 = &(Pn[N+2]) ;
  Cmph = &(Pnp1[N+3]) ; Smph = &(Cmph[N+3]) ;

  WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(x0, xf, &r, &th, &ph) ;
  Cth = COS(th) ; Sth = SIN(th) ; 

  /*initialize recursions*/
  WBFMM_FUNCTION_NAME(wbfmm_legendre_init)(Cth, Sth,
					   &(Pn[0]), &(Pnp1[0]), &(Pnp1[1])) ;
  Cmph[0] = 1.0 ; Smph[0] = 0.0 ;
  Cmph[1] = COS(ph) ; Smph[1] = SIN(ph) ;

  /*first two terms by hand*/
  n = 0 ; 
  m = 0 ;
  rn = 1.0/r/r ;
  idx = n*n ;
  anm = SQRT((WBFMM_REAL)(2*n+1)/(2*n+3)*(n+m+1)*(n-m+1)) ;
  b1  = SQRT((WBFMM_REAL)(2*n+1)/(2*n+3)*(n+m+1)*(n+m+2)) ;
  Snm   = rn*Pnp1[m+0]*anm ;
  Snmp1 = rn*Pnp1[m+1]*b1 ;
  for ( i = 0 ; i < nq ; i ++ ) {
    cr = cfft[cstr*idx+i] ;
    field[fstr*i+0] -= cr*Snmp1*Cmph[1] ;
    field[fstr*i+1] -= cr*Snmp1*Smph[1] ;
    field[fstr*i+2] -= cr*Snm ;
  }

  if ( N == 0 ) return 0 ;
  
  n = 1 ; 
  m = 0 ; 
  rn /= r ;
  Cmph[n+1] = Cmph[n]*Cmph[1] - Smph[n]*Smph[1] ;
  Smph[n+1] = Smph[n]*Cmph[1] + Cmph[n]*Smph[1] ;
  WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pn, &Pnp1,
						      n, Cth, Sth) ;
  idx = n*n ;
  anm = SQRT((WBFMM_REAL)(2*n+1)/(2*n+3)*(n+m+1)*(n-m+1)) ;
  b1  = SQRT((WBFMM_REAL)(2*n+1)/(2*n+3)*(n+m+1)*(n+m+2)) ;
  Snm   = rn*Pnp1[m+0]*anm ;
  Snmp1 = rn*Pnp1[m+1]*b1 ;
  for ( i = 0 ; i < nq ; i ++ ) {
    cr = cfft[cstr*idx+i] ;
    field[fstr*i+0] -= cr*Snmp1*Cmph[1] ;
    field[fstr*i+1] -= cr*Snmp1*Smph[1] ;
    field[fstr*i+2] -= cr*Snm ;
  }

  m = 1 ; 
  idx = wbfmm_index_laplace_nm(n,m) ;
  anm = SQRT((WBFMM_REAL)(2*n+1)/(2*n+3)*(n+m+1)*(n-m+1)) ;
  b1  = SQRT((WBFMM_REAL)(2*n+1)/(2*n+3)*(n+m+1)*(n+m+2)) ;
  b2  = SQRT((WBFMM_REAL)(2*n+1)/(2*n+3)*(n-m+1)*(n-m+2)) ;
  Snmm1 = rn*Pnp1[m-1]*b2 ;
  Snm   = rn*Pnp1[m+0]*anm*2.0 ;
  Snmp1 = rn*Pnp1[m+1]*b1 ;
  for ( i = 0 ; i < nq ; i ++ ) {
    cr = cfft[cstr*(idx+0)+i] ; ci = cfft[cstr*(idx+1)+i] ;
    
    field[fstr*i+0] += Snmm1*(cr*Cmph[m-1] - ci*Smph[m-1]) ;      
    field[fstr*i+0] -= Snmp1*(cr*Cmph[m+1] - ci*Smph[m+1]) ;      

    field[fstr*i+1] -= Snmm1*(ci*Cmph[m-1] + cr*Smph[m-1]) ;
    field[fstr*i+1] -= Snmp1*(ci*Cmph[m+1] + cr*Smph[m+1]) ;
    
    field[fstr*i+2] -= Snm  *(cr*Cmph[m+0] - ci*Smph[m+0]) ;
  }

  for ( n = 2 ; n <= N ; n ++ ) {
    rn /= r ;
    WBFMM_FUNCTION_NAME(wbfmm_legendre_recursion_array)(&Pn, &Pnp1,
							n, Cth, Sth) ;
    Cmph[n+1] = Cmph[n]*Cmph[1] - Smph[n]*Smph[1] ;
    Smph[n+1] = Smph[n]*Cmph[1] + Cmph[n]*Smph[1] ;

    m = 0 ; 
    idx = n*n ;
    anm = SQRT((WBFMM_REAL)(2*n+1)/(2*n+3)*(n+m+1)*(n-m+1)) ;
    b1  = SQRT((WBFMM_REAL)(2*n+1)/(2*n+3)*(n+m+1)*(n+m+2)) ;
    Snm   = rn*Pnp1[m+0]*anm ;
    Snmp1 = rn*Pnp1[m+1]*b1 ;
    for ( i = 0 ; i < nq ; i ++ ) {
      cr = cfft[cstr*idx+i] ;
      field[fstr*i+0] -= cr*Snmp1*Cmph[1] ;
      field[fstr*i+1] -= cr*Snmp1*Smph[1] ;    
      field[fstr*i+2] -= cr*Snm ;
    }
    
    for ( m = 1 ; m <= n ; m ++ ) {
      idx = wbfmm_index_laplace_nm(n,m) ;
      anm = SQRT((WBFMM_REAL)(2*n+1)/(2*n+3)*(n+m+1)*(n-m+1)) ;
      b1  = SQRT((WBFMM_REAL)(2*n+1)/(2*n+3)*(n+m+1)*(n+m+2)) ;
      b2  = SQRT((WBFMM_REAL)(2*n+1)/(2*n+3)*(n-m+1)*(n-m+2)) ;
      Snmm1 = rn*Pnp1[m-1]*b2 ;
      Snm   = rn*Pnp1[m+0]*anm*2.0 ;
      Snmp1 = rn*Pnp1[m+1]*b1 ;
      for ( i = 0 ; i < nq ; i ++ ) {
	cr = cfft[cstr*(idx+0)+i] ; ci = cfft[cstr*(idx+1)+i] ;

	field[fstr*i+0] += Snmm1*(cr*Cmph[m-1] - ci*Smph[m-1]) ;      
	field[fstr*i+0] -= Snmp1*(cr*Cmph[m+1] - ci*Smph[m+1]) ;      
    
	field[fstr*i+1] -= Snmm1*(ci*Cmph[m-1] + cr*Smph[m-1]) ;
	field[fstr*i+1] -= Snmp1*(ci*Cmph[m+1] + cr*Smph[m+1]) ;
	
	field[fstr*i+2] -= Snm  *(cr*Cmph[m+0] - ci*Smph[m+0]) ;
      }
    }
  }
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_laplace_field_curl)(WBFMM_REAL *xs,
						   gint xstride,
						   WBFMM_REAL *src,
						   gint sstride,
						   gint nq,
						   WBFMM_REAL *normals,
						   gint nstr,
						   WBFMM_REAL *dipoles,
						   gint dstr,
						   gint nsrc,
						   WBFMM_REAL *xf,
						   WBFMM_REAL *field,
						   gint fstr)

{
  gint i ;
  WBFMM_REAL r, th, ph, nR[3] ;

  if ( src == NULL && normals == NULL && dipoles == NULL ) return 0 ;

  if ( nq < 3 )
    g_error("%s: not enough source components (%d) for curl calculation",
	    __FUNCTION__, nq) ;
  
  if ( fstr < 3 )
    g_error("%s: field data stride (%d) must be greater than two",
	    __FUNCTION__, fstr) ;

  g_assert(sstride >= nq) ;
  
  if ( normals != NULL && dipoles == NULL )
    g_error("%s: normals specified but no dipole strengths (dipoles == NULL)",
	    __FUNCTION__) ;

  if ( normals == NULL && dipoles == NULL ) {
    for ( i = 0 ; i < nsrc ; i ++ ) {
      WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(&(xs[i*xstride]), xf, 
							&r, &th, &ph) ;
      if ( r > WBFMM_LOCAL_CUTOFF_RADIUS ) {
	nR[0] = (xf[0] - xs[i*xstride+0])/r/r/r*0.25*M_1_PI ;
	nR[1] = (xf[1] - xs[i*xstride+1])/r/r/r*0.25*M_1_PI ;
	nR[2] = (xf[2] - xs[i*xstride+2])/r/r/r*0.25*M_1_PI ;

	field[0] -= src[i*sstride+2]*nR[1] - src[i*sstride+1]*nR[2] ;
	field[1] -= src[i*sstride+0]*nR[2] - src[i*sstride+2]*nR[0] ;
	field[2] -= src[i*sstride+1]*nR[0] - src[i*sstride+0]*nR[1] ;
      }
    }

    return 0 ;
  }

  g_assert_not_reached() ;
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_tree_laplace_box_local_curl)(wbfmm_tree_t *t,
							    guint level,
							    guint b,
							    WBFMM_REAL *x,
							    WBFMM_REAL *f,
							    gint fstr,
							    WBFMM_REAL *src,
							    gint sstr,
							    WBFMM_REAL
							    *normals,
							    gint nstr,
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
  
  WBFMM_FUNCTION_NAME(wbfmm_laplace_expansion_local_curl_evaluate)(xb, C, 8*nq,
  								   t->order_r[level],
  								   nq, x, f,
  								   fstr,
  								   work) ;
  
  if ( !eval_neighbours ) return 0 ;

  if ( src == NULL && normals == NULL && d == NULL ) return 0 ;
  
  if ( normals != NULL && d == NULL ) {
    g_error("%s: normals specified but no dipole strengths (d == NULL)",
	    __FUNCTION__) ;
  }

  /*add the contribution from sources in neighbour boxes*/
  nnbr = wbfmm_box_neighbours(level, b, neighbours) ;
  g_assert(nnbr >= 0 && nnbr < 28) ;

  if ( normals == NULL && d == NULL ) {
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

	  f[0] -= src[idx*sstr+2]*nR[1] - src[idx*sstr+1]*nR[2] ;
	  f[1] -= src[idx*sstr+0]*nR[2] - src[idx*sstr+2]*nR[0] ;
	  f[2] -= src[idx*sstr+1]*nR[0] - src[idx*sstr+0]*nR[1] ;
	}
      }
    }
    
    return 0 ;
  }

  g_assert_not_reached() ;
  
  if ( src == NULL && normals != NULL ) {
    /*dipoles only*/
    /* g_assert_not_reached() ; */
    WBFMM_REAL th, ph, nr ;
    
    for ( i = 0 ; i < nnbr ; i ++ ) {
      box = boxes[neighbours[i]] ;
      for ( j = 0 ; j < box.n ; j ++ ) {
	idx = t->ip[box.i+j] ;
	xs = wbfmm_tree_point_index(t, idx) ;

	WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(xs, x, &r, &th, &ph) ;
	if ( r > WBFMM_LOCAL_CUTOFF_RADIUS ) {
	  nr =
	    (x[0] - xs[0])*normals[idx*nstr+0] +
	    (x[1] - xs[1])*normals[idx*nstr+1] + 
	    (x[2] - xs[2])*normals[idx*nstr+2] ;
	  nr /= 4.0*M_PI*r*r*r ;
	  for ( k = 0 ; k < nq ; k ++ ) f[k] += d[idx*dstr+k]*nr ;
	}
      }
      
    } 

    return 0 ;
  }
  
  if ( src != NULL && normals != NULL ) {
    /*sources and dipoles*/
    WBFMM_REAL th, ph, nr, g ;
    
    for ( i = 0 ; i < nnbr ; i ++ ) {
      box = boxes[neighbours[i]] ;
      for ( j = 0 ; j < box.n ; j ++ ) {
	idx = t->ip[box.i+j] ;
	xs = wbfmm_tree_point_index(t, idx) ;

	WBFMM_FUNCTION_NAME(wbfmm_cartesian_to_spherical)(xs, x, &r, &th, &ph) ;
	if ( r > WBFMM_LOCAL_CUTOFF_RADIUS ) {
	  nr =
	    (x[0] - xs[0])*normals[idx*nstr+0] +
	    (x[1] - xs[1])*normals[idx*nstr+1] + 
	    (x[2] - xs[2])*normals[idx*nstr+2] ;
	  g = 0.25*M_1_PI/r ;
	  /* nr /= 4.0*M_PI*r*r*r ; */
	  /* nr *= g/r/r ; /\* 4.0*M_PI*r*r*r ; *\/ */
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

