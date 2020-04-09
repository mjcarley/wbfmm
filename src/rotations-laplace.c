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

static gint _wbfmm_rotate_H_laplace_ref_1(WBFMM_REAL *Co, gint cstro,
					  WBFMM_REAL *Ci, gint cstri,
					  gint N,
					  WBFMM_REAL *H,
					  WBFMM_REAL ph, WBFMM_REAL ch,
					  WBFMM_REAL sc)

{
  gint n, m, nu, idxi, idxo ;
  WBFMM_REAL Cmch, Smch, Cnph, Snph, Cch, Sch, Cph, Sph ;
  WBFMM_REAL Hp, Hm ;

  /*initialize recursions*/
  Cph = COS(ph) ; Sph = SIN(ph) ;
  Cch = COS(ch) ; Sch = SIN(ch) ;

  for ( n = 0 ; n <= N ; n ++ ) {
    WBFMM_REAL tr ;

    nu = 0 ; idxo = n*n ;
    m = 0 ; idxi = n*n ;
    Hp = H[wbfmm_rotation_index_numn(nu,m,n)] ;

    tr = Hp*Ci[cstri*idxi] ;
    
    Cmch = 1.0 ; Smch = 0.0 ;

    /* idxi -- ; */
    for ( m = 1 ; m <= n ; m ++ ) {
      wbfmm_cos_sin_recursion(Cmch,Smch,Cch,Sch) ;
      /* idxi += 2 ; */
      idxi = wbfmm_index_laplace_nm(n,m) ;
      
      Hp = H[wbfmm_rotation_index_numn( nu,m,n)] ;
      tr += 2.0*(Ci[cstri*(idxi+0)]*Cmch - Ci[cstri*(idxi+1)]*Smch)*Hp ;
    }
    Co[cstro*idxo] += tr ;
  }
  
  for ( n = 0 ; n <= N ; n ++ ) {
    Cnph = 1.0 ; Snph = 0.0 ;
    idxo = n*n ;

    for ( nu = 1 ; nu <= n ; nu ++ ) {
      WBFMM_REAL ti = 0.0, tr ;

      wbfmm_cos_sin_recursion(Cnph,Snph,Cph,Sph) ;
      idxo += 2 ;

      m = 0 ; idxi = n*n ;
      Hp = H[wbfmm_rotation_index_numn(nu,m,n)] ;

      tr = Hp*Ci[cstri*idxi] ;
      
      Cmch = 1.0 ; Smch = 0.0 ;
      /* idxi -- ; */
      for ( m = 1 ; m <= n ; m ++ ) {
	wbfmm_cos_sin_recursion(Cmch,Smch,Cch,Sch) ;
	/* idxi += 2 ; */
	idxi = wbfmm_index_laplace_nm(n,m) ;
	
	Hp = H[wbfmm_rotation_index_numn( nu,m,n)] ;
	Hm = H[wbfmm_rotation_index_numn(-nu,m,n)] ;
	tr += (Ci[cstri*(idxi+0)]*Cmch - Ci[cstri*(idxi+1)]*Smch)*(Hp+Hm) ;
	ti += (Ci[cstri*(idxi+0)]*Smch + Ci[cstri*(idxi+1)]*Cmch)*(Hp-Hm) ;
      }
      Co[cstro*(idxo+0)] += Cnph*tr + Snph*ti ;
      Co[cstro*(idxo+1)] += Cnph*ti - Snph*tr ;
    }
  }
  
  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_laplace_rotate_H)(WBFMM_REAL *Co, gint cstro,
						 WBFMM_REAL *Ci, gint cstri,
						 gint N, gint nq,
						 WBFMM_REAL *H,
						 WBFMM_REAL ph, WBFMM_REAL ch,
						 WBFMM_REAL sc)

{
  if ( nq == 1 ) 
    return
      _wbfmm_rotate_H_laplace_ref_nq(Co, cstro, Ci, cstri, N, 1, H, ph, ch, sc) ;

  return
    _wbfmm_rotate_H_laplace_ref_nq(Co, cstro, Ci, cstri, N, nq, H, ph, ch, sc) ;
  
  return 0 ;
}
