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

/*
  functions for child-parent, parent-child, and other shifts in tree
  calculations
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <math.h>
#include <string.h>

#include <glib.h>

#include <wbfmm.h>

#include "wbfmm-private.h"

#include <stdio.h>

/* /\*table of \cos m\pi/4 for rotations on upward pass*\/ */
/* extern WBFMM_REAL CmPI_4[] ; */
/* /\*table of \cos n\pi/2 for rotations*\/ */
/* extern WBFMM_REAL CnPI_2[] ; */

/* #define cos_n_PI_4(_n) (CmPI_4[(_n)%8]) */
/* #define sin_n_PI_4(_n) (CmPI_4[((_n)+6)%8]) */
/* #define cos_n_PI_2(_n) (CnPI_2[(_n)%4]) */
/* #define sin_n_PI_2(_n) (CnPI_2[((_n)+3)%4]) */

static inline void increment_cfft_cp_real_bw(WBFMM_REAL *E0,
					     WBFMM_REAL *Ci,
					     WBFMM_REAL H03p, WBFMM_REAL H03m,
					     gint nq,
					     WBFMM_REAL *Cp)

{
  WBFMM_REAL H03 ;
  gint i ;
  
  H03 = H03p + H03m ;
  for ( i = 0 ; i < nq ; i ++ ) {
    Cp[i] +=
      ((Ci[0*nq+0*nq+i] + Ci[0*nq+3*nq+i] +
	Ci[0*nq+1*nq+i] + Ci[0*nq+2*nq+i])*E0[0] +
       (Ci[8*nq+0*nq+i] + Ci[8*nq+3*nq+i] -
	Ci[8*nq+1*nq+i] + Ci[8*nq+2*nq+i])*E0[1])*H03 ;
  }

  return ;
}

static inline void increment_cfft_cp_complex_bw(WBFMM_REAL *E0,
						WBFMM_REAL *E0ph,
						WBFMM_REAL *E1ph,
						WBFMM_REAL *Ci,
						WBFMM_REAL H03p,
						WBFMM_REAL H03m,
						gint nq,
						WBFMM_REAL *Cp)

{
  WBFMM_REAL sH03, dH03, t1r, t1i, t2r, t2i ;
  WBFMM_REAL E0s03[2], E0d03[2] ;
  gint i ;
  
  sH03 = H03p + H03m ; dH03 = H03p - H03m ;

  E0s03[0] = E0[0]*sH03 ; E0s03[1] = E0[1]*sH03 ;
  E0d03[0] = E0[0]*dH03 ; E0d03[1] = E0[1]*dH03 ;
  
  for ( i = 0 ; i < nq ; i ++ ) {
    t1r = Ci[0*nq+0*nq+i]*E0s03[0] + Ci[8*nq+0*nq+i]*E0s03[1] ;
    t1i = Ci[8*nq+0*nq+i]*E0d03[0] - Ci[0*nq+0*nq+i]*E0d03[1] ;

    t2r = Ci[0*nq+2*nq+i]*E0s03[0] + Ci[8*nq+2*nq+i]*E0s03[1] ;
    t2i = Ci[8*nq+2*nq+i]*E0d03[0] - Ci[0*nq+2*nq+i]*E0d03[1] ;

    Cp[0*nq+i] += E0ph[0]*(t2r + t1r) - E0ph[1]*(t2i - t1i) ;
    Cp[8*nq+i] += E0ph[0]*(t2i + t1i) + E0ph[1]*(t2r - t1r) ;

    t1r = Ci[0*nq+1*nq+i]*E0s03[0] - Ci[8*nq+1*nq+i]*E0s03[1] ;
    t1i = Ci[8*nq+1*nq+i]*E0d03[0] + Ci[0*nq+1*nq+i]*E0d03[1] ;

    t2r = Ci[0*nq+3*nq+i]*E0s03[0] + Ci[8*nq+3*nq+i]*E0s03[1] ;
    t2i = Ci[8*nq+3*nq+i]*E0d03[0] - Ci[0*nq+3*nq+i]*E0d03[1] ;

    Cp[0*nq+i] += E1ph[0]*(t2r + t1r) - E1ph[1]*(t2i - t1i) ;
    Cp[8*nq+i] += E1ph[0]*(t2i + t1i) + E1ph[1]*(t2r - t1r) ;
  }

  return ;
}

static inline void increment_buf_cp_real_bw(WBFMM_REAL *E0, WBFMM_REAL *E1,
					    WBFMM_REAL *C,
					    WBFMM_REAL H03p, WBFMM_REAL H03m,
					    gint nq,
					    WBFMM_REAL *buf)

{
  WBFMM_REAL H03 ;
  gint i ;
  
  H03 = H03p + H03m ;

  for ( i = 0 ; i < nq ; i ++ ) {
    buf[0*nq+i] += (C[0*nq+i]*E0[0] - C[8*nq+0*nq+i]*E0[1])*H03 ;
    buf[7*nq+i] += (C[7*nq+i]*E0[0] - C[8*nq+7*nq+i]*E0[1])*H03 ;

    buf[1*nq+i] += (C[1*nq+i]*E1[0] - C[8*nq+1*nq+i]*E1[1])*H03 ;
    buf[6*nq+i] += (C[6*nq+i]*E1[0] - C[8*nq+6*nq+i]*E1[1])*H03 ;

    buf[2*nq+i] += (C[2*nq+i]*E0[0] + C[8*nq+2*nq+i]*E0[1])*H03 ;
    buf[5*nq+i] += (C[5*nq+i]*E0[0] + C[8*nq+5*nq+i]*E0[1])*H03 ;

    buf[3*nq+i] += (C[3*nq+i]*E1[0] + C[8*nq+3*nq+i]*E1[1])*H03 ;
    buf[4*nq+i] += (C[4*nq+i]*E1[0] + C[8*nq+4*nq+i]*E1[1])*H03 ;
  }

  return ;
}

static inline void increment_buf_cp_complex_bw(WBFMM_REAL *E0, WBFMM_REAL *E1,
					       WBFMM_REAL *C,
					       WBFMM_REAL H03p, WBFMM_REAL H03m,
					       gint nq,
					       WBFMM_REAL *buf)

{
  WBFMM_REAL sH03, dH03 ;
  gint i ;
  
  sH03 = H03p + H03m ; dH03 = H03p - H03m ;
  for ( i = 0 ; i < nq ; i ++ ) {
    buf[0*2*nq+2*i+0] += (C[0*nq+i]*E0[0] - C[8*nq+0*nq+i]*E0[1])*sH03 ;
    buf[0*2*nq+2*i+1] += (C[0*nq+i]*E0[1] + C[8*nq+0*nq+i]*E0[0])*dH03 ;

    buf[7*2*nq+2*i+0] += (C[7*nq+i]*E0[0] - C[8*nq+7*nq+i]*E0[1])*sH03 ;
    buf[7*2*nq+2*i+1] += (C[7*nq+i]*E0[1] + C[8*nq+7*nq+i]*E0[0])*dH03 ;

    buf[1*2*nq+2*i+0] += (C[1*nq+i]*E1[0] - C[8*nq+1*nq+i]*E1[1])*sH03 ;
    buf[1*2*nq+2*i+1] += (C[1*nq+i]*E1[1] + C[8*nq+1*nq+i]*E1[0])*dH03 ;

    buf[6*2*nq+2*i+0] += (C[6*nq+i]*E1[0] - C[8*nq+6*nq+i]*E1[1])*sH03 ;
    buf[6*2*nq+2*i+1] += (C[6*nq+i]*E1[1] + C[8*nq+6*nq+i]*E1[0])*dH03 ;

    buf[2*2*nq+2*i+0] += (C[2*nq+i]*E0[0] + C[8*nq+2*nq+i]*E0[1])*sH03 ;
    buf[2*2*nq+2*i+1] += (C[2*nq+i]*E0[1] - C[8*nq+2*nq+i]*E0[0])*dH03 ;

    buf[5*2*nq+2*i+0] += (C[5*nq+i]*E0[0] + C[8*nq+5*nq+i]*E0[1])*sH03 ;
    buf[5*2*nq+2*i+1] += (C[5*nq+i]*E0[1] - C[8*nq+5*nq+i]*E0[0])*dH03 ;

    buf[3*2*nq+2*i+0] += (C[3*nq+i]*E1[0] + C[8*nq+3*nq+i]*E1[1])*sH03 ;
    buf[3*2*nq+2*i+1] -= (C[3*nq+i]*E1[1] - C[8*nq+3*nq+i]*E1[0])*dH03 ;

    buf[4*2*nq+2*i+0] += (C[4*nq+i]*E1[0] + C[8*nq+4*nq+i]*E1[1])*sH03 ;
    buf[4*2*nq+2*i+1] -= (C[4*nq+i]*E1[1] - C[8*nq+4*nq+i]*E1[0])*dH03 ;
  }

  return ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_laplace_child_parent_shift_bw)(WBFMM_REAL *Cp,
							      gint Np,
							      WBFMM_REAL *Cc,
							      gint Nc,
							      gint nq,
							      WBFMM_REAL *H03, 
							      gint Lh,
							      WBFMM_REAL wb,
							      WBFMM_REAL *work)

/*
  child-parent shift for Laplace problem, based on Helmholtz version
  in shift.c, but allowing for vector inputs, and with indexing
  changed to Laplace convention

  03: "lower" boxes (think of Morton index)
  47: "upper" boxes

  assumes all coefficients are densely packed in groups of 8*nq
  coefficients of the same index, in order of child box Morton index
*/

{
  WBFMM_REAL *Cr, *cr, *ci, buf[128], H, H1, c, tn[64] = {0.0} ;
  WBFMM_REAL E0[2], E1[2], E0ph[2], E1ph[2], t ;
  gint nu, n, m, ic, ip, str, i, j, sgn, sgnn ;

  g_assert(nq <= 8) ;
  t = SQRT(3.0)*0.5*wb ;
  
  /*stride in number of elements per entry*/
  str = 8*nq ;

  /*used to store the rotated child coefficients*/
  Cr = &(work[str*(Np+1)*(Np+1)]) ;
  memset(Cr, 0, str*(Nc+1)*(Nc+1)*sizeof(WBFMM_REAL)) ;

  /*rotate child box coefficients using Cr as temporary storage*/
  for ( n = 0 ; n <= Nc ; n ++ ) {
    memset(buf, 0, 16*nq*sizeof(WBFMM_REAL)) ;
    nu = 0 ; ip = n*n ;
    m  = 0 ; ic = n*n ;
    
    H = H03[wbfmm_rotation_index_numn(nu,m,n)] ;
    for ( i = 0 ; i < 8*nq ; i ++ ) buf[i] = H*Cc[str*ic+i] ;

    for ( m = 1 ; m <= n ; m ++ ) {
      ic = wbfmm_index_laplace_nm(n,m) ;

      E0[0] = cos_n_PI_4(m)   ; E0[1] = sin_n_PI_4(m) ;
      E1[0] = cos_n_PI_4(3*m) ; E1[1] = sin_n_PI_4(3*m) ;

      increment_buf_cp_real_bw(E0, E1, &(Cc[str*ic]),
			       H03[wbfmm_rotation_index_numn( nu,m,n)],
			       H03[wbfmm_rotation_index_numn(-nu,m,n)],
			       nq, buf) ;
    }

    for ( i = 0 ; i < 8*nq ; i ++ ) Cr[str*ip+i] += buf[i] ;

    for ( nu = 1 ; nu <= n ; nu ++ ) {
      memset(buf, 0, 16*nq*sizeof(WBFMM_REAL)) ;

      E0ph[0] = cos_n_PI_2(nu) ; E0ph[1] = sin_n_PI_2(nu) ;
      
      ip = wbfmm_index_laplace_nm(n,nu) ;
      m = 0 ; ic = n*n ;
      H = H03[wbfmm_rotation_index_numn(nu,m,n)] ;
      for ( i = 0 ; i < 8*nq ; i ++ ) buf[2*i+0] = H*Cc[str*ic+i] ;
      
      for ( m = 1 ; m <= n ; m ++ ) {
	E0[0] = cos_n_PI_4(m)   ; E0[1] = sin_n_PI_4(m) ;
	E1[0] = cos_n_PI_4(3*m) ; E1[1] = sin_n_PI_4(3*m) ;
      
	ic = wbfmm_index_laplace_nm(n,m) ;
	increment_buf_cp_complex_bw(E0, E1, &(Cc[str*ic]),
				    H03[wbfmm_rotation_index_numn( nu,m,n)],
				    H03[wbfmm_rotation_index_numn(-nu,m,n)],
				    nq, buf) ;
      }

      cr = &(Cr[str*(ip+0)]) ; ci = &(Cr[str*(ip+1)]) ;
      for ( i = 0 ; i < nq ; i ++ ) {
	cr[0*nq+i] += buf[0*2*nq+2*i+0]*E0ph[0] - buf[0*2*nq+2*i+1]*E0ph[1] ;
	ci[0*nq+i] += buf[0*2*nq+2*i+1]*E0ph[0] + buf[0*2*nq+2*i+0]*E0ph[1] ;

	cr[7*nq+i] += buf[7*2*nq+2*i+0]*E0ph[0] - buf[7*2*nq+2*i+1]*E0ph[1] ;
	ci[7*nq+i] += buf[7*2*nq+2*i+1]*E0ph[0] + buf[7*2*nq+2*i+0]*E0ph[1] ;

	cr[1*nq+i] += buf[1*2*nq+2*i+0]*E0ph[0] + buf[1*2*nq+2*i+1]*E0ph[1] ;
	ci[1*nq+i] += buf[1*2*nq+2*i+1]*E0ph[0] - buf[1*2*nq+2*i+0]*E0ph[1] ;

	cr[6*nq+i] += buf[6*2*nq+2*i+0]*E0ph[0] + buf[6*2*nq+2*i+1]*E0ph[1] ;
	ci[6*nq+i] += buf[6*2*nq+2*i+1]*E0ph[0] - buf[6*2*nq+2*i+0]*E0ph[1] ;

	cr[2*nq+i] += buf[2*2*nq+2*i+0]*E0ph[0] + buf[2*2*nq+2*i+1]*E0ph[1] ;
	ci[2*nq+i] -= buf[2*2*nq+2*i+1]*E0ph[0] - buf[2*2*nq+2*i+0]*E0ph[1] ;

	cr[5*nq+i] += buf[5*2*nq+2*i+0]*E0ph[0] + buf[5*2*nq+2*i+1]*E0ph[1] ;
	ci[5*nq+i] -= buf[5*2*nq+2*i+1]*E0ph[0] - buf[5*2*nq+2*i+0]*E0ph[1] ;

	cr[3*nq+i] += buf[3*2*nq+2*i+0]*E0ph[0] - buf[3*2*nq+2*i+1]*E0ph[1] ;
	ci[3*nq+i] += buf[3*2*nq+2*i+1]*E0ph[0] + buf[3*2*nq+2*i+0]*E0ph[1] ;

	cr[4*nq+i] += buf[4*2*nq+2*i+0]*E0ph[0] - buf[4*2*nq+2*i+1]*E0ph[1] ;
	ci[4*nq+i] += buf[4*2*nq+2*i+1]*E0ph[0] + buf[4*2*nq+2*i+0]*E0ph[1] ;
      }
    }
  }

  /*Cr now contains rotated child coefficients, translate and store in
    work*/
  memset(work, 0, str*(Np+1)*(Np+1)*sizeof(WBFMM_REAL)) ;
  tn[0] = 1.0 ;
  m  = 0 ; sgnn = 1 ;
  for ( ip = n = 0 ; n <= Np ; (n ++), (ip = n*n) ) {
    sgn = sgnn ;
    for ( ic = nu = 0 ; nu <= MIN(n, Nc) ; (nu ++), (ic = nu*nu) ) {
      c = wbfmm_coaxial_translation_SS_cfft(n, nu, m)*tn[n-nu] ;
      g_assert(tn[n-nu] != 0.0) ;
      for ( j = 0 ; j < nq ; j ++ ) {
	for ( i = 0 ; i < 4 ; i ++ ) {
	  work[str*ip+0*nq+i*nq+j] +=
	    c*(Cr[str*ic+i*nq+j] + sgn*Cr[str*ic+(7-i)*nq+j]) ;
	}
      }
      sgn = -sgn ;
    }
    tn[n+1] = -tn[n]*t ;
    sgnn = -sgnn ;
  }

  for ( m = 1 ; m <= Np ; m ++ ) {
    sgnn = 1 ;
    for ( n = m ; n <= Np ; n ++ ) {
      ip = wbfmm_index_laplace_nm(n,m) ;
      sgn = sgnn ;
      for ( nu = m ; nu <= MIN(n,Nc) ; nu ++ ) {
  	ic = wbfmm_index_laplace_nm(nu,m) ;
  	c = wbfmm_coaxial_translation_SS_cfft(n, nu, m)*tn[n-nu] ;
	for ( j = 0 ; j < nq ; j ++ ) {
	  for ( i = 0 ; i < 4 ; i ++ ) {
	    work[str*(ip+0)+0*nq+i*nq+j] +=
	      c*(Cr[str*(ic+0)+i*nq+j] + sgn*Cr[str*(ic+0)+(7-i)*nq+j]) ;
	    work[str*(ip+1)+0*nq+i*nq+j] +=
	      c*(Cr[str*(ic+1)+i*nq+j] + sgn*Cr[str*(ic+1)+(7-i)*nq+j]) ;
	  }
	}
	sgn = -sgn ;
      }
      sgnn = -sgnn ;
    }
  }
  
  /*work now contains rotated and shifted coefficients, summed in
    pairs, perform reverse rotation into parent coefficient array*/
  for ( n = 0 ; n <= Np ; n ++ ) {
    nu = 0 ; ip = n*n ;
    m  = 0 ; ic = n*n ;

    H1 = H03[wbfmm_rotation_index_numn(nu,m,n)] ;
    for ( i = 0 ; i < nq ; i ++ ) {
      Cp[str*ip+i] +=
	H1*(work[str*ic+0*nq+i] + work[str*ic+1*nq+i] +
	    work[str*ic+2*nq+i] + work[str*ic+3*nq+i]) ;
    }

    E0ph[0] = cos_n_PI_4(nu) ; E0ph[1] = sin_n_PI_4(nu) ;
    for ( m = 1 ; m <= n ; m ++ ) {
      E0[0] = cos_n_PI_2(m) ; E0[1] = sin_n_PI_2(m) ;

      ic = wbfmm_index_laplace_nm(n,m) ;
      increment_cfft_cp_real_bw(E0, &(work[str*ic]),
				H03[wbfmm_rotation_index_numn( nu,m,n)],
				H03[wbfmm_rotation_index_numn(-nu,m,n)],
				nq, &(Cp[str*ip])) ;
    }

    for ( nu = 1 ; nu <= n ; nu ++ ) {
      E0ph[0] = cos_n_PI_4(  nu) ; E0ph[1] = sin_n_PI_4(  nu) ;
      E1ph[0] = cos_n_PI_4(3*nu) ; E1ph[1] = sin_n_PI_4(3*nu) ;

      ip = wbfmm_index_laplace_nm(n,nu) ;
      m = 0 ; ic = n*n ;

      H1 = H03[wbfmm_rotation_index_numn(nu,m,n)] ;
      for ( i = 0 ; i < nq ; i ++ ) {
      	Cp[str*(ip+0)+i] +=
	  ((work[str*ic+2*nq+i] + work[str*ic+0*nq+i])*E0ph[0] +
	   (work[str*ic+1*nq+i] + work[str*ic+3*nq+i])*E1ph[0])*H1 ;
      	Cp[str*(ip+1)+i] +=
	  ((work[str*ic+2*nq+i] - work[str*ic+0*nq+i])*E0ph[1] +
	   (work[str*ic+3*nq+i] - work[str*ic+1*nq+i])*E1ph[1])*H1 ;
      }

      for ( m = 1 ; m <= n ; m ++ ) {
	E0[0] = cos_n_PI_2(m) ; E0[1] = sin_n_PI_2(m) ;

	ic = wbfmm_index_laplace_nm(n,m) ;
	increment_cfft_cp_complex_bw(E0, E0ph, E1ph, &(work[str*ic]),
				     H03[wbfmm_rotation_index_numn( nu,m,n)],
				     H03[wbfmm_rotation_index_numn(-nu,m,n)],
				     nq, &(Cp[str*ip])) ;
      }
    }
  }
  
  return 0 ;
}

#if 0
static inline void increment_buf_pc_real(WBFMM_REAL *E, WBFMM_REAL *C,
					 WBFMM_REAL H03p, WBFMM_REAL H03m,
					 WBFMM_REAL H47p, WBFMM_REAL H47m,
					 gint nq,
					 WBFMM_REAL *buf)

{
  WBFMM_REAL H03, H47 ;
  gint i ;
  
  H03 = H03p + H03m ; H47 = H47p + H47m ;

  for ( i = 0 ; i < nq ; i ++ ) {
    buf[7*nq+i] += (C[0*nq+i]*E[2] + C[8*nq+i]*E[3])*H47 ;
    buf[6*nq+i] += (C[0*nq+i]*E[0] + C[8*nq+i]*E[1])*H47 ;
    buf[5*nq+i] += (C[0*nq+i]*E[2] - C[8*nq+i]*E[3])*H47 ;
    buf[4*nq+i] += (C[0*nq+i]*E[0] - C[8*nq+i]*E[1])*H47 ;

    buf[3*nq+i] += (C[0*nq+i]*E[2] + C[8*nq+i]*E[3])*H03 ;
    buf[2*nq+i] += (C[0*nq+i]*E[0] + C[8*nq+i]*E[1])*H03 ;
    buf[1*nq+i] += (C[0*nq+i]*E[2] - C[8*nq+i]*E[3])*H03 ;
    buf[0*nq+i] += (C[0*nq+i]*E[0] - C[8*nq+i]*E[1])*H03 ;
  }

  return ;
}

static inline void increment_buf_pc_complex(WBFMM_REAL *E,
					    WBFMM_REAL *C,
					    WBFMM_REAL H03p, WBFMM_REAL H03m,
					    WBFMM_REAL H47p, WBFMM_REAL H47m,
					    gint nq,
					    WBFMM_REAL *buf)

{
  WBFMM_REAL sH03, sH47, dH03, dH47, Ci, Cr ;
  WBFMM_REAL CrE[4], CiE[4] ;
  gint i ;
  
  sH03 = H03p + H03m ; sH47 = H47p + H47m ;
  dH03 = H03p - H03m ; dH47 = H47p - H47m ;
  for ( i = 0 ; i < nq ; i ++ ) {
    Cr = C[0*nq+i] ; Ci = C[8*nq+i] ;

    CrE[0] = Cr*E[0] ; CrE[1] = Cr*E[1] ;
    CrE[2] = Cr*E[2] ; CrE[3] = Cr*E[3] ;
    CiE[0] = Ci*E[0] ; CiE[1] = Ci*E[1] ;
    CiE[2] = Ci*E[2] ; CiE[3] = Ci*E[3] ;
    
    buf[7*2*nq+2*i+0] += (CrE[2] + CiE[3])*sH47 ;
    buf[7*2*nq+2*i+1] += (CiE[2] - CrE[3])*dH47 ;

    buf[6*2*nq+2*i+0] += (CrE[0] + CiE[1])*sH47 ;
    buf[6*2*nq+2*i+1] += (CiE[0] - CrE[1])*dH47 ;

    buf[5*2*nq+2*i+0] += (CrE[2] - CiE[3])*sH47 ;
    buf[5*2*nq+2*i+1] += (CiE[2] + CrE[3])*dH47 ;

    buf[4*2*nq+2*i+0] += (CrE[0] - CiE[1])*sH47 ;
    buf[4*2*nq+2*i+1] += (CiE[0] + CrE[1])*dH47 ;

    buf[3*2*nq+2*i+0] += (CrE[2] + CiE[3])*sH03 ;
    buf[3*2*nq+2*i+1] += (CiE[2] - CrE[3])*dH03 ;

    buf[2*2*nq+2*i+0] += (CrE[0] + CiE[1])*sH03 ;
    buf[2*2*nq+2*i+1] += (CiE[0] - CrE[1])*dH03 ;

    buf[1*2*nq+2*i+0] += (CrE[2] - CiE[3])*sH03 ;
    buf[1*2*nq+2*i+1] += (CiE[2] + CrE[3])*dH03 ;

    buf[0*2*nq+2*i+0] += (CrE[0] - CiE[1])*sH03 ;
    buf[0*2*nq+2*i+1] += (CiE[0] + CrE[1])*dH03 ;
  }

  return ;
}

static inline void increment_cfft_pc_real(WBFMM_REAL *E,  WBFMM_REAL *C,
					  WBFMM_REAL H03p, WBFMM_REAL H03m,
					  WBFMM_REAL H47p, WBFMM_REAL H47m,
					  gint nq, WBFMM_REAL *buf)

{
  WBFMM_REAL H03, H47 ;
  gint i ;
  
  H03 = H03p + H03m ; H47 = H47p + H47m ;

  for ( i = 0 ; i < 4*nq ; i ++ ) {
    buf[0*nq+i] += (C[0*nq+0*nq+i]*E[0] - C[8*nq+0*nq+i]*E[1])*H03 ; 
    buf[4*nq+i] += (C[0*nq+4*nq+i]*E[0] - C[8*nq+4*nq+i]*E[1])*H47 ; 
  }

  return ;
}

static inline void increment_cfft_pc_complex(WBFMM_REAL *E,   WBFMM_REAL *C,
					     WBFMM_REAL H03p, WBFMM_REAL H03m,
					     WBFMM_REAL H47p, WBFMM_REAL H47m,
					     gint nq, WBFMM_REAL *buf)

{
  WBFMM_REAL sH03, sH47, dH03, dH47 ;
  gint i ;
  
  sH03 = H03p + H03m ; sH47 = H47p + H47m ;
  dH03 = H03p - H03m ; dH47 = H47p - H47m ;

  for ( i = 0 ; i < 4*nq ; i ++ ) {
    buf[0*nq+2*i+0] += (C[0*nq+0*nq+i]*E[0] - C[8*nq+0*nq+i]*E[1])*sH03 ; 
    buf[0*nq+2*i+1] += (C[0*nq+0*nq+i]*E[1] + C[8*nq+0*nq+i]*E[0])*dH03 ; 
    buf[8*nq+2*i+0] += (C[0*nq+4*nq+i]*E[0] - C[8*nq+4*nq+i]*E[1])*sH47 ; 
    buf[8*nq+2*i+1] += (C[0*nq+4*nq+i]*E[1] + C[8*nq+4*nq+i]*E[0])*dH47 ; 
  }

  return ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_laplace_parent_child_shift)(WBFMM_REAL *Cc,
							   gint Nc,
							   WBFMM_REAL *Cp,
							   gint Np,
							   gint nq,
							   WBFMM_REAL *H03, 
							   WBFMM_REAL *H47,
							   gint Lh,
							   WBFMM_REAL wb,
							   WBFMM_REAL *work)

/*
  parent-child shift for Laplace problem, based on Helmholtz version
  in shift.c, but allowing for vector inputs, and with indexing
  changed to Laplace convention

  03: "lower" boxes (think of Morton index)
  47: "upper" boxes

  assumes all coefficients are densely packed in groups of 8*nq
  coefficients of the same index, in order of child box Morton index

  Cp points to the first entry of the parent box
*/

{
  WBFMM_REAL *Cr, *Ct, *cr, *ci, H, c, tn[64] = {1.0} ;
  WBFMM_REAL E0ph[2], E1ph[2], t, E[4] ;
  WBFMM_REAL buf[1024] ;
  gint nu, n, m, ic, ip, str, i ;
  
  g_assert(nq <= 64) ;
  
  /*translation length*/
  /* t = SQRT(3.0)*0.25*wb ; */
  t = SQRT(0.1875)*wb ;
  /*stride in number of elements per entry*/
  str = 8*nq ;

  /*used to store the rotated parent coefficients*/
  Cr = &(work[str*(Nc+1)*(Nc+1)]) ;
  memset(Cr, 0, str*(Np+1)*(Np+1)*sizeof(WBFMM_REAL)) ;

  /* buf = work ; */
  /*rotate parent box coefficients using Cr as temporary storage*/
  for ( n = 0 ; n <= Np ; n ++ ) {
    nu = 0 ; ic = n*n ;
    m  = 0 ; ip = n*n ;

    /*single set of parent coefficients must be rotated eight ways*/
    H = H03[wbfmm_rotation_index_numn(nu,m,n)] ;
    for ( i = 0 ; i < nq ; i ++ ) {
      buf[0*nq+i] = buf[1*nq+i] = buf[2*nq+i] = buf[3*nq+i] = H*Cp[str*ip+i] ;
    }
    
    H = H47[wbfmm_rotation_index_numn(nu,m,n)] ;
    for ( i = 0 ; i < nq ; i ++ ) {
      buf[4*nq+i] = buf[5*nq+i] = buf[6*nq+i] = buf[7*nq+i] = H*Cp[str*ip+i] ;
    }
    
    for ( m = 1 ; m <= n ; m ++ ) {
      ip = wbfmm_index_laplace_nm(n,m) ;

      E[0] = cos_n_PI_4(3*m) ; E[1] = -sin_n_PI_4(3*m) ;
      E[2] = cos_n_PI_4(m)   ; E[3] = -sin_n_PI_4(m) ;
      
      increment_buf_pc_real(E, &(Cp[str*ip]),
      			    H03[wbfmm_rotation_index_numn( nu,m,n)],
      			    H03[wbfmm_rotation_index_numn(-nu,m,n)],
      			    H47[wbfmm_rotation_index_numn( nu,m,n)],
      			    H47[wbfmm_rotation_index_numn(-nu,m,n)],
      			    nq, buf) ;
    }

    for ( i = 0 ; i < 8*nq ; i ++ ) Cr[str*ic+i] += buf[i] ;

    for ( nu = 1 ; nu <= n ; nu ++ ) {
      memset(buf, 0, 16*nq*sizeof(WBFMM_REAL)) ;

      E0ph[0] = cos_n_PI_2(nu) ; E0ph[1] = sin_n_PI_2(nu) ;
      
      ic = wbfmm_index_laplace_nm(n,nu) ;
      m = 0 ; ip = n*n ;
      H = H03[wbfmm_rotation_index_numn(nu,m,n)] ;
      for ( i = 0 ; i < nq ; i ++ ) {
	buf[0*2*nq+2*i+0] = buf[1*2*nq+2*i+0] =
	  buf[2*2*nq+2*i+0] = buf[3*2*nq+2*i+0] = H*Cp[str*ip+i] ;
      }
      H = H47[wbfmm_rotation_index_numn(nu,m,n)] ;
      for ( i = 0 ; i < nq ; i ++ ) {
	buf[4*2*nq+2*i+0] = buf[5*2*nq+2*i+0] =
	  buf[6*2*nq+2*i+0] = buf[7*2*nq+2*i+0] = H*Cp[str*ip+i] ;
      }
      
      for ( m = 1 ; m <= n ; m ++ ) {
	E[0] = cos_n_PI_4(3*m) ; E[1] = -sin_n_PI_4(3*m) ;
	E[2] = cos_n_PI_4(m)   ; E[3] = -sin_n_PI_4(m) ;
      
	ip = wbfmm_index_laplace_nm(n,m) ;
	increment_buf_pc_complex(E, &(Cp[str*ip]),
				 H03[wbfmm_rotation_index_numn( nu,m,n)],
				 H03[wbfmm_rotation_index_numn(-nu,m,n)],
				 H47[wbfmm_rotation_index_numn( nu,m,n)],
				 H47[wbfmm_rotation_index_numn(-nu,m,n)],
				 nq, buf) ;
      }

      cr = &(Cr[str*(ic+0)]) ; ci = &(Cr[str*(ic+1)]) ;
      for ( i = 0 ; i < 8*nq ; i ++ ) {
      	cr[i] += buf[2*i+0]*E0ph[0] - buf[2*i+1]*E0ph[1] ;
      	ci[i] += buf[2*i+1]*E0ph[0] + buf[2*i+0]*E0ph[1] ;
      }
    }
  }
  
  /*Cr now contains rotated parent coefficients, translate and store in
    Ct*/
  Ct = work ;
  memset(Ct, 0, str*(Nc+1)*(Nc+1)*sizeof(WBFMM_REAL)) ;
  m  = 0 ;
  for ( ip = nu = 0 ; nu <= Np ;
	(nu ++), (ip = nu*nu), (tn[nu] = tn[nu-1]*t) ) {
    for ( ic = n = 0 ; n <= MIN(nu,Nc) ; (n ++), (ic = n*n) ) {
      c = wbfmm_coaxial_translation_RR_cfft(n, nu, m)*tn[nu-n] ;
      /* g_assert(tn[nu-n] != 0.0) ; */
      for ( i = 0 ; i < 8*nq ; i ++ ) Ct[str*ic+i] += c*Cr[str*ip+i] ;
    }
  }
  
  for ( m = 1 ; m <= Nc ; m ++ ) {
    for ( nu = m ; nu <= Np ; nu ++ ) {
      ip = wbfmm_index_laplace_nm(nu,m) ;
      for ( n = m ; n <= MIN(nu,Nc) ; n ++ ) {
	ic = wbfmm_index_laplace_nm(n,m) ;
	c = wbfmm_coaxial_translation_RR_cfft(n, nu, m)*tn[nu-n] ;
	/* g_assert(tn[n-nu] != 0.0) ; */
	for ( i = 0 ; i < 16*nq ; i ++ ) {
	  Ct[str*(ic+0)+i] += c*Cr[str*(ip+0)+i] ;
	}
      }
    }
  }

  /*Ct now contains rotated and shifted coefficients, perform
    reverse rotation into child coefficient array*/

  for ( n = 0 ; n <= Nc ; n ++ ) {
    nu = 0 ; ic = n*n ;
    m  = 0 ; ip = n*n ;

    H = H03[wbfmm_rotation_index_numn(nu,m,n)] ;
    for ( i = 0 ; i < 4*nq ; i ++ ) buf[i] = H*Ct[str*ip+i] ;
    
    H = H47[wbfmm_rotation_index_numn(nu,m,n)] ;
    for ( i = 4*nq ; i < 8*nq ; i ++ ) buf[i] = H*Ct[str*ip+i] ;

    for ( m = 1 ; m <= n ; m ++ ) {
      ip = wbfmm_index_laplace_nm(n,m) ;

      E[0] = cos_n_PI_2(m) ; E[1] = -sin_n_PI_2(m) ;

      increment_cfft_pc_real(E, &(Ct[str*ip]),
			     H03[wbfmm_rotation_index_numn( nu,m,n)],
			     H03[wbfmm_rotation_index_numn(-nu,m,n)],
			     H47[wbfmm_rotation_index_numn( nu,m,n)],
			     H47[wbfmm_rotation_index_numn(-nu,m,n)],
			     nq, buf) ;
    }

    for ( i = 0 ; i < 8*nq ; i ++ ) Cc[str*ic+i] += buf[i] ;

    for ( nu = 1 ; nu <= n ; nu ++ ) {
      memset(buf, 0, 16*nq*sizeof(WBFMM_REAL)) ;

      E0ph[0] = cos_n_PI_4(3*nu) ; E0ph[1] = sin_n_PI_4(3*nu) ;
      E1ph[0] = cos_n_PI_4(nu)   ; E1ph[1] = sin_n_PI_4(nu) ;
      
      ic = wbfmm_index_laplace_nm(n,nu) ;
      m = 0 ; ip = n*n ;
      H = H03[wbfmm_rotation_index_numn(nu,m,n)] ;
      for ( i = 0 ; i < 4*nq ; i ++ ) buf[2*i+0] = H*Ct[str*ip+i] ;
      H = H47[wbfmm_rotation_index_numn(nu,m,n)] ;
      for ( i = 4*nq ; i < 8*nq ; i ++ ) buf[2*i+0] = H*Ct[str*ip+i] ;
      
      for ( m = 1 ; m <= n ; m ++ ) {
	E[0] = cos_n_PI_2(m)   ; E[1] = -sin_n_PI_2(m) ;
      
	ip = wbfmm_index_laplace_nm(n,m) ;
	increment_cfft_pc_complex(E, &(Ct[str*ip]),
				  H03[wbfmm_rotation_index_numn( nu,m,n)],
				  H03[wbfmm_rotation_index_numn(-nu,m,n)],
				  H47[wbfmm_rotation_index_numn( nu,m,n)],
				  H47[wbfmm_rotation_index_numn(-nu,m,n)],
				  nq, buf) ;
      }

      cr = &(Cc[str*(ic+0)]) ; ci = &(Cc[str*(ic+1)]) ;
      for ( i = 0 ; i < nq ; i ++ ) {
	cr[7*nq+i] += buf[7*2*nq+2*i+0]*E1ph[0] + buf[7*2*nq+2*i+1]*E1ph[1] ;
	ci[7*nq+i] += buf[7*2*nq+2*i+1]*E1ph[0] - buf[7*2*nq+2*i+0]*E1ph[1] ;

	cr[6*nq+i] += buf[6*2*nq+2*i+0]*E0ph[0] + buf[6*2*nq+2*i+1]*E0ph[1] ;
	ci[6*nq+i] += buf[6*2*nq+2*i+1]*E0ph[0] - buf[6*2*nq+2*i+0]*E0ph[1] ;

	cr[5*nq+i] += buf[5*2*nq+2*i+0]*E1ph[0] - buf[5*2*nq+2*i+1]*E1ph[1] ;
	ci[5*nq+i] += buf[5*2*nq+2*i+1]*E1ph[0] + buf[5*2*nq+2*i+0]*E1ph[1] ;

	cr[4*nq+i] += buf[4*2*nq+2*i+0]*E0ph[0] - buf[4*2*nq+2*i+1]*E0ph[1] ;
	ci[4*nq+i] += buf[4*2*nq+2*i+1]*E0ph[0] + buf[4*2*nq+2*i+0]*E0ph[1] ;

	cr[3*nq+i] += buf[3*2*nq+2*i+0]*E1ph[0] + buf[3*2*nq+2*i+1]*E1ph[1] ;
	ci[3*nq+i] += buf[3*2*nq+2*i+1]*E1ph[0] - buf[3*2*nq+2*i+0]*E1ph[1] ;

	cr[2*nq+i] += buf[2*2*nq+2*i+0]*E0ph[0] + buf[2*2*nq+2*i+1]*E0ph[1] ;
	ci[2*nq+i] += buf[2*2*nq+2*i+1]*E0ph[0] - buf[2*2*nq+2*i+0]*E0ph[1] ;

	cr[1*nq+i] += buf[1*2*nq+2*i+0]*E1ph[0] - buf[1*2*nq+2*i+1]*E1ph[1] ;
	ci[1*nq+i] += buf[1*2*nq+2*i+1]*E1ph[0] + buf[1*2*nq+2*i+0]*E1ph[1] ;

	cr[0*nq+i] += buf[0*2*nq+2*i+0]*E0ph[0] - buf[0*2*nq+2*i+1]*E0ph[1] ;
	ci[0*nq+i] += buf[0*2*nq+2*i+1]*E0ph[0] + buf[0*2*nq+2*i+0]*E0ph[1] ;
      }
    }
  }

  return 0 ;
}
#endif
