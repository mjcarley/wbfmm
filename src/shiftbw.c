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

static inline void sincos_recursion(WBFMM_REAL Epr[], WBFMM_REAL Epi[],
				    WBFMM_REAL Enr[], WBFMM_REAL Eni[],
				    WBFMM_REAL Cch[], WBFMM_REAL Sch[])

{
  WBFMM_REAL tmp ;

  /*Ep increment argument by \chi, En decrement by \chi*/
  tmp = Epr[0] ; 
  Epr[0] = Epr[0]*Cch[0] - Epi[0]*Sch[0] ;
  Epi[0] = tmp*Sch[0] + Epi[0]*Cch[0] ;
  tmp = Epr[1] ; 
  Epr[1] = Epr[1]*Cch[1] - Epi[1]*Sch[1] ;
  Epi[1] = tmp*Sch[1] + Epi[1]*Cch[1] ;

  tmp = Enr[0] ; 
  Enr[0] = Enr[0]*Cch[0] + Eni[0]*Sch[0] ;
  Eni[0] = -tmp*Sch[0] + Eni[0]*Cch[0] ;
  tmp = Enr[1] ; 
  Enr[1] = Enr[1]*Cch[1] + Eni[1]*Sch[1] ;
  Eni[1] = -tmp*Sch[1] + Eni[1]*Cch[1] ;

  return ;
}

static inline void increment_buf_cp_bw(WBFMM_REAL Er[], WBFMM_REAL Ei[],
				       WBFMM_REAL *Cc, WBFMM_REAL buf[])

/*
  increment buffer entries for child-parent shift operations

  Er: \cos(m\chi_{1}-\nu\phi_{1}) \cos(m\chi_{2}-\nu\phi_{2})
  Ei: \sin(m\chi_{1}-\nu\phi_{1}) \sin(m\chi_{2}-\nu\phi_{2})
*/

{
  buf[2*0+0] += (Er[0]*Cc[2*0+0] - Ei[0]*Cc[2*0+1]) ;
  buf[2*0+1] += (Er[0]*Cc[2*0+1] + Ei[0]*Cc[2*0+0]) ;

  buf[2*7+0] += (Er[0]*Cc[2*7+0] - Ei[0]*Cc[2*7+1]) ;
  buf[2*7+1] += (Er[0]*Cc[2*7+1] + Ei[0]*Cc[2*7+0]) ;

  buf[2*1+0] += (Er[1]*Cc[2*1+0] - Ei[1]*Cc[2*1+1]) ;
  buf[2*1+1] += (Er[1]*Cc[2*1+1] + Ei[1]*Cc[2*1+0]) ;

  buf[2*6+0] += (Er[1]*Cc[2*6+0] - Ei[1]*Cc[2*6+1]) ;
  buf[2*6+1] += (Er[1]*Cc[2*6+1] + Ei[1]*Cc[2*6+0]) ;
    
  buf[2*2+0] += (Er[0]*Cc[2*2+0] + Ei[0]*Cc[2*2+1]) ;
  buf[2*2+1] += (Er[0]*Cc[2*2+1] - Ei[0]*Cc[2*2+0]) ;
    
  buf[2*5+0] += (Er[0]*Cc[2*5+0] + Ei[0]*Cc[2*5+1]) ;
  buf[2*5+1] += (Er[0]*Cc[2*5+1] - Ei[0]*Cc[2*5+0]) ;

  buf[2*3+0] += (Er[1]*Cc[2*3+0] + Ei[1]*Cc[2*3+1]) ;
  buf[2*3+1] += (Er[1]*Cc[2*3+1] - Ei[1]*Cc[2*3+0]) ;
    
  buf[2*4+0] += (Er[1]*Cc[2*4+0] + Ei[1]*Cc[2*4+1]) ;
  buf[2*4+1] += (Er[1]*Cc[2*4+1] - Ei[1]*Cc[2*4+0]) ;

  return ;
}

static inline void increment_buf_pc(WBFMM_REAL Er[], WBFMM_REAL Ei[],
				    WBFMM_REAL *Cp, 
				    WBFMM_REAL H03, WBFMM_REAL H47,
				    WBFMM_REAL buf[])

/*
  increment buffer entries for parent-child shift operations

  this is the same as increment_buf_cp except that rotation angles are
  taken in the opposite order and Cp contains only one coefficient
  (the parent box coefficient which is rotated eight ways into the
  child entries)
*/

{
  buf[2*7+0] += H47*(Er[0]*Cp[0] - Ei[0]*Cp[1]) ;
  buf[2*7+1] += H47*(Er[0]*Cp[1] + Ei[0]*Cp[0]) ;

  buf[2*6+0] += H47*(Er[1]*Cp[0] - Ei[1]*Cp[1]) ;
  buf[2*6+1] += H47*(Er[1]*Cp[1] + Ei[1]*Cp[0]) ;

  buf[2*5+0] += H47*(Er[0]*Cp[0] + Ei[0]*Cp[1]) ;
  buf[2*5+1] += H47*(Er[0]*Cp[1] - Ei[0]*Cp[0]) ;
    
  buf[2*4+0] += H47*(Er[1]*Cp[0] + Ei[1]*Cp[1]) ;
  buf[2*4+1] += H47*(Er[1]*Cp[1] - Ei[1]*Cp[0]) ;
    
  buf[2*3+0] += H03*(Er[0]*Cp[0] - Ei[0]*Cp[1]) ;
  buf[2*3+1] += H03*(Er[0]*Cp[1] + Ei[0]*Cp[0]) ;

  buf[2*2+0] += H03*(Er[1]*Cp[0] - Ei[1]*Cp[1]) ;
  buf[2*2+1] += H03*(Er[1]*Cp[1] + Ei[1]*Cp[0]) ;

  buf[2*1+0] += H03*(Er[0]*Cp[0] + Ei[0]*Cp[1]) ;
  buf[2*1+1] += H03*(Er[0]*Cp[1] - Ei[0]*Cp[0]) ;
    
  buf[2*0+0] += H03*(Er[1]*Cp[0] + Ei[1]*Cp[1]) ;
  buf[2*0+1] += H03*(Er[1]*Cp[1] - Ei[1]*Cp[0]) ;

  return ;
}

static inline void increment_cfft_cp_bw(WBFMM_REAL Er[], WBFMM_REAL Ei[],
					WBFMM_REAL *Cr, WBFMM_REAL H03,
					WBFMM_REAL *Cp)

{
  Cp[0] += 
    H03*(Er[0]*(Cr[0] + Cr[4]) -  Ei[0]*(Cr[1]-Cr[5]) + 
	 Er[1]*(Cr[2] + Cr[6]) -  Ei[1]*(Cr[3]-Cr[7])) ;

  Cp[1] += 
    H03*(Er[0]*(Cr[1]+Cr[5]) + Ei[0]*(Cr[0]-Cr[4]) +
	 Er[1]*(Cr[3]+Cr[7]) + Ei[1]*(Cr[2]-Cr[6])) ;

  return ;
}

static inline void increment_cfft_pc(WBFMM_REAL Er[], WBFMM_REAL Ei[],
				     WBFMM_REAL *Cr, 
				     WBFMM_REAL H03, WBFMM_REAL H47,
				     WBFMM_REAL *buf)

{
  buf[2*7+0] += H47*(Er[0]*Cr[2*7+0] - Ei[0]*Cr[2*7+1]) ;
  buf[2*7+1] += H47*(Er[0]*Cr[2*7+1] + Ei[0]*Cr[2*7+0]) ;

  buf[2*6+0] += H47*(Er[1]*Cr[2*6+0] - Ei[1]*Cr[2*6+1]) ;
  buf[2*6+1] += H47*(Er[1]*Cr[2*6+1] + Ei[1]*Cr[2*6+0]) ;

  buf[2*5+0] += H47*(Er[0]*Cr[2*5+0] + Ei[0]*Cr[2*5+1]) ;
  buf[2*5+1] += H47*(Er[0]*Cr[2*5+1] - Ei[0]*Cr[2*5+0]) ;
    
  buf[2*4+0] += H47*(Er[1]*Cr[2*4+0] + Ei[1]*Cr[2*4+1]) ;
  buf[2*4+1] += H47*(Er[1]*Cr[2*4+1] - Ei[1]*Cr[2*4+0]) ;
    
  buf[2*3+0] += H03*(Er[0]*Cr[2*3+0] - Ei[0]*Cr[2*3+1]) ;
  buf[2*3+1] += H03*(Er[0]*Cr[2*3+1] + Ei[0]*Cr[2*3+0]) ;

  buf[2*2+0] += H03*(Er[1]*Cr[2*2+0] - Ei[1]*Cr[2*2+1]) ;
  buf[2*2+1] += H03*(Er[1]*Cr[2*2+1] + Ei[1]*Cr[2*2+0]) ;

  buf[2*1+0] += H03*(Er[0]*Cr[2*1+0] + Ei[0]*Cr[2*1+1]) ;
  buf[2*1+1] += H03*(Er[0]*Cr[2*1+1] - Ei[0]*Cr[2*1+0]) ;
    
  buf[2*0+0] += H03*(Er[1]*Cr[2*0+0] + Ei[1]*Cr[2*0+1]) ;
  buf[2*0+1] += H03*(Er[1]*Cr[2*0+1] - Ei[1]*Cr[2*0+0]) ;

  return ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_child_parent_shift_bw)(WBFMM_REAL *Cp, gint Np,
						      WBFMM_REAL *Cc, gint Nc,
						      WBFMM_REAL *H03, gint Lh,
						      WBFMM_REAL *transf,
						      WBFMM_REAL *transb,
						      gint Ls,
						      WBFMM_REAL *work)

/*
  03: "lower" boxes (think of Morton index)

  assumes all coefficients are densely packed in groups of eight
  coefficients of the same index, in order of child box Morton index
*/

{
  WBFMM_REAL *Cr, tmp ;
  WBFMM_REAL Epr[2], Epi[2], Enr[2], Eni[2], Cch[2], Sch[2] ;
  WBFMM_REAL Cnph, Snph, S0, Cnph0, Snph0 ;
  WBFMM_REAL Cnch, Snch, Cn3ch, Sn3ch ;
  WBFMM_REAL Cnch0, Snch0, Cn3ch0, Sn3ch0 ;
  gint nu, l, n, m, ic, ip, ih, str, i, sgn, off, offc ;

  /*stride in number of complex elements per entry*/
  str = 8 ;
  /*used for trigonometric recursions*/
  Cch[0] = Sch[0] = Sch[1] = M_SQRT1_2 ; Cch[1] = -M_SQRT1_2 ;

  Cnph0 = 1.0 ; Snph0 = 0.0 ; S0 = 1.0 ;
  /*used to store the rotated child coefficients*/
  /* Cr = &(work[8*2*wbfmm_coefficient_index_nm(Np+1,0)]) ; */
  Cr = &(work[8*2*wbfmm_coefficient_number(Np)]) ;

  /*rotate child box coefficients using Cr as temporary storage*/
  for ( n = 0 ; n <= Nc ; n ++ ) {
    Cnph = Cnph0 ; Snph = Snph0 ;
    for ( nu = -n ; nu <= n ; nu ++ ) {
      WBFMM_REAL buf[16] = {0.0} ;
      WBFMM_REAL HEr[2], HEi[2] ;
      m = 0 ; 
      Epr[0] = Enr[0] = Epr[1] = Enr[1] = Cnph ;
      Epi[0] = Eni[0] = Snph ;
      Epi[1] = Eni[1] = -Epi[0] ;

      ic = wbfmm_coefficient_index_nm(n,m) ;
      ih = wbfmm_rotation_index_numn(nu,m,n) ;

      offc = 2*str*ic ;
      HEr[0] = Epr[0]*H03[ih] ; HEr[1] = Epr[1]*H03[ih] ;
      HEi[0] = Epi[0]*H03[ih] ; HEi[1] = Epi[1]*H03[ih] ;
      increment_buf_cp_bw(HEr, HEi, &(Cc[offc]), buf) ;

      for ( m = 1 ; m <= n ; m ++ ) {
	sincos_recursion(Epr, Epi, Enr, Eni, Cch, Sch) ;

	ic = wbfmm_coefficient_index_nm(n,m) ;
	ih = wbfmm_rotation_index_numn(nu,m,n) ;
	offc = 2*str*ic ;

	HEr[0] = Epr[0]*H03[ih] ; HEr[1] = Epr[1]*H03[ih] ;
	HEi[0] = Epi[0]*H03[ih] ; HEi[1] = Epi[1]*H03[ih] ;
	increment_buf_cp_bw(HEr, HEi, &(Cc[offc]), buf) ;

	ic = wbfmm_coefficient_index_nm(n,-m) ;
	ih = wbfmm_rotation_index_numn(-nu,m,n) ;
	offc = 2*str*ic ;

	HEr[0] = Enr[0]*H03[ih] ; HEr[1] = Enr[1]*H03[ih] ;
	HEi[0] = Eni[0]*H03[ih] ; HEi[1] = Eni[1]*H03[ih] ;
	increment_buf_cp_bw(HEr, HEi, &(Cc[offc]), buf) ;
      }

      ip = wbfmm_coefficient_index_nm(n,nu) ; off = 2*str*ip ;      
      for ( i = 0 ; i < 16 ; i ++ ) Cr[off+i] = buf[i] ;

      /*trigonometric recursions taking note that C0 = 0*/
      tmp = Cnph ; Cnph = -Snph ; Snph = tmp ;
    }
    tmp = Cnph0 ; Cnph0 = Snph0 ; Snph0 = -tmp ;
  }

  /*Cr now contains the rotated Cc coefficients, apply coaxial
    translation to all coefficients to shift to centre of parent box,
    using BW translations on "upper" child boxes (4-7) to combine them
    with their opposite numbers in the lower boxes (0-3)*/
  for ( n = 0 ; n <= Np ; n ++ ) {
    for ( m = -n ; m <= n ; m ++ ) {
      WBFMM_REAL buf[8] = {0.0} ;
      /*loop on input and coefficients*/
      for ( l = ABS(m) ; l <= Nc ; l ++ ) {
	ic = wbfmm_coefficient_index_nm(l, m) ;
	offc = 2*str*ic ;
	ih = wbfmm_coaxial_index(l, m, n) ;
	sgn = wbfmm_coaxial_index_sgn(l, m, n) ;
	for ( i = 0 ; i < 4 ; i ++ ) {
	  /*forward translation on lower child boxes*/
	  buf[2*i+0] += transf[ih]*sgn*Cr[offc+2*i+0] ;
	  buf[2*i+1] += transf[ih]*sgn*Cr[offc+2*i+1] ;
	  /*bumwise translation on upper child boxes*/
	  buf[2*i+0] += transb[ih]*sgn*Cr[offc+2*(7-i)+0] ;
	  buf[2*i+1] += transb[ih]*sgn*Cr[offc+2*(7-i)+1] ;
	}
      }
      ip = wbfmm_coefficient_index_nm(n, m) ;
      /*offset into work is halved because of summation of upper and lower
       child contributions*/
      off = 2*str*ip/2 ;
      for ( i = 0 ; i < 8 ; i ++ ) work[off+i] = buf[i] ;
    }
  }

  /*work now contains the rotated coefficients shifted to the centre
    of the parent box, which must be reverse rotated and summed into Cp*/
  Cch[0] = Cch[1] = 0.0 ; Sch[0] = -1 ; Sch[1] = 1.0 ;
  S0 = M_SQRT1_2 ;
  Cnch0 = Cn3ch0 = 1.0 ; Snch0 = Sn3ch0 = 0.0 ;
  for ( n = 0 ; n <= Np ; n ++ ) {
    Cnch  = Cnch0  ; Snch  = Snch0 ;
    Cn3ch = Cn3ch0 ; Sn3ch = Sn3ch0 ;
    for ( nu = -n ; nu <= n ; nu ++ ) {
      ip = wbfmm_coefficient_index_nm(n,nu) ; off = 2*str*ip ;
      Cp[off+0] = Cp[off+1] = 0.0 ;

      m = 0 ; 
      Epr[0] = Enr[0] =  Cnch  ;
      Epi[0] = Eni[0] = -Snch  ;
      Epr[1] = Enr[1] =  Cn3ch ;
      Epi[1] = Eni[1] = -Sn3ch ;

      ic = wbfmm_coefficient_index_nm(n,m) ;
      ih = wbfmm_rotation_index_numn(nu,m,n) ;
      /* offc = 2*str*ic ; */
      offc = 2*str*ic/2 ;

      increment_cfft_cp_bw(Epr, Epi, &(work[offc]), H03[ih], &(Cp[off])) ;

      for ( m = 1 ; m <= n ; m ++ ) {
	sincos_recursion(Epr, Epi, Enr, Eni, Cch, Sch) ;

	ic = wbfmm_coefficient_index_nm(n,m) ;
	ih = wbfmm_rotation_index_numn(nu,m,n) ;
        /* offc = 2*str*ic ; */
        offc = 2*str*ic/2 ;

	increment_cfft_cp_bw(Epr, Epi, &(work[offc]), H03[ih], &(Cp[off])) ;

	ic = wbfmm_coefficient_index_nm(n,-m) ;
	ih = wbfmm_rotation_index_numn(-nu,m,n) ;
        /* offc = 2*str*ic ; */
        offc = 2*str*ic/2 ;

	increment_cfft_cp_bw(Enr, Eni, &(work[offc]), H03[ih], &(Cp[off])) ;
      }
      tmp = Cnch ; Cnch = S0*(Cnch - Snch) ; Snch = S0*(tmp + Snch) ;
      tmp = Cn3ch ; Cn3ch = -S0*(Cn3ch + Sn3ch) ; Sn3ch = S0*(tmp - Sn3ch) ;
    }
    tmp = Cnch0 ; 
    Cnch0 = S0*(Cnch0 + Snch0) ; Snch0 = S0*(-tmp + Snch0) ;
    tmp = Cn3ch0 ;
    Cn3ch0 = S0*(-Cn3ch0 + Sn3ch0) ; Sn3ch0 = -S0*(tmp + Sn3ch0) ;
  }
  
  return 0 ;
}
