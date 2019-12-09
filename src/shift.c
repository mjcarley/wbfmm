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

WBFMM_REAL _wbfmm_shifts_th[49] = {0.0} ;
WBFMM_REAL _wbfmm_shifts_ph[17] = {0.0} ;
/* WBFMM_REAL _wbfmm_shifts_r[19] = {0.0} ; */
WBFMM_REAL _wbfmm_shifts_r[WBFMM_SHIFTS_R_NUMBER] = {0.0} ;

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

static inline void increment_buf_cp(WBFMM_REAL Er[], WBFMM_REAL Ei[],
				    WBFMM_REAL *Cc, 
				    WBFMM_REAL H03, WBFMM_REAL H47,
				    WBFMM_REAL buf[])

/*
  increment buffer entries for child-parent shift operations

  Er: \cos(m\chi_{1}-\nu\phi_{1}) \cos(m\chi_{2}-\nu\phi_{2})
  Ei: \sin(m\chi_{1}-\nu\phi_{1}) \sin(m\chi_{2}-\nu\phi_{2})
*/

{
  buf[2*0+0] += H03*(Er[0]*Cc[2*0+0] - Ei[0]*Cc[2*0+1]) ;
  buf[2*0+1] += H03*(Er[0]*Cc[2*0+1] + Ei[0]*Cc[2*0+0]) ;

  buf[2*1+0] += H03*(Er[1]*Cc[2*1+0] - Ei[1]*Cc[2*1+1]) ;
  buf[2*1+1] += H03*(Er[1]*Cc[2*1+1] + Ei[1]*Cc[2*1+0]) ;

  buf[2*2+0] += H03*(Er[0]*Cc[2*2+0] + Ei[0]*Cc[2*2+1]) ;
  buf[2*2+1] += H03*(Er[0]*Cc[2*2+1] - Ei[0]*Cc[2*2+0]) ;
    
  buf[2*3+0] += H03*(Er[1]*Cc[2*3+0] + Ei[1]*Cc[2*3+1]) ;
  buf[2*3+1] += H03*(Er[1]*Cc[2*3+1] - Ei[1]*Cc[2*3+0]) ;
    
  buf[2*4+0] += H47*(Er[0]*Cc[2*4+0] - Ei[0]*Cc[2*4+1]) ;
  buf[2*4+1] += H47*(Er[0]*Cc[2*4+1] + Ei[0]*Cc[2*4+0]) ;

  buf[2*5+0] += H47*(Er[1]*Cc[2*5+0] - Ei[1]*Cc[2*5+1]) ;
  buf[2*5+1] += H47*(Er[1]*Cc[2*5+1] + Ei[1]*Cc[2*5+0]) ;

  buf[2*6+0] += H47*(Er[0]*Cc[2*6+0] + Ei[0]*Cc[2*6+1]) ;
  buf[2*6+1] += H47*(Er[0]*Cc[2*6+1] - Ei[0]*Cc[2*6+0]) ;
    
  buf[2*7+0] += H47*(Er[1]*Cc[2*7+0] + Ei[1]*Cc[2*7+1]) ;
  buf[2*7+1] += H47*(Er[1]*Cc[2*7+1] - Ei[1]*Cc[2*7+0]) ;

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

static inline void increment_cfft_cp(WBFMM_REAL Er[], WBFMM_REAL Ei[],
				     WBFMM_REAL *Cr, 
				     WBFMM_REAL H03, WBFMM_REAL H47,
				     WBFMM_REAL *Cp)

{
  Cp[0] += 
    H03*(Er[0]*(Cr[0] + Cr[4]) -  Ei[0]*(Cr[1]-Cr[5]) + 
	 Er[1]*(Cr[2] + Cr[6]) -  Ei[1]*(Cr[3]-Cr[7])) +
    H47*(Er[0]*(Cr[8] + Cr[12]) - Ei[0]*(Cr[9]-Cr[13]) +
	 Er[1]*(Cr[10]+ Cr[14]) - Ei[1]*(Cr[11]-Cr[15])) ;

  Cp[1] += 
    H03*(Er[0]*(Cr[1]+Cr[5]) + Ei[0]*(Cr[0]-Cr[4]) +
	 Er[1]*(Cr[3]+Cr[7]) + Ei[1]*(Cr[2]-Cr[6])) +
    H47*(Er[0]*(Cr[9]+Cr[13]) + Ei[0]*(Cr[8]-Cr[12]) +
	 Er[1]*(Cr[11]+Cr[15]) + Ei[1]*(Cr[10]-Cr[14])) ;

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

gint WBFMM_FUNCTION_NAME(wbfmm_child_parent_shift)(WBFMM_REAL *Cp, gint Np,
						   WBFMM_REAL *Cc, gint Nc,
						   WBFMM_REAL *H03, 
						   WBFMM_REAL *H47, gint Lh,
						   WBFMM_REAL *trans, gint Ls,
						   WBFMM_REAL *work)

/*
  03: "lower" boxes (think of Morton index)
  47: "upper" boxes

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
  Cr = &(work[8*2*wbfmm_coefficient_index_nm(Np+1,0)]) ;

  /*rotate child box coefficients using Cr as temporary storage*/
  for ( n = 0 ; n <= Nc ; n ++ ) {
    Cnph = Cnph0 ; Snph = Snph0 ;
    for ( nu = -n ; nu <= n ; nu ++ ) {
      WBFMM_REAL buf[16] = {0.0} ;
      m = 0 ; 
      Epr[0] = Enr[0] = Epr[1] = Enr[1] = Cnph ;
      Epi[0] = Eni[0] = Snph ;
      Epi[1] = Eni[1] = -Epi[0] ;

      ic = wbfmm_coefficient_index_nm(n,m) ;
      ih = wbfmm_rotation_index_numn(nu,m,n) ;

      offc = 2*str*ic ;
      increment_buf_cp(Epr, Epi, &(Cc[offc]), H03[ih], H47[ih], buf) ;

      for ( m = 1 ; m <= n ; m ++ ) {
	sincos_recursion(Epr, Epi, Enr, Eni, Cch, Sch) ;

	ic = wbfmm_coefficient_index_nm(n,m) ;
	ih = wbfmm_rotation_index_numn(nu,m,n) ;
	offc = 2*str*ic ;

	increment_buf_cp(Epr, Epi, &(Cc[offc]), H03[ih], H47[ih], buf) ;

	ic = wbfmm_coefficient_index_nm(n,-m) ;
	ih = wbfmm_rotation_index_numn(-nu,m,n) ;
	offc = 2*str*ic ;

	increment_buf_cp(Enr, Eni, &(Cc[offc]), H03[ih], H47[ih], buf) ;
      }

      ip = wbfmm_coefficient_index_nm(n,nu) ; off = 2*str*ip ;      
      for ( i = 0 ; i < 16 ; i ++ ) Cr[off+i] = buf[i] ;

      /*trigonometric recursions taking note that C0 = 0*/
      tmp = Cnph ; Cnph = -Snph ; Snph = tmp ;
    }
    tmp = Cnph0 ; Cnph0 = Snph0 ; Snph0 = -tmp ;
  }

  /*Cr now contains the rotated Cc coefficients, apply coaxial
    translation to all coefficients to shift to centre of parent box*/
  for ( n = 0 ; n <= Np ; n ++ ) {
    for ( m = -n ; m <= n ; m ++ ) {
      WBFMM_REAL buf[16] = {0.0} ;
      /*loop on input and coefficients*/
      for ( l = ABS(m) ; l <= Nc ; l ++ ) {
	ic = wbfmm_coefficient_index_nm(l, m) ;
	offc = 2*str*ic ;
	ih = wbfmm_coaxial_index(l, m, n) ;
	sgn = wbfmm_coaxial_index_sgn(l, m, n) ;
	for ( i = 0 ; i < 16 ; i ++ ) buf[i] += trans[ih]*sgn*Cr[offc+i] ;
      }
      ip = wbfmm_coefficient_index_nm(n, m) ;
      off = 2*str*ip ;
      for ( i = 0 ; i < 16 ; i ++ ) work[off+i] = buf[i] ;
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
      offc = 2*str*ic ;

      increment_cfft_cp(Epr, Epi, &(work[offc]), H03[ih], H47[ih], &(Cp[off])) ;

      for ( m = 1 ; m <= n ; m ++ ) {
	sincos_recursion(Epr, Epi, Enr, Eni, Cch, Sch) ;

	ic = wbfmm_coefficient_index_nm(n,m) ;
	ih = wbfmm_rotation_index_numn(nu,m,n) ;
        offc = 2*str*ic ;

	increment_cfft_cp(Epr, Epi, &(work[offc]), H03[ih], H47[ih], 
			  &(Cp[off])) ;

	ic = wbfmm_coefficient_index_nm(n,-m) ;
	ih = wbfmm_rotation_index_numn(-nu,m,n) ;
        offc = 2*str*ic ;

	increment_cfft_cp(Enr, Eni, &(work[offc]), H03[ih], H47[ih], 
			  &(Cp[off])) ;
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

gint WBFMM_FUNCTION_NAME(wbfmm_parent_child_shift)(WBFMM_REAL *Cc, gint Nc,
						   WBFMM_REAL *Cp, gint Np,
						   WBFMM_REAL *H03, 
						   WBFMM_REAL *H47, gint Lh,
						   WBFMM_REAL *trans, gint Ls,
						   WBFMM_REAL *work)

/*
  03: "lower" boxes (think of Morton index)
  47: "upper" boxes

  assumes all coefficients are densely packed in groups of eight
  coefficients of the same index, in order of child box Morton index

  contents of child coefficients Cc are overwritten
*/

{
  WBFMM_REAL *Cr, tmp ;
  WBFMM_REAL Epr[2], Epi[2], Enr[2], Eni[2], Cch[2], Sch[2] ;
  WBFMM_REAL Cnph, Snph, S0, Cnph0, Snph0 ;
  WBFMM_REAL Cnch, Snch, Cn3ch, Sn3ch ;
  WBFMM_REAL Cnch0, Snch0, Cn3ch0, Sn3ch0 ;
  gint nu, l, n, m, ic, ip, ih, str, i, sgn, offc, offp ;

  /*stride in number of complex elements per entry*/
  str = 8 ;
  /*used for trigonometric recursions*/
  Cch[0] = Sch[0] = Sch[1] = M_SQRT1_2 ; Cch[1] = -M_SQRT1_2 ;

  Cnph0 = 1.0 ; Snph0 = 0.0 ; S0 = 1.0 ;
  /*used to store the rotated parent coefficients*/
  Cr = &(work[8*2*wbfmm_coefficient_index_nm(Nc+1,0)]) ;

  /*rotate parent box coefficients using Cr as temporary storage*/
  for ( n = 0 ; n <= Np ; n ++ ) {
    Cnph = Cnph0 ; Snph = Snph0 ;
    for ( nu = -n ; nu <= n ; nu ++ ) {
      WBFMM_REAL buf[16] = {0.0} ;
      m = 0 ; 
      Epr[0] = Enr[0] = Epr[1] = Enr[1] = Cnph ;
      Epi[0] = Eni[0] = Snph ;
      Epi[1] = Eni[1] = -Epi[0] ;

      ip = wbfmm_coefficient_index_nm(n,m) ;
      ih = wbfmm_rotation_index_numn(nu,m,n) ;

      offp = 2*str*ip ;
      increment_buf_pc(Epr, Epi, &(Cp[offp]), H03[ih], H47[ih], buf) ;

      for ( m = 1 ; m <= n ; m ++ ) {
	sincos_recursion(Epr, Epi, Enr, Eni, Cch, Sch) ;

	ip = wbfmm_coefficient_index_nm(n,m) ;
	ih = wbfmm_rotation_index_numn(nu,m,n) ;
	offp = 2*str*ip ;

	increment_buf_pc(Epr, Epi, &(Cp[offp]), H03[ih], H47[ih], buf) ;

	ip = wbfmm_coefficient_index_nm(n,-m) ;
	ih = wbfmm_rotation_index_numn(-nu,m,n) ;
	offp = 2*str*ip ;

	increment_buf_pc(Enr, Eni, &(Cp[offp]), H03[ih], H47[ih], buf) ;
      }

      ip = wbfmm_coefficient_index_nm(n,nu) ; offp = 2*str*ip ;      
      for ( i = 0 ; i < 16 ; i ++ ) Cr[offp+i] = buf[i] ;

      /*trigonometric recursions taking note that C0 = 0*/
      tmp = Cnph ; Cnph = -Snph ; Snph = tmp ;
    }
    tmp = Cnph0 ; Cnph0 = Snph0 ; Snph0 = -tmp ;
  }

  /*Cr now contains the rotated Cp coefficients, apply coaxial
    translation to all coefficients to shift to centre of child boxes*/
  for ( n = 0 ; n <= Nc ; n ++ ) {
    for ( m = -n ; m <= n ; m ++ ) {
      ic = wbfmm_coefficient_index_nm(n,m) ; offc = 2*str*ic ;
      WBFMM_REAL *buf = &(work[offc]) ;
      memset(buf, 0, 16*sizeof(WBFMM_REAL)) ;
      /*loop on input and coefficients*/
      for ( l = ABS(m) ; l <= Np ; l ++ ) {
	ip = wbfmm_coefficient_index_nm(l,m) ; offp = 2*str*ip ;
	ih = wbfmm_coaxial_index(l, m, n) ;
	sgn = wbfmm_coaxial_index_sgn(l, m, n) ;
	for ( i = 0 ; i < 16 ; i ++ ) buf[i] += trans[ih]*sgn*Cr[offp+i] ;
      }
    }
  }

  /*work now contains the rotated coefficients shifted to the centres
    of the child boxes, which must be reverse rotated and summed into
    Cc, eight at a time*/
  Cch[0] = Cch[1] = 0.0 ; Sch[0] = -1 ; Sch[1] = 1.0 ;
  S0 = M_SQRT1_2 ;
  Cnch0 = Cn3ch0 = 1.0 ; Snch0 = Sn3ch0 = 0.0 ;
  for ( n = 0 ; n <= Nc ; n ++ ) {
    Cnch  = Cnch0  ; Snch  = Snch0 ;
    Cn3ch = Cn3ch0 ; Sn3ch = Sn3ch0 ;
    for ( nu = -n ; nu <= n ; nu ++ ) {
      WBFMM_REAL buf[16] = {0.0} ;

      m = 0 ; 
      Epr[0] = Enr[0] =  Cnch  ;
      Epi[0] = Eni[0] = -Snch  ;
      Epr[1] = Enr[1] =  Cn3ch ;
      Epi[1] = Eni[1] = -Sn3ch ;

      ic = wbfmm_coefficient_index_nm(n,m) ;
      ih = wbfmm_rotation_index_numn(nu,m,n) ;
      offc = 2*str*ic ;

      increment_cfft_pc(Epr, Epi, &(work[offc]), H03[ih], H47[ih], buf) ;

      for ( m = 1 ; m <= n ; m ++ ) {
	sincos_recursion(Epr, Epi, Enr, Eni, Cch, Sch) ;

	ic = wbfmm_coefficient_index_nm(n,m) ;
	ih = wbfmm_rotation_index_numn(nu,m,n) ;
        offc = 2*str*ic ;

	increment_cfft_pc(Epr, Epi, &(work[offc]), H03[ih], H47[ih], buf) ;

	ic = wbfmm_coefficient_index_nm(n,-m) ;
	ih = wbfmm_rotation_index_numn(-nu,m,n) ;
        offc = 2*str*ic ;

	increment_cfft_pc(Enr, Eni, &(work[offc]), H03[ih], H47[ih], buf) ;
      }
      tmp = Cnch ; Cnch = S0*(Cnch - Snch) ; Snch = S0*(tmp + Snch) ;
      tmp = Cn3ch ; Cn3ch = -S0*(Cn3ch + Sn3ch) ; Sn3ch = S0*(tmp - Sn3ch) ;

      /*shift buf into child boxes*/
      ic = wbfmm_coefficient_index_nm(n,nu) ; offc = 2*str*ic ;
      for ( i = 0 ; i < 16 ; i ++ ) Cc[offc+i] = buf[i] ;
    }

    tmp = Cnch0 ; 
    Cnch0 = S0*(Cnch0 + Snch0) ; Snch0 = S0*(-tmp + Snch0) ;
    tmp = Cn3ch0 ;
    Cn3ch0 = S0*(-Cn3ch0 + Sn3ch0) ; Sn3ch0 = -S0*(tmp + Sn3ch0) ;
  }

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_shift_angles_list4)(gint i, gint j, gint k,
						   WBFMM_REAL *th,
						   WBFMM_REAL *ph,
						   WBFMM_REAL *ch,
						   WBFMM_REAL *rs)

{
  gint idx = (i+3)*49+(j+3)*7+k+3 ;
  gint ip, str = 4 ;

  /*_wbfmm_shift_angles holds indices of entries in tables of pre-computed
    rotations and translations; 

    ch and ph indices are 1-based, to allow negative values for
    particular transformations, which compacts the lookup tables by
    about half

    translation distance rs must be multiplied by the box width at the
    appropriate level in the tree
*/

  ip = _wbfmm_shift_angles[str*idx+0] ;
  /*64 used as placeholder for i==j==k==0, to keep indexing consistent*/
  g_assert(ip < 49) ; 
  *th = _wbfmm_shifts_th[ip] ;

  ip = _wbfmm_shift_angles[str*idx+1] ;
  g_assert(ip != 64) ;
  *ph = (ip >= 0 ? _wbfmm_shifts_ph[ip-1] : -_wbfmm_shifts_ph[-1-ip]) ;

  ip = _wbfmm_shift_angles[str*idx+2] ;
  g_assert(ip != 64) ;
  *ch = (ip >= 0 ? _wbfmm_shifts_ph[ip-1] : -_wbfmm_shifts_ph[-1-ip]) ;

  ip = _wbfmm_shift_angles[str*idx+3] ;
  g_assert(ip >= 0 && ip < 19) ;
  *rs = _wbfmm_shifts_r[ip] ;

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_shift_angle_table_init)(void)

{
  /*
    note th[36] is rotation angle for child-parent upward shifts from
    lower half of parent box, and th[12] for shifts from the upper
    half
   */

  _wbfmm_shifts_th[ 0] = ACOS(-SQRT(9.0/9)) ;
  _wbfmm_shifts_th[ 1] = ACOS(-SQRT(9.0/10)) ;
  _wbfmm_shifts_th[ 2] = ACOS(-SQRT(9.0/11)) ;
  _wbfmm_shifts_th[ 3] = ACOS(-SQRT(4.0/5)) ;
  _wbfmm_shifts_th[ 4] = ACOS(-SQRT(9.0/13)) ;
  _wbfmm_shifts_th[ 5] = ACOS(-SQRT(4.0/6)) ;
  _wbfmm_shifts_th[ 6] = ACOS(-SQRT(9.0/14)) ;
  _wbfmm_shifts_th[ 7] = ACOS(-SQRT(9.0/17)) ;
  _wbfmm_shifts_th[ 8] = ACOS(-SQRT(9.0/18)) ;
  _wbfmm_shifts_th[ 9] = ACOS(-SQRT(9.0/19)) ;
  _wbfmm_shifts_th[10] = ACOS(-SQRT(4.0/9)) ;
  _wbfmm_shifts_th[11] = ACOS(-SQRT(9.0/22)) ;
  _wbfmm_shifts_th[12] = ACOS(-SQRT(9.0/27)) ;
  _wbfmm_shifts_th[13] = ACOS(-SQRT(4.0/13)) ;
  _wbfmm_shifts_th[14] = ACOS(-SQRT(4.0/14)) ;
  _wbfmm_shifts_th[15] = ACOS(-SQRT(4.0/17)) ;
  _wbfmm_shifts_th[16] = ACOS(-SQRT(1.0/5)) ;
  _wbfmm_shifts_th[17] = ACOS(-SQRT(4.0/22)) ;
  _wbfmm_shifts_th[18] = ACOS(-SQRT(1.0/6)) ;
  _wbfmm_shifts_th[19] = ACOS(-SQRT(1.0/9)) ;
  _wbfmm_shifts_th[20] = ACOS(-SQRT(1.0/10)) ;
  _wbfmm_shifts_th[21] = ACOS(-SQRT(1.0/11)) ;
  _wbfmm_shifts_th[22] = ACOS(-SQRT(1.0/14)) ;
  _wbfmm_shifts_th[23] = ACOS(-SQRT(1.0/19)) ;
  _wbfmm_shifts_th[24] = ACOS( SQRT(0.0/18)) ;

  _wbfmm_shifts_th[25] = ACOS(SQRT(1.0/19)) ;
  _wbfmm_shifts_th[26] = ACOS(SQRT(1.0/14)) ;
  _wbfmm_shifts_th[27] = ACOS(SQRT(1.0/11)) ;
  _wbfmm_shifts_th[28] = ACOS(SQRT(1.0/10)) ;
  _wbfmm_shifts_th[29] = ACOS(SQRT(1.0/9)) ;
  _wbfmm_shifts_th[30] = ACOS(SQRT(1.0/6)) ;
  _wbfmm_shifts_th[31] = ACOS(SQRT(4.0/22)) ;
  _wbfmm_shifts_th[32] = ACOS(SQRT(1.0/5)) ;
  _wbfmm_shifts_th[33] = ACOS(SQRT(4.0/17)) ;
  _wbfmm_shifts_th[34] = ACOS(SQRT(4.0/14)) ;
  _wbfmm_shifts_th[35] = ACOS(SQRT(4.0/13)) ;
  _wbfmm_shifts_th[36] = ACOS(SQRT(9.0/27)) ;
  _wbfmm_shifts_th[37] = ACOS(SQRT(9.0/22)) ;
  _wbfmm_shifts_th[38] = ACOS(SQRT(4.0/9)) ;
  _wbfmm_shifts_th[39] = ACOS(SQRT(9.0/19)) ;
  _wbfmm_shifts_th[40] = ACOS(SQRT(9.0/18)) ;
  _wbfmm_shifts_th[41] = ACOS(SQRT(9.0/17)) ;
  _wbfmm_shifts_th[42] = ACOS(SQRT(9.0/14)) ;
  _wbfmm_shifts_th[43] = ACOS(SQRT(4.0/6)) ;
  _wbfmm_shifts_th[44] = ACOS(SQRT(9.0/13)) ;
  _wbfmm_shifts_th[45] = ACOS(SQRT(4.0/5)) ;
  _wbfmm_shifts_th[46] = ACOS(SQRT(9.0/11)) ;
  _wbfmm_shifts_th[47] = ACOS(SQRT(9.0/10)) ;
  _wbfmm_shifts_th[48] = ACOS(SQRT(9.0/9)) ;

  _wbfmm_shifts_ph[ 0] = 0.0 ;
  _wbfmm_shifts_ph[ 1] = ATAN2(1,  3) ;
  _wbfmm_shifts_ph[ 2] = ATAN2(1,  2) ;
  _wbfmm_shifts_ph[ 3] = ATAN2(2,  3) ;
  _wbfmm_shifts_ph[ 4] = ATAN2(1,  1) ;
  _wbfmm_shifts_ph[ 5] = ATAN2(3,  2) ;
  _wbfmm_shifts_ph[ 6] = ATAN2(2,  1) ;
  _wbfmm_shifts_ph[ 7] = ATAN2(3,  1) ;
  _wbfmm_shifts_ph[ 8] = ATAN2(1,  0) ;
  _wbfmm_shifts_ph[ 9] = ATAN2(3, -1) ;
  _wbfmm_shifts_ph[10] = ATAN2(2, -1) ;
  _wbfmm_shifts_ph[11] = ATAN2(3, -2) ;
  _wbfmm_shifts_ph[12] = ATAN2(1, -1) ;
  _wbfmm_shifts_ph[13] = ATAN2(2, -3) ;
  _wbfmm_shifts_ph[14] = ATAN2(1, -2) ;
  _wbfmm_shifts_ph[15] = ATAN2(1, -3) ;
  _wbfmm_shifts_ph[16] = M_PI ;

  _wbfmm_shifts_r[ 0] = SQRT( 4.0) ;
  _wbfmm_shifts_r[ 1] = SQRT( 5.0) ;
  _wbfmm_shifts_r[ 2] = SQRT( 6.0) ;
  _wbfmm_shifts_r[ 3] = SQRT( 8.0) ;
  _wbfmm_shifts_r[ 4] = SQRT( 9.0) ;
  _wbfmm_shifts_r[ 5] = SQRT(10.0) ;
  _wbfmm_shifts_r[ 6] = SQRT(11.0) ;
  _wbfmm_shifts_r[ 7] = SQRT(12.0) ;
  _wbfmm_shifts_r[ 8] = SQRT(13.0) ;
  _wbfmm_shifts_r[ 9] = SQRT(14.0) ;
  _wbfmm_shifts_r[10] = SQRT(17.0) ;
  _wbfmm_shifts_r[11] = SQRT(18.0) ;
  _wbfmm_shifts_r[12] = SQRT(19.0) ;
  _wbfmm_shifts_r[13] = SQRT(22.0) ;
  _wbfmm_shifts_r[14] = SQRT(27.0) ;

  return 0 ;
}

wbfmm_shift_operators_t 
WBFMM_FUNCTION_NAME(*wbfmm_shift_operators_new)(guint L,
						gboolean bw,
						WBFMM_REAL *work)

{
  wbfmm_shift_operators_t *op = NULL ;
  WBFMM_REAL *rotations ;
  guint nerot ;
  gint i ;

  /*
    workspace size = wbfmm_element_number_rotation(2*L) elements
  */

  op = (wbfmm_shift_operators_t *)g_malloc0(sizeof(wbfmm_shift_operators_t)) ;

  op->size = sizeof(WBFMM_REAL) ;

  op->Lmax = L ;

  op->bw = bw ;
  
  if ( _wbfmm_shifts_th[1] == 0.0 ) 
    g_error("%s: rotation table not initiated; call "
	    "wbfmm_shift_angle_table_init()", 
	    __FUNCTION__) ;

  /*number of elements in rotation operators*/
  nerot = op->nerot = wbfmm_element_number_rotation(L) ;

  /*rotations are the same at all levels, so allocate them here*/
  rotations = op->rotations = g_malloc(49*nerot*sizeof(WBFMM_REAL)) ;

  /*first (0,\pi) rotations, as listed in _wbfmm_shifts_th*/
  for ( i = 0 ; i < 49 ; i ++ ) {
    WBFMM_FUNCTION_NAME(wbfmm_coefficients_H_rotation)(&(rotations[i*nerot]),
						       L, _wbfmm_shifts_th[i],
						       work) ;
  }

  for ( i = 0 ; i <= WBFMM_TREE_MAX_DEPTH ; i ++ ) {
    op->SR[i] = op->SS[i] = NULL ; op->L[i] = 0 ;
  }

  return op ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_shift_operators_coaxial_SR_init)
  (wbfmm_shift_operators_t *w, WBFMM_REAL D, guint level, guint L,
   WBFMM_REAL k, WBFMM_REAL *work)

{
  WBFMM_REAL wb, *SR ;
  guint nb, ne ;
  gint i ;

  /*
    initialize the (complex) singular to regular translation operators
    needed for the list 4 shifts at level `level' for a domain of total
    width D
  */

  nb = 1 << level ;
  wb = D/nb ;  

  ne = wbfmm_element_number_coaxial(L) ;
  if ( w->bw )
    SR = w->SR[level] = g_malloc0(2*WBFMM_SHIFTS_R_NUMBER*2*ne*
				  sizeof(WBFMM_REAL)) ;
  else
    SR = w->SR[level] = g_malloc0(  WBFMM_SHIFTS_R_NUMBER*2*ne*
				    sizeof(WBFMM_REAL)) ;
  w->L[level] = L ;

  for ( i = 0 ; i < WBFMM_SHIFTS_R_NUMBER ; i ++ ) {
    WBFMM_FUNCTION_NAME(wbfmm_coefficients_SR_coaxial)
      (&(SR[i*2*ne]), L, k*wb*_wbfmm_shifts_r[i], work) ;    
  }

  if ( !(w->bw) ) return 0 ;

  /*if we are using the backward translation algorithm on list 4 boxes,
    the negative translation operators follow the positive*/
  for ( i = 0 ; i < WBFMM_SHIFTS_R_NUMBER ; i ++ ) {
    WBFMM_FUNCTION_NAME(wbfmm_coefficients_SR_coaxial)
      (&(SR[(WBFMM_SHIFTS_R_NUMBER+i)*2*ne]), L, -k*wb*_wbfmm_shifts_r[i],
       work) ;    
  }

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_shift_operators_coaxial_SS_init)
     (wbfmm_shift_operators_t *w, WBFMM_REAL D, guint level, 
      guint L, WBFMM_REAL k, WBFMM_REAL *work)

{
  WBFMM_REAL wb, *SS ;
  guint nb, ne ;

  /*
    initialize the (real) singular to singular (also regular to
    regular) translation operators needed for the upward pass from
    level `level', to level-1, with expansion order L
  */

  nb = 1 << level ;
  wb = D/nb ;

  ne = wbfmm_element_number_coaxial(L) ;
  if ( w->bw )
    SS = w->SS[level] = g_malloc0(2*ne*sizeof(WBFMM_REAL)) ;
  else
    SS = w->SS[level] = g_malloc0(  ne*sizeof(WBFMM_REAL)) ;

  WBFMM_FUNCTION_NAME(wbfmm_coefficients_RR_coaxial)(SS, L, k*wb*0.5*SQRT(3.0), 
						     work) ;    
  if ( !(w->bw) ) return 0 ;

  WBFMM_FUNCTION_NAME(wbfmm_coefficients_RR_coaxial)(&(SS[ne]), L,
						     -k*wb*0.5*SQRT(3.0), 
						     work) ;      
  
  return 0 ;
}
