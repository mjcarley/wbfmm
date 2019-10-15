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
  calculations for the Laplace problem
*/

#ifdef _HAVE_CONFIG_H_
#include <config.h>
#endif /*_HAVE_CONFIG_H_*/

#include <math.h>
#include <string.h>

#include <glib.h>

#include <wbfmm.h>

#include "wbfmm-private.h"

extern gint _wbfmm_shift_angles[] ;
extern WBFMM_REAL _wbfmm_shifts_ph[], _wbfmm_shifts_ch[], _wbfmm_shifts_r[] ;

#if 0
static gboolean coefficient_check_laplace(WBFMM_REAL *C, gint cstr,
					  gint N, gint nq)

{
  gint n, m, idx ;

  n = 0 ; m = 0 ; idx = n*n ;
  if ( C[cstr*idx+0] != C[cstr*idx+1] ) {
    fprintf(stderr, "%s: mismatch at n=%d, m=%d, idx=%d, (%lg,%lg)\n",
	    __FUNCTION__, n, m, idx, C[cstr*idx+0], C[cstr*idx+1]) ;
    return FALSE ;
  }

  n = 1 ; m = 0 ; idx = n*n ;
  if ( C[cstr*idx+0] != C[cstr*idx+1] ) {
    fprintf(stderr, "%s: mismatch at n=%d, m=%d, idx=%d, (%lg,%lg)\n",
	    __FUNCTION__, n, m, idx, C[cstr*idx+0], C[cstr*idx+1]) ;
    return FALSE ;
  }
  
  m = 1 ;
  idx = wbfmm_index_laplace_nm(n,m) ;
  if ( C[cstr*idx+0] != C[cstr*idx+1] ) {
    fprintf(stderr, "%s: mismatch at n=%d, m=%d, idx=%d, (%lg,%lg)\n",
	    __FUNCTION__, n, m, idx, C[cstr*idx+0], C[cstr*idx+1]) ;
    return FALSE ;
  }
  
  for ( n = 2 ; n <= N ; n ++ ) {
    m = 0 ; idx = n*n ;
    if ( C[cstr*idx+0] != C[cstr*idx+1] ) {
    fprintf(stderr, "%s: mismatch at n=%d, m=%d, idx=%d, (%lg,%lg)\n",
	    __FUNCTION__, n, m, idx, C[cstr*idx+0], C[cstr*idx+1]) ;
      return FALSE ;
    }
    for ( m = 1 ; m <= n ; m ++ ) {
      idx = wbfmm_index_laplace_nm(n,m) ;
      if ( C[cstr*idx+0] != C[cstr*idx+1] ) {
	fprintf(stderr, "%s: mismatch at n=%d, m=%d, idx=%d, (%lg,%lg)\n",
		__FUNCTION__, n, m, idx, C[cstr*idx+0], C[cstr*idx+1]) ;
	return FALSE ;
      }
    }
  }
  
  return TRUE ;
}
#endif

gint WBFMM_FUNCTION_NAME(wbfmm_laplace_downward_pass)(wbfmm_tree_t *t,
						      wbfmm_shift_operators_t
						      *op,
						      guint level,
						      WBFMM_REAL *work)

{
  guint nb, Ns, Nr, ni, nerot, ncs, ncr, idx4, j, Nc, Np, nq ;
  gint ith, iph ;
  guint64 ip, ic, ilist[378] ;
  wbfmm_box_t *bp, *bc ;
  WBFMM_REAL *rotations, ph, ch, *H, *wks, *wkr ;
  WBFMM_REAL *H03, *H47, wb, r ;
  /* gint i ; */
  
  g_assert(level > 1) ;
  g_assert(wbfmm_tree_problem(t) == WBFMM_PROBLEM_LAPLACE) ;

  nq = wbfmm_tree_source_size(t) ;
  /*number of boxes at this level*/
  nb = 1 << 3*(level) ;

  /*parent box width*/
  wb = wbfmm_tree_width(t)/(1 << (level)) ;
    
  /*singular and regular expansion orders at this level*/
  Ns = t->order_s[level] ; Nr = t->order_r[level] ;
  ncs = (Ns+1)*(Ns+1) ;
  ncr = (Nr+1)*(Nr+1) ;
  wks = work ; wkr = &(wks[nq*ncs]) ;

  /*boxes at this level (parent)*/
  bp = t->boxes[level] ;

  nerot = op->nerot ;
  /*rotation operators*/
  rotations = (WBFMM_REAL *)(op->rotations) ;

  /*interaction list 4, loop on boxes at this level*/
  for ( ip = 0 ; ip < nb ; ip ++ ) {
    /*locate boxes in interaction list*/
    ni = wbfmm_box_interaction_list_4(level, ip, ilist, TRUE) ;
    /*loop on interaction list and compute SR-shifted fields*/
    for ( j = 0 ; j < ni ; j ++ ) {
      /*find entry in shift lookup tables*/
      idx4 = ilist[2*j+1] ;
      /*index of rotation operator*/
      ith = _wbfmm_shift_angles[4*idx4+0] ;
      H = &(rotations[ith*nerot]) ;
      
      /*index of translation operator*/
      ith = _wbfmm_shift_angles[4*idx4+3] ;
      r = wb*_wbfmm_shifts_r[ith] ;

      /*rotation angles \phi and \chi*/
      iph =  _wbfmm_shift_angles[4*idx4+1] ;
      ph = (iph >= 0 ? _wbfmm_shifts_ph[iph-1] : -_wbfmm_shifts_ph[-1-iph]) ;
      iph =  _wbfmm_shift_angles[4*idx4+2] ;
      ch = (iph >= 0 ? _wbfmm_shifts_ph[iph-1] : -_wbfmm_shifts_ph[-1-iph]) ;
      
      /*clear workspace*/
      memset(wks, 0, nq*(ncs+ncr)*sizeof(WBFMM_REAL)) ;
      /*rotate singular coefficients into wks*/
      WBFMM_FUNCTION_NAME(wbfmm_rotate_H_laplace)(wks, nq,
						  bp[ilist[2*j+0]].mps, 8*nq,
						  Ns, nq, H, ph, ch) ;
      /*translate into wkr*/
      WBFMM_FUNCTION_NAME(wbfmm_coaxial_translate_SR_laplace)(wkr, nq, Nr,
							      wks, nq, Ns, 
							      nq, r) ;
      /*rotate regular coefficients into mpr*/
      WBFMM_FUNCTION_NAME(wbfmm_rotate_H_laplace)(bp[ip].mpr,
						  8*nq,
						  wkr, nq,
						  Nr, nq,
						  H, ch, ph) ;
      /* g_assert(coefficient_check_laplace(bp[ip].mpr, 8*nq, Nr, nq)) ; */
    }

  }

  /*no downward shift at the deepest level*/
  if ( level == t-> depth ) return 0 ;

  /*rotation operators for parent-child shifts*/
  H03 = &(rotations[12*nerot]) ;
  H47 = &(rotations[36*nerot]) ;
  Np = t->order_r[level  ] ;
  Nc = t->order_r[level+1] ;
  bc = t->boxes[level+1] ;

  for ( ip = 0 ; ip < nb ; ip ++ ) {
    ic = wbfmm_box_first_child(ip) ;
    WBFMM_FUNCTION_NAME(wbfmm_parent_child_shift_laplace)
      ((WBFMM_REAL *)(bc[ic].mpr), Nc,
       (WBFMM_REAL *)(bp[ip].mpr), Np,
       nq, H03, H47, Np,
       wb, work) ;
  }

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_laplace_upward_pass)(wbfmm_tree_t *t,
						    wbfmm_shift_operators_t *op,
						    guint level,
						    WBFMM_REAL *work)

{
  guint np, Np, Nc ;
  guint64 ip, ic;
  wbfmm_box_t *bp, *bc ;
  WBFMM_REAL *H03, *H47, *rotations, wb ;

  g_assert(level > 1) ;
  g_assert(wbfmm_tree_problem(t) == WBFMM_PROBLEM_LAPLACE) ;

  /*number of parent boxes into which to shift child data*/
  np = 1 << 3*(level-1) ;

  /*expansion orders at parent and child levels*/
  Np = t->order_s[level-1] ;
  Nc = t->order_s[level  ] ;

  /*parent and child boxes*/
  bp = t->boxes[level-1] ;
  bc = t->boxes[level  ] ;

  /*child box width*/
  wb = wbfmm_tree_width(t)/(1 << (level)) ;
  
  /*rotation and shift operators*/
  rotations = (WBFMM_REAL *)(op->rotations) ;
  H03 = &(rotations[36*(op->nerot)]) ;
  H47 = &(rotations[12*(op->nerot)]) ;

  g_assert(H03 != NULL) ;
  g_assert(H47 != NULL) ;

  for ( ip = 0 ; ip < np ; ip ++ ) {
    /*locate first child of parent box*/
    ic = wbfmm_box_first_child(ip) ;
    WBFMM_FUNCTION_NAME(wbfmm_child_parent_shift_laplace)((WBFMM_REAL *)
							  (bp[ip].mps),
							  Np,
							  (WBFMM_REAL *)
							  (bc[ic].mps),
							  Nc,
							  t->nq, 
							  H03, H47, Np, 
							  wb, work) ;
  }

  return 0 ;
}
