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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <math.h>
#include <string.h>

#include <glib.h>

#include <wbfmm.h>

#include "wbfmm-private.h"

static inline void _wbfmm_laplace_downward_pass(guint level, guint64 ip,
						wbfmm_box_t *bp,
						WBFMM_REAL wb,
						guint Ns, guint Nr,
						WBFMM_REAL *rotations,
						guint nerot,
						WBFMM_REAL *wks,
						WBFMM_REAL *wkr,
						gint nq)
{
  guint ni, j, idx4 ;
  gint ith, iph ;
  guint64 ilist[378] ;
  WBFMM_REAL *C, *Cn, *H, ph, ch, r ;
  
  /*locate boxes in interaction list*/
  ni = wbfmm_box_interaction_list_4(level, ip, ilist, TRUE) ;
  C = (WBFMM_REAL *)(bp[ip].mpr) ;

  /*loop on interaction list and compute SR-shifted fields*/
  for ( j = 0 ; j < ni ; j ++ ) {
    /*coefficients of neighbour box*/
    Cn = (WBFMM_REAL *)(bp[ilist[2*j+0]].mps) ;
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
    
    /*rotate singular coefficients into wks*/
    WBFMM_FUNCTION_NAME(wbfmm_laplace_rotate_H)(wks, nq,
						Cn, 8*nq,
						Ns, nq, H, ph, ch, 0.0) ;
    /*translate into wkr*/
    WBFMM_FUNCTION_NAME(wbfmm_laplace_coaxial_translate_SR)(wkr, nq, Nr,
							    wks, nq, Ns, 
							    nq, r, 0.0) ;
    /*rotate regular coefficients into mpr*/
    WBFMM_FUNCTION_NAME(wbfmm_laplace_rotate_H)(C, 8*nq,
						wkr, nq,
						Nr, nq,
						H, ch, ph, 1.0) ;
  }

  return ;
}

static inline void _wbfmm_diagonal_shift(guint64 grid[], gint idx4,
					 WBFMM_REAL *rotations, guint nerot,
					 gint Nr, gint Ns,
					 wbfmm_box_t *bp, WBFMM_REAL *target,
					 WBFMM_REAL *wks,
					 WBFMM_REAL *wkr,
					 WBFMM_REAL wb,
					 gint nq)

{
  WBFMM_REAL *H, ch, ph, r ;
  gint ith, ic, ix ;

  if ( grid[idx4] == 0 ) return ;

  /*index of box in list 4 grid*/
  ic = grid[idx4] - 1 ;
  /*mark box as processed*/
  grid[idx4] = 0 ;
  
  /*index of rotation operator*/
  ith = _wbfmm_shift_angles[4*idx4+0] ;
  H = &(rotations[ith*nerot]) ;
  
  /*index of translation operator*/
  ix = _wbfmm_shift_angles[4*idx4+3] ;
  r = wb*_wbfmm_shifts_r[ix] ;
  
  /*rotation angles \phi and \chi*/
  ith =  _wbfmm_shift_angles[4*idx4+1] ;
  ph = (ith >= 0 ? _wbfmm_shifts_ph[ith-1] : -_wbfmm_shifts_ph[-1-ith]) ;
  ith =  _wbfmm_shift_angles[4*idx4+2] ;
  ch = (ith >= 0 ? _wbfmm_shifts_ph[ith-1] : -_wbfmm_shifts_ph[-1-ith]) ;
      
  /*rotate singular coefficients into wks*/
  WBFMM_FUNCTION_NAME(wbfmm_laplace_rotate_H)(wks, nq,
					      bp[ic].mps, 8*nq,
					      Ns, nq, H, ph, ch, 0.0) ;
  /*translate into wkr*/
  WBFMM_FUNCTION_NAME(wbfmm_laplace_coaxial_translate_SR)(wkr, nq, Nr,
							  wks, nq, Ns,
							  nq, r, 0.0) ;
  if ( grid[342-idx4] != 0 ) {
    ic = grid[342-idx4] - 1 ;
    WBFMM_FUNCTION_NAME(wbfmm_laplace_rotate_H)(wks, nq,
						bp[ic].mps, 8*nq,
						Ns, nq, H, ph, ch, 0.0) ;
    WBFMM_FUNCTION_NAME(wbfmm_laplace_coaxial_translate_SR)(wkr, nq, Nr,
							    wks, nq, Ns,
							    nq, -r, 1.0) ;
    grid[342-idx4] = 0 ;
  }
  /*rotate regular coefficients into target*/
  WBFMM_FUNCTION_NAME(wbfmm_laplace_rotate_H)(target, 8*nq,
					      wkr, nq, Nr, nq, H, ch, ph, 1.0) ;
  
  return ;
}

static inline void _wbfmm_shift_up(guint64 grid[], gint idx4,
				   gint Nr, gint Ns,
				   wbfmm_box_t *bp, WBFMM_REAL *target,
				   WBFMM_REAL wb, gint nq)
{
  WBFMM_REAL r ;
  gint ith, ic ;

  if ( grid[idx4] == 0 ) return ;
  
  ic = grid[idx4] - 1 ;
  ith = _wbfmm_shift_angles[4*idx4+3] ;
  r = wb*_wbfmm_shifts_r[ith] ;
  WBFMM_FUNCTION_NAME(wbfmm_laplace_coaxial_translate_SR)(target, 8*nq, Nr,
							  bp[ic].mps, 8*nq, Ns, 
							  nq, r, 1.0) ;
  grid[idx4] = 0 ;

  return ;
}

static inline void _wbfmm_shift_down(guint64 grid[], gint idx4,
				     gint Nr, gint Ns,
				     wbfmm_box_t *bp, WBFMM_REAL *target,
				     WBFMM_REAL wb, gint nq)
{
  WBFMM_REAL r ;
  gint ith, ic ;

  if ( grid[idx4] == 0 ) return ;
  
  ic = grid[idx4] - 1 ;
  ith = _wbfmm_shift_angles[4*idx4+3] ;
  r = wb*_wbfmm_shifts_r[ith] ;
  WBFMM_FUNCTION_NAME(wbfmm_laplace_coaxial_translate_SR)(target, 8*nq, Nr,
							  bp[ic].mps, 8*nq, Ns, 
							  nq, -r, 1.0) ;
  grid[idx4] = 0 ;

  return ;
}

static inline void _wbfmm_diagonal_shift_3(guint64 grid[],
					   guint idx4f[], guint idx4b[],
					   WBFMM_REAL *rotations, guint nerot,
					   gint Nr, gint Ns,
					   wbfmm_box_t *bp, WBFMM_REAL *target,
					   WBFMM_REAL *wks, WBFMM_REAL *wkr,
					   WBFMM_REAL wb, gint nq)

{
  WBFMM_REAL *H, ch, ph, sc, r ;
  gint ith, ic, ix, i ;

  if ( grid[idx4f[0]] == 0 && grid[idx4f[1]] == 0 &&
       grid[idx4b[0]] == 0 && grid[idx4b[1]] == 0 )
    return ;

  /*index of rotation operator, assuming indices are all compatible*/
  ith = _wbfmm_shift_angles[4*idx4f[0]+0] ;
  H = &(rotations[ith*nerot]) ;
  /*rotation angles \phi and \chi*/
  ith =  _wbfmm_shift_angles[4*idx4f[0]+1] ;
  ph = (ith >= 0 ? _wbfmm_shifts_ph[ith-1] : -_wbfmm_shifts_ph[-1-ith]) ;
  ith =  _wbfmm_shift_angles[4*idx4f[0]+2] ;
  ch = (ith >= 0 ? _wbfmm_shifts_ph[ith-1] : -_wbfmm_shifts_ph[-1-ith]) ;

  sc = 0.0 ;
  for ( i = 0 ; i < 2 ; i ++ ) {
    /*index of translation operator*/
    ix = _wbfmm_shift_angles[4*idx4f[i]+3] ;
    r = wb*_wbfmm_shifts_r[ix] ;
    if ( grid[idx4f[i]] != 0 ) {
      ic = grid[idx4f[i]] - 1 ;

      /*rotate singular coefficients into wks*/
      WBFMM_FUNCTION_NAME(wbfmm_laplace_rotate_H)(wks, nq,
						  bp[ic].mps, 8*nq,
						  Ns, nq, H, ph, ch, 0.0) ;
      /*translate into wkr*/
      WBFMM_FUNCTION_NAME(wbfmm_laplace_coaxial_translate_SR)(wkr, nq, Nr,
							      wks, nq, Ns,
							      nq, r, sc) ;
      sc = 1.0 ;
      grid[idx4f[i]] = 0 ;
    }
    if ( grid[idx4b[i]] != 0 ) {
      ic = grid[idx4b[i]] - 1 ;
      /*rotate singular coefficients into wks*/
      WBFMM_FUNCTION_NAME(wbfmm_laplace_rotate_H)(wks, nq,
						  bp[ic].mps, 8*nq,
						  Ns, nq, H, ph, ch, 0.0) ;
      /*translate into wkr*/
      WBFMM_FUNCTION_NAME(wbfmm_laplace_coaxial_translate_SR)(wkr, nq, Nr,
							      wks, nq, Ns,
							      nq, -r, sc) ;
      sc = 1.0 ;
      grid[idx4b[i]] = 0 ;
    }
  }

  /*rotate regular coefficients into mpr*/
  WBFMM_FUNCTION_NAME(wbfmm_laplace_rotate_H)(target, 8*nq,
					      wkr, nq, Nr, nq, H, ch, ph, 1.0) ;

  return ;
}

static inline void _wbfmm_laplace_downward_pass_bw(guint level, guint64 ip,
						   wbfmm_box_t *bp,
						   WBFMM_REAL wb,
						   guint Ns, guint Nr,
						   WBFMM_REAL *rotations,
						   guint nerot,
						   WBFMM_REAL *wks,
						   WBFMM_REAL *wkr,
						   gint nq)
{
  guint64 grid[343] = {0} ;
  guint idx4, idx4f[2], idx4b[2] ;
  
  /*locate boxes in interaction list*/
  wbfmm_box_interaction_grid_4(level, ip, grid) ;

  /*deal with special cases combining rotations*/
  idx4f[0] = (-3+3)*49 + (-3+3)*7 + (-3+3) ;
  idx4f[1] = (-2+3)*49 + (-2+3)*7 + (-2+3) ;
  idx4b[1] = ( 2+3)*49 + ( 2+3)*7 + ( 2+3) ;
  idx4b[0] = ( 3+3)*49 + ( 3+3)*7 + ( 3+3) ;
  _wbfmm_diagonal_shift_3(grid, idx4f, idx4b, rotations, nerot,
			  Nr, Ns, bp, bp[ip].mpr, wks, wkr, wb, nq) ;
  idx4f[0] = ( 3+3)*49 + (-3+3)*7 + (-3+3) ;
  idx4f[1] = ( 2+3)*49 + (-2+3)*7 + (-2+3) ;
  idx4b[1] = (-2+3)*49 + ( 2+3)*7 + ( 2+3) ;
  idx4b[0] = (-3+3)*49 + ( 3+3)*7 + ( 3+3) ;
  _wbfmm_diagonal_shift_3(grid, idx4f, idx4b, rotations, nerot,
			  Nr, Ns, bp, bp[ip].mpr, wks, wkr, wb, nq) ;
  idx4f[0] = (-3+3)*49 + ( 3+3)*7 + (-3+3) ;
  idx4f[1] = (-2+3)*49 + ( 2+3)*7 + (-2+3) ;
  idx4b[1] = ( 2+3)*49 + (-2+3)*7 + ( 2+3) ;
  idx4b[0] = ( 3+3)*49 + (-3+3)*7 + ( 3+3) ;
  _wbfmm_diagonal_shift_3(grid, idx4f, idx4b, rotations, nerot,
			  Nr, Ns, bp, bp[ip].mpr, wks, wkr, wb, nq) ;
  idx4b[0] = (-3+3)*49 + (-3+3)*7 + ( 3+3) ;
  idx4b[1] = (-2+3)*49 + (-2+3)*7 + ( 2+3) ;
  idx4f[1] = ( 2+3)*49 + ( 2+3)*7 + (-2+3) ;
  idx4f[0] = ( 3+3)*49 + ( 3+3)*7 + (-3+3) ;
  _wbfmm_diagonal_shift_3(grid, idx4f, idx4b, rotations, nerot,
			  Nr, Ns, bp, bp[ip].mpr, wks, wkr, wb, nq) ;
  idx4f[0] = (-3+3)*49 + ( 0+3)*7 + ( 0+3) ;
  idx4f[1] = (-2+3)*49 + ( 0+3)*7 + ( 0+3) ;
  idx4b[1] = ( 2+3)*49 + ( 0+3)*7 + ( 0+3) ;
  idx4b[0] = ( 3+3)*49 + ( 0+3)*7 + ( 0+3) ;
  _wbfmm_diagonal_shift_3(grid, idx4f, idx4b, rotations, nerot,
			  Nr, Ns, bp, bp[ip].mpr, wks, wkr, wb, nq) ;
  idx4f[0] = ( 0+3)*49 + (-3+3)*7 + ( 0+3) ;
  idx4f[1] = ( 0+3)*49 + (-2+3)*7 + ( 0+3) ;
  idx4b[1] = ( 0+3)*49 + ( 2+3)*7 + ( 0+3) ;
  idx4b[0] = ( 0+3)*49 + ( 3+3)*7 + ( 0+3) ;
  _wbfmm_diagonal_shift_3(grid, idx4f, idx4b, rotations, nerot,
			  Nr, Ns, bp, bp[ip].mpr, wks, wkr, wb, nq) ;

  for ( idx4 = 0 ; idx4 < 168 ; idx4 ++ ) {
    _wbfmm_diagonal_shift(grid, idx4, rotations, nerot,
			  Nr, Ns, bp, bp[ip].mpr, wks, wkr, wb, nq) ;
  }

  /*vertically displaced above and below target box, no rotations
    required*/
  idx4 = (0+3)*49 + (0+3)*7 + (-3+3) ;
  _wbfmm_shift_up(grid, idx4, Nr, Ns, bp, bp[ip].mpr, wb, nq) ;
  idx4 = (0+3)*49 + (0+3)*7 + (-2+3) ;
  _wbfmm_shift_up(grid, idx4, Nr, Ns, bp, bp[ip].mpr, wb, nq) ;
  idx4 = (0+3)*49 + (0+3)*7 + ( 2+3) ;
  _wbfmm_shift_down(grid, idx4, Nr, Ns, bp, bp[ip].mpr, wb, nq) ;
  idx4 = (0+3)*49 + (0+3)*7 + ( 3+3) ;
  _wbfmm_shift_down(grid, idx4, Nr, Ns, bp, bp[ip].mpr, wb, nq) ;
    
  for ( idx4 = 175 ; idx4 < 343 ; idx4 ++ ) {
    _wbfmm_diagonal_shift(grid, idx4, rotations, nerot,
  			  Nr, Ns, bp, bp[ip].mpr, wks, wkr, wb, nq) ;
  }
  
  return ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_laplace_downward_pass)(wbfmm_tree_t *t,
						      wbfmm_shift_operators_t
						      *op,
						      guint level,
						      WBFMM_REAL *work)

{
  guint nb, Ns, Nr, nerot, ncr, Nc, Np, nq ;
  guint64 ip, ic ;
  wbfmm_box_t *bp, *bc ;
  WBFMM_REAL *rotations, *wks, *wkr ;
  WBFMM_REAL *H03, *H47, wb ;
  
  g_assert(level > 1) ;
  g_assert(wbfmm_tree_problem(t) == WBFMM_PROBLEM_LAPLACE) ;

  nq = wbfmm_tree_source_size(t) ;
  /*number of boxes at this level*/
  nb = 1 << 3*(level) ;

  /*parent box width*/
  wb = wbfmm_tree_width(t)/(1 << (level)) ;
    
  /*singular and regular expansion orders at this level*/
  Ns = t->order_s[level] ; Nr = t->order_r[level] ;
  ncr = (Nr+1)*(Nr+1) ;
  wkr = work ; wks = &(wkr[nq*ncr]) ;

  /*boxes at this level (parent)*/
  bp = t->boxes[level] ;

  nerot = op->nerot ;
  /*rotation operators*/
  rotations = (WBFMM_REAL *)(op->rotations) ;

  /*interaction list 4, loop on boxes at this level*/
  if ( op->bw )
    for ( ip = 0 ; ip < nb ; ip ++ ) {
      _wbfmm_laplace_downward_pass_bw(level, ip, bp, wb, Ns, Nr, rotations,
				      nerot, wks, wkr, nq) ;
    }
  else
    for ( ip = 0 ; ip < nb ; ip ++ ) {
      _wbfmm_laplace_downward_pass(level, ip, bp, wb, Ns, Nr, rotations,
				   nerot, wks, wkr, nq) ;
  }

  /*no downward shift at the deepest level*/
  if ( level == t-> depth ) return 0 ;

  /*rotation operators for parent-child shifts*/
  /* H03 = &(rotations[12*nerot]) ; */
  /* H47 = &(rotations[36*nerot]) ; */
  H03 = &(rotations[36*nerot]) ;
  H47 = &(rotations[12*nerot]) ;
  Np = t->order_r[level  ] ;
  Nc = t->order_r[level+1] ;
  bc = t->boxes[level+1] ;

  for ( ip = 0 ; ip < nb ; ip ++ ) {
    ic = wbfmm_box_first_child(ip) ;
    WBFMM_FUNCTION_NAME(wbfmm_laplace_parent_child_shift)
      ((WBFMM_REAL *)(bc[ic].mpr), Nc,
       (WBFMM_REAL *)(bp[ip].mpr), Np,
       nq, H03, H47, Np, wb, work) ;
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
  /* H03 = &(rotations[36*(op->nerot)]) ; */
  /* H47 = &(rotations[12*(op->nerot)]) ; */
  H03 = &(rotations[12*(op->nerot)]) ;
  H47 = &(rotations[36*(op->nerot)]) ;

  g_assert(H03 != NULL) ;
  g_assert(H47 != NULL) ;

  for ( ip = 0 ; ip < np ; ip ++ ) {
    /*locate first child of parent box*/
    ic = wbfmm_box_first_child(ip) ;
    WBFMM_FUNCTION_NAME(wbfmm_laplace_child_parent_shift)((WBFMM_REAL *)
							  (bp[ip].mps),
							  Np,
							  (WBFMM_REAL *)
							  (bc[ic].mps),
							  Nc,
							  wbfmm_tree_source_size(t), 
							  H03, H47, Np, 
							  wb, work) ;
  }

  return 0 ;
}
