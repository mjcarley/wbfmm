/* This file is part of WBFMM, a Wide-Band Fast Multipole Method code
 *
 * Copyright (C) 2019, 2020 Michael Carley
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

/* #define WBFMM_CHECK_ISNAN */

#ifdef WBFMM_CHECK_ISNAN
#include <stdlib.h>

static gint check_isnan(char *name, WBFMM_REAL *f, gint n)

{
  gint i ;

  for ( i = 0 ; i < n ; i ++ ) {
    if ( isnan(f[i]) ) {
      fprintf(stderr, "%s: NaN at element %d of %d\n", name, i, n) ;
      exit(1) ;
    }
  }
  
  return 0 ;
}
#endif /*WBFMM_CHECK_ISNAN*/

static inline void _wbfmm_downward_pass_box(guint level, guint64 ip,
					    wbfmm_box_t *b,
					    guint Ns, guint Nr,
					    WBFMM_REAL *rotations,
					    guint nerot,
					    WBFMM_REAL *shifts,
					    guint necx,
					    WBFMM_REAL *wks, guint ncs,
					    WBFMM_REAL *wkr, guint ncr,
					    gint nq)

{
  guint idx4, ni ;
  guint64 ilist[378] ;
  gint j, iph, ith ;
  WBFMM_REAL ph, ch, *H, *Cx ;

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
    g_assert(ith >= 0) ;
    Cx = &(shifts[ith*necx]) ;

    /*rotation angles \phi and \chi*/
    iph =  _wbfmm_shift_angles[4*idx4+1] ;
    ph = (iph >= 0 ? _wbfmm_shifts_ph[iph-1] : -_wbfmm_shifts_ph[-1-iph]) ;
    iph =  _wbfmm_shift_angles[4*idx4+2] ;
    ch = (iph >= 0 ? _wbfmm_shifts_ph[iph-1] : -_wbfmm_shifts_ph[-1-iph]) ;
      
    /*rotate singular coefficients into wks*/
    WBFMM_FUNCTION_NAME(wbfmm_rotate_H)(wks, nq, b[ilist[2*j+0]].mps, 8*nq,
					Ns, nq, H, ph, ch, 0.0) ;
    /*translate into wkr*/
    WBFMM_FUNCTION_NAME(wbfmm_coaxial_translate)(wkr, nq, Nr, wks, nq, Ns, nq,
						 Cx, Nr, TRUE, 0.0) ;
    /*rotate regular coefficients into mpr*/
    WBFMM_FUNCTION_NAME(wbfmm_rotate_H)(b[ip].mpr, 8*nq, wkr, nq, Nr, nq,
					H, ch, ph, 1.0) ;
  }
  
  return ;
}

static inline void _wbfmm_shift_up(guint64 grid[], gint idx4,
				   WBFMM_REAL *shifts, guint necx,
				   gint Nr, gint Ns,
				   wbfmm_box_t *bp, WBFMM_REAL *target,
				   gint nq)
{
  WBFMM_REAL *Cx ;
  gint ith, ic ;

  if ( grid[idx4] == 0 ) return ;
  
  ic = grid[idx4] - 1 ;
  grid[idx4] = 0 ;

  if ( bp[ic].n == 0 ) return ;  
  
  ith = _wbfmm_shift_angles[4*idx4+3] ;
  Cx = &(shifts[(ith*2+0)*necx]) ;
  WBFMM_FUNCTION_NAME(wbfmm_coaxial_translate)(target, 8*nq, Nr,
					       bp[ic].mps, 8*nq, Ns, nq,
					       Cx, Nr, TRUE, 1.0) ;

  return ;
}

static inline void _wbfmm_shift_down(guint64 grid[], gint idx4,
				     WBFMM_REAL *shifts, guint necx,
				     gint Nr, gint Ns,
				     wbfmm_box_t *bp, WBFMM_REAL *target,
				     gint nq)
{
  WBFMM_REAL *Cx ;
  gint ith, ic ;

  if ( grid[idx4] == 0 ) return ;
  
  ic = grid[idx4] - 1 ;
  grid[idx4] = 0 ;

  if ( bp[ic].n == 0 ) return ;

  ith = _wbfmm_shift_angles[4*idx4+3] ;
  Cx = &(shifts[(ith*2+1)*necx]) ;
  WBFMM_FUNCTION_NAME(wbfmm_coaxial_translate)(target, 8*nq, Nr,
					       bp[ic].mps, 8*nq, Ns, nq,
					       Cx, Nr, TRUE, 1.0) ;
  return ;
}

static inline void _wbfmm_diagonal_shift(guint64 grid[], gint idx4,
					 WBFMM_REAL *rotations, guint nerot,
					 WBFMM_REAL *shifts, guint necx,
					 gint Nr, gint Ns,
					 wbfmm_box_t *bp, WBFMM_REAL *target,
					 WBFMM_REAL *wks, gint ncs,
					 WBFMM_REAL *wkr, gint ncr,
					 gint nq)

{
  WBFMM_REAL *H, *Cx, ch, ph, wt ;
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
  Cx = &(shifts[(2*ix+0)*necx]) ;
  
  /*rotation angles \phi and \chi*/
  ith =  _wbfmm_shift_angles[4*idx4+1] ;
  ph = (ith >= 0 ? _wbfmm_shifts_ph[ith-1] : -_wbfmm_shifts_ph[-1-ith]) ;
  ith =  _wbfmm_shift_angles[4*idx4+2] ;
  ch = (ith >= 0 ? _wbfmm_shifts_ph[ith-1] : -_wbfmm_shifts_ph[-1-ith]) ;
      
  wt = 0.0 ;
  if ( bp[ic].n != 0 ) {  
    /*rotate singular coefficients into wks*/
    WBFMM_FUNCTION_NAME(wbfmm_rotate_H)(wks, nq, bp[ic].mps, 8*nq, Ns, nq,
					H, ph, ch, 0.0) ;
    /*translate into wkr*/
    WBFMM_FUNCTION_NAME(wbfmm_coaxial_translate)(wkr, nq, Nr,
						 wks, nq, Ns, nq,
						 Cx, Nr, TRUE, wt) ;
    /* Cx, Nr, TRUE, 0.0) ; */
    wt = 1.0 ;
  }
  if ( grid[342-idx4] != 0 ) {
    ic = grid[342-idx4] - 1 ;
    grid[342-idx4] = 0 ;
    if ( bp[ic].n != 0 ) {  
      WBFMM_FUNCTION_NAME(wbfmm_rotate_H)(wks, nq, bp[ic].mps, 8*nq, Ns, nq,
					  H, ph, ch, 0.0) ;
      Cx = &(shifts[(2*ix+1)*necx]) ;
      WBFMM_FUNCTION_NAME(wbfmm_coaxial_translate)(wkr, nq, Nr,
						   wks, nq, Ns,
						   nq, Cx, Nr, TRUE, wt) ;
						   /* nq, Cx, Nr, TRUE, 1.0) ; */
      wt = 1.0 ;
    }
  }

  if ( wt == 0.0 ) return ;

  /*rotate regular coefficients into mpr*/
  WBFMM_FUNCTION_NAME(wbfmm_rotate_H)(target, 8*nq, wkr, nq, Nr, nq, H, ch, ph,
				      1.0) ;
  
  return ;
}

static inline void _wbfmm_diagonal_shift_3(guint64 grid[],
					   guint idx4f[], guint idx4b[],
					   WBFMM_REAL *rotations, guint nerot,
					   WBFMM_REAL *shifts, guint necx,
					   gint Nr, gint Ns,
					   wbfmm_box_t *bp, WBFMM_REAL *target,
					   WBFMM_REAL *wks, gint ncs,
					   WBFMM_REAL *wkr, gint ncr,
					   gint nq)

{
  WBFMM_REAL *H, *Cx, ch, ph, sc ;
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
    if ( grid[idx4f[i]] != 0 ) {
      ic = grid[idx4f[i]] - 1 ;
      grid[idx4f[i]] = 0 ;
      /*index of translation operator*/
      ix = _wbfmm_shift_angles[4*idx4f[i]+3] ;
      Cx = &(shifts[(2*ix+0)*necx]) ;
      if ( bp[ic].n != 0 ) {
	/*rotate singular coefficients into wks*/
	WBFMM_FUNCTION_NAME(wbfmm_rotate_H)(wks, nq, bp[ic].mps, 8*nq,
					    Ns, nq, H, ph, ch, 0.0) ;
	/*translate into wkr*/
	WBFMM_FUNCTION_NAME(wbfmm_coaxial_translate)(wkr, nq, Nr,
						     wks, nq, Ns,
						     nq, Cx, Nr, TRUE, sc) ;
	sc = 1.0 ;
      }
    }
    if ( grid[idx4b[i]] != 0 ) {
      ic = grid[idx4b[i]] - 1 ;
      grid[idx4b[i]] = 0 ;
      /*index of translation operator*/
      ix = _wbfmm_shift_angles[4*idx4f[i]+3] ;
      Cx = &(shifts[(2*ix+1)*necx]) ;
      if ( bp[ic].n != 0 ) {
	/*rotate singular coefficients into wks*/
	WBFMM_FUNCTION_NAME(wbfmm_rotate_H)(wks, nq, bp[ic].mps, 8*nq,
					    Ns, nq, H, ph, ch, 0.0) ;
	/*translate into wkr*/
	WBFMM_FUNCTION_NAME(wbfmm_coaxial_translate)(wkr, nq, Nr,
						     wks, nq, Ns, nq,
						     Cx, Nr, TRUE, sc) ;
	sc = 1.0 ;
      }
    }
  }

  if ( sc == 0.0 ) return ;

  /*rotate regular coefficients into mpr*/
  WBFMM_FUNCTION_NAME(wbfmm_rotate_H)(target, 8*nq, wkr, nq, Nr, nq, H, ch, ph,
				      1.0) ;
  
  return ;
}

static inline void _wbfmm_downward_pass_box_bw(guint level, guint64 ip,
					       wbfmm_box_t *b,
					       guint Ns, guint Nr,
					       WBFMM_REAL *rotations,
					       guint nerot,
					       WBFMM_REAL *shifts,
					       guint necx,
					       WBFMM_REAL *wks, guint ncs,
					       WBFMM_REAL *wkr, guint ncr,
					       gint nq)

{
  guint idx4, idx4f[2], idx4b[2] ;
  guint64 grid[343] = {0} ;

  /* if ( nq != 1 ) */
  /*   g_error("%s: not checked for nq (%d) > 1", __FUNCTION__, nq) ; */
  
  /*locate boxes in interaction list*/
  wbfmm_box_interaction_grid_4(level, ip, grid) ;
  /*deal with special cases combining rotations*/
  idx4f[0] = (-3+3)*49 + (-3+3)*7 + (-3+3) ;
  idx4f[1] = (-2+3)*49 + (-2+3)*7 + (-2+3) ;
  idx4b[1] = ( 2+3)*49 + ( 2+3)*7 + ( 2+3) ;
  idx4b[0] = ( 3+3)*49 + ( 3+3)*7 + ( 3+3) ;
  _wbfmm_diagonal_shift_3(grid, idx4f, idx4b,
			  rotations, nerot, shifts, necx,
			  Nr, Ns, b, b[ip].mpr, wks, ncs, wkr, ncr, nq) ;
  idx4f[0] = ( 3+3)*49 + (-3+3)*7 + (-3+3) ;
  idx4f[1] = ( 2+3)*49 + (-2+3)*7 + (-2+3) ;
  idx4b[1] = (-2+3)*49 + ( 2+3)*7 + ( 2+3) ;
  idx4b[0] = (-3+3)*49 + ( 3+3)*7 + ( 3+3) ;
  _wbfmm_diagonal_shift_3(grid, idx4f, idx4b,
			  rotations, nerot, shifts, necx,
			  Nr, Ns, b, b[ip].mpr, wks, ncs, wkr, ncr, nq) ;
  idx4f[0] = (-3+3)*49 + ( 3+3)*7 + (-3+3) ;
  idx4f[1] = (-2+3)*49 + ( 2+3)*7 + (-2+3) ;
  idx4b[1] = ( 2+3)*49 + (-2+3)*7 + ( 2+3) ;
  idx4b[0] = ( 3+3)*49 + (-3+3)*7 + ( 3+3) ;
  _wbfmm_diagonal_shift_3(grid, idx4f, idx4b,
			  rotations, nerot, shifts, necx,
			  Nr, Ns, b, b[ip].mpr, wks, ncs, wkr, ncr, nq) ;
  idx4b[0] = (-3+3)*49 + (-3+3)*7 + ( 3+3) ;
  idx4b[1] = (-2+3)*49 + (-2+3)*7 + ( 2+3) ;
  idx4f[1] = ( 2+3)*49 + ( 2+3)*7 + (-2+3) ;
  idx4f[0] = ( 3+3)*49 + ( 3+3)*7 + (-3+3) ;
  _wbfmm_diagonal_shift_3(grid, idx4f, idx4b,
			  rotations, nerot, shifts, necx,
			  Nr, Ns, b, b[ip].mpr, wks, ncs, wkr, ncr, nq) ;

  idx4f[0] = (-3+3)*49 + ( 0+3)*7 + ( 0+3) ;
  idx4f[1] = (-2+3)*49 + ( 0+3)*7 + ( 0+3) ;
  idx4b[1] = ( 2+3)*49 + ( 0+3)*7 + ( 0+3) ;
  idx4b[0] = ( 3+3)*49 + ( 0+3)*7 + ( 0+3) ;
  _wbfmm_diagonal_shift_3(grid, idx4f, idx4b,
			  rotations, nerot, shifts, necx,
			  Nr, Ns, b, b[ip].mpr, wks, ncs, wkr, ncr, nq) ;
  idx4f[0] = ( 0+3)*49 + (-3+3)*7 + ( 0+3) ;
  idx4f[1] = ( 0+3)*49 + (-2+3)*7 + ( 0+3) ;
  idx4b[1] = ( 0+3)*49 + ( 2+3)*7 + ( 0+3) ;
  idx4b[0] = ( 0+3)*49 + ( 3+3)*7 + ( 0+3) ;
  _wbfmm_diagonal_shift_3(grid, idx4f, idx4b,
			  rotations, nerot, shifts, necx,
			  Nr, Ns, b, b[ip].mpr, wks, ncs, wkr, ncr, nq) ;
    
  for ( idx4 = 0 ; idx4 < 168 ; idx4 ++ ) {
    _wbfmm_diagonal_shift(grid, idx4, rotations, nerot, shifts, necx,
			  Nr, Ns, b, b[ip].mpr, wks, ncs, wkr, ncr, nq) ;
  }
  /*vertically displaced above and below target box, no rotations
    required*/
  idx4 = (0+3)*49 + (0+3)*7 + (-3+3) ;
  _wbfmm_shift_up(grid, idx4, shifts, necx, Nr, Ns, b, b[ip].mpr, nq) ;
  idx4 = (0+3)*49 + (0+3)*7 + (-2+3) ;
  _wbfmm_shift_up(grid, idx4, shifts, necx, Nr, Ns, b, b[ip].mpr, nq) ;
  idx4 = (0+3)*49 + (0+3)*7 + ( 2+3) ;
  _wbfmm_shift_down(grid, idx4, shifts, necx, Nr, Ns, b, b[ip].mpr, nq) ;
  idx4 = (0+3)*49 + (0+3)*7 + ( 3+3) ;
  _wbfmm_shift_down(grid, idx4, shifts, necx, Nr, Ns, b, b[ip].mpr, nq) ;

  for ( idx4 = 175 ; idx4 < 343 ; idx4 ++ ) {
    _wbfmm_diagonal_shift(grid, idx4, rotations, nerot, shifts, necx,
			  Nr, Ns, b, b[ip].mpr, wks, ncs, wkr, ncr, nq) ;
  }

  return ;
}

static gpointer downward_pass_thread(gpointer idata)

{
  gpointer *data = idata ;
  gpointer *pdata = data[1] ;
  guint nb, Ns, Nr, nerot, necx, ncs, ncr, level ;
  WBFMM_REAL *rotations, *shifts, *work, *wkr, *wks ;
  wbfmm_box_t *bp ;
  wbfmm_tree_t *t ;
  wbfmm_shift_operators_t *op ;
  gint i, nq, nth ;
  guint64 ip ;

  i = GPOINTER_TO_INT(data[0]) ;

  t     = pdata[WBFMM_DOWNWARD_PASS_TREE] ;
  op    = pdata[WBFMM_DOWNWARD_PASS_OP] ;
  level = *((guint *)pdata[WBFMM_DOWNWARD_PASS_LEVEL]) ;
  work  = (WBFMM_REAL *)pdata[WBFMM_DOWNWARD_PASS_WORK] ;
  nq    = *((gint *)pdata[WBFMM_DOWNWARD_PASS_NQ]) ;
  nth   = *((gint *)pdata[WBFMM_DOWNWARD_PASS_NTHREAD]) ;

  /*number of boxes at this level*/
  nb = 1 << 3*(level) ;

  /*singular and regular expansion orders at this level*/
  Ns = t->order_s[level] ; Nr = t->order_r[level] ;
  ncs = wbfmm_coefficient_number(Ns) ;
  ncr = wbfmm_coefficient_number(Nr) ;
  wks = work ; wkr = &(wks[2*ncs*nq]) ;

  /*boxes at this level (parent)*/
  bp = t->boxes[level] ;

  nerot = op->nerot ;
  /*rotation and shift operators*/
  rotations = (WBFMM_REAL *)(op->rotations) ;
  shifts    = (WBFMM_REAL *)(op->SR[level]) ;

  /*number of elements in translation operators*/
  necx  = 2*wbfmm_element_number_coaxial(op->L[level]) ;

  /* fprintf(stderr, "level %u; thread %d\n", level, i) ; */

  wks = &(work[i*2*(ncr+ncs)*nq]) ; wkr = &(wks[2*ncs*nq]) ;
  
  if ( op->bw )
    for ( ip = i ; ip < nb ; ip += nth ) {
      _wbfmm_downward_pass_box_bw(level, ip, bp, Ns, Nr, rotations, nerot,
  				  shifts, necx, wks, ncs, wkr, ncr, nq) ;
    }
  else
    for ( ip = i ; ip < nb ; ip += nth ) {
      _wbfmm_downward_pass_box(level, ip, bp, Ns, Nr, rotations, nerot,
  			       shifts, necx, wks, ncs, wkr, ncr, nq) ;
    }

  return NULL ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_downward_pass_ref)(wbfmm_tree_t *t,
						  wbfmm_shift_operators_t *op,
						  guint level,
						  WBFMM_REAL *work,
						  gint nthreads)

{
  guint nb, Ns, Nr, nerot, necx, ncs, ncr, Nc, Np ;
  guint64 ic, ip ;
  gint nq ;
  wbfmm_box_t *bp, *bc ;
  WBFMM_REAL *rotations, *shifts, *wks, *wkr, *H03, *H47, *trans ;

  nq = wbfmm_tree_source_size(t) ;
  
  g_assert(level > 1) ;
  g_assert(wbfmm_tree_problem(t) == WBFMM_PROBLEM_HELMHOLTZ) ;

  /*number of boxes at this level*/
  nb = 1 << 3*(level) ;

  /*singular and regular expansion orders at this level*/
  Ns = t->order_s[level] ; Nr = t->order_r[level] ;
  ncs = wbfmm_coefficient_number(Ns) ;
  ncr = wbfmm_coefficient_number(Nr) ;
  wks = work ; wkr = &(wks[2*ncs*nq]) ;

  /*boxes at this level (parent)*/
  bp = t->boxes[level] ;

  nerot = op->nerot ;
  /*rotation and shift operators*/
  rotations = (WBFMM_REAL *)(op->rotations) ;
  shifts    = (WBFMM_REAL *)(op->SR[level]) ;

  /*number of elements in translation operators*/
  necx  = 2*wbfmm_element_number_coaxial(op->L[level]) ;

#ifdef _OPENMP
  if ( nthreads == 0 ) {
    if ( op->bw )
      for ( ip = 0 ; ip < nb ; ip ++ ) {
	_wbfmm_downward_pass_box_bw(level, ip, bp, Ns, Nr, rotations, nerot,
				    shifts, necx, wks, ncs, wkr, ncr, nq) ;
      }
    else
      for ( ip = 0 ; ip < nb ; ip ++ ) {
	_wbfmm_downward_pass_box(level, ip, bp, Ns, Nr, rotations, nerot,
				 shifts, necx, wks, ncs, wkr, ncr, nq) ;
      }
  } else {
    GThread *threads[WBFMM_THREAD_NUMBER_MAX] ;
    gint nth, i ;
    guint nproc ;
    gpointer data[WBFMM_DOWNWARD_PASS_DATA_SIZE],
      main_data[2*WBFMM_THREAD_NUMBER_MAX] ;

    g_assert(nthreads <= WBFMM_THREAD_NUMBER_MAX) ;
    data[WBFMM_DOWNWARD_PASS_WORK] = work ;
    data[WBFMM_DOWNWARD_PASS_LEVEL] = &level ;
    data[WBFMM_DOWNWARD_PASS_NQ] = &nq ;
    data[WBFMM_DOWNWARD_PASS_NTHREAD] = &nth ;
    data[WBFMM_DOWNWARD_PASS_TREE] = t ;
    data[WBFMM_DOWNWARD_PASS_OP] = op ;

    nproc = g_get_num_processors() ;
    if ( nthreads < 0 ) nth = nproc ; else nth = nthreads ;
    if ( nth > nproc )
      g_error("%s: not enough processes (%u) for requested number of "
	      "threads (%d)", __FUNCTION__, nproc, nth) ;
    for ( i = 0 ; i < nth ; i ++ ) {
      main_data[2*i+0] = GINT_TO_POINTER(i) ;
      main_data[2*i+1] = data ;
      threads[i] = g_thread_new(NULL, downward_pass_thread, &(main_data[2*i])) ;
    }
    /*make sure all threads complete before we move on*/
    for ( i = 0 ; i < nth ; i ++ ) g_thread_join(threads[i]) ;
  }
#else /*_OPENMP*/
  if ( op->bw )
    for ( ip = 0 ; ip < nb ; ip ++ ) {
      _wbfmm_downward_pass_box_bw(level, ip, bp, Ns, Nr, rotations, nerot,
				  shifts, necx, wks, ncs, wkr, ncr, nq) ;
    }
  else
    for ( ip = 0 ; ip < nb ; ip ++ ) {
      _wbfmm_downward_pass_box(level, ip, bp, Ns, Nr, rotations, nerot,
			       shifts, necx, wks, ncs, wkr, ncr, nq) ;
    }
#endif /*_OPENMP*/
  /* /\*interaction list 4, loop on boxes at this level*\/ */
  /* if ( op->bw ) */
  /*   for ( ip = 0 ; ip < nb ; ip ++ ) { */
  /*     _wbfmm_downward_pass_box_bw(level, ip, bp, Ns, Nr, rotations, nerot, */
  /* 				  shifts, necx, wks, ncs, wkr, ncr, nq) ; */
  /*   } */
  /* else */
  /*   for ( ip = 0 ; ip < nb ; ip ++ ) { */
  /*     _wbfmm_downward_pass_box(level, ip, bp, Ns, Nr, rotations, nerot, */
  /* 			       shifts, necx, wks, ncs, wkr, ncr, nq) ; */
  /*   } */

  /*no downward shift at the deepest level*/
  if ( level == t-> depth ) return 0 ;

  /*rotation operators for parent-child shifts*/
  H03 = &(rotations[36*nerot]) ;
  H47 = &(rotations[12*nerot]) ;
  Np = t->order_r[level  ] ;
  Nc = t->order_r[level+1] ;
  bc = t->boxes[level+1] ;
  trans = (WBFMM_REAL *)(op->SS[level+1]) ;
  for ( ip = 0 ; ip < nb ; ip ++ ) {
    ic = wbfmm_box_first_child(ip) ;
    WBFMM_FUNCTION_NAME(wbfmm_parent_child_shift)((WBFMM_REAL *)(bc[ic].mpr),
						  Nc,
						  (WBFMM_REAL *)(bp[ip].mpr),
						  Np,
						  H03, H47, Np,
						  trans, Np,
						  wbfmm_tree_source_size(t),
						  work) ;
#ifdef WBFMM_CHECK_ISNAN
    g_assert_not_reached() ; /*check no longer needed*/
    check_isnan("downward pass, parent-child shift",
		bc[ic].mpr, 2*wbfmm_coefficient_index_nm(Nc+1,0)) ;
#endif /*WBFMM_CHECK_ISNAN*/
    
  }

  return 0 ;
}

gint WBFMM_FUNCTION_NAME(wbfmm_upward_pass)(wbfmm_tree_t *t,
					    wbfmm_shift_operators_t *op,
					    guint level, WBFMM_REAL *work)

{
  guint np, Np, Nc ;
  guint64 ip, ic;
  wbfmm_box_t *bp, *bc ;
  WBFMM_REAL *H03, *H47, *trans, *rotations ;

  g_assert(level > 1) ;
  g_assert(wbfmm_tree_problem(t) == WBFMM_PROBLEM_HELMHOLTZ) ;

  /*number of parent boxes into which to shift child data*/
  np = 1 << 3*(level-1) ;

  /*expansion orders at parent and child levels*/
  Np = t->order_s[level-1] ;
  Nc = t->order_s[level  ] ;

  /*parent and child boxes*/
  bp = t->boxes[level-1] ;
  bc = t->boxes[level  ] ;

  /*rotation and shift operators*/
  rotations = (WBFMM_REAL *)(op->rotations) ;
  /* H03 = &(rotations[36*(op->nerot)]) ; */
  /* H47 = &(rotations[12*(op->nerot)]) ; */
  H03 = &(rotations[12*(op->nerot)]) ;
  H47 = &(rotations[36*(op->nerot)]) ;
  trans = (WBFMM_REAL *)(op->SS[level]) ;

  g_assert(H03 != NULL) ;
  g_assert(H47 != NULL) ;
  g_assert(trans != NULL) ;
  for ( ip = 0 ; ip < np ; ip ++ ) {
    /* if ( bp[ip].n != 0 ) { */
      /*locate first child of parent box*/
      ic = wbfmm_box_first_child(ip) ;
      WBFMM_FUNCTION_NAME(wbfmm_child_parent_shift)((WBFMM_REAL *)(bp[ip].mps),
						    Np,
						    (WBFMM_REAL *)(bc[ic].mps),
						    Nc,
						    H03, H47, Np, 
						    trans, Np,
						    wbfmm_tree_source_size(t),
						    work) ;
    /* } */
  }

  return 0 ;
}
