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
#include <glib.h>
#include <string.h>
#include <stdio.h>
#include <inttypes.h>

#include <wbfmm.h>

#include "wbfmm-private.h"

#define _DATA_WIDTH       8
#define _DATA_TREE        0
/* #define _DATA_ */


/**
 * @file   tree.c
 * @author Michael Carley <ensmjc@rpc-ensmjc.bath.ac.uk>
 * @date   Mon Jun 24 14:51:21 2019
 * 
 * @brief  
 * 
 * 
 */

/** 
 * Add a new level to an existing octree. The function assigns memory
 * for, and initializes, a new layer of boxes of type ::wbfmm_box_t
 *
 * @ingroup boxes
 * 
 * @param t an existing ::wbfmm_tree_t
 * 
 * @return 0 on success.
 */

gint wbfmm_tree_add_level(wbfmm_tree_t *t)

{
  guint nb ;

  if ( t->depth >= WBFMM_TREE_MAX_DEPTH ) 
    g_error("%s: tree already at maximum depth %u", 
	    __FUNCTION__, WBFMM_TREE_MAX_DEPTH) ;

  t->depth ++ ;
  nb = 1 << (3*(t->depth)) ;
  t->boxes[t->depth] = (wbfmm_box_t *)g_malloc0(nb*sizeof(wbfmm_box_t)) ;

  return 0 ;
}


/** 
 *
 * @ingroup targets
 * 
 * @param l a ::wbfmm_target_list_t generated using ::wbfmm_target_list_new
 * or ::wbfmm_target_list_new_f;
 * @param field a field definition made up of ::wbfmm_field_t
 * 
 * @return 0 on success.
 */

gint wbfmm_target_list_coefficients_init(wbfmm_target_list_t *l,
					 guint field)

{
  wbfmm_tree_t *t = wbfmm_target_list_tree(l) ;
  gint nr, nc ;

  /*order of regular expansions in leaf boxes*/
  nr = t->order_r[t->depth] ;
  l->field = field ;
  
  switch ( wbfmm_tree_problem(t) ) {
  default:
    g_error("%s: unrecognized problem %u",
	    __FUNCTION__, wbfmm_tree_problem(t)) ;
    break ;
  case WBFMM_PROBLEM_LAPLACE:
    nc = (nr+1)*(nr+1) ;
    break ;
  case WBFMM_PROBLEM_HELMHOLTZ:
    nc = 2*wbfmm_coefficient_index_nm(nr+1, 0) ;
    break ;
  }

  switch ( field ) {
  default:
    g_error("%s: unrecognized field definition (%u)", __FUNCTION__, field) ;
    break ;
  case WBFMM_FIELD_SCALAR:   nc *= 1 ; break ;
  case WBFMM_FIELD_GRADIENT: nc *= 3 ; break ; 
  }
  
  l->nc = nc ;
  l->cfft = g_malloc0(nc*wbfmm_target_list_point_number_max(l)*(l->size)) ;

  return 0 ;
}

/** 
 *
 * @ingroup boxes
 * 
 * Set the regular expansion coefficients at a level to zero, for when
 * a tree is being reused.
 *
 * @param t a ::wbfmm_tree_t
 * @param level a level of \a t
 * 
 * @return 0 on success.
 */

gint wbfmm_tree_coefficients_zero(wbfmm_tree_t *t, guint level)

{
  gint nb, nc ;

  nb = 1 << (3*level) ;
  nc = t->order_r[level] ;
  /* nc = (nc+1)*(nc+1) ; */
  nc = wbfmm_coefficient_number(nc) ;
  if ( wbfmm_tree_problem(t) == WBFMM_PROBLEM_HELMHOLTZ ) nc *= 2 ;
  memset(t->mpr[level], 0, nb*nc*(t->nq)*(t->size)) ;
  /* nc = t->order_s[level] ; */
  /* nc = (nc+1)*(nc+1) ; */
  /* memset(t->mps[level], 0, nb*nc*(t->nq)*(t->size)) ; */

  return 0 ;
}
