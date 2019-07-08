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

#ifdef _HAVE_CONFIG_H_
#include <config.h>
#endif /*_HAVE_CONFIG_H_*/

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
