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

/**
 * @defgroup indexing Indexing and lookup operations
 *
 * @brief Indexing functions for accessing tree data structures
 *
 * Functions for indexing and lookup in tree data structures,
 * including finding neighbours and interaction lists, based on the
 * methods of Gumerov, Duraiswami, and Borovikov, Data Structures,
 * Optimal Choice of Parameters, and Complexity Results for
 * Generalized Multilevel Fast Multipole Methods in d Dimensions, 2003
 *
 * http://users.umiacs.umd.edu/~gumerov/PDFs/cs-tr-4458.pdf
 * 
 * Code for Morton indexing operations is taken from:
 * https://www.forceflow.be/2013/10/07/morton-encodingdecoding-through-bit-interleaving-implementations/
 *
 */

/* @{ */


#define _is_neighbour(_i,_j,_k,_i0,_j0,_k0)		   \
  ((ABSDIFF((_i),(_i0))<=1) && (ABSDIFF((_j),(_j0))<=1) && \
   (ABSDIFF((_k),(_k0))<=1))

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>

#include <glib.h>

#include <wbfmm.h>

#include "wbfmm-private.h"

/* 
   Based on public domain code from
   https://stackoverflow.com/questions/49748864/morton-reverse-encoding-for-a-3d-grid
*/

/* 
   Morton encoding in binary (components 21-bit: 0..2097151)
   0zyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyx 
*/


static inline void morton_decode(guint64 m, 
				 guint32 *xto, guint32 *yto, guint32 *zto)

{
    const guint64  
      mask0 = 
      BITMASK_0000000001000001000001000001000001000001000001000001000001000001,
      mask1 = 
      BITMASK_0000001000001000001000001000001000001000001000001000001000001000,
      mask2 = 
      BITMASK_0001000000000000000000000000000000000000000000000000000000000000,
      mask3 = 
      BITMASK_0000000000000011000000000011000000000011000000000011000000000011,
      mask4 = 
      BITMASK_0000000111000000000011000000000011000000000011000000000011000000,
      mask5 = 
      BITMASK_0000000000000000000000000000000000001111000000000000000000001111,
      mask6 = 
      BITMASK_0000000000000000000000001111000000000000000000001111000000000000,
      mask7 = 
      BITMASK_0000000000011111000000000000000000000000000000000000000000000000,
      mask8 = 
      BITMASK_0000000000000000000000000000000000000000000000000000000011111111,
      mask9 = 
      BITMASK_0000000000000000000000000001111111111111000000000000000000000000;
    guint64  
      x = m,
      y = m >> 1,
      z = m >> 2;

    /* 000c00c00c00c00c00c00c00c00c00c00c00c00c00c00c00c00c00c00c00c00c */
    x = (x & mask0) | ((x & mask1) >> 2) | ((x & mask2) >> 4);
    y = (y & mask0) | ((y & mask1) >> 2) | ((y & mask2) >> 4);
    z = (z & mask0) | ((z & mask1) >> 2) | ((z & mask2) >> 4);
    /* 0000000ccc0000cc0000cc0000cc0000cc0000cc0000cc0000cc0000cc0000cc */
    x = (x & mask3) | ((x & mask4) >> 4);
    y = (y & mask3) | ((y & mask4) >> 4);
    z = (z & mask3) | ((z & mask4) >> 4);
    /* 00000000000ccccc00000000cccc00000000cccc00000000cccc00000000cccc */
    x = (x & mask5) | ((x & mask6) >> 8) | ((x & mask7) >> 16);
    y = (y & mask5) | ((y & mask6) >> 8) | ((y & mask7) >> 16);
    z = (z & mask5) | ((z & mask6) >> 8) | ((z & mask7) >> 16);
    /* 000000000000000000000000000ccccccccccccc0000000000000000cccccccc */
    x = (x & mask8) | ((x & mask9) >> 16);
    y = (y & mask8) | ((y & mask9) >> 16);
    z = (z & mask8) | ((z & mask9) >> 16);
    /* 0000000000000000000000000000000000000000000ccccccccccccccccccccc */
    if (xto) *xto = x;
    if (yto) *yto = y;
    if (zto) *zto = z;
}

static inline guint64 morton_encode(guint32 xsrc, guint32 ysrc, guint32 zsrc)

{
    const guint64  
      mask0 = 
      BITMASK_0000000001000001000001000001000001000001000001000001000001000001,
      mask1 =
      BITMASK_0000001000001000001000001000001000001000001000001000001000001000,
      mask2 =
      BITMASK_0001000000000000000000000000000000000000000000000000000000000000,
      mask3 =
      BITMASK_0000000000000011000000000011000000000011000000000011000000000011,
      mask4 =
      BITMASK_0000000111000000000011000000000011000000000011000000000011000000,
      mask5 =
      BITMASK_0000000000000000000000000000000000001111000000000000000000001111,
      mask6 =
      BITMASK_0000000000000000000000001111000000000000000000001111000000000000,
      mask7 =
      BITMASK_0000000000011111000000000000000000000000000000000000000000000000,
      mask8 =
      BITMASK_0000000000000000000000000000000000000000000000000000000011111111,
      mask9 =
      BITMASK_0000000000000000000000000001111111111111000000000000000000000000;
    guint64  x = xsrc, y = ysrc, z = zsrc ;
    /* 0000000000000000000000000000000000000000000ccccccccccccccccccccc */
    x = (x & mask8) | ((x << 16) & mask9) ;
    y = (y & mask8) | ((y << 16) & mask9) ;
    z = (z & mask8) | ((z << 16) & mask9) ;
    /* 000000000000000000000000000ccccccccccccc0000000000000000cccccccc */
    x = (x & mask5) | ((x << 8) & mask6) | ((x << 16) & mask7) ;
    y = (y & mask5) | ((y << 8) & mask6) | ((y << 16) & mask7) ;
    z = (z & mask5) | ((z << 8) & mask6) | ((z << 16) & mask7) ;
    /* 00000000000ccccc00000000cccc00000000cccc00000000cccc00000000cccc */
    x = (x & mask3) | ((x << 4) & mask4) ;
    y = (y & mask3) | ((y << 4) & mask4) ;
    z = (z & mask3) | ((z << 4) & mask4) ;
    /* 0000000ccc0000cc0000cc0000cc0000cc0000cc0000cc0000cc0000cc0000cc */
    x = (x & mask0) | ((x << 2) & mask1) | ((x << 4) & mask2) ;
    y = (y & mask0) | ((y << 2) & mask1) | ((y << 4) & mask2) ;
    z = (z & mask0) | ((z << 2) & mask1) | ((z << 4) & mask2) ;
    /* 000c00c00c00c00c00c00c00c00c00c00c00c00c00c00c00c00c00c00c00c00c */
    return x | (y << 1) | (z << 2) ;
}

/*
  find the coordinates of a box at a given level containing the point x,
  coordinates generated using wbfmm_point_index3d(...)
*/

guint64 wbfmm_point_locate_box(guint64 x, guint level)

{
  guint64 i ;

  i = x >> 3*(20-level) ;

  return i ;
}

gint wbfmm_point_from_index(guint64 i, guint32 *x, guint32 *y, guint32 *z)

{
  morton_decode(i, x, y, z) ;

  return 0 ;
}

/** 
 * Generate a Morton index for a box with corner at integer
 * coordinates (i,j,k).
 * 
 * @param i x index of bottom left hand corner;
 * @param j y index of bottom left hand corner;
 * @param k z index of bottom left hand corner.
 * 
 * @return Morton index for (\a i, \a j, \a k).
 */

guint64 wbfmm_box_index(guint32 i, guint32 j, guint32 k)

{
  guint64 idx ;

  idx = morton_encode(i, j, k) ;

  return idx ;
}

/** 
 * Compute indices for bottom left hand corner of box defined by its
 * Morton index, as generated by ::wbfmm_box_index
 * 
 * @param idx index of box corner;
 * @param i on output, x index of bottom left hand corner of box;
 * @param j on output, y index of bottom left hand corner of box;
 * @param k on output, z index of bottom left hand corner of box.
 * 
 * @return 0 on success.
 */

gint wbfmm_box_location(guint64 idx, guint32 *i, guint32 *j, guint32 *k)

{
  morton_decode(idx, i, j, k) ;

  return 0 ;
}

guint64 wbfmm_box_parent(guint64 idx)

{
  guint64 p ;

  p = idx >> 3 ;

  return p ;
}

guint64 wbfmm_box_first_child(guint64 idx)

{
  guint64 c ;

  c = idx << 3 ;

  return c ;
}


gint wbfmm_box_neighbours(guint level, guint64 idx, guint64 *neighbours)

/*
  indices of (potential) neighbours at the same level as box with
  index idx,

  return number of neighbours in list on exit
*/

{
  gint n ;
  guint32 i, j, k, i0, i1, j0, j1, k0, k1, ii, jj, kk, nbox ;

  morton_decode(idx, &i, &j, &k) ;
  /*number of boxes per side on this level*/
  nbox = 1 << level ;
  i0 = j0 = k0 = 0 ; i1 = j1 = k1 = nbox - 1 ;

  if ( i > 1 ) i0 = i - 1 ; 
  if ( j > 1 ) j0 = j - 1 ; 
  if ( k > 1 ) k0 = k - 1 ; 

  if ( i < nbox - 2 ) i1 = i + 1 ;
  if ( j < nbox - 2 ) j1 = j + 1 ;
  if ( k < nbox - 2 ) k1 = k + 1 ;

  n = 0 ;
  for ( ii = i0 ; ii <= i1 ; ii ++ ) {
    for ( jj = j0 ; jj <= j1 ; jj ++ ) {
      for ( kk = k0 ; kk <= k1 ; kk ++ ) {
	neighbours[n] = wbfmm_box_index(ii, jj, kk) ; n ++ ;
      }
    }
  }

  return n ;
}

static gint compare_wbfmm_shift_angles(gconstpointer i1, gconstpointer i2)

{
  const guint64 *ii1 = i1, *ii2 = i2 ;
  gint s1, s2 ;

  s1 = ii1[1] ; s2 = ii2[1] ;

  /*compare based on entries in _wbfmm_shift_angles*/
  if ( _wbfmm_shift_angles[4*s1+0] < _wbfmm_shift_angles[4*s2+0] ) return -1 ;
  if ( _wbfmm_shift_angles[4*s1+0] > _wbfmm_shift_angles[4*s2+0] ) return  1 ;

  if ( _wbfmm_shift_angles[4*s1+1] < _wbfmm_shift_angles[4*s2+1] ) return -1 ;
  if ( _wbfmm_shift_angles[4*s1+1] > _wbfmm_shift_angles[4*s2+1] ) return  1 ;

  if ( _wbfmm_shift_angles[4*s1+2] < _wbfmm_shift_angles[4*s2+2] ) return -1 ;
  if ( _wbfmm_shift_angles[4*s1+2] > _wbfmm_shift_angles[4*s2+2] ) return  1 ;

  return 0 ;
}

/** 
 * @brief Find the local interaction list for a specified box
 *
 * Find the indices of boxes on a given level of a tree which interact
 * directly with a specified box (list 4 in Gumerov and Duraiswami's
 * notation). These are boxes which are children of neighbours of the
 * parent of the specified box, and separated from it by at least one
 * box. On exit \a list contains entries made up of two integers, a
 * box index and the index for looking up rotation and translation
 * operations.
 * 
 * @param level tree level for list;
 * @param idx index of box whose interaction list is to be found;
 * @param list on exit contains list of interacting boxes and specification
 * of rotation required;
 * @param sort if TRUE, sort \a list so that boxes with the same rotation
 * angles are grouped together
 * 
 * @return number of entries in \a list
 */

gint wbfmm_box_interaction_list_4(guint level, guint64 idx, 
				  guint64 *list, gboolean sort)

{
  gint n, dx, dy, dz ;
  guint32 i, j, k, i0, i1, j0, j1, k0, k1, ii, jj, kk, nbox ;
  guint64 ishift ;
  guint32 ic, jc, kc ;
  
  morton_decode(idx, &i, &j, &k) ;
  /*number of boxes per side on this level*/
  nbox = 1 << level ;

  /*reduce to even part (coordinates of bottom corner of parent box)*/
  i0 = i - i%2 ; j0 = j - j%2 ; k0 = k - k%2 ; 

  i1 = MIN(i0+3,nbox-1) ; j1 = MIN(j0+3,nbox-1) ; k1 = MIN(k0+3,nbox-1) ;

  if ( i0 > 0 ) i0 -= 2 ;
  if ( j0 > 0 ) j0 -= 2 ;
  if ( k0 > 0 ) k0 -= 2 ;

  n = 0 ;
  for ( ii = i0 ; ii <= i1 ; ii ++ ) {
    for ( jj = j0 ; jj <= j1 ; jj ++ ) {
      for ( kk = k0 ; kk <= k1 ; kk ++ ) {
	if ( !_is_neighbour(i,j,k,ii,jj,kk) ) {
	  list[2*n+0] = wbfmm_box_index(ii, jj, kk) ;
	  wbfmm_box_location(list[2*n+0], &ic, &jc, &kc) ;
	  g_assert((ic == ii) && (jc == jj) && (kc == kk)) ;
	  /*index into shift angle table*/
	  dx = ( ii > i ? ii - i : -(gint)(i-ii)) ;
	  dy = ( jj > j ? jj - j : -(gint)(j-jj)) ;
	  dz = ( kk > k ? kk - k : -(gint)(k-kk)) ;
	  ishift = (guint64)((dx+3)*49+(dy+3)*7+dz+3) ;
	  list[2*n+1] = ishift ;
	  /* fprintf(stderr, "%lu %lu (%d %d %d)\n", */
	  /* 	  list[2*n+0], list[2*n+1], dx, dy, dz) ; */
	  n ++ ;
	}
      }
    }
  }

  g_assert(n < 190) ;

  if ( !sort ) return n ;

  /*sort list to group boxes with the same rotations together*/
  qsort(list, n, 2*sizeof(guint64), compare_wbfmm_shift_angles) ;

  return n ;
}

gint wbfmm_box_interaction_index(gint i, gint j, gint k)

/*
  index of rotation angles for a given list 4 box, defined by
  -3 <= i,j,k <= 3, 
*/

{
  gint idx = (i+3)*49+(j+3)*7+k+3 ;
  
  g_assert(idx >= 0 && idx < 343) ;

  return _wbfmm_shift_angles[4*idx+0] ;
}

gint wbfmm_box_interaction_grid_4(guint level, guint64 idx, guint64 grid[])

{
  guint32 i, j, k, i0, i1, j0, j1, k0, k1, ii, jj, kk, nbox, idxg ;
    
  memset(grid, 0, 343*sizeof(guint64)) ;
  
  morton_decode(idx, &i, &j, &k) ;
  /*number of boxes per side on this level*/
  nbox = 1 << level ;

  /*reduce to even part (coordinates of bottom corner of parent box)*/
  i0 = i - i%2 ; j0 = j - j%2 ; k0 = k - k%2 ; 
  
  i1 = MIN(i0+3,nbox-1) ; j1 = MIN(j0+3,nbox-1) ; k1 = MIN(k0+3,nbox-1) ;

  if ( i0 > 0 ) i0 -= 2 ;
  if ( j0 > 0 ) j0 -= 2 ;
  if ( k0 > 0 ) k0 -= 2 ;

  for ( ii = i0 ; ii <= i1 ; ii ++ ) {
    for ( jj = j0 ; jj <= j1 ; jj ++ ) {
      for ( kk = k0 ; kk <= k1 ; kk ++ ) {
	if ( !_is_neighbour(i,j,k,ii,jj,kk) ) {
	  idxg = (ii-i+3)*49 + (jj-j+3)*7 + kk - k + 3 ;
	  g_assert(idxg < 343 && idxg >=0) ;
	  grid[idxg] = wbfmm_box_index(ii, jj, kk) + 1 ;
	}
      }
    }
  }

  return 0 ;
}

/* @} */
