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
 * @file   wbfmm.h
 * @author Michael Carley <ensmjc@rpc-ensmjc.bath.ac.uk>
 * @date   Mon Jun 24 10:28:46 2019
 * 
 * @brief  Header for Wide Band FMM library
 * 
 * 
 */


#ifndef _WBFMM_H_INCLUDED_
#define _WBFMM_H_INCLUDED_

#include <stdio.h>

#define WBFMM_INDEX_SCALE (1LU << 63)
#define WBFMM_INDEX_SHIFT (1U << 20)

/**
 * @struct wbfmm_box_t
 * @ingroup boxes
 *
 * Data type for octree boxes
 *
 */

typedef struct {
  guint32 i, /**< index of first source point in box */
      n ;  /**< number of points in box */
  gpointer mps, /**< pointer to singular multipole expansion data */
    mpr ;   /**< pointer to regular multipole expansion data */

} wbfmm_box_t ;

/*maximum depth of an octree*/
#define WBFMM_TREE_MAX_DEPTH 16

/**
 * @struct wbfmm_tree_t
 * @ingroup boxes
 *
 * Data type for octrees
 *
 */

typedef struct {
  wbfmm_box_t *boxes[WBFMM_TREE_MAX_DEPTH+1] ; 
  /**< arrays of boxes at each level */
  guint 
  maxpoints, /**< maximum number of points in tree */
    npoints, /**< number of points in tree */
    *ip, /**< indices of points, sorted by Morton index */
    depth, /**< depth of tree */
    order_s[WBFMM_TREE_MAX_DEPTH+1], 
  /**< order of singular expansions at each level */
    order_r[WBFMM_TREE_MAX_DEPTH+1] ; 
  /**< order of regular expansions at each level */
  gchar 
  x[24], /**< origin of tree domain cube */
    *points ; /**< point coordinates */
  gpointer 
  *mps[WBFMM_TREE_MAX_DEPTH+1], 
  /**< singular expansion data at each level */
    *mpr[WBFMM_TREE_MAX_DEPTH+1] ; 
  /**< regular expansion data at each level */
  gsize pstr ; /**< stride in point data */
  gdouble D ; /**< width of domain cube */
} wbfmm_tree_t ;

/**
 * @struct wbfmm_shift_operators_t
 * @ingroup boxes
 *
 * Data type holding operators for upward and downward passes and
 * interaction calculations at each level
 *
 */

typedef struct {
  gsize size ; /**< size of data type, i.e. float or double */
  guint Lmax, /** < maximum order of expansions */
    nlevels, /** <  number of levels in tree */
    L[WBFMM_TREE_MAX_DEPTH+1], /** < maximum order of expansion per level */
    nerot ; /** <  number of elements in rotation operators */
  gpointer 
  SR[WBFMM_TREE_MAX_DEPTH+1],   
  /** < singular-to-regular coaxial translations */
    SS[WBFMM_TREE_MAX_DEPTH+1],
  /** < singular-to-singular (regular-to-regular) coaxial translations */
    rotations ; /** < rotation operations (H) */
} wbfmm_shift_operators_t ;

#define wbfmm_tree_point_number(_t) ((_t)->npoints)
#define wbfmm_tree_point_number_max(_t) ((_t)->maxpoints)
#define wbfmm_tree_depth(_t) ((_t)->depth)
#define wbfmm_tree_boxes_level(_t,_l) ((_t)->boxes[(_l)])
#define wbfmm_tree_width(_t)  ((_t)->D)
#define wbfmm_tree_origin(_t) ((gpointer)(&((_t)->x[0]))) 

#define wbfmm_element_number_coaxial(_N)			\
  ( ((_N)+1)*((_N)+2)*((_N)+3)/6 + ((_N)+1)*((_N)+2)/2 + (_N)+1)
#define wbfmm_element_number_rotation(_N)			\
  ((((_N)+1)*(4*(_N)+3)*((_N)+2)/6 + (((_N)+1)+((_N+1)))*((_N+1)+1) + (_N)))

/*index of parent in previous level*/
#define wbfmm_box_index_parent(_i) ((_i)/8)
/*index of first child in next level, add 0, ..., 7 to get full list*/
#define wbfmm_box_index_first_child(_i) ((_i)*8)

#define wbfmm_coefficient_index_nm(_n,_m) ((_n)*((_n)+1)+(_m))

gint wbfmm_cartesian_to_spherical(gdouble *x0, gdouble *x,
				  gdouble *r, gdouble *th, gdouble *ph) ;
gint wbfmm_shift_coordinates(gdouble *x, gdouble *y,
			     gdouble *ix, gdouble *iy, gdouble *iz,
			     gdouble *r) ;

gint wbfmm_legendre_recursion_array(gdouble **Pnm1, gdouble **Pn, gint n,
				     gdouble C, gdouble S) ;
gint wbfmm_bessel_j_recursion(gdouble *jnm1, gdouble *jn, 
			       gdouble x, gint n) ;

gint wbfmm_bessel_h_recursion(gdouble *hnm1, gdouble *hn, 
			       gdouble x, gint n) ;
gint wbfmm_bessel_j_init(gdouble x, gdouble *j0, gdouble *j1) ;
gint wbfmm_bessel_h_init(gdouble x, gdouble *h0, gdouble *h1) ;
gint wbfmm_legendre_init(gdouble C, gdouble S, 
			  gdouble *P0, gdouble *P10, gdouble *P11) ;


gint wbfmm_expansion_h_cfft(gdouble k, gint N, 
			    gdouble *x0,
			    gdouble *xs,
			    gdouble *q,
			    gdouble *cfft, gint cstr, gdouble *work) ;
gint wbfmm_expansion_dipole_h_cfft(gdouble k, gint N, 
				   gdouble *x0,
				   gdouble *xs,
				   gdouble *fx,
				   gdouble *fy,
				   gdouble *fz,
				   gdouble *cfft, gint cstr,
				   gdouble *work) ;

gint wbfmm_expansion_h_evaluate(gdouble k, gdouble *x0,
				gdouble *cfft, 
				gint cstr,
				gint N, 
				gdouble *xf, gdouble *field,
				gdouble *work) ;
gint wbfmm_expansion_j_evaluate(gdouble k, gdouble *x0,
				gdouble *cfft, 
				gint cstr,
				gint N, 
				gdouble *xf, gdouble *field,
				gdouble *work) ;

gint wbfmm_total_field(gdouble k,
		       gdouble *xs, gint xstride,
		       gdouble *src, gint sstride,
		       gint nsrc,
		       gdouble *xf, gdouble *field) ;
gint wbfmm_total_dipole_field(gdouble k,
			      gdouble *xs, gint xstride,
			      gdouble *src, gint sstride,
			      gint nsrc,
			      gdouble *xf, gdouble *field) ;
gint wbfmm_coordinate_transform(gdouble *x, 
				gdouble *ix, gdouble *iy, gdouble *iz,
				gdouble *y) ;

gint wbfmm_coefficients_RR_coaxial(gdouble *cfftRR, gint L,
				   gdouble kr, gdouble *work) ;
gint wbfmm_coefficients_SR_coaxial(gdouble *cfftSR, gint L,
				   gdouble kr, gdouble *work) ;
gint wbfmm_coaxial_translate(gdouble *Co, gint cstro, gint No,
			     gdouble *Ci, gint cstri, gint Ni,
			     gdouble *cfft, gint L,
			     gboolean complex) ;

gint wbfmm_rotation_angles(gdouble *ix, gdouble *iy, gdouble *iz, 
			   gdouble *jx, gdouble *jy, gdouble *jz, 
			   gdouble *th, gdouble *ph, gdouble *ch) ;
gint wbfmm_coefficients_H_rotation(gdouble *H, gint N, gdouble th,
				   gdouble *work) ;
gint wbfmm_rotate_H(gdouble *Co, gint cstro, gint N, gdouble *Ci, 
		    gint cstri, gdouble *H,
		    gdouble ph, gdouble ch) ;

gint wbfmm_child_parent_shift(gdouble *Cp, gint Np,
			      gdouble *Cc, gint Nc,
			      gdouble *H03, 
			      gdouble *H47, gint Lh,
			      gdouble *shift, gint Ls,
			      gdouble *work) ;
gint wbfmm_parent_child_shift(gdouble *Cc, gint Nc,
			      gdouble *Cp, gint Np,
			      gdouble *H03, 
			      gdouble *H47, gint Lh,
			      gdouble *shift, gint Ls,
			      gdouble *work) ;
gint wbfmm_shift_angles_list4(gint i, gint j, gint k,
			      gdouble *th, gdouble *ph,
			      gdouble *ch, gdouble *rs) ;
gint wbfmm_shift_angle_table_init(void) ;
wbfmm_shift_operators_t *wbfmm_shift_operators_new(guint L, gdouble *work) ;
gint wbfmm_shift_operators_coaxial_SR_init(wbfmm_shift_operators_t *w, 
					   gdouble D, guint level, guint L,
					   gdouble k, gdouble *work) ;
gint wbfmm_shift_operators_coaxial_SS_init(wbfmm_shift_operators_t *w, 
					   gdouble D, guint level, 
					   guint L, gdouble k, gdouble *work) ;
gint wbfmm_upward_pass(wbfmm_tree_t *t,
		       wbfmm_shift_operators_t *op,
		       guint level, gdouble *work) ;
gint wbfmm_downward_pass(wbfmm_tree_t *t,
			 wbfmm_shift_operators_t *op,
			 guint level, gdouble *work) ;
gint wbfmm_tree_box_field(wbfmm_tree_t *t, guint level,
			  guint b, gdouble k,
			  gdouble *x, gdouble *f, gdouble *work) ;
gint wbfmm_tree_box_local_field(wbfmm_tree_t *t, guint level,
				guint b, gdouble k,
				gdouble *x, gdouble *f, 
				gdouble *src, gint sstr,
				gboolean eval_neighbours,
				gdouble *work) ;
guint64 wbfmm_point_box(wbfmm_tree_t *t, guint level, gdouble *x) ;

gint wbfmm_tree_refine(wbfmm_tree_t *t) ;
gint wbfmm_tree_add_level(wbfmm_tree_t *tree) ;
gint wbfmm_tree_add_points(wbfmm_tree_t *t, 
			   gpointer pts, guint npts, gsize stride) ;

/*indexing and octrees*/
guint64 wbfmm_point_index_3d(gdouble *x, gdouble *c, gdouble D) ;

wbfmm_tree_t *wbfmm_tree_new(gdouble *x, gdouble D, guint maxpoints) ;
gint wbfmm_tree_coefficient_init(wbfmm_tree_t *t,
				 guint l, guint nr, guint ns) ;
gint wbfmm_tree_leaf_expansions(wbfmm_tree_t *t, gdouble k,
				gdouble *src, gint sstr,
				gdouble *work) ;

/*assorted utilities, handy for debugging*/
gint wbfmm_box_location_from_index(guint64 i, guint32 level,
				   gdouble *x0, gdouble D,
				   gdouble *x, gdouble *wb) ;
gint wbfmm_tree_box_centre(wbfmm_tree_t *t, guint level,
			   guint64 b, gdouble *xb,
			   gdouble *wb) ;
gint wbfmm_points_origin_width(gdouble *x, gint str, gint n,
			      gdouble *xmin, gdouble *D) ;
gint wbfmm_shift_angles(gdouble *xi, gdouble *xj,
			gdouble *th, gdouble *ph,
			gdouble *ch, gdouble *r) ;
gint wbfmm_tree_write_sources(wbfmm_tree_t *t, 
			      gdouble *q, gint stride,
			      FILE *f) ;
gint wbfmm_rotation_write_coefficients(gdouble *H, gint N, FILE *f) ;

/*single-precision functions*/
gint wbfmm_cartesian_to_spherical_f(gfloat *x0, gfloat *x,
				  gfloat *r, gfloat *th, gfloat *ph) ;
gint wbfmm_shift_coordinates_f(gfloat *x, gfloat *y,
			     gfloat *ix, gfloat *iy, gfloat *iz,
			     gfloat *r) ;

gint wbfmm_legendre_recursion_array_f(gfloat **Pnm1, gfloat **Pn, gint n,
				     gfloat C, gfloat S) ;
gint wbfmm_bessel_j_recursion_f(gfloat *jnm1, gfloat *jn, 
			       gfloat x, gint n) ;

gint wbfmm_bessel_h_recursion_f(gfloat *hnm1, gfloat *hn, 
			       gfloat x, gint n) ;
gint wbfmm_bessel_j_init_f(gfloat x, gfloat *j0, gfloat *j1) ;
gint wbfmm_bessel_h_init_f(gfloat x, gfloat *h0, gfloat *h1) ;
gint wbfmm_legendre_init_f(gfloat C, gfloat S, 
			  gfloat *P0, gfloat *P10, gfloat *P11) ;


gint wbfmm_expansion_h_cfft_f(gfloat k, gint N, 
			      gfloat *x0, gfloat *xs, gfloat *q,
			      gfloat *cfft, gint cstr, gfloat *work) ;
gint wbfmm_expansion_dipole_h_cfft_f(gfloat k, gint N, 
				     gfloat *x0,
				     gfloat *xs,
				     gfloat *fx,
				     gfloat *fy,
				     gfloat *fz,
				     gfloat *cfft, gint cstr,
				     gfloat *work) ;

gint wbfmm_expansion_h_evaluate_f(gfloat k, gfloat *x0,
				  gfloat *cfft, 
				  gint cstr,
				  gint N, 
				  gfloat *xf, gfloat *field,
				  gfloat *work) ;
gint wbfmm_expansion_j_evaluate_f(gfloat k, gfloat *x0,
				  gfloat *cfft, 
				  gint cstr,
				  gint N, 
				  gfloat *xf, gfloat *field,
				  gfloat *work) ;

gint wbfmm_total_field_f(gfloat k,
			 gfloat *xs, gint xstride,
			 gfloat *src, gint sstride,
			 gint nsrc,
			 gfloat *xf, gfloat *field) ;
gint wbfmm_total_dipole_field_f(gfloat k,
				gfloat *xs, gint xstride,
				gfloat *src, gint sstride,
				gint nsrc,
				gfloat *xf, gfloat *field) ;

gint wbfmm_coordinate_transform_f(gfloat *x, 
				  gfloat *ix, gfloat *iy, gfloat *iz,
				  gfloat *y) ;

gint wbfmm_coefficients_RR_coaxial_f(gfloat *cfftRR, gint L,
				     gfloat kr, gfloat *work) ;
gint wbfmm_coefficients_SR_coaxial_f(gfloat *cfftSR, gint L,
				     gfloat kr, gfloat *work) ;
gint wbfmm_coaxial_translate_f(gfloat *Co, gint cstro, gint No,
			       gfloat *Ci, gint cstri, gint Ni,
			       gfloat *cfft, gint L,
			       gboolean complex) ;

gint wbfmm_rotation_angles_f(gfloat *ix, gfloat *iy, gfloat *iz, 
			   gfloat *jx, gfloat *jy, gfloat *jz, 
			   gfloat *th, gfloat *ph, gfloat *ch) ;
gint wbfmm_coefficients_H_rotation_f(gfloat *H, gint N, gfloat th,
				     gfloat *work) ;
gint wbfmm_rotate_H_f(gfloat *Co, gint cstro, gint N, gfloat *Ci, gint cstri,
		      gfloat *H, gfloat ph, gfloat ch) ;

/*indexing and octrees*/
guint64 wbfmm_point_index_3d_f(gfloat *x, gfloat *c, gfloat D) ;

wbfmm_tree_t *wbfmm_tree_new_f(gfloat *x, gfloat D, guint maxpoints) ;
gint wbfmm_tree_coefficient_init_f(wbfmm_tree_t *t,
				   guint l, guint nr, guint ns) ;
gint wbfmm_tree_leaf_expansions_f(wbfmm_tree_t *t, gfloat k,
				  gfloat *src, gint sstr,
				  gfloat *work) ;
gint wbfmm_tree_refine_f(wbfmm_tree_t *t) ;
gint wbfmm_tree_add_points_f(wbfmm_tree_t *t, 
			     gpointer pts, guint npts, gsize stride) ;

/*assorted utilities, handy for debugging*/
gint wbfmm_box_location_from_index_f(guint64 i, guint32 level,
				     gfloat *x0, gfloat D,
				     gfloat *x, gfloat *wb) ;
gint wbfmm_tree_box_centre_f(wbfmm_tree_t *t, guint level,
			     guint64 b, gfloat *xb,
			     gfloat *wb) ;
gint wbfmm_child_parent_shift_f(gfloat *Cp, gint Np,
				gfloat *Cc, gint Nc,
				gfloat *H03, 
				gfloat *H47, gint Lh,
				gfloat *shift, gint Ls,
				gfloat *work) ;
gint wbfmm_parent_child_shift_f(gfloat *Cc, gint Nc,
				gfloat *Cp, gint Np,
				gfloat *H03, 
				gfloat *H47, gint Lh,
				gfloat *shift, gint Ls,
				gfloat *work) ;
gint wbfmm_points_origin_width_f(gfloat *x, gint str, gint n,
				gfloat *xmin, gfloat *D) ;

gint wbfmm_shift_angles_list4_f(gint i, gint j, gint k,
				gfloat *th, gfloat *ph,
				gfloat *ch, gfloat *rs) ;
gint wbfmm_shift_angles_f(gfloat *xi, gfloat *xj,
			  gfloat *th, gfloat *ph,
			  gfloat *ch, gfloat *r) ;

gint wbfmm_tree_write_sources_f(wbfmm_tree_t *t, 
				gfloat *q, gint stride,
				FILE *f) ;

gint wbfmm_rotation_write_coefficients_f(gfloat *H, gint N, FILE *f) ;

gint wbfmm_shift_angle_table_init_f(void) ;
wbfmm_shift_operators_t *wbfmm_shift_operators_new_f(guint L,
						     gfloat *work) ;
gint wbfmm_shift_operators_coaxial_SR_init_f(wbfmm_shift_operators_t *w, 
					     gfloat D, guint level, guint L,
					     gfloat k, gfloat *work) ;
gint wbfmm_shift_operators_coaxial_SS_init_f(wbfmm_shift_operators_t *w, 
					     gfloat D, guint level, 
					     guint L, gfloat k, gfloat *work) ;
gint wbfmm_upward_pass_f(wbfmm_tree_t *t,
			 wbfmm_shift_operators_t *op,
			 guint level, gfloat *work) ;
gint wbfmm_downward_pass_f(wbfmm_tree_t *t,
			   wbfmm_shift_operators_t *op,
			   guint level, gfloat *work) ;
gint wbfmm_tree_box_field_f(wbfmm_tree_t *t, guint level,
			    guint b, gfloat k,
			    gfloat *x, gfloat *f, gfloat *work) ;
gint wbfmm_tree_box_local_field_f(wbfmm_tree_t *t, guint level,
				  guint b, gfloat k,
				  gfloat *x, gfloat *f, 
				  gfloat *src, gint sstr,
				  gboolean eval_neighbours,
				  gfloat *work) ;
guint64 wbfmm_point_box_f(wbfmm_tree_t *t, guint level, gfloat *x) ;

/*precision independent functions*/
guint64 wbfmm_point_locate_box(guint64 x, guint level) ;
gint wbfmm_point_from_index(guint64 i, guint32 *x, guint32 *y, guint32 *z) ;
guint64 wbfmm_box_index(guint32 i, guint32 j, guint32 k) ;
gint wbfmm_box_location(guint64 idx, guint32 *i, guint32 *j, guint32 *k) ;
guint64 wbfmm_box_parent(guint64 idx) ;
guint64 wbfmm_box_first_child(guint64 idx) ;
gint wbfmm_tree_print(FILE *f, wbfmm_tree_t *t, guint level,
		      gboolean print_empty) ;
gint wbfmm_logging_init(FILE *f, gchar *p, 
			GLogLevelFlags log_level,
			gpointer exit_func, gboolean timed) ;
gint wbfmm_box_neighbours(guint level, guint64 idx, guint64 *neighbours) ;
gint wbfmm_box_interaction_list_4(guint level, guint64 idx, 
				  guint64 *list, gboolean sort) ;
gint wbfmm_box_interaction_index(gint i, gint j, gint k) ;

#endif /*_WBFMM_H_INCLUDED_*/
