/* This file is part of WBFMM, a Wide-Band Fast Multipole Method code
 *
 * Copyright (C) 2019, 2021 Michael Carley
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


#ifndef WBFMM_H_INCLUDED
#define WBFMM_H_INCLUDED

#include <stdio.h>

#define WBFMM_INDEX_SCALE (1LU << 63)
#define WBFMM_INDEX_SHIFT (1U << 20)

#define WBFMM_LOCAL_CUTOFF_RADIUS 1e-6
#ifndef WBFMM_THREAD_NUMBER_MAX
#define WBFMM_THREAD_NUMBER_MAX 16
#endif /*WBFMM_THREAD_NUMBER_MAX*/

/**
 * @struct wbfmm_library_config_t
 * @ingroup util
 *
 * Data type to report library compilation settings and other
 * configuration information
 *
 */

typedef struct {
  gboolean
  avx,
    avx2,
    fma,
    openmp ;
  char *switches ;
  gsize
  real_size ;
} wbfmm_library_config_t ;

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
 * Selection of physical problem to be handled by a ::wbfmm_tree_t
 *
 * @ingroup boxes
 *
 */

typedef enum
  {
   WBFMM_PROBLEM_LAPLACE = 1, /**< Laplace equation */
   WBFMM_PROBLEM_HELMHOLTZ = 2 /**< Helmholtz equation */	      
  } wbfmm_problem_t ; 

/** 
 * Selection of field to be calculated
 *
 * @ingroup boxes
 *
 */

typedef enum
  {
   WBFMM_FIELD_POTENTIAL     = 1 << 0, /**< scalar field */
   WBFMM_FIELD_GRADIENT      = 1 << 1, /**< gradient of scalar field */
   WBFMM_FIELD_CURL          = 1 << 2,  /**< curl of vector field */
   WBFMM_FIELD_CURL_GRADIENT = 1 << 3  /**< curl of vector field and its 
					  gradient */
  } wbfmm_field_t ;

/**
 * @struct wbfmm_tree_t
 * @ingroup boxes
 *
 * Data type for octrees
 *
 */

typedef struct {
  wbfmm_problem_t problem ;
  wbfmm_box_t *boxes[WBFMM_TREE_MAX_DEPTH+1] ; 
  /**< arrays of boxes at each level */
  gboolean sorted ;
  guint 
  maxpoints, /**< maximum number of points in tree */
    npoints, /**< number of points in tree */
    *ip, /**< indices of points, sorted by Morton index */
    nq, /**< number of source components*/
    depth, /**< depth of tree */
    order_s[WBFMM_TREE_MAX_DEPTH+1], 
  /**< order of singular expansions at each level */
    order_r[WBFMM_TREE_MAX_DEPTH+1] ; 
  /**< order of regular expansions at each level */
  char 
  x[24], /**< origin of tree domain cube */
    *normals, /**< normal coordinates */
    *points ; /**< point coordinates */
  gpointer 
  mps[WBFMM_TREE_MAX_DEPTH+1], 
  /**< singular expansion data at each level */
    mpr[WBFMM_TREE_MAX_DEPTH+1] ; 
  /**< regular expansion data at each level */
  gsize
  size,  /**< size of floating point type in data (float, double, etc) */
    nstr, /**< stride in normal data */
    pstr ; /**< stride in point data */
  gdouble D ; /**< width of domain cube */
} wbfmm_tree_t ;

typedef enum
  {
   WBFMM_SOURCE_MONOPOLE = 1 << 0,
   WBFMM_SOURCE_DIPOLE   = 1 << 1
  } wbfmm_source_t ;

/**
 * @struct wbfmm_target_list_t
 * @ingroup boxes
 *
 * Data type for target point lists
 *
 */

typedef struct {
  wbfmm_tree_t *t ; /**< tree containing source data */
  guint
  field,     /**< field specifier (see ::wbfmm_tree_box_local_field) */
    source,  /**< source specifier */
    maxpoints, /**< maximum number of points in target list */
    npoints, /**< number of points in target list */
    *ip,     /**< indices of points, sorted by Morton index */
    nc ;     /**< number of coefficients (size of blocks of coefficients) */
  guint32
  *boxes ;   /**< box indices of points */
  char 
  *points ; /**< field point coordinates */
  gsize
  size,      /**< size of floating point type in data (float, double, etc) */
    pstr ;   /**< stride in point data */
  gint
  *ibox,    /**< start and end of source index lists for each box */
    *isrc,  /**< source index lists for each box */
    *ics ;  /**< start of near-field coefficients for each target */
  gpointer
  cfft,     /**< coefficients of regular expansions in boxes */
    csrc ;  /**< coefficients of near-field (direct) interactions, 
	       point-by-point */
  gboolean
  complex ; /**< complex-valued field or real */
} wbfmm_target_list_t ;

#define wbfmm_target_list_field(_l) ((_l)->field)
#define wbfmm_target_point_box(_l,_i) ((_l)->boxes[(_i)])

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
  gboolean
  bw ; /** < operators allocated for backward translation algorithm */
} wbfmm_shift_operators_t ;

#define wbfmm_tree_point_number(_t) ((_t)->npoints)
#define wbfmm_tree_point_number_max(_t) ((_t)->maxpoints)
#define wbfmm_tree_depth(_t) ((_t)->depth)
#define wbfmm_tree_boxes_level(_t,_l) ((_t)->boxes[(_l)])
#define wbfmm_tree_width(_t)  ((_t)->D)
#define wbfmm_tree_origin(_t) ((gpointer)(&((_t)->x[0]))) 
#define wbfmm_tree_problem(_t) ((_t)->problem)
#define wbfmm_tree_source_size(_t) ((_t)->nq)

#define wbfmm_target_list_tree(_t) ((_t)->t)
#define wbfmm_target_list_point_number(_t) ((_t)->npoints)
#define wbfmm_target_list_point_number_max(_t) ((_t)->maxpoints)
#define wbfmm_target_list_gradient(_t) ((_t)->grad)

#define wbfmm_element_number_coaxial(_N)			\
  ( ((_N)+1)*((_N)+2)*((_N)+3)/6 + ((_N)+1)*((_N)+2)/2 + (_N)+1)
#define wbfmm_element_number_rotation(_N)			\
  ((((_N)+1)*(4*(_N)+3)*((_N)+2)/6 + (((_N)+1)+((_N+1)))*((_N+1)+1) + (_N)))
#define wbfmm_T_rotation_matrix_size(_n) (((_n)+1)*(2*(_n)+1)*(2*(_n)+3)/3)

/*index of parent in previous level*/
#define wbfmm_box_index_parent(_i) ((_i)/8)
/*index of first child in next level, add 0, ..., 7 to get full list*/
#define wbfmm_box_index_first_child(_i) ((_i)*8)

#define wbfmm_coefficient_index_nm(_n,_m) ((_n)*((_n)+1)+(_m))
#define wbfmm_coefficient_number(_n) \
  wbfmm_coefficient_index_nm(((_n)+1), (-((_n)+1)))
#define wbfmm_conjugate_index_nm(_n,_m) (((_n)*((_n)+1)/2)+(_m))

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
			    gdouble *q, gint nq,
			    gdouble *cfft, gint cstr, gdouble *work) ;
gint wbfmm_expansion_dipole_h_cfft(gdouble k, gint N, 
				   gdouble *x0,
				   gdouble *xs,
				   gdouble *fx,
				   gdouble *fy,
				   gdouble *fz,
				   gint nq,				   
				   gdouble *cfft, gint cstr,
				   gdouble *work) ;
gint wbfmm_expansion_normal_h_cfft(gdouble k, gint N, 
				   gdouble *x0,
				   gdouble *xs,
				   gdouble *normal,
				   gdouble *q, gint nq,
				   gdouble *cfft, gint cstr,
				   gdouble *work) ;

gint wbfmm_expansion_h_evaluate(gdouble k, gdouble *x0,
				gdouble *cfft, 
				gint cstr,
				gint N, gint nq,
				gdouble *xf, gdouble *field,
				gint fstr,
				gdouble *work) ;
gint wbfmm_expansion_h_grad_evaluate(gdouble k, gdouble *x0,
				     gdouble *cfft, 
				     gint cstr,
				     gint N, gint nq,
				     gdouble *xf,
				     gdouble *field, gint fstr,
				     gdouble *work) ;
gint wbfmm_expansion_j_evaluate(gdouble k, gdouble *x0,
				gdouble *cfft, 
				gint cstr,
				gint N,  gint nq,
				gdouble *xf, gdouble *field,
				gint fstr, 
				gdouble *work) ;
gint wbfmm_expansion_j_grad_evaluate(gdouble k, gdouble *x0,
				     gdouble *cfft, 
				     gint cstr,
				     gint N, gint nq,
				     gdouble *xf,
				     gdouble *field, gint fstr,
				     gdouble *work) ;

gint wbfmm_total_field(gdouble k,
		       gdouble *xs, gint xstride,
		       gdouble *src, gint sstride,
		       gdouble *normals, gint nstr,
		       gdouble *dipoles, gint dstr,
		       gint nq, gint nsrc,
		       gdouble *xf, gdouble *field, gint fstr) ;
gint wbfmm_total_field_grad(gdouble k,
			    gdouble *xs, gint xstride,
			    gdouble *src, gint sstride,
			    gdouble *normals, gint nstr,
			    gdouble *dipoles, gint dstr,
			    gint nq, gint nsrc,
			    gdouble *xf, gdouble *field, gint fstr) ;
gint wbfmm_total_dipole_field(gdouble k,
			      gdouble *xs, gint xstride,
			      gdouble *src, gint sstride,
			      gint nsrc,
			      gdouble *xf, gdouble *field) ;
gint wbfmm_total_normal_field(gdouble k,
			      gdouble *xs, gint xstride,
			      gdouble *ns, gint nstride,
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
gint wbfmm_coaxial_translate_ref(gdouble *Co, gint cstro, gint No,
				 gdouble *Ci, gint cstri, gint Ni,
				 gint nq,
				 gdouble *cfft, gint L,
				 gboolean complex, gdouble sc) ;

gint wbfmm_rotation_angles(gdouble *ix, gdouble *iy, gdouble *iz, 
			   gdouble *jx, gdouble *jy, gdouble *jz, 
			   gdouble *th, gdouble *ph, gdouble *ch) ;
gint wbfmm_coefficients_H_rotation(gdouble *H, gint N, gdouble th,
				   gdouble *work) ;
gint wbfmm_coefficients_H_to_T(gdouble *H, gint N,
			       gdouble th, 
			       gdouble ph, 
			       gdouble ch,
			       gdouble *T) ;
gint wbfmm_coefficients_H_to_T_f(gfloat *H, gint N,
				 gfloat th, 
				 gfloat ph, 
				 gfloat ch,
				 gfloat *T) ;
  
gint wbfmm_rotate_H_ref(gdouble *Co, gint cstro, 
			gdouble *Ci, gint cstri,
			gint N,
			gint nq,
			gdouble *H,
			gdouble ph, gdouble ch, gdouble sc) ;
gint wbfmm_rotate_H_avx(gdouble *Co, gint cstro, 
			gdouble *Ci, gint cstri,
			gint N,
			gint nq,
			gdouble *H,
			gdouble ph, gdouble ch, gdouble sc) ;
gint wbfmm_rotate_T(gdouble *Co, gint cstro, gdouble *Ci, gint cstri,
		    gint N, gint nq, gdouble *T, gdouble *sc) ;
gint wbfmm_rotate_T_f(gfloat *Co, gint cstro, gfloat *Ci, gint cstri,
		      gint N, gint nq, gfloat *T, gfloat *sc) ;

gint wbfmm_local_coefficients(gdouble k, gdouble *x, gint N,
			      guint field, gdouble *cfft, gdouble *work) ;
gint wbfmm_local_coefficients_f(gfloat k, gfloat *x, gint N,
				guint field, gfloat *cfft, gfloat *work) ;

gint wbfmm_expansion_apply(gdouble *C, gint cstr, gint nq,
			   gdouble *ec, gint N, guint field,
			   gdouble *f, gint fstr) ;
gint wbfmm_expansion_apply_f(gfloat *C, gint cstr, gint nq,
			     gfloat *ec, gint N, guint field,
			     gfloat *f, gint fstr) ;

gint wbfmm_laplace_expansion_cfft(gint N, gdouble *x0,
				  gdouble *xs, gdouble *q, gint nq,
				  gdouble *cfft, gint cstr,
				  gdouble *work) ;
gint wbfmm_laplace_expansion_normal_cfft(gint N, gdouble *x0,
					 gdouble *xs,
					 gdouble *normal,
					 gdouble *q,
					 gint nq,
					 gdouble *cfft,
					 gint cstr,
					 gdouble *work) ;
gint wbfmm_laplace_expansion_dipole_cfft(gint N,
					 gdouble *x0,
					 gdouble *xs,
					 gdouble *fx,
					 gdouble *fy,
					 gdouble *fz,
					 gint nq,
					 gdouble *cfft,
					 gint cstr,
					 gdouble *work) ;
gint wbfmm_laplace_field_direct(gdouble *xs, gint xstride,
				gdouble *n, gint nstr,
				gint nsrc,
				gdouble *src, gint sstr,
				gdouble *d, gint dstr,
				gint nq,
				guint field,
				gdouble *xf, gdouble *f, gint fstr) ;
gint wbfmm_laplace_field_direct_f(gfloat *xs, gint xstride,
				  gfloat *n, gint nstr,
				  gint nsrc,
				  gfloat *src, gint sstr,
				  gfloat *d, gint dstr,
				  gint nq,
				  guint field,
				  gfloat *xf, gfloat *f, gint fstr) ;


gint wbfmm_laplace_field(gdouble *xs, gint xstride,
			 gdouble *src, gint sstride,
			 gint nq,
			 gdouble *normals, gint nstr,
			 gdouble *dipoles, gint dstr,
			 gint nsrc,
			 gdouble *xf, gdouble *field,
			 gint fstr) ;

gint wbfmm_laplace_field_grad(gdouble *xs, gint xstride,
			      gdouble *src, gint sstride,
			      gint nq,
			      gdouble *normals, gint nstr,
			      gdouble *dipoles, gint dstr,
			      gint nsrc,
			      gdouble *xf, gdouble *field, gint fstr) ;
gint wbfmm_laplace_field_laplacian(gdouble *xs,
				   gint xstride,
				   gdouble *src,
				   gint sstride,
				   gint nq,
				   gdouble *normals,
				   gint nstr,
				   gdouble *dipoles,
				   gint dstr,
				   gint nsrc,
				   gdouble *xf,
				   gdouble *field,
				   gint fstr) ;

gint wbfmm_laplace_field_curl(gdouble *xs, gint xstride,
			      gdouble *src, gint sstride,
			      gint nq,
			      gdouble *normals, gint nstr,
			      gdouble *dipoles, gint dstr,
			      gint nsrc,
			      gdouble *xf, gdouble *field, gint fstr) ;
gint wbfmm_laplace_expansion_grad_evaluate(gdouble *x0, gdouble *cfft,
					   gint cstr, gint N, gint nq,
					   gdouble *xf, gdouble *field,
					   gint fstr,
					   gdouble *work) ;
gint wbfmm_laplace_expansion_curl_evaluate(gdouble *x0, gdouble *cfft,
					   gint cstr, gint N, gint nq,
					   gdouble *xf, gdouble *field,
					   gint fstr, gdouble *work) ;
gint wbfmm_laplace_expansion_curl_evaluate_f(gfloat *x0, gfloat *cfft,
					     gint cstr, gint N, gint nq,
					     gfloat *xf, gfloat *field,
					     gint fstr, gfloat *work) ;

gint wbfmm_laplace_field_grad_f(gfloat *xs, gint xstride,
				gfloat *src, gint sstride,
				gint nq,
				gfloat *normals, gint nstr,
				gfloat *dipoles, gint dstr,
				gint nsrc,
				gfloat *xf, gfloat *field, gint fstr) ;
gint wbfmm_laplace_field_laplacian_f(gfloat *xs,
				     gint xstride,
				     gfloat *src,
				     gint sstride,
				     gint nq,
				     gfloat *normals,
				     gint nstr,
				     gfloat *dipoles,
				     gint dstr,
				     gint nsrc,
				     gfloat *xf,
				     gfloat *field,
				     gint fstr) ;
gint wbfmm_laplace_field_curl_f(gfloat *xs, gint xstride,
				gfloat *src, gint sstride,
				gint nq,
				gfloat *normals, gint nstr,
				gfloat *dipoles, gint dstr,
				gint nsrc,
				gfloat *xf, gfloat *field, gint fstr) ;
gint wbfmm_laplace_expansion_grad_evaluate_f(gfloat *x0, gfloat *cfft,
					     gint cstr, gint N, gint nq,
					     gfloat *xf, gfloat *field,
					     gint fstr,
					     gfloat *work) ;
gint wbfmm_laplace_expansion_evaluate(gdouble *x0, gdouble *cfft,
				      gint cstr, gint N, gint nq,
				      gdouble *xf, gdouble *field,
				      gdouble *work) ;
gint wbfmm_laplace_expansion_evaluate_f(gfloat *x0, gfloat *cfft,
					gint cstr, gint N, gint nq,
					gfloat *xf, gfloat *field,
					gfloat *work) ;
gint wbfmm_laplace_expansion_local_eval(gdouble *x0,
					gdouble *cfft, gint cstr, 
					gint N,	gint nq,
					guint field,
					gdouble *xf,
					gdouble *f,gint fstr,
					gdouble *work) ;
gint wbfmm_laplace_expansion_local_eval_f(gfloat *x0,
					  gfloat *cfft, gint cstr, 
					  gint N, gint nq,
					  guint field,
					  gfloat *xf,
					  gfloat *f,gint fstr,
					  gfloat *work) ;

gint wbfmm_laplace_expansion_local_laplacian_evaluate(gdouble *x0,
						      gdouble *cfft,
						      gint cstr, 
						      gint N,
						      gint nq,
						      gdouble *xf,
						      gdouble *field,
						      gint fstr,
						      gdouble *work) ;
gint wbfmm_laplace_expansion_local_laplacian_evaluate_f(gfloat *x0,
							gfloat *cfft,
							gint cstr, 
							gint N,
							gint nq,
							gfloat *xf,
							gfloat *field,
							gint fstr,
							gfloat *work) ;
gint wbfmm_laplace_expansion_local_curl_evaluate(gdouble *x0, gdouble *cfft,
						 gint cstr, gint N,
						 gint nq, gdouble *xf,
						 gdouble *field,
						 gint fstr, gdouble *work) ;
gint wbfmm_laplace_expansion_local_curl_evaluate_f(gfloat *x0, gfloat *cfft,
						   gint cstr, gint N,
						   gint nq, gfloat *xf,
						   gfloat *field,
						   gint fstr, gfloat *work) ;

gint wbfmm_laplace_coaxial_translate_init(gint N) ;
gint wbfmm_laplace_coaxial_translate_init_f(gint N) ;

gint wbfmm_laplace_expansion_cfft_f(gint N, gfloat *x0,
				    gfloat *xs, gfloat *q, gint nq,
				    gfloat *cfft, gint cstr,
				    gfloat *work) ;
gint wbfmm_laplace_expansion_normal_cfft_f(gint N, gfloat *x0,
					   gfloat *xs,
					   gfloat *normal,
					   gfloat *q,
					   gint nq,
					   gfloat *cfft,
					   gint cstr,
					   gfloat *work) ;
gint wbfmm_laplace_expansion_dipole_cfft_f(gint N,
					   gfloat *x0,
					   gfloat *xs,
					   gfloat *fx,
					   gfloat *fy,
					   gfloat *fz,
					   gint nq,
					   gfloat *cfft,
					   gint cstr,
					   gfloat *work) ;
gint wbfmm_laplace_field_f(gfloat *xs, gint xstride,
			   gfloat *src, gint sstride,
			   gint nq,
			   gfloat *normals, gint nstr,
			   gfloat *dipoles, gint dstr,
			   gint nsrc,
			   gfloat *xf, gfloat *field,
			   gint fstr) ;
			   
gint wbfmm_laplace_coaxial_translate_SS(gdouble *Co, gint cstro, gint No,
					gdouble *Ci, gint cstri, gint Ni,
					gint nq, gdouble t, gdouble sc) ;
gint wbfmm_laplace_coaxial_translate_SS_f(gfloat *Co, gint cstro, gint No,
					  gfloat *Ci, gint cstri, gint Ni,
					  gint nq, gfloat t, gfloat sc) ;
gint wbfmm_laplace_coaxial_translate_RR(gdouble *Co, gint cstro, gint No,
					gdouble *Ci, gint cstri, gint Ni,
					gint nq, gdouble t, gdouble sc) ;
gint wbfmm_laplace_coaxial_translate_RR_f(gfloat *Co, gint cstro, gint No,
					  gfloat *Ci, gint cstri, gint Ni,
					  gint nq, gfloat t, gfloat sc) ;
gint wbfmm_laplace_coaxial_translate_SR(gdouble *Co, gint cstro, gint No,
					gdouble *Ci, gint cstri, gint Ni,
					gint nq, gdouble t, gdouble sc) ;
gint wbfmm_laplace_coaxial_translate_SR_f(gfloat *Co, gint cstro, gint No,
					  gfloat *Ci, gint cstri, gint Ni,
					  gint nq, gfloat t, gfloat sc) ;
gint wbfmm_laplace_rotate_H_ref(gdouble *Co, gint cstro,
				gdouble *Ci, gint cstri,
				gint N, gint nq,
				gdouble *H,
				gdouble ph, gdouble ch, gdouble sc) ;
gint wbfmm_laplace_rotate_H_ref_f(gfloat *Co, gint cstro,
				  gfloat *Ci, gint cstri,
				  gint N, gint nq,
				  gfloat *H,
				  gfloat ph, gfloat ch, gfloat sc) ;
gint wbfmm_laplace_rotate_H_avx(gdouble *Co, gint cstro,
				gdouble *Ci, gint cstri,
				gint N, gint nq,
				gdouble *H,
				gdouble ph, gdouble ch, gdouble sc) ;

gint wbfmm_laplace_child_parent_shift(gdouble *Cp, gint Np,
				      gdouble *Cc, gint Nc,
				      gint nq,
				      gdouble *H03, gdouble *H47,
				      gint Lh,
				      gdouble t,
				      gdouble *work) ;
gint wbfmm_laplace_child_parent_shift_f(gfloat *Cp, gint Np,
					gfloat *Cc, gint Nc,
					gint nq,
					gfloat *H03, gfloat *H47,
					gint Lh,
					gfloat t,
					gfloat *work) ;
gint wbfmm_laplace_parent_child_shift(gdouble *Cc, gint Nc,
				      gdouble *Cp, gint Np,
				      gint nq,
				      gdouble *H03, gdouble *H47,
				      gint Lh,
				      gdouble t,
				      gdouble *work) ;
gint wbfmm_laplace_parent_child_shift_f(gfloat *Cc, gint Nc,
					gfloat *Cp, gint Np,
					gint nq,
					gfloat *H03, gfloat *H47,
					gint Lh,
					gfloat t,
					gfloat *work) ;
gint wbfmm_laplace_child_parent_shift_bw(gdouble *Cp, gint Np,
					 gdouble *Cc, gint Nc,
					 gint nq,
					 gdouble *H03, gint Lh,
					 gdouble wb,
					 gdouble *work) ;
gint wbfmm_laplace_child_parent_shift_bw_f(gfloat *Cp, gint Np,
					   gfloat *Cc, gint Nc,
					   gint nq,
					   gfloat *H03, gint Lh,
					   gfloat wb,
					   gfloat *work) ;

  gint wbfmm_tree_laplace_coefficient_init(wbfmm_tree_t *t,
					 guint l, 
					 guint nr,
					 guint ns) ;
gint wbfmm_tree_laplace_coefficient_init_f(wbfmm_tree_t *t,
					   guint l, 
					   guint nr,
					   guint ns) ;
gint wbfmm_tree_laplace_leaf_expansions(wbfmm_tree_t *t,
					gdouble *src,
					gint sstr,
					gdouble *dipoles,
					gint dstr,
					gboolean zero_expansions,
					gdouble *work) ;
gint wbfmm_tree_laplace_leaf_expansions_f(wbfmm_tree_t *t,
					  gfloat *src,
					  gint sstr,
					  gfloat *dipoles,
					  gint dstr,
					  gboolean zero_expansions,
					  gfloat *work) ;
gint wbfmm_laplace_downward_pass_ref(wbfmm_tree_t *t,
				     wbfmm_shift_operators_t *op,
				     guint level,
				     gdouble *work, gint nthreads) ;
gint wbfmm_laplace_downward_pass_avx(wbfmm_tree_t *t,
				     wbfmm_shift_operators_t *op,
				     guint level,
				     gdouble *work, gint nthreads) ;
gint wbfmm_laplace_downward_pass_ref_f(wbfmm_tree_t *t,
				       wbfmm_shift_operators_t *op,
				       guint level,
				       gfloat *work, gint nthreads) ;
gint wbfmm_laplace_upward_pass(wbfmm_tree_t *t,
				 wbfmm_shift_operators_t *op,
				 guint level,
				 gdouble *work) ;
gint wbfmm_laplace_upward_pass_f(wbfmm_tree_t *t,
				   wbfmm_shift_operators_t *op,
				   guint level,
				   gfloat *work) ;

gint wbfmm_laplace_box_field(wbfmm_tree_t *t, guint level, guint b,
			     gdouble *src, gint sstr,
			     gdouble *d, gint dstr,
			     guint field,
			     gboolean eval_neighbours,
			     gdouble *x, gdouble *f, gint fstr, gdouble *work) ;
gint wbfmm_laplace_box_field_f(wbfmm_tree_t *t, guint level, guint b,
			       gfloat *src, gint sstr,
			       gfloat *d, gint dstr,
			       guint field,
			       gboolean eval_neighbours,
			       gfloat *x, gfloat *f, gint fstr, gfloat *work) ;

gint wbfmm_laplace_local_coefficients(gdouble *x, gint N,
				      guint field, gdouble *cfft,
				      gdouble *work) ;
gint wbfmm_laplace_local_coefficients_f(gfloat *x, gint N,
					guint field, gfloat *cfft,
					gfloat *work) ;
gint wbfmm_laplace_field_coefficients(gdouble *x, gint N,
				      guint field, gdouble *cfft,
				      gdouble *work) ;
gint wbfmm_laplace_field_coefficients_f(gfloat *x, gint N,
					guint field, gfloat *cfft,
					gfloat *work) ;
gint wbfmm_laplace_expansion_apply(gdouble *C, gint cstr, gint nq,
				   gdouble *ec, gint N, guint field,
				   gdouble *f, gint fstr) ;
gint wbfmm_laplace_expansion_apply_f(gfloat *C, gint cstr, gint nq,
				     gfloat *ec, gint N, guint field,
				     gfloat *f, gint fstr) ;


gint wbfmm_child_parent_shift_bw(gdouble *Cp, gint Np, gdouble *Cc, gint Nc,
				 gdouble *H03, gint Lh,
				 gdouble *transf, gdouble *transb, gint Ls,
				 gdouble *work) ;
gint wbfmm_child_parent_shift_bw_f(gfloat *Cp, gint Np, gfloat *Cc, gint Nc,
				   gfloat *H03, gint Lh,
				   gfloat *transf, gfloat *transb, gint Ls,
				   gfloat *work) ;

gint wbfmm_child_parent_shift(gdouble *Cp, gint Np,
			      gdouble *Cc, gint Nc,
			      gdouble *H03, 
			      gdouble *H47, gint Lh,
			      gdouble *shift, gint Ls,
			      gint nq,
			      gdouble *work) ;
gint wbfmm_parent_child_shift(gdouble *Cc, gint Nc,
			      gdouble *Cp, gint Np,
			      gdouble *H03, 
			      gdouble *H47, gint Lh,
			      gdouble *shift, gint Ls,
			      gint nq,
			      gdouble *work) ;
gint wbfmm_shift_angles_list4(gint i, gint j, gint k,
			      gdouble *th, gdouble *ph,
			      gdouble *ch, gdouble *rs) ;
gint wbfmm_shift_angle_table_init(void) ;
wbfmm_shift_operators_t *wbfmm_shift_operators_new(guint L,
						   gboolean bw,
						   gdouble *work) ;
gint wbfmm_shift_operators_coaxial_SR_init(wbfmm_shift_operators_t *w, 
					   gdouble D, guint level, guint L,
					   gdouble k, gdouble *work) ;
gint wbfmm_shift_operators_coaxial_SS_init(wbfmm_shift_operators_t *w, 
					   gdouble D, guint level, 
					   guint L, gdouble k, gdouble *work) ;
gint wbfmm_upward_pass(wbfmm_tree_t *t,
		       wbfmm_shift_operators_t *op,
		       guint level, gdouble *work) ;
gint wbfmm_downward_pass_ref(wbfmm_tree_t *t,
			     wbfmm_shift_operators_t *op,
			     guint level, gdouble *work, gint nthreads) ;
gint wbfmm_downward_pass_avx(wbfmm_tree_t *t,
			     wbfmm_shift_operators_t *op,
			     guint level, gdouble *work, gint nthreads) ;
gint wbfmm_tree_box_field(wbfmm_tree_t *t, guint level,
			  guint b, gdouble k,
			  gdouble *x, gdouble *f, gint fstr, gdouble *work) ;
gint wbfmm_tree_box_local_field(wbfmm_tree_t *t, guint level,
				guint b, gdouble k,
				gdouble *x, gdouble *f, gint fstr,
				gdouble *src, gint sstr,
				gdouble *d, gint dstr,
				gboolean eval_neighbours,
				guint field,
				gdouble *work) ;
guint64 wbfmm_point_box(wbfmm_tree_t *t, guint level, gdouble *x) ;

gint wbfmm_tree_refine(wbfmm_tree_t *t) ;
gint wbfmm_tree_add_level(wbfmm_tree_t *tree) ;
gint wbfmm_tree_add_points(wbfmm_tree_t *t, 
			   gpointer pts, gsize pstr,
			   gpointer normals, gsize nstr,
			   guint npts, gboolean sorted) ;
gint wbfmm_tree_sort_points(wbfmm_tree_t *t, 
			    gpointer pts, gsize psize,
			    guint npts) ;
gint wbfmm_target_list_coefficients_init(wbfmm_target_list_t *l,
					 guint field) ;

/*indexing and octrees*/
guint64 wbfmm_point_index_3d(gdouble *x, gdouble *c, gdouble D) ;

wbfmm_tree_t *wbfmm_tree_new(gdouble *x, gdouble D, guint maxpoints) ;
gint wbfmm_tree_coefficient_init(wbfmm_tree_t *t,
				 guint l, guint nr, guint ns) ;
gint wbfmm_tree_coefficient_clear(wbfmm_tree_t *t, guint l) ;
gint wbfmm_tree_leaf_expansions(wbfmm_tree_t *t, gdouble k,
				gdouble *src, gint sstr,
				gdouble *dipoles, gint dstr,
				gboolean zero_expansions,
				gdouble *work) ;

gint wbfmm_truncation_number(wbfmm_tree_t *t, gdouble k, guint level,
			     gdouble tol) ;

/*assorted utilities, handy for debugging*/
gint wbfmm_box_location_from_index(guint64 i, guint32 level,
				   gdouble *x0, gdouble D,
				   gdouble *x, gdouble *wb) ;
gint wbfmm_tree_box_centre(wbfmm_tree_t *t, guint level,
			   guint64 b, gdouble *xb,
			   gdouble *wb) ;
gint wbfmm_points_origin_width(gdouble *x, gint str, gint n,
			       gdouble *xmin,
			       gdouble *xmax,
			       gdouble *D,
			       gboolean init_limits) ;
			       
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
			      gfloat *x0, gfloat *xs, gfloat *q,  gint nq,
			      gfloat *cfft, gint cstr, gfloat *work) ;
gint wbfmm_expansion_dipole_h_cfft_f(gfloat k, gint N, 
				     gfloat *x0,
				     gfloat *xs,
				     gfloat *fx,
				     gfloat *fy,
				     gfloat *fz,
				     gint nq,				     
				     gfloat *cfft, gint cstr,
				     gfloat *work) ;
gint wbfmm_expansion_normal_h_cfft_f(gfloat k, gint N, 
				     gfloat *x0,
				     gfloat *xs,
				     gfloat *normal,
				     gfloat *q, gint nq,
				     gfloat *cfft, gint cstr,
				     gfloat *work) ;

gint wbfmm_expansion_h_evaluate_f(gfloat k, gfloat *x0,
				  gfloat *cfft, 
				  gint cstr,
				  gint N, gint nq,
				  gfloat *xf, gfloat *field,
				  gint fstr,
				  gfloat *work) ;
gint wbfmm_expansion_h_grad_evaluate_f(gfloat k, gfloat *x0,
				       gfloat *cfft, 
				       gint cstr,
				       gint N, gint nq,
				       gfloat *xf,
				       gfloat *field, gint fstr,
				       gfloat *work) ;
gint wbfmm_expansion_j_evaluate_f(gfloat k, gfloat *x0,
				  gfloat *cfft, 
				  gint cstr,
				  gint N, gint nq,
				  gfloat *xf, gfloat *field,
				  gint fstr, 
				  gfloat *work) ;
gint wbfmm_expansion_j_grad_evaluate_f(gfloat k, gfloat *x0,
				       gfloat *cfft, 
				       gint cstr,
				       gint N, gint nq,
				       gfloat *xf,
				       gfloat *field, gint fstr,
				       gfloat *work) ;

gint wbfmm_total_field_f(gfloat k,
			 gfloat *xs, gint xstride,
			 gfloat *src, gint sstride,
			 gfloat *normals, gint nstr,
			 gfloat *dipoles, gint dstr,
			 gint nq, gint nsrc,
			 gfloat *xf, gfloat *field, gint fstr) ;
gint wbfmm_total_field_grad_f(gfloat k,
			      gfloat *xs, gint xstride,
			      gfloat *src, gint sstride,
			      gfloat *normals, gint nstr,
			      gfloat *dipoles, gint dstr,
			      gint nq, gint nsrc,
			      gfloat *xf, gfloat *field, gint fstr) ;
gint wbfmm_total_dipole_field_f(gfloat k,
				gfloat *xs, gint xstride,
				gfloat *src, gint sstride,
				gint nsrc,
				gfloat *xf, gfloat *field) ;
gint wbfmm_total_normal_field_f(gfloat k,
				gfloat *xs, gint xstride,
				gfloat *ns, gint nstride,
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
gint wbfmm_coaxial_translate_ref_f(gfloat *Co, gint cstro, gint No,
				   gfloat *Ci, gint cstri, gint Ni,
				   gint nq,
				   gfloat *cfft, gint L,
				   gboolean complex, gfloat sc) ;

gint wbfmm_rotation_angles_f(gfloat *ix, gfloat *iy, gfloat *iz, 
			   gfloat *jx, gfloat *jy, gfloat *jz, 
			   gfloat *th, gfloat *ph, gfloat *ch) ;
gint wbfmm_coefficients_H_rotation_f(gfloat *H, gint N, gfloat th,
				     gfloat *work) ;

gint wbfmm_rotate_H_ref_f(gfloat *Co, gint cstro, 
			  gfloat *Ci, gint cstri,
			  gint N, gint nq, gfloat *H,
			  gfloat ph, gfloat ch, gfloat sc) ;
gint wbfmm_rotate_H_avx_f(gfloat *Co, gint cstro, 
			  gfloat *Ci, gint cstri,
			  gint N, gint nq, gfloat *H,
			  gfloat ph, gfloat ch, gfloat sc) ;

/*indexing and octrees*/
guint64 wbfmm_point_index_3d_f(gfloat *x, gfloat *c, gfloat D) ;

wbfmm_tree_t *wbfmm_tree_new_f(gfloat *x, gfloat D, guint maxpoints) ;
gint wbfmm_tree_coefficient_init_f(wbfmm_tree_t *t,
				   guint l, guint nr, guint ns) ;
gint wbfmm_tree_coefficient_clear_f(wbfmm_tree_t *t, guint l) ;
gint wbfmm_tree_leaf_expansions_f(wbfmm_tree_t *t, gfloat k,
				  gfloat *src, gint sstr,
				  gfloat *dipoles, gint dstr,
				  gboolean zero_expansions,
				  gfloat *work) ;

gint wbfmm_tree_refine_f(wbfmm_tree_t *t) ;
gint wbfmm_tree_add_points_f(wbfmm_tree_t *t, 
			     gpointer pts, gsize pstr,
			     gpointer normals, gsize nstr,
			     guint npts, gboolean sorted) ;
gint wbfmm_tree_sort_points_f(wbfmm_tree_t *t, 
			      gpointer pts, gsize psize,
			      guint npts) ;
gint wbfmm_truncation_number_f(wbfmm_tree_t *t, gfloat k, guint level,
			       gfloat tol) ;

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
				gint nq,
				gfloat *work) ;
gint wbfmm_parent_child_shift_f(gfloat *Cc, gint Nc,
				gfloat *Cp, gint Np,
				gfloat *H03, 
				gfloat *H47, gint Lh,
				gfloat *shift, gint Ls,
				gint nq,
				gfloat *work) ;
gint wbfmm_points_origin_width_f(gfloat *x, gint str, gint n,
				 gfloat *xmin, gfloat *xmax, gfloat *D,
				 gboolean init_limits) ;

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
						     gboolean bw,
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
gint wbfmm_downward_pass_ref_f(wbfmm_tree_t *t,
			       wbfmm_shift_operators_t *op,
			       guint level, gfloat *work, gint nthreads) ;
gint wbfmm_tree_box_field_f(wbfmm_tree_t *t, guint level,
			    guint b, gfloat k,
			    gfloat *x, gfloat *f, gint fstr, gfloat *work) ;
gint wbfmm_tree_box_local_field_f(wbfmm_tree_t *t, guint level,
				  guint b, gfloat k,
				  gfloat *x, gfloat *f, gint fstr,
				  gfloat *src, gint sstr,
				  /* gfloat *normals, gint nstr, */
				  gfloat *d, gint dstr,
				  gboolean eval_neighbours,
				  guint field,
				  gfloat *work) ;
guint64 wbfmm_point_box_f(wbfmm_tree_t *t, guint level, gfloat *x) ;

wbfmm_target_list_t *wbfmm_target_list_new(wbfmm_tree_t *t, guint npts) ;
wbfmm_target_list_t *wbfmm_target_list_new_f(wbfmm_tree_t *t, guint npts) ;

gint wbfmm_target_list_add_points(wbfmm_target_list_t *l,
				  gpointer pts, gsize pstr,
				  guint npts) ;				  
gint wbfmm_target_list_add_points_f(wbfmm_target_list_t *l,
				    gpointer pts, gsize pstr,
				    guint npts) ;
gint wbfmm_laplace_target_list_local_coefficients(wbfmm_target_list_t *l,
						  guint source, gdouble *work) ;
gint wbfmm_laplace_target_list_local_coefficients_f(wbfmm_target_list_t *l,
						    guint source,
						    gfloat *work) ;
gint wbfmm_target_list_local_coefficients(wbfmm_target_list_t *l,
					  gdouble k,
					  gdouble *work) ;
gint wbfmm_target_list_local_coefficients_f(wbfmm_target_list_t *l,
					    gfloat k,
					    gfloat *work) ;
gint wbfmm_target_list_local_field(wbfmm_target_list_t *l,
				   gdouble *src, gint sstr,
				   gdouble *nsrc, gint nstr,
				   gdouble *f, gint fstr) ;
gint wbfmm_target_list_local_field_f(wbfmm_target_list_t *l,
				     gfloat *src, gint sstr,
				     gfloat *nsrc, gint nstr,
				     gfloat *f, gint fstr) ;

/*precision independent functions*/
guint64 wbfmm_point_locate_box(guint64 x, guint level) ;
gint wbfmm_point_from_index(guint64 i, guint32 *x, guint32 *y, guint32 *z) ;
guint64 wbfmm_box_index(guint32 i, guint32 j, guint32 k) ;
gint wbfmm_box_location(guint64 idx, guint32 *i, guint32 *j, guint32 *k) ;
guint64 wbfmm_box_parent(guint64 idx) ;
guint64 wbfmm_box_first_child(guint64 idx) ;
gint wbfmm_tree_print(FILE *f, wbfmm_tree_t *t, guint level,
		      gboolean print_empty) ;
gint wbfmm_tree_print_f(FILE *f, wbfmm_tree_t *t, guint level,
			gboolean print_empty) ;
gint wbfmm_logging_init(FILE *f, char *p, 
			GLogLevelFlags log_level,
			gpointer exit_func, gboolean timed) ;
gint wbfmm_box_neighbours(guint level, guint64 idx, guint64 *neighbours) ;
gint wbfmm_box_interaction_list_4(guint level, guint64 idx, 
				  guint64 *list, gboolean sort) ;
gint wbfmm_box_interaction_grid_4(guint level, guint64 idx, guint64 list[]) ;
gint wbfmm_box_interaction_index(gint i, gint j, gint k) ;
gint wbfmm_library_config(wbfmm_library_config_t *c) ;
gint wbfmm_library_config_print(wbfmm_library_config_t *c, FILE *f) ;
gint wbfmm_tree_coefficients_zero(wbfmm_tree_t *t, guint level) ;

/*compile time switches for compiler options*/
#ifdef WBFMM_USE_AVX

#define wbfmm_rotate_H(_Co,_cstro,_Ci,_cstri,_N,_nq,_H,_ph,_ch,_sc)	\
  wbfmm_rotate_H_avx(_Co,_cstro,_Ci,_cstri,_N,_nq,_H,_ph,_ch,_sc)

#define wbfmm_rotate_H_f(_Co,_cstro,_Ci,_cstri,_N,_nq,_H,_ph,_ch,_sc)	\
  wbfmm_rotate_H_ref_f(_Co,_cstro,_Ci,_cstri,_N,_nq,_H,_ph,_ch,_sc)

#define wbfmm_coaxial_translate(_Co,_cstro,_No,_Ci,_cstri,_Ni,_nq,_cft,_L,_c,_sc) \
  wbfmm_coaxial_translate_ref(_Co,_cstro,_No,_Ci,_cstri,_Ni,_nq,_cft,_L,_c,_sc)
#define wbfmm_coaxial_translate_f(_Co,_cstro,_No,_Ci,_cstri,_Ni,_nq,_cft,_L,_c,_sc) \
  wbfmm_coaxial_translate_ref_f(_Co,_cstro,_No,_Ci,_cstri,_Ni,_nq,_cft,_L,_c,_sc)

#define wbfmm_downward_pass(_t,_op,_level,_work,_nth)	\
  wbfmm_downward_pass_avx(_t,_op,_level,_work,_nth)
#define wbfmm_downward_pass_f(_t,_op,_level,_work,_nth)		\
  wbfmm_downward_pass_ref_f(_t,_op,_level,_work,_nth)

#define wbfmm_laplace_rotate_H(_Co,_cstro,_Ci,_cstri,_N,_nq,_H,_ph,_ch,_sc) \
  wbfmm_laplace_rotate_H_avx(_Co,_cstro,_Ci,_cstri,_N,_nq,_H,_ph,_ch,_sc)
#define wbfmm_laplace_rotate_H_f(_Co,_cstro,_Ci,_cstri,_N,_nq,_H,_ph,_ch,_sc) \
  wbfmm_laplace_rotate_H_ref_f(_Co,_cstro,_Ci,_cstri,_N,_nq,_H,_ph,_ch,_sc)

#define wbfmm_laplace_downward_pass(_t,_op,_level,_work,_nth)	\
  wbfmm_laplace_downward_pass_avx(_t,_op,_level,_work,_nth)
#define wbfmm_laplace_downward_pass_f(_t,_op,_level,_work,_nth)	\
  wbfmm_laplace_downward_pass_ref_f(_t,_op,_level,_work,_nth)

#else

#define wbfmm_rotate_H(_Co,_cstro,_Ci,_cstri,_N,_nq,_H,_ph,_ch,_sc)	\
  wbfmm_rotate_H_ref(_Co,_cstro,_Ci,_cstri,_N,_nq,_H,_ph,_ch,_sc)

#define wbfmm_rotate_H_f(_Co,_cstro,_Ci,_cstri,_N,_nq,_H,_ph,_ch,_sc)	\
  wbfmm_rotate_H_ref_f(_Co,_cstro,_Ci,_cstri,_N,_nq,_H,_ph,_ch,_sc)

#define wbfmm_coaxial_translate(_Co,_cstro,_No,_Ci,_cstri,_Ni,_nq,_cft,_L,_c,_sc) \
  wbfmm_coaxial_translate_ref(_Co,_cstro,_No,_Ci,_cstri,_Ni,_nq,_cft,_L,_c,_sc)
#define wbfmm_coaxial_translate_f(_Co,_cstro,_No,_Ci,_cstri,_Ni,_nq,_cft,_L,_c,_sc) \
  wbfmm_coaxial_translate_ref_f(_Co,_cstro,_No,_Ci,_cstri,_Ni,_nq,_cft,_L,_c,_sc)

#define wbfmm_downward_pass(_t,_op,_level,_work,_nth)	\
  wbfmm_downward_pass_ref(_t,_op,_level,_work,_nth)
#define wbfmm_downward_pass_f(_t,_op,_level,_work,_nth)	\
  wbfmm_downward_pass_ref_f(_t,_op,_level,_work,_nth)

#define wbfmm_laplace_rotate_H(_Co,_cstro,_Ci,_cstri,_N,_nq,_H,_ph,_ch,_sc) \
  wbfmm_laplace_rotate_H_ref(_Co,_cstro,_Ci,_cstri,_N,_nq,_H,_ph,_ch,_sc)
#define wbfmm_laplace_rotate_H_f(_Co,_cstro,_Ci,_cstri,_N,_nq,_H,_ph,_ch,_sc) \
  wbfmm_laplace_rotate_H_ref_f(_Co,_cstro,_Ci,_cstri,_N,_nq,_H,_ph,_ch,_sc)

#define wbfmm_laplace_downward_pass(_t,_op,_level,_work,_nth)	\
  wbfmm_laplace_downward_pass_ref(_t,_op,_level,_work,_nth)
#define wbfmm_laplace_downward_pass_f(_t,_op,_level,_work,_nth)	\
  wbfmm_laplace_downward_pass_ref_f(_t,_op,_level,_work,_nth)

#endif

#endif /*WBFMM_H_INCLUDED*/
