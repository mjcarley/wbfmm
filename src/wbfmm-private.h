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

#ifndef WBFMM_PRIVATE_H_INCLUDED 
#define WBFMM_PRIVATE_H_INCLUDED

#include <stdio.h>

#ifdef WBFMM_SINGLE_PRECISION

#define WBFMM_REAL gfloat

#define WBFMM_FUNCTION_NAME(_func) _func##_f

#define SQRT(_x) sqrtf((_x))
#define CBRT(_x) cbrtf((_x))
#define SIN(_x) sinf((_x))
#define COS(_x) cosf((_x))
#define ACOS(_x) acosf((_x))
#define ATAN(_x) atanf((_x))
#define ATAN2(_y,_x) atan2f((_y),(_x))
#define LOG(_x) logf((_x))

#else

#define WBFMM_REAL gdouble

#define WBFMM_FUNCTION_NAME(_func) _func

#define SQRT(_x) sqrt((_x))
#define CBRT(_x) cbrt((_x))
#define SIN(_x) sin((_x))
#define COS(_x) cos((_x))
#define ACOS(_x) acos((_x))
#define ATAN(_x) atan((_x))
#define ATAN2(_y,_x) atan2((_y),(_x))
#define LOG(_x) log((_x))

#endif /*WBFMM_SINGLE_PRECISION*/

#define wbfmm_cos_sin_recursion(_Cn,_Sn,_C,_S)	\
  do { WBFMM_REAL _tmp = (_Cn) ;		\
  (_Cn) = (_Cn)*(_C) - (_Sn)*(_S) ;		\
  (_Sn) = (_Sn)*(_C) + (_tmp)*(_S) ;		\
  } while (0)

#define wbfmm_tree_point_index(_t,_i)		\
  ((WBFMM_REAL *)(&((_t)->points[(_i)*((_t)->pstr)])))
#define wbfmm_target_list_point_index(_t,_i)		\
  ((WBFMM_REAL *)(&((_t)->points[(_i)*((_t)->pstr)])))
/* #define wbfmm_tree_origin(_t) ((WBFMM_REAL *)(&((_t)->x[0])))  */

#define wbfmm_rotation_index_numn(_nu,_m,_n)			\
  ((_n)*(4*(_n)-1)*((_n)+1)/6 + ((_nu)+(_n))*((_n)+1) + (_m))

#define wbfmm_index_laplace_nm(_n,_m) ((_n)*(_n)+(2*(_m))-1)

#define _wbfmm_SS_coefficient_index_nmnu(_n,_m,_nu)	\
  ((_n)*((_n)+2)*((_n)+1)/6 + (_m)*(2*(_n)-(_m)+1)/2 + (_nu))

#define wbfmm_coaxial_translation_SS_cfft(_n, _nd, _m)	\
  (_wbfmm_SS_coefficients_laplace[_wbfmm_SS_coefficient_index_nmnu((_n),(_m),(_nd))])

#define _wbfmm_RR_coefficient_index_nmnu(_n,_m,_nu)		    \
  ((_wbfmm_translation_Nmax+1)*(_wbfmm_translation_Nmax+2)/2*(_m) + \
   (_nu)*((_nu)+1)/2 + (_n))

#define wbfmm_coaxial_translation_RR_cfft(_n, _nd, _m)	\
  (_wbfmm_RR_coefficients_laplace[_wbfmm_RR_coefficient_index_nmnu((_n),(_m),(_nd))])

#define _wbfmm_SR_coefficient_index_nmnu(_n,_m,_nu)			\
  ((6*(_wbfmm_translation_Nmax+1)*					\
    (_wbfmm_translation_Nmax+2)+1-2*(_m)*(2*_wbfmm_translation_Nmax+3)+	\
    2*(_m)*(_m))*(_m)/6 +						\
   ((_n)-(_m))*(_wbfmm_translation_Nmax+1-(_m)) - (_m) + (_nu))
#define wbfmm_coaxial_translation_SR_cfft(_n, _nd, _m)	\
  (_wbfmm_SR_coefficients_laplace[_wbfmm_SR_coefficient_index_nmnu((_n),(_m),(_nd))])
  
/* extern WBFMM_REAL *_wbfmm_SS_coefficients_laplace, */
/*   *_wbfmm_RR_coefficients_laplace, *_wbfmm_SR_coefficients_laplace ; */
/* extern gint _wbfmm_translation_Nmax ; */
/* extern gint _wbfmm_shift_angles[] ; */
/* extern WBFMM_REAL _wbfmm_shifts_ph[], _wbfmm_shifts_ch[], _wbfmm_shifts_r[] ; */

#define WBFMM_SHIFTS_R_NUMBER 15

#define ABSDIFF(_i,_j) ((_i) >= (_j) ? ((_i)-(_j)) : ((_j)-(_i)))
#define IS_EVEN(_i) (((_i)%2==0)?1:0)

#define BITMASK_0000000001000001000001000001000001000001000001000001000001000001 UINT64_C(18300341342965825)
#define BITMASK_0000001000001000001000001000001000001000001000001000001000001000 UINT64_C(146402730743726600)
#define BITMASK_0001000000000000000000000000000000000000000000000000000000000000 UINT64_C(1152921504606846976)
/* 0000000ccc0000cc0000cc0000cc0000cc0000cc0000cc0000cc0000cc0000cc */
#define BITMASK_0000000000000011000000000011000000000011000000000011000000000011 UINT64_C(844631138906115)
#define BITMASK_0000000111000000000011000000000011000000000011000000000011000000 UINT64_C(126113986927919296)
/* 00000000000ccccc00000000cccc00000000cccc00000000cccc00000000cccc */
#define BITMASK_0000000000000000000000000000000000001111000000000000000000001111 UINT64_C(251658255)
#define BITMASK_0000000000000000000000001111000000000000000000001111000000000000 UINT64_C(1030792212480)
#define BITMASK_0000000000011111000000000000000000000000000000000000000000000000 UINT64_C(8725724278030336)
/* 000000000000000000000000000ccccccccccccc0000000000000000cccccccc */
#define BITMASK_0000000000000000000000000000000000000000000000000000000011111111 UINT64_C(255)
#define BITMASK_0000000000000000000000000001111111111111000000000000000000000000 UINT64_C(137422176256)
/* ccccccccccccccccccccc */
#define BITMASK_21BITS  UINT64_C(2097151)

extern const gdouble WBFMM_FACTORIALS[] ;
extern const gfloat WBFMM_FACTORIALS_F[] ;

#ifdef WBFMM_SINGLE_PRECISION
#define wbfmm_factorial(_n) ((WBFMM_FACTORIALS_F[(_n)]))
#else  /*WBFMM_SINGLE_PRECISION*/
#define wbfmm_factorial(_n) ((WBFMM_FACTORIALS[(_n)]))
#endif /*WBFMM_SINGLE_PRECISION*/

#define wbfmm_coaxial_index_lmn(_l,_m,_n)		\
  ((_l)*((_l)+1)*((_l)+2)/6 + (_n)*((_n)+1)/2 + (_m))

/*G&D (4.81) and (4.82)*/
#define wbfmm_coaxial_index(_l,_m,_n)			     \
  ((_l)>=(_n) ? wbfmm_coaxial_index_lmn((_l),ABS(_m),(_n)) : \
   wbfmm_coaxial_index_lmn((_n),ABS(_m),(_l)))
#define wbfmm_coaxial_index_sgn(_l,_m,_n)	\
  ((_l)>=(_n) ? 1 : minus_one_pow((_n)+(_l)))

#define minus_one_pow(_n) ((2*((_n)/2) == (_n) ? 1 : -1))

#define yes_if_true(_t)  ((_t) == TRUE ? "yes" : "no") 

/* gint print_bits_uint(FILE *f, guint x) ; */

gint _wbfmm_bessel_j_scaled_init(gdouble x, gdouble *j0, gdouble *j1) ;

gdouble recursion_anm(gint n, gint m) ;
gdouble recursion_bnm(gint n, gint m) ;

gfloat recursion_anm_f(gint n, gint m) ;
gfloat recursion_bnm_f(gint n, gint m) ;

#define WBFMM_DOWNWARD_PASS_DATA_SIZE  8
#define WBFMM_DOWNWARD_PASS_LEVEL      0
#define WBFMM_DOWNWARD_PASS_WORK       1
#define WBFMM_DOWNWARD_PASS_NQ         2
#define WBFMM_DOWNWARD_PASS_NTHREAD    3
#define WBFMM_DOWNWARD_PASS_TREE       4
#define WBFMM_DOWNWARD_PASS_OP         5

extern gdouble CmPI_4[], CnPI_2[] ;
extern gfloat CmPI_4f[], CnPI_2f[] ;
extern gdouble *_wbfmm_SS_coefficients_laplace,
  *_wbfmm_RR_coefficients_laplace,
  *_wbfmm_SR_coefficients_laplace ;
extern gfloat *_wbfmm_SS_coefficients_laplace_f,
  *_wbfmm_RR_coefficients_laplace_f,
  *_wbfmm_SR_coefficients_laplace_f ;
extern gint _wbfmm_translation_Nmax ;
extern gdouble _wbfmm_shifts_th[], _wbfmm_shifts_ph[], _wbfmm_shifts_r[] ;
extern gfloat _wbfmm_shifts_th_f[], _wbfmm_shifts_ph_f[], _wbfmm_shifts_r_f[] ;
extern gint _wbfmm_shift_angles[] ;

#ifdef WBFMM_SINGLE_PRECISION
#define cos_n_PI_4(_n) (CmPI_4f[(_n)%8])
#define sin_n_PI_4(_n) (CmPI_4f[((_n)+6)%8])
#define cos_n_PI_2(_n) (CnPI_2f[(_n)%4])
#define sin_n_PI_2(_n) (CnPI_2f[((_n)+3)%4])
#else /*WBFMM_SINGLE_PRECISION*/
#define cos_n_PI_4(_n) (CmPI_4[(_n)%8])
#define sin_n_PI_4(_n) (CmPI_4[((_n)+6)%8])
#define cos_n_PI_2(_n) (CnPI_2[(_n)%4])
#define sin_n_PI_2(_n) (CnPI_2[((_n)+3)%4])
#endif

#endif /*WBFMM_PRIVATE_H_INCLUDED*/
