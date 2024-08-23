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

#ifndef WBFMM_DERIVATIVES_H_INCLUDED 
#define WBFMM_DERIVATIVES_H_INCLUDED

#define WBFMM_DERIVATIVE_0_R    0 
#define WBFMM_DERIVATIVE_0_I    1 
#define WBFMM_DERIVATIVE_X_R    0
#define WBFMM_DERIVATIVE_X_I    1
#define WBFMM_DERIVATIVE_Y_R    2
#define WBFMM_DERIVATIVE_Y_I    3
#define WBFMM_DERIVATIVE_Z_R    4
#define WBFMM_DERIVATIVE_Z_I    5

#define WBFMM_DERIVATIVE_XX_R   0
#define WBFMM_DERIVATIVE_XX_I   1
#define WBFMM_DERIVATIVE_YY_R   2
#define WBFMM_DERIVATIVE_YY_I   3
#define WBFMM_DERIVATIVE_ZZ_R   4
#define WBFMM_DERIVATIVE_ZZ_I   5
#define WBFMM_DERIVATIVE_XY_R   6
#define WBFMM_DERIVATIVE_XY_I   7
#define WBFMM_DERIVATIVE_YZ_R   8
#define WBFMM_DERIVATIVE_YZ_I   9
#define WBFMM_DERIVATIVE_ZX_R  10
#define WBFMM_DERIVATIVE_ZX_I  11

/*
 * second partial derivatives of Rnm, for m >= 2
 *
 * where appropriate, derivatives are multiplied by two to account for 
 * symmetry in Fourier coefficients
 */

#define Rnm_derivatives_2(_n,_m,_rnm2,_Pnm2,_Cmph,_Smph,_dRnm)		\
  do  {									\
  gint _off = 2 ;							\
  WBFMM_REAL _anm2[5], _Rnm2[5] ;					\
    g_assert((_m) >= 2) ;						\
									\
    _anm2[_off-2] = SQRT((WBFMM_REAL)(2*(_n)+1)/(2*(_n)-3)*((_n)+(_m))*	\
			 ((_n)+(_m)-1)*((_n)+(_m)-2)*((_n)+(_m)-3)) ;	\
    _anm2[_off-1] = SQRT((WBFMM_REAL)(2*(_n)+1)/(2*(_n)-3)*		\
			 (n*(_n)-(_m)*m)*((_n)+(_m)-1)*((_n)+(_m)-2)) ;	\
    _anm2[_off+0] = SQRT((WBFMM_REAL)(2*(_n)+1)/(2*(_n)-3)*		\
			 (n*(_n)-(_m)*m)*((n-1)*(n-1)-m*m)) ;		\
    _anm2[_off+1] = SQRT((WBFMM_REAL)(2*(_n)+1)/(2*(_n)-3)*		\
			 (n*(_n)-(_m)*m)*((_n)-(_m)-1)*((_n)-(_m)-2)) ;	\
    _anm2[_off+2] = SQRT((WBFMM_REAL)(2*(_n)+1)/(2*(_n)-3)*((_n)-(_m))*	\
			 ((_n)-(_m)-1)*((_n)-(_m)-2)*((_n)-(_m)-3)) ;	\
    _Rnm2[_off-2] = (_rnm2)*(_Pnm2)[(_m)-2] ;				\
    _Rnm2[_off-1] = (_rnm2)*(_Pnm2)[(_m)-1] ;				\
    _Rnm2[_off+0] = (_rnm2)*(_Pnm2)[(_m)+0] ;				\
    _Rnm2[_off+1] = (_rnm2)*(_Pnm2)[(_m)+1] ;				\
    _Rnm2[_off+2] = (_rnm2)*(_Pnm2)[(_m)+2] ;				\
    									\
    (_dRnm)[WBFMM_DERIVATIVE_XX_R] =					\
      2.0*(_anm2[_off+2]*_Rnm2[_off+2]*(_Cmph)[(_m)+2] -		\
	   2.0*_anm2[_off+0]*_Rnm2[_off+0]*(_Cmph)[(_m)+0] +		\
	   _anm2[_off-2]*_Rnm2[_off-2]*(_Cmph)[(_m)-2])*0.25 ;		\
    (_dRnm)[WBFMM_DERIVATIVE_XX_I] =					\
      2.0*(_anm2[_off+2]*_Rnm2[_off+2]*(_Smph)[(_m)+2] -		\
	   2.0*_anm2[_off+0]*_Rnm2[_off+0]*(_Smph)[(_m)+0] +		\
	   _anm2[_off-2]*_Rnm2[_off-2]*(_Smph)[(_m)-2])*0.25 ;		\
    (_dRnm)[WBFMM_DERIVATIVE_YY_R] =					\
      -2.0*(_anm2[_off+2]*_Rnm2[_off+2]*(_Cmph)[(_m)+2] +		\
	    2.0*_anm2[_off+0]*_Rnm2[_off+0]*(_Cmph)[(_m)+0] +		\
	    _anm2[_off-2]*_Rnm2[_off-2]*(_Cmph)[(_m)-2])*0.25 ;		\
    (_dRnm)[WBFMM_DERIVATIVE_YY_I] =					\
      -2.0*(_anm2[_off+2]*_Rnm2[_off+2]*(_Smph)[(_m)+2] +		\
	    2.0*_anm2[_off+0]*_Rnm2[_off+0]*(_Smph)[(_m)+0] +		\
	    _anm2[_off-2]*_Rnm2[_off-2]*(_Smph)[(_m)-2])*0.25 ;		\
    									\
    (_dRnm)[WBFMM_DERIVATIVE_ZZ_R] =					\
      _anm2[_off+0]*_Rnm2[_off+0]*(_Cmph)[(_m)+0]*2 ;			\
    (_dRnm)[WBFMM_DERIVATIVE_ZZ_I] = \
      _anm2[_off+0]*_Rnm2[_off+0]*(_Smph)[(_m)+0]*2 ;			\
    									\
    (_dRnm)[WBFMM_DERIVATIVE_XY_R] =					\
      2.0*(_anm2[_off+2]*_Rnm2[_off+2]*(_Smph)[(_m)+2] -		\
	   _anm2[_off-2]*_Rnm2[_off-2]*(_Smph)[(_m)-2])*0.25 ;		\
    (_dRnm)[WBFMM_DERIVATIVE_XY_I] =					\
      2.0*(-_anm2[_off+2]*_Rnm2[_off+2]*(_Cmph)[(_m)+2] +		\
	   _anm2[_off-2]*_Rnm2[_off-2]*(_Cmph)[(_m)-2])*0.25 ;		\
    (_dRnm)[WBFMM_DERIVATIVE_YZ_R] =					\
      -_anm2[_off+1]*_Rnm2[_off+1]*(_Smph)[(_m)+1] -			\
      _anm2[_off-1]*_Rnm2[_off-1]*(_Smph)[(_m)-1] ;			\
    (_dRnm)[WBFMM_DERIVATIVE_YZ_I] =					\
      _anm2[_off+1]*_Rnm2[_off+1]*(_Cmph)[(_m)+1] +			\
      _anm2[_off-1]*_Rnm2[_off-1]*(_Cmph)[(_m)-1] ;			\
    									\
    (_dRnm)[WBFMM_DERIVATIVE_ZX_R] =					\
      -_anm2[_off+1]*_Rnm2[_off+1]*(_Cmph)[(_m)+1] +			\
      _anm2[_off-1]*_Rnm2[_off-1]*(_Cmph)[(_m)-1] ;			\
    (_dRnm)[WBFMM_DERIVATIVE_ZX_I] =					\
      -_anm2[_off+1]*_Rnm2[_off+1]*(_Smph)[(_m)+1] +			\
      _anm2[_off-1]*_Rnm2[_off-1]*(_Smph)[(_m)-1] ;			\
  } while (0)

/*
 * second partial derivatives of Rnm, for m == 1 (this requires
 * dealing with special cases arising from negative m)
 *
 * where appropriate, derivatives are multiplied by two to account for 
 * symmetry in Fourier coefficients
 */

#define Rnm_derivatives_2m1(_n,_m,_rnm2,_Pnm2,_Cmph,_Smph,_dRnm)	\
  do  {									\
  WBFMM_REAL _Rnm2[5], _anm2[5] ;					\
  gint _off = 2 ;							\
									\
  g_assert((_m) == 1) ;							\
  _anm2[_off-2] = SQRT((WBFMM_REAL)(2*(_n)+1)/(2*(_n)-3)*((_n)+(_m))*	\
		       ((_n)+(_m)-1)*((_n)+(_m)-2)*((_n)+(_m)-3)) ;	\
  _anm2[_off-1] = SQRT((WBFMM_REAL)(2*(_n)+1)/(2*(_n)-3)*		\
		       (n*(_n)-(_m)*m)*((_n)+(_m)-1)*((_n)+(_m)-2)) ;	\
  _anm2[_off+0] = SQRT((WBFMM_REAL)(2*(_n)+1)/(2*(_n)-3)*		\
		       (n*(_n)-(_m)*m)*((n-1)*(n-1)-m*m)) ;		\
  _anm2[_off+1] = SQRT((WBFMM_REAL)(2*(_n)+1)/(2*(_n)-3)*		\
		       (n*(_n)-(_m)*m)*((_n)-(_m)-1)*((_n)-(_m)-2)) ;	\
  _anm2[_off+2] = SQRT((WBFMM_REAL)(2*(_n)+1)/(2*(_n)-3)*((_n)-(_m))*	\
		       ((_n)-(_m)-1)*((_n)-(_m)-2)*((_n)-(_m)-3)) ;	\
  /*DO NOT CHANGE THIS: it comes from the negative order mode*/		\
  _Rnm2[_off-2] = (_rnm2)*(_Pnm2)[1] ;					\
  _Rnm2[_off-1] = (_rnm2)*(_Pnm2)[(_m)-1] ;				\
  _Rnm2[_off+0] = (_rnm2)*(_Pnm2)[(_m)+0] ;				\
  _Rnm2[_off+1] = (_rnm2)*(_Pnm2)[(_m)+1] ;				\
  _Rnm2[_off+2] = (_rnm2)*(_Pnm2)[(_m)+2] ;				\
									\
  /*OR THIS (note sign of third terms in brackets)*/			\
  (_dRnm)[WBFMM_DERIVATIVE_XX_R] =					\
    2.0*(_anm2[_off+2]*_Rnm2[_off+2]*(_Cmph)[(_m)+2] -			\
	 2.0*_anm2[_off+0]*_Rnm2[_off+0]*(_Cmph)[(_m)+0] -		\
	 _anm2[_off-2]*_Rnm2[_off-2]*(_Cmph)[1])*0.25 ;			\
  (_dRnm)[WBFMM_DERIVATIVE_XX_I] =					\
    2.0*(_anm2[_off+2]*_Rnm2[_off+2]*(_Smph)[(_m)+2] -			\
	 2.0*_anm2[_off+0]*_Rnm2[_off+0]*(_Smph)[(_m)+0] +		\
	 _anm2[_off-2]*_Rnm2[_off-2]*(_Smph)[1])*0.25 ;			\
  (_dRnm)[WBFMM_DERIVATIVE_YY_R] =					\
    -2.0*(_anm2[_off+2]*_Rnm2[_off+2]*(_Cmph)[(_m)+2] +			\
	  2.0*_anm2[_off+0]*_Rnm2[_off+0]*(_Cmph)[(_m)+0] -		\
	  _anm2[_off-2]*_Rnm2[_off-2]*(_Cmph)[1])*0.25 ;		\
  (_dRnm)[WBFMM_DERIVATIVE_YY_I] =					\
    -2.0*(_anm2[_off+2]*_Rnm2[_off+2]*(_Smph)[(_m)+2] +			\
	  2.0*_anm2[_off+0]*_Rnm2[_off+0]*(_Smph)[(_m)+0] +		\
	  _anm2[_off-2]*_Rnm2[_off-2]*(_Smph)[1])*0.25 ;		\
  (_dRnm)[WBFMM_DERIVATIVE_ZZ_R] =					\
    _anm2[_off+0]*_Rnm2[_off+0]*(_Cmph)[(_m)+0]*2 ;			\
  (_dRnm)[WBFMM_DERIVATIVE_ZZ_I] =					\
    _anm2[_off+0]*_Rnm2[_off+0]*(_Smph)[(_m)+0]*2 ;			\
									\
  /*LEAVE THIS ALONE TOO (note the signs on the nm2mm2 terms)*/		\
  (_dRnm)[WBFMM_DERIVATIVE_XY_R] =					\
    2.0*( _anm2[_off+2]*_Rnm2[_off+2]*(_Smph)[(_m)+2] -			\
	  _anm2[_off-2]*_Rnm2[_off-2]*(_Smph)[1])*0.25 ;		\
  (_dRnm)[WBFMM_DERIVATIVE_XY_I] =					\
    2.0*(-_anm2[_off+2]*_Rnm2[_off+2]*(_Cmph)[(_m)+2] -			\
	 _anm2[_off-2]*_Rnm2[_off-2]*(_Cmph)[1])*0.25 ;			\
									\
  (_dRnm)[WBFMM_DERIVATIVE_YZ_R] =					\
    -_anm2[_off+1]*_Rnm2[_off+1]*(_Smph)[(_m)+1] -			\
    _anm2[_off-1]*_Rnm2[_off-1]*(_Smph)[(_m)-1] ;			\
  (_dRnm)[WBFMM_DERIVATIVE_YZ_I] =					\
    _anm2[_off+1]*_Rnm2[_off+1]*(_Cmph)[(_m)+1] +			\
    _anm2[_off-1]*_Rnm2[_off-1]*(_Cmph)[(_m)-1] ;			\
  									\
  (_dRnm)[WBFMM_DERIVATIVE_ZX_R] =					\
    -_anm2[_off+1]*_Rnm2[_off+1]*(_Cmph)[(_m)+1] +			\
    _anm2[_off-1]*_Rnm2[_off-1]*(_Cmph)[(_m)-1] ;			\
  (_dRnm)[WBFMM_DERIVATIVE_ZX_I] =					\
    -_anm2[_off+1]*_Rnm2[_off+1]*(_Smph)[(_m)+1] +			\
    _anm2[_off-1]*_Rnm2[_off-1]*(_Smph)[(_m)-1] ;			\
  } while (0)

/*
 * second partial derivatives of Rnm, for m == 0 (this requires
 * dealing with special cases arising from negative m)
 *
 * where appropriate, derivatives are multiplied by two to account for 
 * symmetry in Fourier coefficients
 */

#define Rnm_derivatives_2m0(_n,_m,_rnm2,_Pnm2,_Cmph,_Smph,_dRnm)	\
  do  {									\
  WBFMM_REAL _anm2[5], _Rnm2[5] ;					\
  gint _off = 5 ;							\
  g_assert((_m) == 0) ;							\
									\
  _anm2[_off-2] = SQRT((WBFMM_REAL)(2*n+1)/(2*n-3)*((_n)+(_m))*		\
		       ((_n)+(_m)-1)*((_n)+(_m)-2)*((_n)+(_m)-3)) ;	\
  _anm2[_off+0] = SQRT((WBFMM_REAL)(2*n+1)/(2*n-3)*(n*(_n)-(_m)*m)*	\
			 ((n-1)*(n-1)-m*m)) ;				\
  _anm2[_off+1] = SQRT((WBFMM_REAL)(2*n+1)/(2*n-3)*(n*(_n)-(_m)*m)*	\
		       ((_n)-(_m)-1)*((_n)-(_m)-2)) ;			\
  _anm2[_off+2] = SQRT((WBFMM_REAL)(2*n+1)/(2*n-3)*((_n)-(_m))*		\
		       ((_n)-(_m)-1)*((_n)-(_m)-2)*((_n)-(_m)-3)) ;	\
  									\
  _Rnm2[_off-2] = rnm2*Pnm2[2] ;					\
  _Rnm2[_off+0] = rnm2*Pnm2[(_m)+0] ;					\
  _Rnm2[_off+1] = rnm2*Pnm2[(_m)+1] ;					\
  _Rnm2[_off+2] = rnm2*Pnm2[(_m)+2] ;					\
  (_dRnm)[WBFMM_DERIVATIVE_XX_I] = 0.0 ;				\
  (_dRnm)[WBFMM_DERIVATIVE_YY_I] = 0.0 ;				\
  (_dRnm)[WBFMM_DERIVATIVE_ZZ_I] = 0.0 ;				\
  (_dRnm)[WBFMM_DERIVATIVE_XY_I] = 0.0 ;				\
  (_dRnm)[WBFMM_DERIVATIVE_YZ_I] = 0.0 ;				\
  (_dRnm)[WBFMM_DERIVATIVE_ZZ_I] = 0.0 ;				\
									\
  (_dRnm)[WBFMM_DERIVATIVE_XX_R] =					\
    (_anm2[_off+2]*_Rnm2[_off+2]*(_Cmph)[(_m)+2] -			\
     2.0*_anm2[_off+0]*_Rnm2[_off+0]*(_Cmph)[(_m)+0] +			\
     _anm2[_off-2]*_Rnm2[_off-2]*(_Cmph)[2])*0.25 ;			\
  (_dRnm)[WBFMM_DERIVATIVE_YY_R] =					\
    -(_anm2[_off+2]*_Rnm2[_off+2]*(_Cmph)[(_m)+2] +			\
      2.0*_anm2[_off+0]*_Rnm2[_off+0]*(_Cmph)[(_m)+0] +			\
      _anm2[_off-2]*_Rnm2[_off-2]*(_Cmph)[2])*0.25 ;			\
  (_dRnm)[WBFMM_DERIVATIVE_ZZ_R] =					\
    _anm2[_off+0]*_Rnm2[_off+0]*(_Cmph)[(_m)+0] ;			\
									\
  (_dRnm)[WBFMM_DERIVATIVE_XY_R] =					\
    2.0*( _anm2[_off+2]*_Rnm2[_off+2]*(_Smph)[(_m)+2] -			\
	  _anm2[_off-2]*_Rnm2[_off-2]*(_Smph)[(_m)-2])*0.25 ;		\
  (_dRnm)[WBFMM_DERIVATIVE_XY_I] =					\
    2.0*(-_anm2[_off+2]*_Rnm2[_off+2]*(_Cmph)[(_m)+2] +			\
	 _anm2[_off-2]*_Rnm2[_off-2]*(_Cmph)[(_m)-2])*0.25 ;		\
									\
  (_dRnm)[WBFMM_DERIVATIVE_YZ_R] =					\
    -_anm2[_off+1]*_Rnm2[_off+1]*(_Smph)[(_m)+1] ;			\
  (_dRnm)[WBFMM_DERIVATIVE_ZX_R] =					\
    -_anm2[_off+1]*_Rnm2[_off+1]*(_Cmph)[(_m)+1] ;			\
  } while (0) 

/*
 * partial derivatives of Rnm, for m >= 1
 *
 * where appropriate, derivatives are multiplied by two to account for 
 * symmetry in Fourier coefficients
 */

#define Rnm_derivatives_1(_n,_m,_rnm1,_Pnm1,_Cmph,_Smph,_dRnm)	\
  do  {								\
    WBFMM_REAL _anm[3], _Rnm[3] ;				\
    gint _off = 1 ;						\
    								\
    g_assert((_m) >= 1) ;						\
    									\
    _anm[_off+0] = SQRT((WBFMM_REAL)(2*(_n)+1)/(2*(_n)-1)*		\
			((_n)-(_m))*((_n)+(_m))) ;			\
    _anm[_off+1] = SQRT((WBFMM_REAL)(2*(_n)+1)/(2*(_n)-1)*		\
			((_n)-(_m))*((_n)-(_m)-1)) ;			\
    _anm[_off-1] = SQRT((WBFMM_REAL)(2*(_n)+1)/(2*(_n)-1)*		\
			((_n)+(_m))*((_n)+(_m)-1)) ;			\
    _Rnm[_off-1] = (_rnm1)*(_Pnm1)[(_m)-1]*_anm[_off-1] ;		\
    _Rnm[_off+0] = (_rnm1)*(_Pnm1)[(_m)+0]*_anm[_off+0]*2.0 ;		\
    _Rnm[_off+1] = (_rnm1)*Pnm1[(_m)+1]*_anm[_off+1] ;			\
    									\
    (_dRnm)[WBFMM_DERIVATIVE_X_R] =					\
      -_Rnm[_off+1]*(_Cmph)[(_m)+1] + _Rnm[_off-1]*(_Cmph)[(_m)-1] ;	\
    (_dRnm)[WBFMM_DERIVATIVE_X_I] =					\
      -_Rnm[_off+1]*(_Smph)[(_m)+1] + _Rnm[_off-1]*(_Smph)[(_m)-1] ;	\
    									\
    (_dRnm)[WBFMM_DERIVATIVE_Y_R] =					\
      -_Rnm[_off+1]*(_Smph)[(_m)+1] - _Rnm[_off-1]*(_Smph)[(_m)-1] ;	\
    (_dRnm)[WBFMM_DERIVATIVE_Y_I] =					\
      +_Rnm[_off+1]*(_Cmph)[(_m)+1] + _Rnm[_off-1]*(_Cmph)[(_m)-1] ;	\
    									\
    (_dRnm)[WBFMM_DERIVATIVE_Z_R] = +_Rnm[_off+0]*(_Cmph)[(_m)+0] ;	\
    (_dRnm)[WBFMM_DERIVATIVE_Z_I] = +_Rnm[_off+0]*(_Smph)[(_m)+0] ;	\
  } while (0)

/*
 * partial derivatives of Rnm, for m == 0
 *
 * where appropriate, derivatives are multiplied by two to account for 
 * symmetry in Fourier coefficients
 */

#define Rnm_derivatives_1m0(_n,_m,_rnm1,_Pnm1,_Cmph,_Smph,_dRnm)	\
  do {									\
    WBFMM_REAL _anm[2], _Rnm[2] ;					\
    gint _off = 0 ;							\
    									\
    g_assert((_m) == 0) ;						\
    									\
    (_dRnm)[WBFMM_DERIVATIVE_X_I] = 0.0 ;				\
    (_dRnm)[WBFMM_DERIVATIVE_Y_I] = 0.0 ;				\
    (_dRnm)[WBFMM_DERIVATIVE_Z_I] = 0.0 ;				\
    _anm[_off+0] = SQRT((WBFMM_REAL)(2*(_n)+1)/(2*(_n)-1)*		\
			((_n)-(_m))*((_n)+(_m))) ;			\
    _anm[_off+1] = SQRT((WBFMM_REAL)(2*(_n)+1)/(2*(_n)-1)*		\
			((_n)-(_m))*((_n)-(_m)-1)) ;			\
    _Rnm[_off+0] = (_rnm1)*(_Pnm1)[(_m)]*_anm[_off+0] ;			\
    _Rnm[_off+1] = (_rnm1)*(_Pnm1)[(_m)+1]*_anm[_off+1] ;		\
    (_dRnm)[WBFMM_DERIVATIVE_X_R] = -_Rnm[_off+1]*(_Cmph)[(_m)+1] ;	\
    (_dRnm)[WBFMM_DERIVATIVE_Y_R] = -_Rnm[_off+1]*(_Smph)[(_m)+1] ;	\
    (_dRnm)[WBFMM_DERIVATIVE_Z_R] = _Rnm[_off+0] ;			\
  } while (0)

#endif /*WBFMM_DERIVATIVES_H_INCLUDED*/
