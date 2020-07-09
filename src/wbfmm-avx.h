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

#ifndef WBFMM_AVX_H_INCLUDED
#define WBFMM_AVX_H_INCLUDED

#ifdef HAVE_AVX_INSTRUCTIONS

#ifdef HAVE_FMA_INSTRUCTIONS
#define wbfmm_cos_sin_recursion_avx(_Cn,_Sn,_C,_S)	\
  do { __m256d _rtmp = (_Cn), _op ;			\
    (_Cn) = _mm256_mul_pd((_Cn), (_C)) ;		\
    (_op) = _mm256_mul_pd((_Sn), (_S)) ;		\
    (_Cn) = _mm256_sub_pd((_Cn), _op) ;			\
    (_Sn) = _mm256_mul_pd((_Sn), (_C)) ;		\
    (_Sn) = _mm256_fmadd_pd(_rtmp, (_S), (_Sn)) ;	\
  } while (0)
#else  /*HAVE_FMA_INSTRUCTIONS*/
#define wbfmm_cos_sin_recursion_avx(_Cn,_Sn,_C,_S)	\
  do { __m256d _rtmp = (_Cn), _op ;			\
  _Cn = _mm256_mul_pd((_Cn), (_C)) ;			\
  _op = _mm256_mul_pd((_Sn), (_S)) ;			\
  (_Cn) = _mm256_sub_pd(_Cn, _op) ;			\
  _Sn = _mm256_mul_pd((_Sn), (_C)) ;			\
  _op = _mm256_mul_pd((_rtmp), (_S)) ;			\
  (_Sn) = _mm256_add_pd(_Sn, _op) ;			\
  } while (0)
#endif  /*HAVE_FMA_INSTRUCTIONS*/

#endif /*HAVE_AVX_INSTRUCTIONS*/

#endif /*WBFMM_AVX_H_INCLUDED*/
