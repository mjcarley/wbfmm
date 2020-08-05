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

#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include <glib.h>

#include <wbfmm.h>

#define BUFSIZE 131072

GTimer *timer ;

gchar *tests[] = {"expansion",
		  "translation",
		  "rotation",
		  "translationSR",
		  "translationRR",
		  "shift",
		  "child_parent",
		  "parent_child",
		  "expansion_gradient",
		  "local_gradient",
		  "expansion_normal",
		  ""} ;

#define wbfmm_index_laplace_nm(_n,_m) ((_n)*(_n)+(2*(_m))-1)

gint expansion_test(gint N, gdouble *x0, gdouble *xs,
		    gint xstride, gdouble *src, gint sstride,
		    gint nsrc, gdouble *xf, gint nfld) ;
gint expansion_normal_test(gint N, gdouble *x0, gdouble *xs,
			   gint xstride, gdouble *src, gint sstride,
			   gint nsrc, gdouble *xf, gint nfld) ;
gint expansion_gradient_test(gint N, gdouble *x0, gdouble *xs,
			     gint xstride, gdouble *src, gint sstride,
			     gint nsrc, gdouble *xf, gint nfld) ;
gint translation_test(gint N, gdouble *x0, gdouble *xs,
		      gint xstride, gdouble *src, gint sstride,
		      gint nsrc, gdouble *xf, gint nfld, gdouble t) ;
gint shift_test(gint N, gdouble *xc, gdouble *x0, gdouble *xs,
		gint xstride, gdouble *src, gint sstride,
		gint nsrc, gdouble *xf, gint nfld) ;
gint rotation_test(gint N, 
		   gdouble *x0,
		   gdouble *xs, gint xstride,
		   gdouble *src, gint sstride,
		   gint nsrc,
		   gdouble ix[], gdouble iy[], 
		   gdouble iz[],			  
		   gdouble *xf, gint nfld) ;
gint child_parent_test(gint N, gdouble *x0, gdouble wb, gint quad,
		       gdouble *xs, gint xstride,
		       gdouble *src, gint sstride,
		       gint nsrc, gdouble *xf, gint nfld) ;
gint parent_child_test(gint N, gdouble *x0, gdouble wb, gint quad,
		       gdouble *xs, gint xstride,
		       gdouble *src, gint sstride,
		       gint nsrc, gdouble *xp) ;
gint translation_RR_test(gint N, gdouble *x0, gdouble *xs,
			 gint xstride, gdouble *src, gint sstride,
			 gint nsrc, gdouble *xf, gint nfld, gdouble t) ;
gint translation_SR_test(gint N, gdouble *x0, gdouble *xs,
			 gint xstride, gdouble *src, gint sstride,
			 gint nsrc, gdouble *xf, gint nfld, gdouble t) ;
gint local_gradient_test(gint N, gdouble *x0, gdouble *xs,
			 gint xstride, gdouble *src, gint sstride,
			 gint nsrc, gdouble *xf, gint nfld, gdouble t) ;

static gint parse_test(gchar *arg)

{
  gint i = 0 ;
  
  while ( strlen(tests[i]) != 0) {
    if ( !strcmp(tests[i], arg) ) return i ;
    i ++ ;
  }

  return -1 ;
}

static gint read_data(gchar *ipfile, gdouble *xc, gdouble *x0, 
		      gdouble **xs, gdouble **src, 
		      gint *nsrc, gint *xstride, gint *sstride,
		      gdouble **xf, gint *fstride, gint *nfld,
		      gdouble *ix,gdouble *iy, gdouble *iz)

{
  FILE *f ;
  gchar line[1024] ;
  gint i ;

  f = fopen(ipfile, "r") ;
  if ( f == NULL ) {
    fprintf(stderr, "cannot open %s\n", ipfile) ;
    exit (1) ;
  }

  fprintf(stderr, "reading %s\n", ipfile) ;
  
  while ( fscanf(f, "%[^\n]c", line) != EOF ) {
    if ( strncmp(line, "x0:", 3) == 0 ) {
      sscanf(&(line[3]), "%lg %lg %lg",
	     &(x0[0]), &(x0[1]), &(x0[2])) ;
    }
    if ( strncmp(line, "xc:", 3) == 0 ) {
      sscanf(&(line[3]), "%lg %lg %lg",
	     &(xc[0]), &(xc[1]), &(xc[2])) ;
    }

    if ( strncmp(line, "ix:", 3) == 0 ) {
      sscanf(&(line[3]), "%lg %lg %lg",
	     &(ix[0]), &(ix[1]), &(ix[2])) ;
    }

    if ( strncmp(line, "iy:", 3) == 0 ) {
      sscanf(&(line[3]), "%lg %lg %lg",
	     &(iy[0]), &(iy[1]), &(iy[2])) ;
    }

    if ( strncmp(line, "iz:", 3) == 0 ) {
      sscanf(&(line[3]), "%lg %lg %lg",
	     &(iz[0]), &(iz[1]), &(iz[2])) ;
    }

    if ( strncmp(line, "sources:", 8) == 0 ) {
      sscanf(&(line[8]), "%d", nsrc) ;
      *xstride = 3 ; *sstride = 2 ;
      *xs = (gdouble *)g_malloc((*nsrc)*(*xstride)*sizeof(gdouble)) ;
      *src = (gdouble *)g_malloc((*nsrc)*(*sstride)*sizeof(gdouble)) ;
      for ( i = 0 ; i < *nsrc ; i ++ ) {
	fscanf(f, "%lg %lg %lg %lg %lg",
	       &((*xs)[(*xstride)*i+0]),
	       &((*xs)[(*xstride)*i+1]),
	       &((*xs)[(*xstride)*i+2]),
	       &((*src)[(*sstride)*i+0]),
	       &((*src)[(*sstride)*i+1])) ;
      }
    }
    if ( strncmp(line, "field:", 6) == 0 ) {
      sscanf(&(line[6]), "%d", nfld) ;
      *fstride = 5 ;
      *xf = (gdouble *)g_malloc((*nfld)*(*fstride)*sizeof(gdouble)) ;
      for ( i = 0 ; i < *nfld ; i ++ ) {
	fscanf(f, "%lg %lg %lg",
	       &((*xf)[(*fstride)*i+0]),
	       &((*xf)[(*fstride)*i+1]),
	       &((*xf)[(*fstride)*i+2])) ;
      }
    }

    fscanf(f, "%*c") ;
  }

  fprintf(stderr, "%s read\n", ipfile) ;

  fclose(f) ;

  return 0 ;
}

gint expansion_test(gint N, gdouble *x0, gdouble *xs,
		    gint xstride, gdouble *src, gint sstride,
		    gint nsrc, gdouble *xf, gint nfld)

{
  gint i, nq, cstr ;
  gdouble cfft[BUFSIZE]={0.0}, work[8192]={0.0}, fc[8]={0.0}, ff[8]={0.0} ;
  gdouble eval[BUFSIZE] = {0.0}, fe[8]={0.0} ;
  guint field ;
  
  nq = 2 ;
  fprintf(stderr, "expansion test\n") ;
  fprintf(stderr, "==============\n") ;
  fprintf(stderr,
	  "N = %d\n"
	  "x0 = (%lg, %lg, %lg)\n"
	  "xf = (%lg, %lg, %lg)\n"
	  "nsrc = %d\n",
	  N, x0[0], x0[1], x0[2], xf[0], xf[1], xf[2], nsrc) ;
  field = WBFMM_FIELD_SCALAR ;
  
  /*reference calculation*/
  wbfmm_laplace_field(xs, xstride, src, sstride, nq, NULL, 0, NULL, 0,
			   nsrc, xf, fc) ;

  cstr = 2 ;
  /*multipole expansion*/
  for ( i = 0 ; i < nsrc ; i ++ ) 
    wbfmm_laplace_expansion_cfft(N, x0, &(xs[i*xstride]),
				      &(src[i*sstride]), nq, cfft, cstr,
				      work) ;

  wbfmm_laplace_expansion_evaluate(x0, cfft, cstr, N, nq, xf, ff, work) ;

  /*check pre-computed evaluation method*/
  xf[0] -= x0[0] ; xf[1] -= x0[1] ; xf[2] -= x0[2] ; 
  wbfmm_laplace_field_coefficients(xf, N, field, eval, work) ;
  wbfmm_laplace_expansion_apply(cfft, cstr, nq, eval, N, field, fe, 1) ;
  
  fprintf(stderr, "exact:      ") ;
  for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", fc[i]) ;
  fprintf(stderr, "\n") ;

  fprintf(stderr, "field:      ") ;
  for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", ff[i]) ;
  fprintf(stderr, "\n") ;

  fprintf(stderr, "expansion:  ") ;
  for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", fe[i]) ;
  fprintf(stderr, "\n") ;

  fprintf(stderr, "errors:\n") ;  
  fprintf(stderr, "field:      ") ;
  for ( i = 0 ; i < nq ; i ++ )
    fprintf(stderr, "%1.16e ", fabs(fc[i]-ff[i])) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "expansion:  ") ;
  for ( i = 0 ; i < nq ; i ++ )
    fprintf(stderr, "%1.16e ", fabs(fc[i]-fe[i])) ;
  fprintf(stderr, "\n") ;
  
  return 0 ;
}

gint expansion_normal_test(gint N, gdouble *x0, gdouble *xs,
			   gint xstride, gdouble *src, gint sstride,
			   gint nsrc, gdouble *xf, gint nfld)

{
  gint i, nq, cstr ;
  gdouble cfft[BUFSIZE]={0.0}, work[8192]={0.0}, fc[8]={0.0}, ff[8]={0.0} ;
  gdouble eval[BUFSIZE] = {0.0}, fe[8]={0.0} ;
  gdouble n[3] ;
  guint field ;
  
  nq = 2 ;

  n[0] = 0.3 ; n[1] = -0.2 ; n[2] = 0.6 ;
  /* n[0] = 0.0 ; n[1] = 1.0 ; n[2] = 0.0 ; */
  
  fprintf(stderr, "expansion normal test\n") ;
  fprintf(stderr, "=====================\n") ;
  fprintf(stderr,
	  "N = %d\n"
	  "x0 = (%lg, %lg, %lg)\n"
	  "n  = (%lg, %lg, %lg)\n"	  
	  "xf = (%lg, %lg, %lg)\n"
	  "nsrc = %d\n",
	  N,
	  x0[0], x0[1], x0[2],
	  n[0],  n[1],  n[2], 
	  xf[0], xf[1], xf[2], nsrc) ;
  field = WBFMM_FIELD_SCALAR ;
  
  /*reference calculation*/
  wbfmm_laplace_field(xs, xstride, NULL, sstride, nq, n, 3, src, nq,
			   nsrc, xf, fc) ;

  cstr = 2 ;
  /*multipole expansion*/
  for ( i = 0 ; i < nsrc ; i ++ ) 
    wbfmm_laplace_expansion_normal_cfft(N, x0, &(xs[i*xstride]),
					     n, &(src[i*sstride]), nq,
					     cfft, cstr,
					     work) ;

  wbfmm_laplace_expansion_evaluate(x0, cfft, cstr, N, nq, xf, ff, work) ;

  /*check pre-computed evaluation method*/
  /* xf[0] -= x0[0] ; xf[1] -= x0[1] ; xf[2] -= x0[2] ;  */
  /* wbfmm_laplace_field_coefficients(xf, N, field, eval, work) ; */
  /* wbfmm_laplace_expansion_apply(cfft, cstr, nq, eval, N, field, fe, 1) ; */
  
  fprintf(stderr, "exact:      ") ;
  for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", fc[i]) ;
  fprintf(stderr, "\n") ;

  fprintf(stderr, "field:      ") ;
  for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", ff[i]) ;
  fprintf(stderr, "\n") ;

  /* fprintf(stderr, "expansion:  ") ; */
  /* for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", fe[i]) ; */
  /* fprintf(stderr, "\n") ; */

  fprintf(stderr, "errors:\n") ;  
  fprintf(stderr, "field:      ") ;
  for ( i = 0 ; i < nq ; i ++ )
    fprintf(stderr, "%1.16e ", fabs(fc[i]-ff[i])) ;
  fprintf(stderr, "\n") ;
  /* fprintf(stderr, "expansion:  ") ; */
  /* for ( i = 0 ; i < nq ; i ++ ) */
  /*   fprintf(stderr, "%1.16e ", fabs(fc[i]-fe[i])) ; */
  fprintf(stderr, "\n") ;
  
  return 0 ;
}

gint expansion_gradient_test(gint N, gdouble *x0, gdouble *xs,
			     gint xstride, gdouble *src, gint sstride,
			     gint nsrc, gdouble *xf, gint nfld)

{
  gint i, nq, cstr, fstr ;
  gdouble cfft[BUFSIZE]={0.0}, work[8192]={0.0} ;
  gdouble fc[32]={0.0}, ff[32]={0.0} ;

  nq = 2 ; fstr = 4 ;
  fprintf(stderr, "expansion gradient test\n") ;
  fprintf(stderr, "=======================\n") ;
  fprintf(stderr,
	  "N = %d\n"
	  "x0 = (%lg, %lg, %lg)\n"
	  "xf = (%lg, %lg, %lg)\n"
	  "nsrc = %d\n",
	  N, x0[0], x0[1], x0[2], xf[0], xf[1], xf[2], nsrc) ;

  /*reference calculation*/
  wbfmm_laplace_field_grad(xs, xstride, src, sstride, nq,
				NULL, 0, NULL, 0, nsrc, xf, fc, fstr) ;

  cstr = 3 ;
  /*multipole expansion*/
  for ( i = 0 ; i < nsrc ; i ++ ) 
    wbfmm_laplace_expansion_cfft(N, x0, &(xs[i*xstride]),
				      &(src[i*sstride]), nq, cfft, cstr,
				      work) ;

  wbfmm_laplace_expansion_grad_evaluate(x0, cfft, cstr, N, nq, xf,
					     ff, fstr, work) ;

  /*check pre-computed evaluation method*/
  /* xf[0] -= x0[0] ; xf[1] -= x0[1] ; xf[2] -= x0[2] ;  */
  /* wbfmm_laplace_field_coefficients(xf, N, FALSE, eval, work) ; */
  /* wbfmm_laplace_expansion_apply(cfft, cstr, nq, eval, N, fe) ; */
  
  fprintf(stderr, "exact:     ") ;
  for ( i = 0 ; i < nq ; i ++ )
    fprintf(stderr, "%lg %lg %lg, ",
	    fc[fstr*i+0], fc[fstr*i+1], fc[fstr*i+2]) ;
  fprintf(stderr, "\n") ;

  fprintf(stderr, "expansion: ") ;
  for ( i = 0 ; i < nq ; i ++ )
    fprintf(stderr, "%lg %lg %lg, ",
	    ff[fstr*i+0], ff[fstr*i+1], ff[fstr*i+2]) ;
  fprintf(stderr, "\n") ;

  fprintf(stderr, "error: ") ;
  for ( i = 0 ; i < nq ; i ++ )
    fprintf(stderr, "%lg %lg %lg, ",
	    fabs(ff[fstr*i+0]-fc[fstr*i+0]),
	    fabs(ff[fstr*i+1]-fc[fstr*i+1]),
	    fabs(ff[fstr*i+2]-fc[fstr*i+2])) ;
  fprintf(stderr, "\n") ;

  return 0 ;
}

gint translation_test(gint N, gdouble *x0, gdouble *xs,
		      gint xstride, gdouble *src, gint sstride,
		      gint nsrc, gdouble *xf, gint nfld, gdouble t)

{
  gint i, nq, cstr ;
  gdouble Ci[BUFSIZE]={0.0}, work[8192]={0.0}, Co[BUFSIZE]={0.0} ;
  gdouble fc[8]={0.0}, ff[8]={0.0}, ft[8]={0.0} ;

  nq = 2 ;
  fprintf(stderr, "translation test\n") ;
  fprintf(stderr, "==============\n") ;
  fprintf(stderr,
	  "N = %d\n"
	  "x0 = (%lg, %lg, %lg)\n"
	  "xf = (%lg, %lg, %lg)\n"
	  "t  = %lg\n"
	  "nsrc = %d\n",
	  N, x0[0], x0[1], x0[2], xf[0], xf[1], xf[2], t, nsrc) ;

  /*reference calculation*/
  wbfmm_laplace_field(xs, xstride, src, sstride, nq, NULL, 0, NULL, 0,
			   nsrc, xf, fc) ;

  cstr = 4 ;
  /*multipole expansion*/
  for ( i = 0 ; i < nsrc ; i ++ ) 
    wbfmm_laplace_expansion_cfft(N, x0, &(xs[i*xstride]),
				      &(src[i*sstride]), nq, Ci, cstr,
				      work) ;

  wbfmm_laplace_expansion_evaluate(x0, Ci, cstr, N, nq, xf, ff, work) ;

  /*translate the expansion*/
  for ( i = 0 ; i < 1024 ; i ++ ) Co[i] = -13.0 ;
  wbfmm_laplace_coaxial_translate_SS(Co, cstr, N, Ci, cstr, N, nq, t,
					  0.0) ;

  x0[2] += t ;  
  wbfmm_laplace_expansion_evaluate(x0, Co, cstr, N, nq, xf, ft, work) ;
  
  for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", fc[i]) ;
  fprintf(stderr, "\n") ;
  for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", ff[i]) ;
  fprintf(stderr, "\n") ;
  for ( i = 0 ; i < nq ; i ++ )
    fprintf(stderr, "%1.16e ", fabs(fc[i]-ff[i])) ;
  fprintf(stderr, "\n") ;
  for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", ft[i]) ;
  fprintf(stderr, "\n") ;
  for ( i = 0 ; i < nq ; i ++ )
    fprintf(stderr, "%1.16e ", fabs(fc[i]-ft[i])) ;
  fprintf(stderr, "\n") ;
  
  return 0 ;
}

gint translation_SR_test(gint N, gdouble *x0, gdouble *xs,
			 gint xstride, gdouble *src, gint sstride,
			 gint nsrc, gdouble *xf, gint nfld, gdouble t)

{
  gint i, nq, cstr, Ns, Nr ;
  gdouble Ci[BUFSIZE]={0.0}, work[8192]={0.0}, Co[BUFSIZE]={0.0} ;
  gdouble fc[8]={0.0}, ff[8]={0.0}, ft[8]={0.0} ;

  nq = 2 ;
  Ns = N ; Nr = Ns + 4 ;

  xf[0] = x0[0] + 0.01 ; 
  xf[1] = x0[1] - 0.01 ; 
  xf[2] = x0[2] + t + 0.01 ; 
  
  fprintf(stderr, "singular to regular translation test\n") ;
  fprintf(stderr, "====================================\n") ;
  fprintf(stderr,
	  "Ns = %d\n"
	  "Nr = %d\n"
	  "x0 = (%lg, %lg, %lg)\n"
	  "xf = (%lg, %lg, %lg)\n"
	  "t  = %lg\n"
	  "nsrc = %d\n",
	  Ns, Nr, x0[0], x0[1], x0[2], xf[0], xf[1], xf[2], t, nsrc) ;

  /*reference calculation*/
  wbfmm_laplace_field(xs, xstride, src, sstride, nq, NULL, 0, NULL, 0,
			   nsrc, xf, fc) ;

  cstr = 4 ;
  /*multipole expansion*/
  for ( i = 0 ; i < nsrc ; i ++ ) 
    wbfmm_laplace_expansion_cfft(Ns, x0, &(xs[i*xstride]),
				      &(src[i*sstride]), nq, Ci, cstr,
				      work) ;

  wbfmm_laplace_expansion_evaluate(x0, Ci, cstr, Ns, nq, xf, ff, work) ;

  /*translate the expansion*/
  for ( i = 0 ; i < 1024 ; i ++ ) Co[i] = -13.0 ;
  wbfmm_laplace_coaxial_translate_SR(Co, cstr, Nr, Ci, cstr, Ns, nq, t,
					  0.0) ;

  x0[2] += t ;  
  wbfmm_laplace_expansion_local_evaluate(x0, Co, cstr, Nr, nq,
					      xf, ft, work) ;
  
  for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", fc[i]) ;
  fprintf(stderr, "\n") ;
  for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", ff[i]) ;
  fprintf(stderr, "\n") ;
  for ( i = 0 ; i < nq ; i ++ )
    fprintf(stderr, "%1.16e ", fabs(fc[i]-ff[i])) ;
  fprintf(stderr, "\n") ;
  for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", ft[i]) ;
  fprintf(stderr, "\n") ;
  for ( i = 0 ; i < nq ; i ++ )
    fprintf(stderr, "%1.16e ", fabs(fc[i]-ft[i])) ;
  fprintf(stderr, "\n") ;
  
  return 0 ;
}

gint translation_RR_test(gint N, gdouble *x0, gdouble *xs,
			 gint xstride, gdouble *src, gint sstride,
			 gint nsrc, gdouble *xf, gint nfld, gdouble t)

{
  gint i, nq, cstr, N0, N1 ;
  gdouble Ci[BUFSIZE]={0.0}, Co[BUFSIZE]={0.0}, Cr[BUFSIZE] = {0.0} ;
  gdouble work[BUFSIZE]={0.0}, eval[BUFSIZE] = {0.0} ;
  gdouble fc[8]={0.0}, ff[8]={0.0}, ft[8]={0.0}, fe[8]={0.0} ;
  guint field ;
  
  nq = 2 ;
  N0 = N+4 ; N1 = N ;
  xf[0] = x0[0] + 0.01 ; 
  xf[1] = x0[1] + 0.01 ; 
  xf[2] = x0[2] + t - 0.05 ; 

  field = WBFMM_FIELD_SCALAR ;
  
  fprintf(stderr, "regular to regular translation test\n") ;
  fprintf(stderr, "====================================\n") ;
  fprintf(stderr,
	  "N0 = %d\n"
	  "N1 = %d\n"
	  "x0 = (%lg, %lg, %lg)\n"
	  "xf = (%lg, %lg, %lg)\n"
	  "t  = %lg\n"
	  "nsrc = %d\n",
	  N0, N1, x0[0], x0[1], x0[2], xf[0], xf[1], xf[2], t, nsrc) ;

  /*reference calculation*/
  wbfmm_laplace_field(xs, xstride, src, sstride, nq, NULL, 0, NULL, 0,
			   nsrc, xf, fc) ;

  cstr = 4 ;
  /*multipole expansion*/
  for ( i = 0 ; i < nsrc ; i ++ ) 
    wbfmm_laplace_expansion_cfft(N0, x0, &(xs[i*xstride]),
				      &(src[i*sstride]), nq, Ci, cstr,
				      work) ;

  wbfmm_laplace_expansion_evaluate(x0, Ci, cstr, N0, nq, xf, ff, work) ;

  /*translate the expansion, singular to regular*/
  for ( i = 0 ; i < 1024 ; i ++ ) Co[i] = -13.0 ;
  wbfmm_laplace_coaxial_translate_SR(Co, cstr, N0, Ci, cstr, N0, nq,
					  t*0.9, 0.0) ;
  /*translate the expansion, regular to regular*/
  for ( i = 0 ; i < 1024 ; i ++ ) Cr[i] = -13.0 ;
  wbfmm_laplace_coaxial_translate_RR(Cr, cstr, N1, Co, cstr, N0, nq,
					  t*0.1, 0.0) ;  

  x0[2] += t ;  
  wbfmm_laplace_expansion_local_evaluate(x0, Cr, cstr, N1, nq,
					      xf, ft, work) ;

  /*check pre-computed evaluation method*/
  xf[0] -= x0[0] ; xf[1] -= x0[1] ; xf[2] -= x0[2] ; 
  wbfmm_laplace_local_coefficients(xf, N1, field, eval, work) ;
  wbfmm_laplace_expansion_apply(Cr, cstr, nq, eval, N1, field, fe, 1) ;
					
  fprintf(stderr, "exact:      ") ;
  for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", fc[i]) ;
  fprintf(stderr, "\n") ;

  fprintf(stderr, "field:      ") ;
  for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", ff[i]) ;
  fprintf(stderr, "\n") ;

  fprintf(stderr, "translated: ") ;
  for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", ft[i]) ;
  fprintf(stderr, "\n") ;

  fprintf(stderr, "expansion:  ") ;
  for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", fe[i]) ;
  fprintf(stderr, "\n") ;

  fprintf(stderr, "errors:\n") ;
  fprintf(stderr, "field:      ") ;
  for ( i = 0 ; i < nq ; i ++ )
    fprintf(stderr, "%1.16e ", fabs(fc[i]-ff[i])) ;
  fprintf(stderr, "\n") ;

  fprintf(stderr, "translated: ") ;
  for ( i = 0 ; i < nq ; i ++ )
    fprintf(stderr, "%1.16e ", fabs(fc[i]-ft[i])) ;
  fprintf(stderr, "\n") ;

  fprintf(stderr, "expansion:  ") ;
  for ( i = 0 ; i < nq ; i ++ )
    fprintf(stderr, "%1.16e ", fabs(fc[i]-fe[i])) ;
  fprintf(stderr, "\n") ;


  
  return 0 ;
}

gint local_gradient_test(gint N, gdouble *x0, gdouble *xs,
			 gint xstride, gdouble *src, gint sstride,
			 gint nsrc, gdouble *xf, gint nfld, gdouble t)

{
  gint i, nq, cstr, Ns, Nr, fstr ;
  gdouble Ci[BUFSIZE]={0.0}, work[8192]={0.0}, Co[BUFSIZE]={0.0} ;
  gdouble eval[BUFSIZE] = {0.0} ;
  gdouble fc[32]={0.0}, ff[32]={0.0}, ft[32]={0.0}, fe[32] = {0.0} ;
  guint field ;
  
  nq = 1 ; fstr = 4 ;
  Ns = N ; Nr = Ns ;

  field = WBFMM_FIELD_GRADIENT ;
  
  xf[0] = x0[0] + 0.35 ; 
  xf[1] = x0[1] - 0.31 ; 
  xf[2] = x0[2] + t + 0.5 ; 
  
  fprintf(stderr, "singular to regular translation gradient test\n") ;
  fprintf(stderr, "=============================================\n") ;
  fprintf(stderr,
	  "Ns = %d\n"
	  "Nr = %d\n"
	  "x0 = (%lg, %lg, %lg)\n"
	  "xf = (%lg, %lg, %lg)\n"
	  "t  = %lg\n"
	  "nsrc = %d\n",
	  Ns, Nr, x0[0], x0[1], x0[2], xf[0], xf[1], xf[2], t, nsrc) ;

  /*reference calculation*/
  wbfmm_laplace_field_grad(xs, xstride, src, sstride, nq, NULL, 0, NULL, 0,
				nsrc, xf, fc, fstr) ;

  cstr = 4 ;
  /*multipole expansion*/
  for ( i = 0 ; i < nsrc ; i ++ ) 
    wbfmm_laplace_expansion_cfft(Ns, x0, &(xs[i*xstride]),
				      &(src[i*sstride]), nq, Ci, cstr,
				      work) ;

  wbfmm_laplace_expansion_grad_evaluate(x0, Ci, cstr, Ns, nq, xf,
					     ff, fstr, work) ;

  /*translate the expansion*/
  wbfmm_laplace_coaxial_translate_SR(Co, cstr, Nr, Ci, cstr, Ns, nq, t,
					  0.0) ;

  x0[2] += t ;
  wbfmm_laplace_expansion_local_grad_evaluate(x0, Co, cstr, Nr, nq,
  						   xf, ft, fstr, work) ;

  xf[0] -= x0[0] ; xf[1] -= x0[1] ; xf[2] -= x0[2] ; 
  wbfmm_laplace_local_coefficients(xf, N, field, eval, work) ;
  wbfmm_laplace_expansion_apply(Co, cstr, nq, eval, N, field, fe, fstr) ;

  fprintf(stderr, "exact:       ") ;
  for ( i = 0 ; i < nq ; i ++ )
    fprintf(stderr, "%lg %lg %lg, ",
	    fc[fstr*i+0], fc[fstr*i+1], fc[fstr*i+2]) ;
  fprintf(stderr, "\n") ;

  fprintf(stderr, "expansion:   ") ;
  for ( i = 0 ; i < nq ; i ++ )
    fprintf(stderr, "%lg %lg %lg, ",
	    ff[fstr*i+0], ff[fstr*i+1], ff[fstr*i+2]) ;
  fprintf(stderr, "\n") ;

  fprintf(stderr, "error: ") ;
  for ( i = 0 ; i < nq ; i ++ )
    fprintf(stderr, "%lg %lg %lg, ",
	    fabs(ff[fstr*i+0]-fc[fstr*i+0]),
	    fabs(ff[fstr*i+1]-fc[fstr*i+1]),
	    fabs(ff[fstr*i+2]-fc[fstr*i+2])) ;
  fprintf(stderr, "\n") ;

  fprintf(stderr, "precomputed: ") ;
  for ( i = 0 ; i < nq ; i ++ )
    fprintf(stderr, "%lg %lg %lg, ",
	    fe[fstr*i+0], fe[fstr*i+1], fe[fstr*i+2]) ;
  fprintf(stderr, "\n") ;
  
  fprintf(stderr, "translated:  ") ;
  for ( i = 0 ; i < nq ; i ++ )
    fprintf(stderr, "%lg %lg %lg, ",
	    ft[fstr*i+0], ft[fstr*i+1], ft[fstr*i+2]) ;
  fprintf(stderr, "\n") ;

  fprintf(stderr, "error: ") ;
  for ( i = 0 ; i < nq ; i ++ )
    fprintf(stderr, "%lg %lg %lg, ",
	    fabs(ft[fstr*i+0]-fc[fstr*i+0]),
	    fabs(ft[fstr*i+1]-fc[fstr*i+1]),
	    fabs(ft[fstr*i+2]-fc[fstr*i+2])) ;
  fprintf(stderr, "\n") ;

  
  return 0 ;
}

gint rotation_test(gint N, 
		   gdouble *x0,
		   gdouble *xs, gint xstride,
		   gdouble *src, gint sstride,
		   gint nsrc,
		   gdouble ix[], gdouble iy[], 
		   gdouble iz[],			  
		   gdouble *xf, gint nfld)

{
  gdouble H[BUFSIZE], work[BUFSIZE], th, ph, ch ;
  gdouble ix0[3], iy0[3], iz0[3], y[3] ;
  gdouble Ci[BUFSIZE] = {0.0}, Co[BUFSIZE] = {0.0}, Cc[BUFSIZE] = {0.0} ;
  gdouble fc[32]={0.0}, ff[32]={0.0}, fr[32]={0.0} ;
  gint i, cstri, cstro, nq ;
  gdouble t0 ;
  
  cstri = 4 ; cstro = 4 ; nq = 2 ;

  ix0[0] = 1.0 ; ix0[1] = 0.0 ; ix0[2] = 0.0 ;
  iy0[0] = 0.0 ; iy0[1] = 1.0 ; iy0[2] = 0.0 ;
  iz0[0] = 0.0 ; iz0[1] = 0.0 ; iz0[2] = 1.0 ;

  fprintf(stderr, "rotation test\n") ;
  fprintf(stderr, "==============\n") ;
  fprintf(stderr,
	  "N = %d\n"
	  "x0 = (%lg, %lg, %lg)\n"
	  "xf = (%lg, %lg, %lg)\n"
	  "nsrc = %d\n"
	  "nq = %d\n",
	  N, x0[0], x0[1], x0[2], xf[0], xf[1], xf[2], nsrc, nq) ;

  wbfmm_rotation_angles(ix0, iy0, iz0, ix, iy, iz, &th, &ph, &ch) ;

  fprintf(stderr, "rotation: (%lg,%lg,%lg)\n", th, ph, ch) ;

  /*reference calculation*/
  wbfmm_laplace_field(xs, xstride, src, sstride, nq, NULL, 0, NULL, 0,
			   nsrc, xf, fc) ;

  /*multipole expansion*/
  for ( i = 0 ; i < nsrc ; i ++ ) 
    wbfmm_laplace_expansion_cfft(N, x0, &(xs[i*xstride]),
				      &(src[i*sstride]), nq, Ci, cstri,
				      work) ;

  /*field from unrotated coefficients*/
  wbfmm_laplace_expansion_evaluate(x0, Ci, cstri, N, nq, xf, ff, work) ;

  wbfmm_coefficients_H_rotation(H, N, th, work) ;

  /*fill Co with rubbish to check pre-scaling*/
  memset(Co, 1, BUFSIZE*sizeof(gdouble)) ;
  for ( i = 0 ; i < 1024 ; i ++ ) Co[i] = -13.0 ;
  /*apply the rotation to the coefficients*/
  t0 = g_timer_elapsed(timer, NULL) ;
  wbfmm_laplace_rotate_H_ref(Co, cstro, Ci, cstri, N, nq, H, ph, ch, 0.0) ;
  t0 = g_timer_elapsed(timer, NULL) - t0 ;
  fprintf(stderr, "rotation: time %lg\n", t0) ;

#ifdef WBFMM_USE_AVX
  memset(Co, 1, BUFSIZE*sizeof(gdouble)) ;
  for ( i = 0 ; i < 1024 ; i ++ ) Co[i] = -13.0 ;
  /*apply the rotation to the coefficients*/
  t0 = g_timer_elapsed(timer, NULL) ;
  wbfmm_laplace_rotate_H(Co, cstro, Ci, cstri, N, nq, H, ph, ch, 0.0) ;
  t0 = g_timer_elapsed(timer, NULL) - t0 ;
  fprintf(stderr, "avx rotation: time %lg\n", t0) ;
#endif
  
  memset(Cc, 0, BUFSIZE*sizeof(gdouble)) ;
  /*reverse rotation as check*/
  /* wbfmm_laplace_rotate_H(Cc, cstro, Co, cstro, N, nq, H, ch, ph) ; */

  /* for ( i = 0 ; i < 32 ; i ++ ) */
  /*   fprintf(stderr, "%lg %lg %1.16e\n", */
  /* 	    Ci[cstri*i], Cc[cstro*i], fabs(Ci[cstri*i]-Cc[cstro*i])) ; */
  
  /*from rotated coefficients*/
  wbfmm_coordinate_transform(xf, ix, iy, iz, y) ;
  wbfmm_laplace_expansion_evaluate(x0, Co, cstro, N, nq, y, fr, work) ;

  for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", fc[i]) ;
  fprintf(stderr, "\n") ;
  for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", ff[i]) ;
  fprintf(stderr, "\n") ;
  for ( i = 0 ; i < nq ; i ++ )
    fprintf(stderr, "%1.16e ", fabs(fc[i]-ff[i])) ;
  fprintf(stderr, "\n") ;
  for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", fr[i]) ;
  fprintf(stderr, "\n") ;
  for ( i = 0 ; i < nq ; i ++ )
    fprintf(stderr, "%1.16e ", fabs(fc[i]-fr[i])) ;
  fprintf(stderr, "\n") ;
  
  return 0 ;
}

gint shift_test(gint N, gdouble *xc, gdouble *x0, gdouble *xs,
		gint xstride, gdouble *src, gint sstride,
		gint nsrc, gdouble *xf, gint nfld)

{
  gint i, nq, cstri, cstro, Ni, No ;
  gdouble Ci[BUFSIZE]={0.0}, Co[BUFSIZE]={0.0} ;
  gdouble Cr1[BUFSIZE] = {0.0}, Cr2[BUFSIZE] = {0.0} ;
  gdouble H1[BUFSIZE] = {0.0} ;
  gdouble work[BUFSIZE]={0.0} ;
  gdouble r, th, ph, ch ;
  gdouble fc[8]={0.0}, ff[8]={0.0}, ft[8]={0.0} ;

  nq = 1 ;
  fprintf(stderr, "shift test\n") ;
  fprintf(stderr, "==========\n") ;
  fprintf(stderr,
	  "N = %d\n"
	  "xc = (%lg, %lg, %lg)\n"
	  "x0 = (%lg, %lg, %lg)\n"
	  "xf = (%lg, %lg, %lg)\n"
	  "nsrc = %d\n",
	  N,
	  xc[0], xc[1], xc[2], x0[0], x0[1], x0[2], xf[0], xf[1], xf[2],
	  nsrc) ;

  Ni = N ; No = Ni + 5 ;
  /*reference calculation*/
  wbfmm_laplace_field(xs, xstride, src, sstride, nq, NULL, 0, NULL, 0,
			   nsrc, xf, fc) ;

  cstri = 1 ; cstro = 1 ;
  /*multipole expansion*/
  for ( i = 0 ; i < nsrc ; i ++ ) 
    wbfmm_laplace_expansion_cfft(Ni, x0, &(xs[i*xstride]),
				      &(src[i*sstride]), nq, Ci, cstri,
				      work) ;

  wbfmm_laplace_expansion_evaluate(x0, Ci, cstri, Ni, nq, xf, ff, work) ;

  /*shift coordinates*/
  wbfmm_shift_angles(x0, xc, &th, &ph, &ch, &r) ;

  fprintf(stderr, "rotation: (%lg, %lg, %lg)\n",
	  th, ph, ch) ;
  
  /*rotation operator*/
  wbfmm_coefficients_H_rotation(H1, MAX(Ni, No), th, work) ;

  /*rotate expansion*/
  wbfmm_laplace_rotate_H(Cr1, cstri, Ci, cstri, Ni, nq, H1, ph, ch, 0.0) ;
  
  /*translate the expansion, singular to singular*/
  wbfmm_laplace_coaxial_translate_SS(Cr2, cstro, No, Cr1, cstri, Ni,
					  nq, r, 0.0) ;

  /*reverse rotation*/
  wbfmm_laplace_rotate_H(Co, cstro, Cr2, cstro, No, nq, H1, ch, ph, 0.0) ;

  wbfmm_laplace_expansion_evaluate(xc, Co, cstro, No, nq,
					xf, ft, work) ;
  
  for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", fc[i]) ;
  fprintf(stderr, "\n") ;
  for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", ff[i]) ;
  fprintf(stderr, "\n") ;
  for ( i = 0 ; i < nq ; i ++ )
    fprintf(stderr, "%1.16e ", fabs(fc[i]-ff[i])) ;
  fprintf(stderr, "\n") ;
  for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", ft[i]) ;
  fprintf(stderr, "\n") ;
  for ( i = 0 ; i < nq ; i ++ )
    fprintf(stderr, "%1.16e ", fabs(fc[i]-ft[i])) ;
  fprintf(stderr, "\n") ;

  return 0 ;
}

gint child_parent_test(gint N, gdouble *x0, gdouble wb, gint quad,
		       gdouble *xs, gint xstride,
		       gdouble *src, gint sstride,
		       gint nsrc, gdouble *xf, gint nfld)

{
  gdouble *Cc, *Cp, *Cpbw, xc[3], xsc[3] ;
  gint i, nq, Nc, Np ;
  gdouble H03[BUFSIZE] = {0.0}, H47[BUFSIZE] = {0.0} ;
  gdouble work[BUFSIZE] = {0.0} ;
  gdouble cr[BUFSIZE] = {0.0}, ct[BUFSIZE] = {0.0}, crr[BUFSIZE] = {0.0} ;
  gdouble r, th, ph, ch ;
  gdouble fc[8] = {0.0}, fr[8] = {0.0}, fp[8] = {0.0}, fs[8] = {0.0},
    fb[8] = {0.0} ;
  gdouble th03, th47 ;
  gdouble t0, t1 ;
  
  nq = 2 ;
  Nc = N ; Np = N ;
  Cc   = (gdouble *)g_malloc0((Nc+1)*(Nc+1)*nq*8*sizeof(gdouble)) ;
  Cp   = (gdouble *)g_malloc0((Np+1)*(Np+1)*nq*8*sizeof(gdouble)) ;
  Cpbw = (gdouble *)g_malloc0((Np+1)*(Np+1)*nq*8*sizeof(gdouble)) ;

  fprintf(stderr, "child-parent shift test\n") ;
  fprintf(stderr, "=======================\n") ;
  fprintf(stderr,
	  "Nc = %d; Np = %d;\n"
	  "x0 = (%lg, %lg, %lg)\n"
	  "xf = (%lg, %lg, %lg)\n"
	  "xs = (%lg, %lg, %lg)\n"
	  "wb = %lg;\n"
	  "quad = %d;\n"
	  "nsrc = %d\n",
	  Nc, Np,
	  x0[0], x0[1], x0[2], xf[0], xf[1], xf[2],
	  xs[0], xs[1], xs[2],
	  wb, quad, nsrc) ;

  /*generate child box coefficients*/
  switch ( quad ) {
  default: g_assert_not_reached() ;
  case 0:
    xc[0] = x0[0] - wb/2 ; xc[1] = x0[1] - wb/2 ; xc[2] = x0[2] - wb/2 ;
    break ;
  case 1:
    xc[0] = x0[0] + wb/2 ; xc[1] = x0[1] - wb/2 ; xc[2] = x0[2] - wb/2 ;
    break ;
  case 2:
    xc[0] = x0[0] - wb/2 ; xc[1] = x0[1] + wb/2 ; xc[2] = x0[2] - wb/2 ;
    break ;
  case 3:
    xc[0] = x0[0] + wb/2 ; xc[1] = x0[1] + wb/2 ; xc[2] = x0[2] - wb/2 ;
    break ;
  case 4:
    xc[0] = x0[0] - wb/2 ; xc[1] = x0[1] - wb/2 ; xc[2] = x0[2] + wb/2 ;
    break ;
  case 5:
    xc[0] = x0[0] + wb/2 ; xc[1] = x0[1] - wb/2 ; xc[2] = x0[2] + wb/2 ;
    break ;
  case 6:
    xc[0] = x0[0] - wb/2 ; xc[1] = x0[1] + wb/2 ; xc[2] = x0[2] + wb/2 ;
    break ;
  case 7:
    xc[0] = x0[0] + wb/2 ; xc[1] = x0[1] + wb/2 ; xc[2] = x0[2] + wb/2 ;
    break ;
  }
  
  xsc[0] = xc[0] + xs[0] ; xsc[1] = xc[1] + xs[1] ; xsc[2] = xc[2] + xs[2] ; 
  /* xsc[0] = x0[0]; xsc[1] = x0[1] ; xsc[2] = x0[2] ; */
  /* xsc[0] = xc[0]; xsc[1] = xc[1] ; xsc[2] = xc[2] ; */

  wbfmm_laplace_expansion_cfft(Nc, xc, xsc, src, nq,
				    &(Cc[quad*nq]), nq*8, work) ;
  wbfmm_laplace_expansion_evaluate(xc, &(Cc[quad*nq]), nq*8,
					Nc, nq, xf, fc, work) ;
  wbfmm_laplace_field(xsc, xstride, src, sstride, nq, NULL, 0, NULL, 0,
			   nsrc, xf, fr) ;

  /*rotations and translations for child-parent shift*/
  th03 = acos(sqrt(1.0/3.0)) ; th47 = M_PI - th03 ;
  wbfmm_coefficients_H_rotation(H03, Np, th03, work) ;
  wbfmm_coefficients_H_rotation(H47, Np, th47, work) ;

  /* Np = Nc = 0 ; */
  fprintf(stderr, "child parent shift:\n") ;
  t0 = g_timer_elapsed(timer, NULL) ;
  for ( i = 0 ; i < 1024 ; i ++ ) Cp[i] = -13.0 ;
  for ( i = 0 ; i < 1024 ; i ++ ) work[8*nq*(Np+1)*(Np+1)] = -13.0 ;
  wbfmm_laplace_child_parent_shift(Cp, Np, Cc, Nc, nq, H03, H47, Np,
  					wb, work) ;
  t1 = g_timer_elapsed(timer, NULL) ;
  fprintf(stderr, "shift completed (%lg);\n", t1 - t0) ;
  memset(work, 0, BUFSIZE*sizeof(gdouble)) ;

  t0 = g_timer_elapsed(timer, NULL) ;  
  wbfmm_laplace_child_parent_shift_bw(Cpbw, Np, Cc, Nc, nq, H03, Np,
					   wb, work) ;
  t1 = g_timer_elapsed(timer, NULL) ;
  fprintf(stderr, "bw shift completed (%lg);\n", t1 - t0) ;
  memset(work, 0, BUFSIZE*sizeof(gdouble)) ;

  /*same thing using a sequence of reference functions*/
  wbfmm_shift_angles(xc, x0, &th, &ph, &ch, &r) ;
  /* ph = 0.5*M_PI ; */
  fprintf(stderr, "th = %lg; ph = %lg; ch = %lg;\n", th, ph, ch) ;
  wbfmm_coefficients_H_rotation(H03, Np, th, work) ;
  wbfmm_laplace_rotate_H(cr, nq, &(Cc[quad*nq]), 8*nq, Nc, nq, H03,
			      ph, ch, 0.0) ;

  wbfmm_laplace_coaxial_translate_SS(ct, nq, Np, cr, nq, Nc, nq, r, 0.0) ;

  wbfmm_laplace_rotate_H(crr, nq, ct, nq, Np, nq, H03, ch, ph, 0.0) ;

  wbfmm_laplace_expansion_evaluate(x0, Cp, nq*8, Np, nq, xf, fp, work) ;
  wbfmm_laplace_expansion_evaluate(x0, Cpbw, nq*8, Np, nq, xf, fb, work) ;
  wbfmm_laplace_expansion_evaluate(x0, crr, nq, Np, nq, xf, fs, work) ;

  fprintf(stderr, "ref:   ") ;
  for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", fr[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "child: ") ;
  for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", fc[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "shift: ") ;
  for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", fp[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "backw: ") ;
  for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", fb[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "chain: ") ;
  for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", fs[i]) ;
  fprintf(stderr, "\n") ;

  fprintf(stderr, "child: ") ;
  for ( i = 0 ; i < nq ; i ++ )
    fprintf(stderr, "%1.16e ", fabs(fc[i]-fr[i])) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "shift: ") ;
  for ( i = 0 ; i < nq ; i ++ )
    fprintf(stderr, "%1.16e ", fabs(fp[i]-fr[i])) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "backw: ") ;
  for ( i = 0 ; i < nq ; i ++ )
    fprintf(stderr, "%1.16e ", fabs(fb[i]-fr[i])) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "chain: ") ;
  for ( i = 0 ; i < nq ; i ++ )
    fprintf(stderr, "%1.16e ", fabs(fs[i]-fr[i])) ;
  fprintf(stderr, "\n") ;
  
  return 0 ;
}

gint parent_child_test(gint N, gdouble *x0, gdouble wb, gint quad,
		       gdouble *xs, gint xstride,
		       gdouble *src, gint sstride,
		       gint nsrc, gdouble *xp)

{
  gdouble *Cc, *Cp, *C0, *C1, xc[3], xe[3] ;
  gint i, nq, Nc, Np ;
  gdouble H03[BUFSIZE] = {0.0}, H47[BUFSIZE] = {0.0}, H0[BUFSIZE] = {0.0} ;
  gdouble work[BUFSIZE]={0.0} ;
  gdouble cr[BUFSIZE] = {0.0}, ct[BUFSIZE] = {0.0}, crr[BUFSIZE] = {0.0} ;
  gdouble r, th, ph, ch ;
  gdouble fc[8] = {0.0}, fr[8] = {0.0}, fp[8] = {0.0}, fs[8] = {0.0} ;
  gdouble th03, th47 ;
  /* gdouble *check ; */
  
  nq = 1 ;
  Nc = N ; Np = Nc+4 ;
  Cc = (gdouble *)g_malloc0((Nc+1)*(Nc+1)*nq*8*sizeof(gdouble)) ;
  Cp = (gdouble *)g_malloc0((Np+1)*(Np+1)*nq*8*sizeof(gdouble)) ;
  C0 = (gdouble *)g_malloc0((Np+1)*(Np+1)*nq*1*sizeof(gdouble)) ;
  C1 = (gdouble *)g_malloc0((Np+1)*(Np+1)*nq*1*sizeof(gdouble)) ;

  fprintf(stderr, "parent-child shift test\n") ;
  fprintf(stderr, "=======================\n") ;
  fprintf(stderr,
	  "Nc = %d; Np = %d;\n"
	  "x0 = (%lg, %lg, %lg)\n"
	  "xs = (%lg, %lg, %lg)\n"
	  "xp = (%lg, %lg, %lg)\n"
	  "wb = %lg;\n"
	  "nsrc = %d\n",
	  Nc, Np,
	  x0[0], x0[1], x0[2],
	  xs[0], xs[1], xs[2],
	  xp[0], xp[1], xp[2],
	  wb, nsrc) ;

  /*generate source expansion*/
  for ( i = 0 ; i < nsrc ; i ++ ) 
    wbfmm_laplace_expansion_cfft(N, x0, &(xs[i*xstride]),
				      &(src[i*sstride]), nq, C0, nq,
				      work) ;
  
  /*shift singular expansion to regular expansion about parent box*/
  fprintf(stderr, "setting shift angles\n") ;
  wbfmm_shift_angles(x0, xp, &th, &ph, &ch, &r) ;
  fprintf(stderr, "generating rotation coefficients\n") ;
  wbfmm_coefficients_H_rotation(H0, Np, th, work) ;
  fprintf(stderr, "performing rotation\n") ;
  wbfmm_laplace_rotate_H(C1, nq, C0, nq, Np, nq, H0, ph, ch, 0.0) ;
  memset(C0, 0, sizeof(gdouble)*(Np+1)*(Np+1)*nq) ;
  fprintf(stderr, "performing translation\n") ;
  wbfmm_laplace_coaxial_translate_SR(C0, nq, Np, C1, nq, Np, nq, r, 0.0) ;
  fprintf(stderr, "performing reverse rotation\n") ;
  wbfmm_laplace_rotate_H(Cp, 8*nq, C0, nq, Np, nq, H0, ch, ph, 0.0) ;

  /*perform RR shift from parent centre to child centres*/
  th47 = acos(sqrt(1.0/3.0)) ; th03 = M_PI - th47 ;
  wbfmm_coefficients_H_rotation(H03, Np, th03, work) ;
  wbfmm_coefficients_H_rotation(H47, Np, th47, work) ;
  fprintf(stderr, "parent child shift:\n") ;
  wbfmm_laplace_parent_child_shift(Cc, Nc, Cp, Np, nq, H03, H47, Np,
  					wb, work) ;

  /*evaluate field in child box, wb is the parent box width so child
    box centres lie at wb/4 from the parent box centre*/
  switch ( quad ) {
  default: g_assert_not_reached() ;
  case 0:
    xc[0] = xp[0] - wb/4 ; xc[1] = xp[1] - wb/4 ; xc[2] = xp[2] - wb/4 ;
    break ;
  case 1:
    xc[0] = xp[0] + wb/4 ; xc[1] = xp[1] - wb/4 ; xc[2] = xp[2] - wb/4 ;
    break ;
  case 2:
    xc[0] = xp[0] - wb/4 ; xc[1] = xp[1] + wb/4 ; xc[2] = xp[2] - wb/4 ;
    break ;
  case 3:
    xc[0] = xp[0] + wb/4 ; xc[1] = xp[1] + wb/4 ; xc[2] = xp[2] - wb/4 ;
    break ;
  case 4:
    xc[0] = xp[0] - wb/4 ; xc[1] = xp[1] - wb/4 ; xc[2] = xp[2] + wb/4 ;
    break ;
  case 5:
    xc[0] = xp[0] + wb/4 ; xc[1] = xp[1] - wb/4 ; xc[2] = xp[2] + wb/4 ;
    break ;
  case 6:
    xc[0] = xp[0] - wb/4 ; xc[1] = xp[1] + wb/4 ; xc[2] = xp[2] + wb/4 ;
    break ;
  case 7:
    xc[0] = xp[0] + wb/4 ; xc[1] = xp[1] + wb/4 ; xc[2] = xp[2] + wb/4 ;
    break ;
  }
  
  xe[0] = xc[0]+0.1*wb ; xe[1] = xc[1]-0.05*wb ; xe[2] = xc[2] + 0.1*wb ;
  /* xe[0] = xc[0] ; xe[1] = xc[1] ; xe[2] = xc[2] ; */
  wbfmm_laplace_expansion_local_evaluate(xp, Cp, 8*nq,
					      Np, nq, xe, fp, work) ;

  /*rotation checks*/
  wbfmm_shift_angles(xp, xc, &th, &ph, &ch, &r) ;  
  fprintf(stderr, "th = %lg; ph = %lg; ch = %lg;\n", th, ph, ch) ;
  wbfmm_coefficients_H_rotation(H03, Np, th, work) ;
  wbfmm_laplace_rotate_H(cr, nq, Cp, 8*nq, Np, nq, H03, ph, ch, 0.0) ;
  wbfmm_laplace_coaxial_translate_RR(ct, nq, Nc, cr, nq, Np, nq, r, 0.0) ;
  wbfmm_laplace_rotate_H(crr, nq, ct, nq, Nc, nq, H03, ch, ph, 0.0) ;

  wbfmm_laplace_expansion_local_evaluate(xc, crr, nq,
					      Nc, nq, xe, fs, work) ;
    
  /* check = crr ; */
  /* for ( n = 0 ; n <= 2 ; n ++ ) { */
  /*   i = n*n ; nu = 0 ; */
  /*   fprintf(stderr, "%d %d %lg %lg %lg\n", n, nu, */
  /* 	    check[i*nq+0], Cc[8*nq*i+nq*quad+0], */
  /* 	    fabs(check[i*nq+0]-Cc[8*nq*i+nq*quad+0])) ; */
  /*   for ( nu = 1 ; nu <= n ; nu ++ ) { */
  /*     i = wbfmm_index_laplace_nm(n,nu) ; */
  /*     fprintf(stderr, "%d %d %lg %lg %lg\n", n, nu, */
  /* 	      check[(i+0)*nq+0], Cc[8*nq*(i+0)+nq*quad+0], */
  /* 	      fabs(check[(i+0)*nq+0]-Cc[8*nq*(i+0)+nq*quad+0])) ; */
  /*     fprintf(stderr, "%d %d %lg %lg %lg\n", n, nu, */
  /* 	      check[(i+1)*nq+0], Cc[8*nq*(i+1)+nq*quad+0], */
  /* 	      fabs(check[(i+1)*nq+0]-Cc[8*nq*(i+1)+nq*quad+0])) ; */
  /*   } */
  /* } */
  
  wbfmm_laplace_expansion_local_evaluate(xc, &(Cc[quad*nq]), nq*8,
					      Nc, nq, xe, fc, work) ;
  wbfmm_laplace_field(xs, xstride, src, sstride, nq, NULL, 0, NULL, 0,
			   nsrc, xe, fr) ;

  fprintf(stderr, "ref:    ") ;
  for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", fr[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "child:  ") ;
  for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", fc[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "parent: ") ;
  for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", fp[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "chain:  ") ;
  for ( i = 0 ; i < nq ; i ++ ) fprintf(stderr, "%lg ", fs[i]) ;
  fprintf(stderr, "\n") ;

  fprintf(stderr, "child:  ") ;
  for ( i = 0 ; i < nq ; i ++ )
    fprintf(stderr, "%1.16e ", fabs(fc[i]-fr[i])) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "parent: ") ;
  for ( i = 0 ; i < nq ; i ++ )
    fprintf(stderr, "%1.16e ", fabs(fp[i]-fr[i])) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "chain:  ") ;
  for ( i = 0 ; i < nq ; i ++ )
    fprintf(stderr, "%1.16e ", fabs(fs[i]-fr[i])) ;
  fprintf(stderr, "\n") ;
  
  return 0 ;
}


gint main(gint argc, gchar **argv)

{
  gdouble *xs, *src, *xf, wb ;
  gdouble x, x0[3] = {0.0}, xc[3] = {0.0}, defaults[128] ;
  gdouble ix[3] = {0}, iy[3] = {0}, iz[3] = {0} ;
  gint nsrc, sstride, xstride, fstride, N, nfld, test, quad ;
  gchar ch, *ipfile ;

  N = 16 ; 
  sstride = 4 ; xstride = 3 ; 
  nsrc = 1 ;
  nfld = 1 ;
  x = 0.025 ;
  test = -1 ;
  ipfile = NULL ;
  ix[0] = iy[1] = iz[2] = 1.0 ; 
  xs = src = xf = NULL ;
  wb = 0.125 ;
  quad = 0 ;

  while ( (ch = getopt(argc, argv, "hf:i:n:N:q:t:w:x:")) != EOF ) {
    switch ( ch ) {
    default: g_assert_not_reached() ; break ;
    case 'h':
      fprintf(stderr,
	      "Usage: %s <options>\n\n"
	      "Options:\n"
	      "  -h print this message and exit\n"
	      "  -i <input file name>\n"
	      "  -n <number of sources>\n"
	      "  -N <expansion order>\n"
	      "  -q <test quadrant>\n"
	      "  -t <test name>\n"
	      "  -w <box width>\n"
	      "  -x <distance, e.g. for translation tests>\n",
	      argv[0]) ;
      return 0 ;
      break ;
    case 'f': nfld = atoi(optarg) ; break ;
    case 'i': ipfile = g_strdup(optarg) ; break ;
    case 'n': nsrc = atoi(optarg) ; break ;
    case 'N': N = atoi(optarg) ; break ;
    case 'q': quad = atoi(optarg) ; break ;
    case 't': test = parse_test(optarg) ; break ;
    case 'w': wb = atof(optarg) ; break ;
    case 'x': x = atof(optarg) ; break ;
    }
  }

  timer = g_timer_new() ;

  wbfmm_laplace_coaxial_translate_init(N+2) ;
  
  if ( ipfile != NULL ) {
    read_data(ipfile, xc, x0, &xs, &src, &nsrc, &xstride, &sstride,
	      &xf, &fstride, &nfld, ix, iy, iz) ;
  } else {
    xs = &(defaults[0]) ;
    xs[0] = 0.1 ; xs[1] = 0.2 ; xs[2] = -0.05 ;
    xstride = 3*sizeof(gdouble) ;
    src = &(xs[3]) ; nsrc = 1 ;
    src[0] = 0.5 ; src[1] = -0.9 ;
    xf = &(src[2]) ;
    xf[0] = 1.5 ; xf[1] = 2.5 ; xf[2] = -1.3 ;

    ix[0] = 0.0 ; ix[1] =  1.0/sqrt(2.0) ; ix[2] = 1.0/sqrt(2.0) ;
    iy[0] = 0.0 ; iy[1] = -1.0/sqrt(2.0) ; iy[2] = 1.0/sqrt(2.0) ;
    iz[0] = 1.0 ; iz[1] =  0.0 ;           iz[2] = 0.0 ;
  }

  if ( test == 0 ) {
    expansion_test(N, x0, xs, xstride, src, sstride, nsrc, xf, nfld) ;

    return 0 ;
  }

  if ( test == 1 ) {
    translation_test(N, x0, xs, xstride, src, sstride, nsrc, xf, nfld, x) ;
    
    return 0 ;
  }

  if ( test == 2 ) {
    rotation_test(N, x0, xs, xstride, src, sstride, nsrc,
		  ix, iy, iz, xf, nfld) ;
    return 0 ;
  }
  
  if ( test == 3 ) {
    translation_SR_test(N, x0, xs, xstride, src, sstride, nsrc, xf, nfld, x) ;
    
    return 0 ;
  }

  if ( test == 4 ) {
    translation_RR_test(N, x0, xs, xstride, src, sstride, nsrc, xf, nfld, x) ;
    
    return 0 ;
  }

  if ( test == 5 ) {
    shift_test(N, xc, x0, xs, xstride, src, sstride, nsrc, xf, nfld) ;
    
    return 0 ;
  }

  if ( test == 6 ) {
    child_parent_test(N, x0, wb, quad, xs, xstride, src, sstride, nsrc,
		      xf, nfld) ;
    
    return 0 ;
  }

  if ( test == 7 ) {
    parent_child_test(N, x0, wb, quad, xs, xstride, src, sstride, nsrc, xf) ;
    
    return 0 ;
  }
  
  if ( test == 8 ) {
    expansion_gradient_test(N, x0, xs, xstride, src, sstride, nsrc, xf, nfld) ;

    return 0 ;
  }

  if ( test == 9 ) {
    local_gradient_test(N, x0, xs, xstride, src, sstride, nsrc, xf, nfld, x) ;
    
    return 0 ;
  }

  if ( test == 10 ) {
    expansion_normal_test(N, x0, xs, xstride, src, sstride, nsrc, xf, nfld) ;

    return 0 ;
  }

  return 0 ;
}
