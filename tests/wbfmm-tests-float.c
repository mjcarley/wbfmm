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
/* 262144 */

GTimer *timer ;

gchar *tests[] = {"legendre",
		  "besselj",
		  "besselh",
		  "expansion",
		  "translation",
		  "rotation",
		  "shift",
		  "location",
		  "box_indices",
		  "box_parents",
		  "box_children",
		  "tree",
		  "child_parent",
		  "shift_local",
		  "translation_local",
		  "neighbours",
		  "interaction_4",
		  "interaction_shifts",
		  "parent_child",
		  "expansion_dipole",
		  "rotations_write",
		  "expansion_normal",
		  "expansion_gradient",
		  "local_gradient",
		  ""} ;

gint rotations_write(gint N, gfloat ix[], gfloat iy[],
		     gfloat iz[]) ;
gint legendre_test(gint N, gfloat C) ;
gint besselj_test(gint N, gfloat x) ;
gint besselh_test(gint N, gfloat x) ;
gint expansion_test(gfloat k, gint N, 
		    gfloat *x0,
		    gfloat *xs, gint xstride,
		    gfloat *src, gint sstride,
		    gint nsrc,
		    gfloat *xf, gint nfld) ;
gint expansion_gradient_test(gfloat k, gint N, 
			     gfloat *x0,
			     gfloat *xs, gint xstride,
			     gfloat *src, gint sstride,
			     gint nsrc,
			     gfloat *xf, gint nfld) ;
gint translation_test(gfloat k, gint N, 
		      gfloat xc[], gfloat x,
		      gfloat *xs, gint xstride,
		      gfloat *src, gint sstride,
		      gint nsrc,
		      gfloat *xf, gint fstride, gint nfld) ;
gint translation_local_test(gfloat k, gint N, 
			    gfloat xc[], gfloat x,
			    gfloat *xs, gint xstride,
			    gfloat *src, gint sstride,
			    gint nsrc,
			    gfloat *xf, gint fstride, gint nfld) ;
gint local_gradient_test(gfloat k, gint N, 
			 gfloat xc[], gfloat x,
			 gfloat *xs, gint xstride,
			 gfloat *src, gint sstride,
			 gint nsrc,
			 gfloat *xf, gint fstride, gint nfld) ;
gint rotation_test(gfloat k, gint N, 
		   gfloat *xc,
		   gfloat *xs, gint xstride,
		   gfloat *src, gint sstride,
		   gint nsrc,
		   gfloat ix[], gfloat iy[], 
		   gfloat iz[],			  
		   gfloat *xf, gint nfld) ;
gint shift_test(gfloat k, gint N, 
		       gfloat *x0, gfloat *x1,
		       gfloat *xs, gint xstride,
		       gfloat *src, gint sstride,
		       gint nsrc,
		gfloat *xf, gint fstride, gint nfld) ;
gint shift_local_test(gfloat k, gint N, 
		       gfloat *x0, gfloat *x1,
		       gfloat *xs, gint xstride,
		       gfloat *src, gint sstride,
		       gint nsrc,
		      gfloat *xf, gint fstride, gint nfld) ;
gint location_test(gfloat *x0, gfloat D) ;
gint box_index_test(gfloat *x0, gfloat D, guint level) ;
gint box_parent_test(gfloat *x0, gfloat D, guint level) ;
gint box_children_test(gfloat *x0, gfloat D, guint level) ;
gint child_parent_test(gfloat k, gint Nc, gfloat *xc, 
		       gfloat wb, gint quad) ;
gint parent_child_test(gfloat k, gint Nc, gfloat *xc, gfloat wb) ;
gint tree_test(gfloat *x0, gfloat D, guint npts) ;
gint interaction_shift_test(gfloat wb) ;
gint check_neighbour_list(guint level, guint32 i, guint32 j, 
			  guint32 k, guint64 *neighbours, gint n) ;
gint neighbour_test(guint level) ;
gint check_interaction_list_4(guint level, guint64 idx,
			      guint32 i, guint32 j, guint32 k, 
			      guint64 *list, gint n) ;
gint ilist4_test(guint level) ;
gint expansion_dipole_test(gfloat k, gint N, 
			   gfloat *x0,
			   gfloat *xs, gint xstride,
			   gfloat *src, gint sstride,
			   gint nsrc,
			   gfloat *xf, gint nfld) ;
gint expansion_normal_test(gfloat k, gint N, 
			   gfloat *x0,
			   gfloat *xs, gint xstride,
			   gfloat *src, gint sstride,
			   gint nsrc,
			   gfloat *xf, gint nfld) ;

static gint parse_test(gchar *arg)

{
  gint i = 0 ;
  
  while ( strlen(tests[i]) != 0) {
    if ( !strcmp(tests[i], arg) ) return i ;
    i ++ ;
  }

  return -1 ;
}

static gint read_data(gchar *ipfile, gfloat *xc, gfloat *x0, 
		      gfloat **xs, gfloat **src, 
		      gint *nsrc, gint *xstride, gint *sstride,
		      gfloat **xf, gint *fstride, gint *nfld,
		      gfloat *ix,gfloat *iy, gfloat *iz)

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
    /* fprintf(stderr, "%s\n", line) ; */
    if ( strncmp(line, "x0:", 3) == 0 ) {
      fprintf(stderr, "x0: ") ;
      sscanf(&(line[3]), "%g %g %g", &(x0[0]), &(x0[1]), &(x0[2])) ;
      fprintf(stderr, "%g %g %g\n", x0[0], x0[1], x0[2]) ;
    }
    if ( strncmp(line, "xc:", 3) == 0 ) {
      fprintf(stderr, "xc: ") ;
      sscanf(&(line[3]), "%g %g %g", &(xc[0]), &(xc[1]), &(xc[2])) ;
      fprintf(stderr, "%g %g %g\n", xc[0], xc[1], xc[2]) ;
    }

    if ( strncmp(line, "ix:", 3) == 0 ) {
      fprintf(stderr, "ix: ") ;
      sscanf(&(line[3]), "%g %g %g", &(ix[0]), &(ix[1]), &(ix[2])) ;
      fprintf(stderr, "%g %g %g\n", ix[0], ix[1], ix[2]) ;
    }

    if ( strncmp(line, "iy:", 3) == 0 ) {
      fprintf(stderr, "iy: ") ;
      sscanf(&(line[3]), "%g %g %g", &(iy[0]), &(iy[1]), &(iy[2])) ;
      fprintf(stderr, "%g %g %g\n", iy[0], iy[1], iy[2]) ;
    }

    if ( strncmp(line, "iz:", 3) == 0 ) {
      fprintf(stderr, "iz: ") ;
      sscanf(&(line[3]), "%g %g %g", &(iz[0]), &(iz[1]), &(iz[2])) ;
      fprintf(stderr, "%g %g %g\n", iz[0], iz[1], iz[2]) ;
    }

    if ( strncmp(line, "sources:", 8) == 0 ) {
      fprintf(stderr, "source list: ") ;
      sscanf(&(line[8]), "%d", nsrc) ;
      fprintf(stderr, "%d source", *nsrc) ;
      fprintf(stderr, (*nsrc == 1 ? "\n" : "s\n")) ;
      *xstride = 3 ; *sstride = 2 ;
      *xs = (gfloat *)g_malloc((*nsrc)*(*xstride)*sizeof(gfloat)) ;
      *src = (gfloat *)g_malloc((*nsrc)*(*sstride)*sizeof(gfloat)) ;
      for ( i = 0 ; i < *nsrc ; i ++ ) {
	fscanf(f, "%g %g %g %g %g",
	       &((*xs)[(*xstride)*i+0]),
	       &((*xs)[(*xstride)*i+1]),
	       &((*xs)[(*xstride)*i+2]),
	       &((*src)[(*sstride)*i+0]),
	       &((*src)[(*sstride)*i+1])) ;
      }
    }
    if ( strncmp(line, "field:", 6) == 0 ) {
      fprintf(stderr, "field points: ") ;
      sscanf(&(line[6]), "%d", nfld) ;
      fprintf(stderr, "%d point", *nfld) ;
      fprintf(stderr, (*nfld == 1 ? "\n" : "s\n")) ;
      *fstride = 5 ;
      *xf = (gfloat *)g_malloc((*nfld)*(*fstride)*sizeof(gfloat)) ;
      for ( i = 0 ; i < *nfld ; i ++ ) {
	fscanf(f, "%g %g %g",
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

gint legendre_test(gint N, gfloat C)

{
  gfloat S, *Pnm1, *Pn ;
  gint n, m ;

  Pnm1 = (gfloat *)g_malloc(sizeof(gfloat)*128) ;
  Pn = (gfloat *)g_malloc(sizeof(gfloat)*128) ;
  
  S = sqrt(1-C*C) ;

  wbfmm_legendre_init_f(C, S, &(Pnm1[0]), &(Pn[0]), &(Pn[1])) ;

  fprintf(stdout, "0 0 %1.16e %1.16e\n", C, Pnm1[0]) ;
  fprintf(stdout, "1 0 %1.16e %1.16e\n", C, Pn[0]) ;
  fprintf(stdout, "1 1 %1.16e %1.16e\n", C, Pn[1]) ;

  for ( n = 2 ; n < N ; n ++ ) {
    wbfmm_legendre_recursion_array_f(&Pnm1, &Pn, n-1, C, S) ;
    for ( m = 0 ; m <= n ; m ++ ) 
      fprintf(stdout, "%d %d %1.16e %1.16e\n", n, m, C, Pn[m]) ;
  }

  return 0 ;
}

gint besselj_test(gint N, gfloat x)

{
  gfloat jn, jnm1 ;
  gint n ;

  wbfmm_bessel_j_init_f(x, &jnm1, &jn) ;
  
  fprintf(stdout, "0 %1.16e %1.16e\n", x, jnm1) ;
  fprintf(stdout, "1 %1.16e %1.16e\n", x, jn) ;

  for ( n = 2 ; n < N ; n ++ ) {
    wbfmm_bessel_j_recursion_f(&jnm1, &jn, x, n-1) ;
    fprintf(stdout, "%d %1.16e %1.16e\n", n, x, jn) ;
  }

  return 0 ;
}

gint besselh_test(gint N, gfloat x)

{
  gfloat hn[2], hnm1[2] ;
  gint n ;

  wbfmm_bessel_h_init_f(x, hnm1, hn) ;

  fprintf(stdout, "0 %1.16e %1.16e %1.16e\n", x, hnm1[0], hnm1[1]) ;
  fprintf(stdout, "1 %1.16e %1.16e %1.16e\n", x, hn[0], hn[1]) ;

  for ( n = 2 ; n < N ; n ++ ) {
    wbfmm_bessel_h_recursion_f(hnm1, hn, x, n-1) ;
    fprintf(stdout, "%d %1.16e %1.16e %1.16e\n", n, x, hn[0], hn[1]) ;
  }

  return 0 ;
}

gint expansion_test(gfloat k, gint N, 
		    gfloat *x0,
		    gfloat *xs, gint xstride,
		    gfloat *src, gint sstride,
		    gint nsrc,
		    gfloat *xf, gint nfld)

{
  gfloat cfft[4096] = {0}, work[1024], field[2] ;
  gint i, cstr ;
  gdouble t0 ;

  cstr = 2 ;
  t0 = g_timer_elapsed(timer, NULL) ;
  fprintf(stderr, "%s start: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  for ( i = 0 ; i < nsrc ; i ++ ) {
    wbfmm_expansion_h_cfft_f(k, N, x0, &(xs[i*xstride]), &(src[i*sstride]),
				cfft, cstr, work) ;
  }

  fprintf(stderr, "%s expansion generated: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  for ( i = 0 ; i < nfld ; i ++ ) {
    field[0] = field[1] = 0.0 ;
    wbfmm_expansion_h_evaluate_f(k, x0, cfft, cstr, N, &(xf[i*xstride]), 
				    field, work) ;

    fprintf(stdout, "%g+j%g ", field[0], field[1]) ;

    work[0] = work[1] = 0.0 ;
    wbfmm_total_field_f(k, xs, xstride, src, sstride,
			   NULL, 0, NULL, 0,			   
			   nsrc, &(xf[i*xstride]), work) ;

    fprintf(stdout, "%g+j%g (%g)\n", work[0], work[1], 
	    sqrt((field[0]-work[0])*(field[0]-work[0]) +
		 (field[1]-work[1])*(field[1]-work[1]))) ;
  }

  fprintf(stderr, "%s end: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  return 0 ;
}

gint expansion_gradient_test(gfloat k, gint N, 
			     gfloat *x0,
			     gfloat *xs, gint xstride,
			     gfloat *src, gint sstride,
			     gint nsrc,
			     gfloat *xf, gint nfld)

{
  gfloat cfft[4096] = {0}, work[1024] ;
  gint i, cstr, fstr ;
  gdouble t0 ;

  cstr = 2 ; fstr = 4 ;
  t0 = g_timer_elapsed(timer, NULL) ;
  fprintf(stderr, "%s start: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  for ( i = 0 ; i < nsrc ; i ++ ) {
    wbfmm_expansion_h_cfft_f(k, N, x0, &(xs[i*xstride]), &(src[i*sstride]),
				cfft, cstr, work) ;
  }

  fprintf(stderr, "%s expansion generated: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  for ( i = 0 ; i < nfld ; i ++ ) {
    gfloat field[6] = {0.0}, dz=1e-6 ;

    wbfmm_expansion_h_grad_evaluate_f(k, x0, cfft, cstr, N,
					 &(xf[i*xstride]), 
					 field, fstr, work) ;

    fprintf(stdout,
	    "%g+j*%g %g+j*%g %g+j*%g\n",
	    field[0], field[1], field[2], field[3], field[4], field[5]) ;

    work[0] = work[1] = work[2] = work[3] = work[4] = work[5] = 0.0 ;
    wbfmm_total_field_grad_f(k, xs, xstride, src, sstride,
				NULL, 0, NULL, 0,			   
				nsrc, &(xf[i*xstride]), work) ;

    fprintf(stdout,
	    "%g+j*%g %g+j*%g %g+j*%g\n",
	    work[0], work[1], work[2], work[3], work[4], work[5]) ;

    fprintf(stdout,
	    "%g %g %g\n",
	    sqrt((work[0]-field[0])*(work[0]-field[0]) +
		 (work[1]-field[1])*(work[1]-field[1])),
	    sqrt((work[2]-field[2])*(work[2]-field[2]) +
		 (work[3]-field[3])*(work[3]-field[3])),
	    sqrt((work[4]-field[4])*(work[4]-field[4]) +
		 (work[5]-field[5])*(work[5]-field[5]))) ;  
    
  }

  fprintf(stderr, "%s end: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  return 0 ;
}

gint expansion_dipole_test(gfloat k, gint N, 
			   gfloat *x0,
			   gfloat *xs, gint xstride,
			   gfloat *src, gint sstride,
			   gint nsrc,
			   gfloat *xf, gint nfld)

{
  gfloat cfft[4096] = {0}, work[1024], field[2], dipole[6] ;
  gfloat *fx, *fy, *fz ;
  gint i, cstr ;
  gdouble t0 ;

  nsrc = 1 ;

  fx = &(dipole[0]) ; fy = &(dipole[2]) ; fz = &(dipole[4]) ;
  fx[0] = 1.3 ; fx[1] = -0.9 ;
  fy[0] = 0.3 ; fy[1] =  0.9 ;
  fz[0] = 0.4 ; fz[1] =  0.1 ;
  
  cstr = 2 ; 
  t0 = g_timer_elapsed(timer, NULL) ;
  fprintf(stderr, "%s start: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  for ( i = 0 ; i < nsrc ; i ++ ) {
    wbfmm_expansion_dipole_h_cfft_f(k, N, x0,
				       &(xs[i*xstride]), fx, fy, fz,
				       cfft, cstr, work) ;
  }

  fprintf(stderr, "%s expansion generated: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  for ( i = 0 ; i < nfld ; i ++ ) {
    field[0] = field[1] = 0.0 ;
    wbfmm_expansion_h_evaluate_f(k, x0, cfft, cstr, N, &(xf[i*xstride]), 
				    field, work) ;

    fprintf(stdout, "%g+j*%g ", field[0], field[1]) ;

    work[0] = work[1] = 0.0 ;
    wbfmm_total_dipole_field_f(k, xs, xstride, dipole, sstride, nsrc,
				  &(xf[i*xstride]), work) ;

    fprintf(stdout, "%g+j*%g (%g)\n", work[0], work[1], 
	    sqrt((field[0]-work[0])*(field[0]-work[0]) +
		 (field[1]-work[1])*(field[1]-work[1]))) ;
  }

  fprintf(stderr, "%s end: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  return 0 ;
}

gint expansion_normal_test(gfloat k, gint N, 
			   gfloat *x0,
			   gfloat *xs, gint xstride,
			   gfloat *src, gint sstride,
			   gint nsrc,
			   gfloat *xf, gint nfld)

{
  gfloat cfft[16384] = {0}, work[16384], field[2], dipole[6] ;
  gfloat n[3], q[2] ;
  gfloat *fx, *fy, *fz ;
  gint i, cstr ;
  gdouble t0 ;

  nsrc = 1 ;

  /* n[0] = 0.3 ; n[1] = -0.2 ; n[2] = 0.5 ; */
  n[0] = 0.1 ; n[1] = -0.3 ; n[2] = 1.0 ;
  q[0] = 1.7 ; q[1] = -0.3 ;
  
  /* fx = &(dipole[0]) ; fy = &(dipole[2]) ; fz = &(dipole[4]) ; */
  /* fx[0] = 1.3 ; fx[1] = -0.9 ; */
  /* fy[0] = 0.3 ; fy[1] =  0.9 ; */
  /* fz[0] = 0.4 ; fz[1] =  0.1 ; */

  
  fx = &(dipole[0]) ; fy = &(dipole[2]) ; fz = &(dipole[4]) ;
  fx[0] = n[0]*q[0] ; fx[1] = n[0]*q[1] ;
  fy[0] = n[1]*q[0] ; fy[1] = n[1]*q[1] ;
  fz[0] = n[2]*q[0] ; fz[1] = n[2]*q[1] ;
  
  cstr = 8 ; 
  t0 = g_timer_elapsed(timer, NULL) ;
  fprintf(stderr, "%s start: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  for ( i = 0 ; i < nsrc ; i ++ ) {
    wbfmm_expansion_normal_h_cfft_f(k, N, x0,
    				       &(xs[i*xstride]), n, q,
    				       cfft, cstr, work) ;
    /* wbfmm_expansion_dipole_h_cfft_f(k, N, x0, */
    /* 				       &(xs[i*xstride]), fx, fy, fz, */
    /* 				       cfft, cstr, work) ; */
  }

  fprintf(stderr, "%s expansion generated: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  for ( i = 0 ; i < nfld ; i ++ ) {
    field[0] = field[1] = 0.0 ;
    wbfmm_expansion_h_evaluate_f(k, x0, cfft, cstr, N, &(xf[i*xstride]), 
				    field, work) ;

    fprintf(stdout, "%g+j*%g ", field[0], field[1]) ;

    work[0] = work[1] = 0.0 ;
    wbfmm_total_normal_field_f(k, xs, xstride, n, 1, q, 1, nsrc,
    				  &(xf[i*xstride]), work) ;
    /* wbfmm_total_dipole_field_f(k, xs, xstride, dipole, 1, nsrc, */
    /* 				  &(xf[i*xstride]), work) ; */

    fprintf(stdout, "%g+j*%g (%g)\n", work[0], work[1], 
	    sqrt((field[0]-work[0])*(field[0]-work[0]) +
		 (field[1]-work[1])*(field[1]-work[1]))) ;
  }

  fprintf(stderr, "%s end: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  return 0 ;
}

gint translation_test(gfloat k, gint N, 
		      gfloat xc[], gfloat x,
		      gfloat *xs, gint xstride,
		      gfloat *src, gint sstride,
		      gint nsrc,
		      gfloat *xf, gint fstride, gint nfld)

/*
  coaxial translation from xc to x0 and check on evaluation of field
*/

{
  gfloat Ci[BUFSIZE] = {0}, Co[BUFSIZE] = {0.0}, shift[BUFSIZE] = {0.0} ;
  gfloat field[2], kr, work[BUFSIZE], x0[3] ;
  gint i, Ni, No, cstri, cstro ;
  gdouble t0 ;

  t0 = g_timer_elapsed(timer, NULL) ;
  fprintf(stderr, "%s start: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  Ni = No = N ;
  if ( N > 12 ) Ni = N - 3 ; 
  cstri = 2 ;
  cstro = 3 ;

  /*generate coaxial shift coefficients*/
  /* kr = -k*(xc[2] - x0[2]) ; */
  /* kr = k*(x0[2] - xc[2]) ; */
  kr = k*x ;
  x0[0] = xc[0] ; x0[1] = xc[1] ; x0[2] = xc[2] + x ; 
  wbfmm_coefficients_RR_coaxial_f(shift, N, kr, work) ;
  fprintf(stderr, "%s coaxial coefficients generated: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  /* memset(work, 0, BUFSIZE*sizeof(gfloat)) ; */

  fprintf(stderr, "initial expansion: %g %g %g\n", 
	  xc[0], xc[1], xc[2]) ;
  fprintf(stderr, "shifted expansion: %g %g %g\n", 
	  x0[0], x0[1], x0[2]) ;
  fprintf(stderr, "coaxial translation: kr = %g\n", kr) ;

  /*expand about origin*/
  for ( i = 0 ; i < nsrc ; i ++ ) {
    wbfmm_expansion_h_cfft_f(k, Ni, xc, &(xs[i*xstride]), &(src[i*sstride]),
				Ci, cstri, work) ;
  }

  fprintf(stderr, "%s initial expansion generated: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  /* memset(work, 0, 1024*sizeof(gfloat)) ; */

  /*apply shift*/
  wbfmm_coaxial_translate_f(Co, cstro, No, Ci, cstri, Ni, shift, N,
			       FALSE) ;
  fprintf(stderr, "%s coefficients translated: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  /*compute field in both frames*/
  for ( i = 0 ; i < nfld ; i ++ ) {
    field[0] = field[1] = 0.0 ;
    wbfmm_expansion_h_evaluate_f(k, xc, Ci, cstri, Ni, &(xf[i*fstride]),
				    field, work) ;

    fprintf(stdout, "%g+j*%g ", field[0], field[1]) ;

    field[0] = field[1] = 0.0 ;
    wbfmm_expansion_h_evaluate_f(k, x0, Co, cstro, No, &(xf[i*fstride]),
				    field, work) ;

    fprintf(stdout, "%g+j*%g ", field[0], field[1]) ;

    work[0] = work[1] = 0.0 ;
    wbfmm_total_field_f(k, xs, xstride, src, sstride,
			   NULL, 0, NULL, 0,			   
			   nsrc, &(xf[i*fstride]), work) ;

    fprintf(stdout, "%g+j*%g ", work[0], work[1]) ;

    fprintf(stdout, "(%g)\n",
  	    sqrt((field[0]-work[0])*(field[0]-work[0]) +
  		 (field[1]-work[1])*(field[1]-work[1]))) ;
  }

  fprintf(stderr, "%s ends: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  return 0 ;
}

gint translation_local_test(gfloat k, gint N, 
			    gfloat xc[], gfloat x,
			    gfloat *xs, gint xstride,
			    gfloat *src, gint sstride,
			    gint nsrc,
			    gfloat *xf, gint fstride, gint nfld)

/*
  coaxial translation from xc to x0 and check on evaluation of field
*/

{
  gfloat Ci[BUFSIZE] = {0}, Co[BUFSIZE] = {0.0}, shift[BUFSIZE] = {0.0} ;
  gfloat field[2], kr, work[BUFSIZE], xr[3], x0[3] ;
  gint i, Ni, No, cstri, cstro ;
  gdouble t0 ;

  t0 = g_timer_elapsed(timer, NULL) ;
  fprintf(stderr, "%s start: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  Ni = No = N ;
  if ( N > 12 ) Ni = N - 3 ; 
  cstri = 1 ;
  cstro = 1 ;

  x0[0] = xc[0] ; x0[1] = xc[1] ; x0[2] = xc[2] + x ; 
  
  xr[0] = x0[0] + 0.01 ;
  xr[1] = x0[1] + 0.01 ;
  xr[2] = x0[2] + 0.01 ;

  /*generate coaxial shift coefficients*/
  /* kr = -k*(xc[2] - x0[2]) ; */
  kr = k*(x0[2] - xc[2]) ;
  wbfmm_coefficients_SR_coaxial_f(shift, N, kr, work) ;
  fprintf(stderr, "%s coaxial coefficients generated: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  /* memset(work, 0, BUFSIZE*sizeof(gfloat)) ; */

  fprintf(stderr, "initial expansion: %g %g %g\n", 
	  xc[0], xc[1], xc[2]) ;
  fprintf(stderr, "shifted expansion: %g %g %g\n", 
	  x0[0], x0[1], x0[2]) ;
  fprintf(stderr, "coaxial translation: kr = %g\n", kr) ;

  /*expand about origin*/
  for ( i = 0 ; i < nsrc ; i ++ ) {
    wbfmm_expansion_h_cfft_f(k, Ni, xc, &(xs[i*xstride]), &(src[i*sstride]),
				Ci, cstri, work) ;
  }

  fprintf(stderr, "%s initial expansion generated: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  /*apply shift*/
  wbfmm_coaxial_translate_f(Co, cstro, No, Ci, cstri, Ni, shift, N,
			       TRUE) ;
  fprintf(stderr, "%s coefficients translated: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  /*compute field in both frames*/
  for ( i = 0 ; i < nfld ; i ++ ) {
    field[0] = field[1] = 0.0 ;
    wbfmm_expansion_h_evaluate_f(k, xc, Ci, cstri, Ni, xr, field, work) ;

    fprintf(stdout, "%g+j*%g ", field[0], field[1]) ;

    field[0] = field[1] = 0.0 ;
    wbfmm_expansion_j_evaluate_f(k, x0, Co, cstro, No, xr, field, work) ;

    fprintf(stdout, "%g+j*%g ", field[0], field[1]) ;

    work[0] = work[1] = 0.0 ;
    wbfmm_total_field_f(k, xs, xstride, src, sstride,
			   NULL, 0, NULL, 0,
			   nsrc, xr, work) ;

    fprintf(stdout, "%g+j*%g ", work[0], work[1]) ;

    fprintf(stdout, "(%g)\n",
  	    sqrt((field[0]-work[0])*(field[0]-work[0]) +
  		 (field[1]-work[1])*(field[1]-work[1]))) ;
  }

  fprintf(stderr, "%s ends: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  return 0 ;
}

gint local_gradient_test(gfloat k, gint N, 
			 gfloat xc[], gfloat x,
			 gfloat *xs, gint xstride,
			 gfloat *src, gint sstride,
			 gint nsrc,
			 gfloat *xf, gint fstride, gint nfld)

/*
  coaxial translation from xc to x0 and check on evaluation of field
*/

{
  gfloat Ci[BUFSIZE] = {0}, Co[BUFSIZE] = {0.0}, shift[BUFSIZE] = {0.0} ;
  gfloat fl[32]={0.0}, fe[32]={0.0}, fc[32]={0.0} ;
  gfloat kr, work[BUFSIZE], xr[3], x0[3] ;
  gint i, Ni, No, cstri, cstro, fstr ;
  gdouble t0 ;

  t0 = g_timer_elapsed(timer, NULL) ;
  fprintf(stderr, "%s start: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  Ni = No = N ; fstr = 3 ;
  if ( N > 12 ) Ni = N - 3 ; 
  cstri = 1 ;
  cstro = 1 ;

  x0[0] = xc[0] ; x0[1] = xc[1] ; x0[2] = xc[2] + x ; 
  
  xr[0] = x0[0] + 0.1 ;
  xr[1] = x0[1] + 0.2 ;
  xr[2] = x0[2] + 0.3 ;

  /*generate coaxial shift coefficients*/
  /* kr = -k*(xc[2] - x0[2]) ; */
  kr = k*(x0[2] - xc[2]) ;
  wbfmm_coefficients_SR_coaxial_f(shift, N, kr, work) ;
  fprintf(stderr, "%s coaxial coefficients generated: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  /* memset(work, 0, BUFSIZE*sizeof(gfloat)) ; */

  fprintf(stderr, "initial expansion: %g %g %g\n", 
	  xc[0], xc[1], xc[2]) ;
  fprintf(stderr, "shifted expansion: %g %g %g\n", 
	  x0[0], x0[1], x0[2]) ;
  fprintf(stderr, "coaxial translation: kr = %g\n", kr) ;

  /*expand about origin*/
  for ( i = 0 ; i < nsrc ; i ++ ) {
    wbfmm_expansion_h_cfft_f(k, Ni, xc, &(xs[i*xstride]), &(src[i*sstride]),
				Ci, cstri, work) ;
  }

  fprintf(stderr, "%s initial expansion generated: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  /*apply shift*/
  wbfmm_coaxial_translate_f(Co, cstro, No, Ci, cstri, Ni, shift, N,
			       TRUE) ;
  fprintf(stderr, "%s coefficients translated: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  /*compute field in both frames*/
  for ( i = 0 ; i < nfld ; i ++ ) {
    memset(fe, 0, 32*sizeof(gfloat)) ;
    memset(fl, 0, 32*sizeof(gfloat)) ;
    memset(fc, 0, 32*sizeof(gfloat)) ;
    wbfmm_total_field_grad_f(k, xs, xstride, src, sstride,
				NULL, 0, NULL, 0,
				nsrc, xr, fc) ;

    fprintf(stdout, "direct:      %g+j*%g "
	    "%g+j*%g %g+j*%g\n",
	    fc[0], fc[1], fc[2], fc[3], fc[4], fc[5]) ;

    wbfmm_expansion_h_grad_evaluate_f(k, xc, Ci, cstri, Ni, xr, fe, fstr,
					 work) ;

    fprintf(stdout, "h expansion: %g+j*%g "
	    "%g+j*%g %g+j*%g\n",
	    fe[0], fe[1], fe[2], fe[3], fe[4], fe[5]) ;

    fprintf(stdout, "             %g %g %g\n",
  	    sqrt((fe[0]-fc[0])*(fe[0]-fc[0]) +
  		 (fe[1]-fc[1])*(fe[1]-fc[1])),
  	    sqrt((fe[2]-fc[2])*(fe[2]-fc[2]) +
  		 (fe[3]-fc[3])*(fe[3]-fc[3])),
  	    sqrt((fe[4]-fc[4])*(fe[4]-fc[4]) +
  		 (fe[5]-fc[5])*(fe[5]-fc[5]))
	    ) ;

    wbfmm_expansion_j_grad_evaluate_f(k, x0, Co, cstro, No, xr, fl, fstr,
					 work) ;

    fprintf(stdout, "j expansion: %g+j*%g "
	    "%g+j*%g %g+j*%g\n",
	    fl[0], fl[1], fl[2], fl[3], fl[4], fl[5]) ;

    fprintf(stdout, "             %g %g %g\n",
  	    sqrt((fl[0]-fc[0])*(fl[0]-fc[0]) +
  		 (fl[1]-fc[1])*(fl[1]-fc[1])),
  	    sqrt((fl[2]-fc[2])*(fl[2]-fc[2]) +
  		 (fl[3]-fc[3])*(fl[3]-fc[3])),
  	    sqrt((fl[4]-fc[4])*(fl[4]-fc[4]) +
  		 (fl[5]-fc[5])*(fl[5]-fc[5]))
	    ) ;
  }

  fprintf(stderr, "%s ends: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  return 0 ;
}

gint rotation_test(gfloat k, gint N, 
		   gfloat *xc,
		   gfloat *xs, gint xstride,
		   gfloat *src, gint sstride,
		   gint nsrc,
		   gfloat ix[], gfloat iy[], 
		   gfloat iz[],			  
		   gfloat *xf, gint nfld)

{
  gfloat H[BUFSIZE*2], work[BUFSIZE], th, ph, ch ;
  gfloat ix0[3], iy0[3], iz0[3] ;
  gfloat Ci[BUFSIZE*2] = {0}, Co[BUFSIZE*2] = {0.0} ;
  gfloat furot[2], frot[2], ref[2], y[3] ;
  gint i, cstri, cstro ;
  gdouble t0, dt ;
  
  cstri = 3 ; cstro = 2 ;

  ix0[0] = 1.0 ; ix0[1] = 0.0 ; ix0[2] = 0.0 ;
  iy0[0] = 0.0 ; iy0[1] = 1.0 ; iy0[2] = 0.0 ;
  iz0[0] = 0.0 ; iz0[1] = 0.0 ; iz0[2] = 1.0 ;

  wbfmm_rotation_angles_f(ix0, iy0, iz0, ix, iy, iz, &th, &ph, &ch) ;

  fprintf(stderr, "rotation: (%g,%g,%g)\n", th, ph, ch) ;
  t0 = g_timer_elapsed(timer, NULL) ;
  fprintf(stderr, "%s start: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  /*expand about origin*/
  for ( i = 0 ; i < nsrc ; i ++ ) {
    wbfmm_expansion_h_cfft_f(k, N, xc, &(xs[i*xstride]), &(src[i*sstride]),
				Ci, cstri, work) ;
  }
  /*fill H with rubbish to make sure entries are being set in the
    function call*/
  /* memset(H, 255, BUFSIZE*sizeof(gfloat)) ; */
  wbfmm_coefficients_H_rotation_f(H, N, th, work) ;

  /*apply the rotation to the coefficients*/
  fprintf(stderr, "%s reference rotation: %lg\n",
	  __FUNCTION__, (dt = g_timer_elapsed(timer, NULL)) - t0) ;
  wbfmm_rotate_H_f(Co, cstro, Ci, cstri, N, H, ph, ch) ;
  dt = g_timer_elapsed(timer, NULL) - dt ;
  fprintf(stderr, "%s rotation complete: %lg (%lg)\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0, dt) ;

#ifdef WBFMM_USE_AVX
  memset(Co, 0, 2*BUFSIZE*sizeof(gfloat)) ;
  fprintf(stderr, "%s AVX rotation: %lg\n",
	  __FUNCTION__, (dt = g_timer_elapsed(timer, NULL)) - t0) ;
  wbfmm_rotate_H_f(Co, cstro, Ci, cstri, N, H, ph, ch) ;
  dt = g_timer_elapsed(timer, NULL) - dt ;
  fprintf(stderr, "%s rotation complete: %lg (%lg)\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0, dt) ;
#endif /*WBFMM_USE_AVX*/

  
  for ( i = 0 ; i < nfld ; i ++ ) {
    ref[0] = ref[1] = 0.0 ;
    wbfmm_total_field_f(k, xs, xstride, src, sstride,
			   NULL, 0, NULL, 0,
			   nsrc, &(xf[i*xstride]), ref) ;
    /*computed field on unrotated coefficients*/
    furot[0] = furot[1] = 0.0 ;
    wbfmm_expansion_h_evaluate_f(k, xc, Ci, cstri, N, &(xf[i*xstride]), 
				    furot, work) ;

    fprintf(stdout, "%g+j*%g ", furot[0], furot[1]) ;

    wbfmm_coordinate_transform_f(&(xf[i*xstride]), ix, iy, iz, y) ;

    frot[0] = frot[1] = 0.0 ;
    wbfmm_expansion_h_evaluate_f(k, xc, Co, cstro, N, y, frot, work) ;

    fprintf(stdout, "%g+j*%g ", frot[0], frot[1]) ;

    fprintf(stdout, "(%g, %g)\n",
	    sqrt((furot[0]-ref[0])*(furot[0]-ref[0]) +
		 (furot[1]-ref[1])*(furot[1]-ref[1])),
	    sqrt((frot[0]-ref[0])*(frot[0]-ref[0]) +
		 (frot[1]-ref[1])*(frot[1]-ref[1]))) ;
		 
  }


  return 0 ;
}

gint shift_test(gfloat k, gint N, 
		gfloat *x0, gfloat *x1,
		gfloat *xs, gint xstride,
		gfloat *src, gint sstride,
		gint nsrc,
		gfloat *xf, gint fstride, gint nfld)
  
{
  gfloat *work ;    
  gfloat th0, ph0, ch0, th1, ph1, ch1, kr ;
  gfloat ix0[3], iy0[3], iz0[3], ix[3], iy[3], iz[3] ;
  gfloat *C0, *C1, *Ca, *Cb, *shift, *H0, *H1 ;
  gfloat f0[2], f1[2], ref[2] ;
  gint i, cstr0, cstr1 ;
  gdouble t0 ;

  cstr0 = 3 ; cstr1 = 8 ;

  fprintf(stderr, 
	  "buffer sizes\n"
	  "shift:        %u\n"
	  "rotations:    %u\n"
	  "coefficients: %u\n"
	  "workspace:    %u\n",
	  wbfmm_element_number_coaxial(N),	  
	  wbfmm_element_number_rotation(N),
	  wbfmm_coefficient_index_nm(N+1,0),
	  wbfmm_element_number_rotation(2*N)) ;

  shift = (gfloat *)g_malloc0(wbfmm_element_number_coaxial(N)*
				  sizeof(gfloat)) ;
  C0 = (gfloat *)g_malloc0(2*wbfmm_coefficient_index_nm(N+1,0)*cstr0*
			       sizeof(gfloat)) ;
  C1 = (gfloat *)g_malloc0(2*wbfmm_coefficient_index_nm(N+1,0)*cstr1*
			       sizeof(gfloat)) ;
  Ca = (gfloat *)g_malloc0(2*wbfmm_coefficient_index_nm(N+1,0)*cstr0*
			       sizeof(gfloat)) ;
  Cb = (gfloat *)g_malloc0(2*wbfmm_coefficient_index_nm(N+1,0)*cstr1*
			       sizeof(gfloat)) ;
  H0 = (gfloat *)g_malloc0(wbfmm_element_number_rotation(N)*
			       sizeof(gfloat)) ;
  H1 = (gfloat *)g_malloc0(wbfmm_element_number_rotation(N)*
			       sizeof(gfloat)) ;
  work = (gfloat *)g_malloc0(wbfmm_element_number_rotation(2*N)*
				 sizeof(gfloat)) ;
  
  t0 = g_timer_elapsed(timer, NULL) ;
  fprintf(stderr, "%s start: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  ix0[0] = 1.0 ; ix0[1] = 0.0 ; ix0[2] = 0.0 ;
  iy0[0] = 0.0 ; iy0[1] = 1.0 ; iy0[2] = 0.0 ;
  iz0[0] = 0.0 ; iz0[1] = 0.0 ; iz0[2] = 1.0 ;

  wbfmm_shift_coordinates_f(x0, x1, ix, iy, iz, &kr) ;
  kr *= k ;

  /*generate rotation coefficients in each direction and translation*/
  wbfmm_rotation_angles_f(ix0, iy0, iz0, ix, iy, iz, &th0, &ph0, &ch0) ;
  wbfmm_coefficients_H_rotation_f(H0, N, th0, work) ;
  fprintf(stderr, "%s rotation coefficients generated: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  wbfmm_coefficients_RR_coaxial_f(shift, N, kr, work) ;
  fprintf(stderr, "%s coaxial coefficients generated: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  wbfmm_rotation_angles_f(ix, iy, iz, ix0, iy0, iz0, &th1, &ph1, &ch1) ;
  wbfmm_coefficients_H_rotation_f(H1, N, th1, work) ;
  fprintf(stderr, "%s inverse rotation coefficients generated: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  fprintf(stderr, "rotation: (%g,%g,%g)\n",
	  th0, ph0, ch0) ;
  fprintf(stderr, "rotation: (%g,%g,%g)\n",
	  th1, ph1, ch1) ;

  /*expand about x0*/
  for ( i = 0 ; i < nsrc ; i ++ ) {
    wbfmm_expansion_h_cfft_f(k, N, x0, &(xs[i*xstride]), &(src[i*sstride]),
				C0, cstr0, work) ;
  }
  fprintf(stderr, "%s expansion coefficients generated: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  /*apply the rotation to the coefficients*/
  wbfmm_rotate_H_f(Ca, cstr0, C0, cstr0, N, H0, ph0, ch0) ;
  /*translate by kr*/
  wbfmm_coaxial_translate_f(Cb, cstr1, N, Ca, cstr0, N, shift, N,
			       FALSE) ;
  /*rotate back*/
  /* wbfmm_rotate_H_f(C1, cstr1, N, Cb, cstr1, H1, ph1, ch1) ; */
  wbfmm_rotate_H_f(C1, cstr1, Cb, cstr1, N, H0, ch0, ph0) ;
  fprintf(stderr, "%s expansion coefficients shifted: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  for ( i = 0 ; i < nfld ; i ++ ) {
    ref[0] = ref[1] = 0.0 ;
    wbfmm_total_field_f(k, xs, xstride, src, sstride,
			   NULL, 0, NULL, 0, nsrc,
			   &(xf[i*xstride]), ref) ;
    /*computed field on unshifted coefficients*/
    f0[0] = f0[1] = 0.0 ;
    wbfmm_expansion_h_evaluate_f(k, x0, C0, cstr0, N, 
				    &(xf[i*xstride]), f0, work) ;

    fprintf(stdout, "%g+j*%g ", f0[0], f0[1]) ;

    /*computed field on shifted coefficients*/
    f1[0] = f1[1] = 0.0 ;
    wbfmm_expansion_h_evaluate_f(k, x1, C1, cstr1, N, 
				    &(xf[i*xstride]), f1, work) ;

    fprintf(stdout, "%g+j*%g ", f1[0], f1[1]) ;

    fprintf(stdout, "%g, %g\n",
	    sqrt((f0[0]-ref[0])*(f0[0]-ref[0]) +
		 (f0[1]-ref[1])*(f0[1]-ref[1])),
	    sqrt((f1[0]-ref[0])*(f1[0]-ref[0]) +
		 (f1[1]-ref[1])*(f1[1]-ref[1]))) ;
		 
  }

  fprintf(stderr, "%s end: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  return 0 ;
}

gint location_test(gfloat *x0, gfloat D)

{
  guint64 xi, bi, box ;
  guint32 xu, yu, zu ;
  gfloat x[3], xb[3], wb ;
  guint level, nb ;
  gint i ;

  level = 1 ;
  nb = 1 << level ;
  wb = D/nb ;

  fprintf(stderr, "%u boxes of width %g at level %u\n", nb, wb, level) ;

  for ( i = 0 ; i < 1 ; i ++ ) {
    /*generate a random box*/
    xb[0] = x0[0] + g_random_int_range(0, nb)*wb ;
    xb[1] = x0[1] + g_random_int_range(0, nb)*wb ;
    xb[2] = x0[2] + g_random_int_range(0, nb)*wb ;
    xb[0] = xb[1] = xb[2] = 0.5 ;
    xb[0] = 0 ;
    bi = wbfmm_point_index_3d_f(xb, x0, D) ;

    x[0] = xb[0] + g_random_double_range(0, wb) ;
    x[1] = xb[1] + g_random_double_range(0, wb) ;
    x[2] = xb[2] + g_random_double_range(0, wb) ;
    xi = wbfmm_point_index_3d_f(x, x0, D) ;
    box = wbfmm_point_locate_box(xi, level) ;
    fprintf(stderr, "(%g,%g,%g) (%g,%g,%g) (%lu,%lu,%lu,",
	    xb[0], xb[1], xb[2],
	    x[0], x[1], x[2], 
	    xi, box, bi) ;
    box = wbfmm_point_locate_box(bi, level) ;
    fprintf(stderr, "%lu) ", box) ;
    wbfmm_point_from_index(box, &xu, &yu, &zu) ;
    fprintf(stderr, "(%u,%u,%u)\n", xu, yu, zu) ;
  }

  return 0 ;
}

gint box_index_test(gfloat *x0, gfloat D, guint level)

{
  guint64 idx ;
  gfloat x[3], wb, w ;
  guint nb ;
  guint32 i, j, k ;

  nb = 1 << level ;
  wb = D/nb ;

  fprintf(stderr, "%u boxes of width %g at level %u\n", nb, wb, level) ;

  for ( i = 0 ; i < nb ; i ++ ) {
    x[0] = x0[0] + (0.5+i)*wb ;
    for ( j = 0 ; j < nb ; j ++ ) {
      x[1] = x0[1] + (0.5+j)*wb ;
      for ( k = 0 ; k < nb ; k ++ ) {
	x[2] = x0[2] + (0.5+k)*wb ;
	idx = wbfmm_box_index(i, j, k) ;
	wbfmm_box_location_from_index_f(idx, level, x0, D, x, &w) ;
	fprintf(stdout, "%lu %u %u %u %u %g %g %g %g\n",
		idx, level, i, j, k, x[0], x[1], x[2], w) ;
      }
    }
  }

  return 0 ;
}

gint box_parent_test(gfloat *x0, gfloat D, guint level)

{
  guint64 idx, p ;
  gfloat x[3], wb, w ;
  guint nb ;
  guint32 i, j, k ;

  g_assert(level > 0) ;

  nb = 1 << level ;
  wb = D/nb ;

  fprintf(stderr, "%u boxes of width %g at level %u\n", nb, wb, level) ;

  for ( i = 0 ; i < nb ; i ++ ) {
    x[0] = x0[0] + (0.5+i)*wb ;
    for ( j = 0 ; j < nb ; j ++ ) {
      x[1] = x0[1] + (0.5+j)*wb ;
      for ( k = 0 ; k < nb ; k ++ ) {
	x[2] = x0[2] + (0.5+k)*wb ;
	idx = wbfmm_box_index(i, j, k) ;
	wbfmm_box_location_from_index_f(idx, level, x0, D, x, &w) ;
	fprintf(stdout, "%lu %u %u %u %u %g %g %g %g ",
		idx, level, i, j, k, x[0], x[1], x[2], w) ;
	p = wbfmm_box_parent(idx) ;
	wbfmm_box_location_from_index_f(p, level-1, x0, D, x, &w) ;
	fprintf(stdout, "%lu %g %g %g %g\n", 
		p, x[0], x[1], x[2], w) ;
      }
    }
  }

  return 0 ;
}

gint box_children_test(gfloat *x0, gfloat D, guint level)

{
  guint64 idx, p ;
  gfloat x[3], xp[3], wb, w, wc ;
  guint nb ;
  guint32 i, j, k, c ;

  nb = 1 << level ;
  wb = D/nb ;

  fprintf(stderr, "%u boxes of width %g at level %u\n", nb, wb, level) ;

  for ( i = 0 ; i < nb ; i ++ ) {
    x[0] = x0[0] + (0.5+i)*wb ;
    for ( j = 0 ; j < nb ; j ++ ) {
      x[1] = x0[1] + (0.5+j)*wb ;
      for ( k = 0 ; k < nb ; k ++ ) {
	x[2] = x0[2] + (0.5+k)*wb ;
	idx = wbfmm_box_index(i, j, k) ;
	wbfmm_box_location_from_index_f(idx, level, x0, D, x, &w) ;
	p = wbfmm_box_first_child(idx) ;
	for ( c = 0 ; c < 8 ; c ++ ) {
	  wbfmm_box_location_from_index_f(p, level+1, x0, D, xp, &wc) ;
	  fprintf(stdout, 
		  "%lu %u %u %u %u %g %g %g %g "
		  "%lu %g %g %g %g\n",
		  idx, level, i, j, k, x[0], x[1], x[2], w,
		  p, xp[0], xp[1], xp[2], wc) ;
	  p ++ ;
	}
      }
    }
  }

  return 0 ;
}

gint tree_test(gfloat *x0, gfloat D, guint npts) 

{
  wbfmm_tree_t *t ;
  gfloat *pts ;
  gsize pstr ;
  gint i, str ;

  t = wbfmm_tree_new_f(x0, D, npts) ;
  str = 7 ;

  /*stride of str to test strided data*/
  pts = (gfloat *)g_malloc(str*npts*sizeof(gfloat)) ;
  pstr = str*sizeof(gfloat) ;

  for ( i = 0 ; i < npts ; i ++ ) {
    pts[i*str+0] = x0[0] + D*g_random_double_range(0,D) ;
    pts[i*str+1] = x0[1] + D*g_random_double_range(0,D) ;
    pts[i*str+2] = x0[2] + D*g_random_double_range(0,D) ;
  }

  wbfmm_tree_add_points_f(t, (gpointer)pts, npts, pstr) ;

  wbfmm_tree_refine_f(t) ;
  wbfmm_tree_refine_f(t) ;

  wbfmm_tree_print(stdout, t, 0, TRUE) ;

  return 0 ;
}

gint child_parent_test(gfloat k, gint Nc, gfloat *xc, 
		       gfloat wb, gint quad) 

{
  gfloat *H03, *H47, *shiftf, *shiftb, *work, *child, *parent ;
  gfloat xsrc[1536], src[1024], xb[3], xf[3] ;
  gfloat fc[2] = {0.0}, fp[2] = {0.0}, fe[2] = {0.0} ;
  gfloat th03, th47 ;
  gfloat kr ;
  gdouble t0, t1 ;
  gint Np, sizew ;

  th03 = acos(sqrt(1.0/3.0)) ; th47 = M_PI - th03 ; 

  Np = Nc + 4 ;
  xf[0] = xc[0] + 5*wb ; 
  xf[1] = xc[1] - 3*wb ; 
  xf[2] = xc[2] + 7*wb ; 

  H03 = (gfloat *)g_malloc0(wbfmm_element_number_rotation(Np)*
				sizeof(gfloat)) ;
  H47 = (gfloat *)g_malloc0(wbfmm_element_number_rotation(Np)*
				sizeof(gfloat)) ;
  shiftf = (gfloat *)g_malloc0(wbfmm_element_number_coaxial(Np)*
				   sizeof(gfloat)) ;
  shiftb = (gfloat *)g_malloc0(wbfmm_element_number_coaxial(Np)*
				   sizeof(gfloat)) ;
  child = (gfloat *)g_malloc0(8*2*wbfmm_coefficient_index_nm(Nc+1,0)*
				  sizeof(gfloat)) ;
  parent = (gfloat *)g_malloc0(8*2*wbfmm_coefficient_index_nm(Np+1,0)*
				   sizeof(gfloat)) ;

  /*calculate and allocate the workspace*/
  sizew = wbfmm_element_number_rotation(2*Np) ;
  sizew = MAX(sizew, 16*(wbfmm_coefficient_index_nm(Np+1,0) +
			 wbfmm_coefficient_index_nm(Nc+1,0))) ;

  work = (gfloat *)g_malloc0(sizew*sizeof(gfloat)) ;

  /*take xc as the centre of the parent box and generate sources*/
  switch ( quad ) {
  default: g_assert_not_reached() ; break ;
  case 0: 
    xb[0] = xc[0] - wb/2 ; xb[1] = xc[1] - wb/2 ; xb[2] = xc[2] - wb/2 ;
    break ;
  case 1: 
    xb[0] = xc[0] + wb/2 ; xb[1] = xc[1] - wb/2 ; xb[2] = xc[2] - wb/2 ;
    break ;
  case 2: 
    xb[0] = xc[0] - wb/2 ; xb[1] = xc[1] + wb/2 ; xb[2] = xc[2] - wb/2 ;
    break ;
  case 3: 
    xb[0] = xc[0] + wb/2 ; xb[1] = xc[1] + wb/2 ; xb[2] = xc[2] - wb/2 ;
    break ;
  case 4: 
    xb[0] = xc[0] - wb/2 ; xb[1] = xc[1] - wb/2 ; xb[2] = xc[2] + wb/2 ;
    break ;
  case 5: 
    xb[0] = xc[0] + wb/2 ; xb[1] = xc[1] - wb/2 ; xb[2] = xc[2] + wb/2 ;
    break ;
  case 6: 
    xb[0] = xc[0] - wb/2 ; xb[1] = xc[1] + wb/2 ; xb[2] = xc[2] + wb/2 ;
    break ;
  case 7: 
    xb[0] = xc[0] + wb/2 ; xb[1] = xc[1] + wb/2 ; xb[2] = xc[2] + wb/2 ;
    break ;  }

  xsrc[0] = xb[0] + 0.125*wb ;
  xsrc[1] = xb[1] + 0.125*wb ;
  xsrc[2] = xb[2] + 0.125*wb ;

  src[0] = 0.5 ; src[1] = -0.3 ;

  /*wb is child box width, or parent box half width*/
  kr = k*sqrt(3.0)*0.5*wb ;

  t0 = g_timer_elapsed(timer, NULL) ;
  fprintf(stderr, "%s start: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  wbfmm_coefficients_H_rotation_f(H03, Np, th03, work) ;
  wbfmm_coefficients_H_rotation_f(H47, Np, th47, work) ;
  fprintf(stderr, "%s rotation coefficients generated: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  wbfmm_coefficients_RR_coaxial_f(shiftf, Np,  kr, work) ;
  wbfmm_coefficients_RR_coaxial_f(shiftb, Np, -kr, work) ;
  fprintf(stderr, "%s coaxial coefficients generated: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  /*generate the child box expansions*/
  wbfmm_expansion_h_cfft_f(k, Nc, xb, xsrc, src, 
			      &(child[2*quad]), 8, work) ;
  fprintf(stderr, "%s child expansion generated: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  t1 = g_timer_elapsed(timer, NULL) ;
  wbfmm_child_parent_shift_f(parent, Np, child, Nc, H03, H47, Np,
  				shiftf, Np, work) ;
  fprintf(stderr, "%s child expansions shifted to parent: %lg (%lg)\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0,
	  g_timer_elapsed(timer, NULL) - t1) ;
  t1 = g_timer_elapsed(timer, NULL) ;
  wbfmm_child_parent_shift_bw_f(parent, Np, child, Nc, H03, Np,
				   shiftf, shiftb, Np, work) ;

  fprintf(stderr, "%s child expansions shifted to parent (BW): %lg (%lg)\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0, 
	  g_timer_elapsed(timer, NULL) - t1) ;

  /*calculate child, parent, and exact fields*/
  wbfmm_expansion_h_evaluate_f(k, xb, &(child[2*quad]), 8, Nc, 
				  xf, fc, work) ;
  wbfmm_expansion_h_evaluate_f(k, xc, &(parent[0]), 8, Np, 
				  xf, fp, work) ;
  wbfmm_total_field_f(k, xsrc, 3, src, 2, NULL, 0, NULL, 0, 1, xf, fe) ;

  fprintf(stderr, "%g+j*%g ", fc[0], fc[1]) ;
  fprintf(stderr, "%g+j*%g ", fp[0], fp[1]) ;
  fprintf(stderr, "(%g, %g)\n", 
	  sqrt((fc[0]-fe[0])*(fc[0]-fe[0]) +
	       (fc[1]-fe[1])*(fc[1]-fe[1])),
	  sqrt((fp[0]-fe[0])*(fp[0]-fe[0]) +
	       (fp[1]-fe[1])*(fp[1]-fe[1]))) ;

  return 0 ;
}

gint parent_child_test(gfloat k, gint Nc, gfloat *xc, 
		       gfloat wb)

{
  gfloat *H03, *H47, *shift, *work, *child, *parent, *SRshift ;
  gfloat xsrc[1536], src[1024], xb[24], xs[3], xf[3] ;
  gfloat fc[2] = {0.0}, fp[2] = {0.0}, fe[2] = {0.0} ;
  /*other way round to child-parent shift rotations because the shifts
    are in the oppposite direction*/
  gfloat th47 = 0.955316618124509, th03 = 2.18627603546528 ;
  gfloat kr, len, dx ;
  gdouble t0 ;
  gint Np, sizew, i, pq ;

  Np = Nc + 4 ;
  /*set a parent quadrant to check indexing*/
  pq = 0 ;
  /*source outside parent box generating field inside*/
  /*xc is centre of parent box*/
  len = 7.5*wb ;
  xs[0] = xc[0] ;
  xs[1] = xc[1] ;
  xs[2] = xc[2] - len ;

  H03 = (gfloat *)g_malloc0(wbfmm_element_number_rotation(Np)*
				sizeof(gfloat)) ;
  H47 = (gfloat *)g_malloc0(wbfmm_element_number_rotation(Np)*
				sizeof(gfloat)) ;
  shift = (gfloat *)g_malloc0(wbfmm_element_number_coaxial(Np)*
				  sizeof(gfloat)) ;
  SRshift = (gfloat *)g_malloc0(wbfmm_element_number_coaxial(Np)*
				    2*sizeof(gfloat)) ;
  child = (gfloat *)g_malloc0(8*2*wbfmm_coefficient_index_nm(Nc+1,0)*
				  sizeof(gfloat)) ;
  parent = (gfloat *)g_malloc0(8*2*wbfmm_coefficient_index_nm(Np+1,0)*
				   sizeof(gfloat)) ;

  /*calculate and allocate the workspace*/
  sizew = wbfmm_element_number_rotation(2*Np) ;
  sizew = MAX(sizew, 32*(wbfmm_coefficient_index_nm(Np+1,0) +
			 wbfmm_coefficient_index_nm(Nc+1,0))) ;

  work = (gfloat *)g_malloc0(sizew*sizeof(gfloat)) ;

  /*take xc as the centre of the parent box*/
  dx = 0.25*wb ;
  xb[0*3+0] = xc[0] - dx ; xb[0*3+1] = xc[1] - dx ; xb[0*3+2] = xc[2] - dx ;
  xb[1*3+0] = xc[0] + dx ; xb[1*3+1] = xc[1] - dx ; xb[1*3+2] = xc[2] - dx ;
  xb[2*3+0] = xc[0] - dx ; xb[2*3+1] = xc[1] + dx ; xb[2*3+2] = xc[2] - dx ;
  xb[3*3+0] = xc[0] + dx ; xb[3*3+1] = xc[1] + dx ; xb[3*3+2] = xc[2] - dx ;
  xb[4*3+0] = xc[0] - dx ; xb[4*3+1] = xc[1] - dx ; xb[4*3+2] = xc[2] + dx ;
  xb[5*3+0] = xc[0] + dx ; xb[5*3+1] = xc[1] - dx ; xb[5*3+2] = xc[2] + dx ;
  xb[6*3+0] = xc[0] - dx ; xb[6*3+1] = xc[1] + dx ; xb[6*3+2] = xc[2] + dx ;
  xb[7*3+0] = xc[0] + dx ; xb[7*3+1] = xc[1] + dx ; xb[7*3+2] = xc[2] + dx ;

  xsrc[0] = xs[0] + 0.*wb ;
  xsrc[1] = xs[1] + 0.*wb ;
  xsrc[2] = xs[2] + 0.*wb ;
  src[0] = 0.5 ; src[1] = -0.3 ;

  /*wb is parent box width*/
  kr = k*sqrt(3.0)*0.25*wb ;

  t0 = g_timer_elapsed(timer, NULL) ;
  fprintf(stderr, "%s start: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  wbfmm_coefficients_H_rotation_f(H03, Np, th03, work) ;
  wbfmm_coefficients_H_rotation_f(H47, Np, th47, work) ;
  fprintf(stderr, "%s rotation coefficients generated: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  wbfmm_coefficients_RR_coaxial_f(shift, Np, kr, work) ;
  fprintf(stderr, "%s coaxial coefficients generated: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  wbfmm_coefficients_SR_coaxial_f(SRshift, Np, k*len, work) ;
  fprintf(stderr, "%s (S|R) coaxial coefficients generated: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  /*generate the source expansion in child for now*/
  wbfmm_expansion_h_cfft_f(k, Np, xs, xsrc, src, child, 1, work) ;
  fprintf(stderr, "%s child expansion generated: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;
  wbfmm_coaxial_translate_f(&(parent[2*pq]), 8, Np, child, 1, Np, 
			       SRshift, Np, TRUE) ;
  
  /*wipe the child data*/
  memset(child, 0, 8*2*wbfmm_coefficient_index_nm(Nc+1,0)*sizeof(gfloat)) ;

  /*the parent box now holds the data for the shift to the children*/
  wbfmm_parent_child_shift_f(child, Nc, parent, Np, H03, H47, Np,
  				shift, Np, work) ;

  fprintf(stderr, "%s parent expansion shifted to children: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;


  /*calculate child, parent, and exact fields in each child box*/
  for ( i = 0 ; i < 8 ; i ++ ) {
    fe[0] = fe[1] = fc[0] = fc[1] = fp[0] = fp[1] = 0.0 ;
    xf[0] = xb[3*i+0] + wb*0.1 ; 
    xf[1] = xb[3*i+1] - wb*0.05 ; 
    xf[2] = xb[3*i+2] + wb*0.1 ; 

    wbfmm_total_field_f(k, xsrc, 3, src, 2, NULL, 0, NULL, 0, 1, xf, fe) ;
    wbfmm_expansion_j_evaluate_f(k, xc, &(parent[2*pq]), 8, Np,
			       xf, fp, work) ;
    wbfmm_expansion_j_evaluate_f(k, &(xb[3*i]), &(child[2*i]), 8, Nc,
			       xf, fc, work) ;
    fprintf(stderr, "%0.4f+j*%0.4f %0.4f+j*%0.4f %0.4f+j*%0.4f (%lg, %lg)\n",
	    fp[0], fp[1], fc[0], fc[1], fe[0], fe[1],
	    sqrt((fp[0]-fe[0])*(fp[0]-fe[0]) + (fp[1]-fe[1])*(fp[1]-fe[1])),
	    sqrt((fc[0]-fe[0])*(fc[0]-fe[0]) + (fc[1]-fe[1])*(fc[1]-fe[1]))) ;
  }

  return 0 ;
}

gint shift_local_test(gfloat k, gint N, 
		      gfloat *x0, gfloat *x1,
		      gfloat *xs, gint xstride,
		      gfloat *src, gint sstride,
		      gint nsrc,
		      gfloat *xf, gint fstride, gint nfld)
  
{
  gfloat *work ;    
  gfloat th0, ph0, ch0, kr ;
  gfloat ix0[3], iy0[3], iz0[3], ix[3], iy[3], iz[3] ;
  gfloat *C0, *C1, *Ca, *Cb, *shift, *H0 ;
  gfloat f0[2], f1[2], ref[2], xr[3] ;
  gint i, cstr0, cstr1, N0, N1 ;
  gdouble t0 ;

  cstr0 = 3 ; cstr1 = 7 ;

  N0 = N ; N1 = N + 4 ;
  fprintf(stderr, 
	  "buffer sizes\n"
	  "shift:        %u\n"
	  "rotations:    %u\n"
	  "coefficients: %u\n"
	  "workspace:    %u\n",
	  wbfmm_element_number_coaxial(N1),
	  wbfmm_element_number_rotation(N1),
	  wbfmm_coefficient_index_nm(N1+1,0),
	  wbfmm_element_number_rotation(2*N1)) ;

  shift = (gfloat *)g_malloc0(2*wbfmm_element_number_coaxial(N1)*
				  sizeof(gfloat)) ;
  C0 = (gfloat *)g_malloc0(2*wbfmm_coefficient_index_nm(N0+1,0)*cstr0*
			       sizeof(gfloat)) ;
  C1 = (gfloat *)g_malloc0(2*wbfmm_coefficient_index_nm(N1+1,0)*cstr1*
			       sizeof(gfloat)) ;
  Ca = (gfloat *)g_malloc0(2*wbfmm_coefficient_index_nm(N0+1,0)*cstr0*
			       sizeof(gfloat)) ;
  Cb = (gfloat *)g_malloc0(2*wbfmm_coefficient_index_nm(N1+1,0)*cstr1*
			       sizeof(gfloat)) ;
  H0 = (gfloat *)g_malloc0(wbfmm_element_number_rotation(N0)*
			       sizeof(gfloat)) ;

  work = (gfloat *)g_malloc0(wbfmm_element_number_rotation(2*N1)*
				 sizeof(gfloat)) ;
  
  t0 = g_timer_elapsed(timer, NULL) ;
  fprintf(stderr, "%s start: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  ix0[0] = 1.0 ; ix0[1] = 0.0 ; ix0[2] = 0.0 ;
  iy0[0] = 0.0 ; iy0[1] = 1.0 ; iy0[2] = 0.0 ;
  iz0[0] = 0.0 ; iz0[1] = 0.0 ; iz0[2] = 1.0 ;

  xr[0] = x1[0] + 0.01 ;
  xr[1] = x1[1] + 0.01 ;
  xr[2] = x1[2] + 0.01 ;

  wbfmm_shift_coordinates_f(x0, x1, ix, iy, iz, &kr) ;
  kr *= k ;

  /*generate rotation coefficients in each direction and translation*/
  wbfmm_rotation_angles_f(ix0, iy0, iz0, ix, iy, iz, &th0, &ph0, &ch0) ;

  wbfmm_shift_angles_f(x0, x1, &th0, &ph0, &ch0, &kr) ;
  kr *= k ;
  wbfmm_coefficients_H_rotation_f(H0, N0, th0, work) ;
  fprintf(stderr, "%s rotation coefficients generated: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  wbfmm_coefficients_SR_coaxial_f(shift, N1, kr, work) ;
  fprintf(stderr, "%s coaxial coefficients generated: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  /* wbfmm_rotation_angles_f(ix, iy, iz, ix0, iy0, iz0, &th1, &ph1, &ch1) ; */
  /* wbfmm_coefficients_H_rotation_f(H1, N1, th1, work) ; */
  fprintf(stderr, "%s inverse rotation coefficients generated: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  fprintf(stderr, "rotation: (%g,%g,%g)\n",
	  th0, ph0, ch0) ;
  /* fprintf(stderr, "rotation: (%g,%g,%g)\n", */
  /* 	  th1, ph1, ch1) ; */

  /*expand about x0*/
  for ( i = 0 ; i < nsrc ; i ++ ) {
    wbfmm_expansion_h_cfft_f(k, N0, x0, &(xs[i*xstride]), &(src[i*sstride]),
				C0, cstr0, work) ;
  }
  fprintf(stderr, "%s expansion coefficients generated: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  /*apply the rotation to the coefficients*/
  wbfmm_rotate_H_f(Ca, cstr0, C0, cstr0, N0, H0, ph0, ch0) ;
  /*translate by kr*/
  wbfmm_coaxial_translate_f(Cb, cstr1, N1, Ca, cstr0, N0, shift, N1,
			       TRUE) ;
  /*rotate back*/
  wbfmm_rotate_H_f(C1, cstr1, Cb, cstr1, N1, H0, ch0, ph0) ;
  fprintf(stderr, "%s expansion coefficients shifted: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  for ( i = 0 ; i < nfld ; i ++ ) {
    ref[0] = ref[1] = 0.0 ;
    wbfmm_total_field_f(k, xs, xstride, src, sstride, NULL, 0, NULL, 0,
			   nsrc, xr, ref) ;
    /*computed field on unshifted coefficients*/
    f0[0] = f0[1] = 0.0 ;
    wbfmm_expansion_h_evaluate_f(k, x0, C0, cstr0, N0, xr, f0, work) ;

    fprintf(stdout, "%g+j*%g ", f0[0], f0[1]) ;

    /*computed field on shifted coefficients*/
    f1[0] = f1[1] = 0.0 ;
    wbfmm_expansion_j_evaluate_f(k, x1, C1, cstr1, N1, xr, f1, work) ;

    fprintf(stdout, "%g+j*%g ", f1[0], f1[1]) ;

    fprintf(stdout, "%g, %g\n",
	    sqrt((f0[0]-ref[0])*(f0[0]-ref[0]) +
		 (f0[1]-ref[1])*(f0[1]-ref[1])),
	    sqrt((f1[0]-ref[0])*(f1[0]-ref[0]) +
		 (f1[1]-ref[1])*(f1[1]-ref[1]))) ;
		 
  }

  fprintf(stderr, "%s end: %lg\n",
	  __FUNCTION__, g_timer_elapsed(timer, NULL) - t0) ;

  return 0 ;
}

gint check_neighbour_list(guint level, guint32 i, guint32 j, 
			  guint32 k, guint64 *neighbours, gint n)

{
  guint32 in, jn, kn ;
  gint ii ;
  guint nb, di, dj, dk ;

  nb = 1 << level ;
  for ( ii = 0 ; ii < n ; ii ++ ) {
    wbfmm_box_location(neighbours[ii], &in, &jn, &kn) ;
    g_assert(in < nb) ; g_assert(jn < nb) ; g_assert(kn < nb) ;
    if ( i >= in ) di = i-in ; else di = in-i ; g_assert(di <= 1) ;
    if ( j >= jn ) dj = j-jn ; else dj = jn-j ; g_assert(dj <= 1) ;
    if ( k >= kn ) dk = k-kn ; else dk = kn-k ; g_assert(dk <= 1) ;
  }

  return 0 ;
}

gint neighbour_test(guint level)

{
  guint nb ;
  guint32 i, j, k ;
  guint64 idx, neighbours[27] ;
  gint n ;

  nb = 1 << level ;

  for ( i = 0 ; i < nb ; i ++ ) {
    for ( j = 0 ; j < nb ; j ++ ) {
      for ( k = 0 ; k < nb ; k ++ ) {
	idx = wbfmm_box_index(i, j, k) ;
	n = wbfmm_box_neighbours(level, idx, neighbours) ;
	check_neighbour_list(level, i, j, k, neighbours, n) ;
	fprintf(stderr, "%u %u %u %d PASS\n", i, j, k, n) ;
      }
    }
  }

  return 0 ;
}

gint check_interaction_list_4(guint level, guint64 idx,
			      guint32 i, guint32 j, guint32 k, 
			      guint64 *list, gint n)

{
  guint32 in, jn, kn, ip, jp, kp ;
  guint64 p, pn, ishift ;
  gint ii ;
  gint nb, di, dj, dk ;

  /*
    list 4 boxes are children of neighbours of the box's parent,
    but are not themselves neighbours of the box
  */
  nb = 1 << level ;

  p = wbfmm_box_parent(idx) ;
  wbfmm_box_location(p, &ip, &jp, &kp) ;

  for ( ii = 0 ; ii < n ; ii ++ ) {
    wbfmm_box_location(list[2*ii+0], &in, &jn, &kn) ;
    g_assert(in < nb) ; g_assert(jn < nb) ; g_assert(kn < nb) ;
    /* if ( i >= in ) di = i-in ; else di = -(gint)(in-i) ; */
    /* if ( j >= jn ) dj = j-jn ; else dj = -(gint)(jn-j) ; */
    /* if ( k >= kn ) dk = k-kn ; else dk = -(gint)(kn-k) ; */
    di = (gint)in - (gint)i ;
    dj = (gint)jn - (gint)j ;
    dk = (gint)kn - (gint)k ;
    /*check boxes are not neighbours*/
    g_assert(!((di*di <= 1) && (dj*dj <= 1) && (dk*dk <= 1))) ;
    
    
    /*check index in shift angle table is correct*/
    ishift = (guint64)((di+3)*49+(dj+3)*7+dk+3) ;

    /* fprintf(stderr, "%lu %lu %lu (%d %d %d)\n", */
    /* 	    list[2*ii+0], list[2*ii+1], ishift, di, dj, dk) ; */
    
    if ( list[2*ii+1] != ishift ) {
      g_error("  %s: ishift %u wrong, %lu should be %lu",
	      __FUNCTION__, ii, list[2*ii+1], ishift) ;
    }
    g_assert(list[2*ii+1] == ishift) ;

    /*check parents ARE neighbours*/
    pn = wbfmm_box_parent(list[2*ii+0]) ;
    wbfmm_box_location(pn, &in, &jn, &kn) ;

    if ( ip >= in ) di = ip-in ; else di = in-ip ; g_assert(di <= 1) ;
    if ( jp >= jn ) dj = jp-jn ; else dj = jn-jp ; g_assert(dj <= 1) ;
    if ( kp >= kn ) dk = kp-kn ; else dk = kn-kp ; g_assert(dk <= 1) ;

  }

  return 0 ;
}

gboolean index_in_list(guint64 *list, gint n, guint64 idx)

{
  gint i ;

  for ( i = 0 ; i < n ; i ++ ) {
    if ( list[2*i+0] == idx ) return TRUE ;
  }
  
  return FALSE ;
}
  
gint check_list_grid(guint level, guint64 idx,
		     guint32 i, guint32 j, guint32 k, 
		     guint64 *list, gint n, guint64 grid[])

{
  guint32 in, jn, kn, ip, jp, kp ;
  guint64 p, pn, ishift ;
  gint nb, di, dj, dk, ii, nn ;

  nb = 1 << level ;

  p = wbfmm_box_parent(idx) ;
  wbfmm_box_location(p, &ip, &jp, &kp) ;

  nn = 0 ;
  for ( di = -3 ; di <= 3 ; di ++ ) {
    for ( dj = -3 ; dj <= 3 ; dj ++ ) {
      for ( dk = -3 ; dk <= 3 ; dk ++ ) {
	ii = (di+3)*49 + (dj+3)*7 + dk + 3 ;
	if ( grid[ii] != 0 ) {
	  pn = grid[ii] - 1 ;
	  g_assert(index_in_list(list, n, pn)) ;
	  nn ++ ;
	}
      }
    }
  }

  /* g_assert(nn == n) ; */
  if ( nn != n )
    g_error("%s: wrong number of entries in grid (%d should be %d)",
	    __FUNCTION__, nn, n) ;    
  
  return 0 ;
}

gint ilist4_test(guint level)

{
  guint nb, p ;
  guint32 i, j, k ;
  guint64 idx, list[378], grid[343] ;
  gint n ;

  nb = 1 << level ;

  fprintf(stderr, "interaction list test\n") ;
  fprintf(stderr, "level = %u\n", level) ;
  
  for ( i = 0 ; i < nb ; i ++ ) {
    for ( j = 0 ; j < nb ; j ++ ) {
      for ( k = 0 ; k < nb ; k ++ ) {
	idx = wbfmm_box_index(i, j, k) ;
	/* fprintf(stderr, "************\n") ;	 */
	n = wbfmm_box_interaction_list_4(level, idx, list, FALSE) ;
	p = wbfmm_box_interaction_grid_4(level, idx, grid) ;
	check_interaction_list_4(level, idx, i, j, k, list, n) ;
	check_list_grid(level, idx, i, j, k, list, n, grid) ;
	fprintf(stderr, "%2u %2u %2u %2u %3d PASS\n", i, j, k, p, n) ;
      }
    }
  }

  return 0 ;
}

gint interaction_shift_test(gfloat wb)

{
  gint i, j, k, idx ;
  gfloat x0[3] = {0.0}, xb[3], th0, ph0, ch0, th1, ph1, ch1, r, rs ;
  gfloat ix0[3], iy0[3], iz0[3], ix[3], iy[3], iz[3] ;
  gdouble emax ;

  ix0[0] = 1.0 ; ix0[1] = 0.0 ; ix0[2] = 0.0 ;
  iy0[0] = 0.0 ; iy0[1] = 1.0 ; iy0[2] = 0.0 ;
  iz0[0] = 0.0 ; iz0[1] = 0.0 ; iz0[2] = 1.0 ;

  wbfmm_shift_angle_table_init_f() ;
  emax = 0.0 ;

  for ( i = -3 ; i <= 3 ; i ++ ) {
    xb[0] = x0[0] + wb*i ;
    for ( j = -3 ; j <= 3 ; j ++ ) {
      xb[1] = x0[1] + wb*j ;
      for ( k = -3 ; k <= 3 ; k ++ ) {
	xb[2] = x0[2] + wb*k ;
	idx = (i+3)*49+(j+3)*7+k+3 ;
	/*shift is from xb to x0*/
	if ( i != 0 || j != 0 || k != 0 ) {
	  wbfmm_shift_coordinates_f(xb, x0, ix, iy, iz, &r) ;
	  wbfmm_rotation_angles_f(ix0, iy0, iz0, ix, iy, iz, 
				     &th0, &ph0, &ch0) ;
	  wbfmm_shift_angles_list4_f(i, j, k, &th1, &ph1, &ch1, &rs) ;
	  /* wbfmm_rotation_angles_f(ix, iy, iz, ix0, iy0, iz0,  */
	  /* 			     &th1, &ph1, &ch1) ; */
	  fprintf(stdout, 
		  "%d %d %d %d %g %g %g %g\n",
		  idx, i, j, k,
		  th0, ph0, ch0, r) ;
	  emax = MAX(emax, fabs(th1-th0)) ;
	  emax = MAX(emax, fabs(ch1-ch0)) ;
	  emax = MAX(emax, fabs(ph1-ph0)) ;
	  emax = MAX(emax, fabs(r - rs*wb)) ;
	  if ( emax > 1e-4 ) {
	    fprintf(stderr, "maximum error %lg, line %d\n",
	  	    emax, idx) ;
	    fprintf(stderr, 
		    "%g %g %g %g\n", 
		    th0, ph0, ch0, r) ;
	    fprintf(stderr, 
		    "%g %g %g %g\n", 
		    th1, ph1, ch1, rs*wb) ;
	    return -1 ;
	  }
	    
	} else {
	  fprintf(stdout, 
		  "%d %d %d %d 8192 8192 8192 0.0\n",
		  idx, i, j, k) ;
		  /* th0, ph0, ch0) ; */

	}
      }
    }
  }

  fprintf(stderr, "maximum error: %lg\n", emax) ;

  return 0 ;
}

gint rotations_write(gint N, gfloat ix[], gfloat iy[],
		     gfloat iz[])

{
  gfloat H[BUFSIZE], work[BUFSIZE], th, ph, ch ;
  gfloat ix0[3], iy0[3], iz0[3] ;

  ix0[0] = 1.0 ; ix0[1] = 0.0 ; ix0[2] = 0.0 ;
  iy0[0] = 0.0 ; iy0[1] = 1.0 ; iy0[2] = 0.0 ;
  iz0[0] = 0.0 ; iz0[1] = 0.0 ; iz0[2] = 1.0 ;

  wbfmm_rotation_angles_f(ix0, iy0, iz0, ix, iy, iz, &th, &ph, &ch) ;

  fprintf(stderr, "rotation: (%g,%g,%g)\n", th, ph, ch) ;

  wbfmm_coefficients_H_rotation_f(H, N, th, work) ;

  wbfmm_rotation_write_coefficients_f(H, N, stdout) ;

  return 0 ;
}


gint main(gint argc, gchar **argv)

{
  gfloat *xs, *src, *xf, wb ;
  gfloat x, k, r, x0[3] = {0.0}, xc[3] = {0.0}, defaults[128] ;
  gfloat ix[3] = {0}, iy[3] = {0}, iz[3] = {0} ;
  gint nsrc, sstride, xstride, fstride, N, i, nfld, test, quad ;
  gchar ch, *ipfile ;
  guint level ;

  k = 2.0 ; N = 16 ; 
  sstride = 2 ; xstride = 3 ; 
  r = 0.25 ;
  nsrc = 1 ;
  nfld = 1 ;
  x = 0.025 ;
  test = -1 ;
  ipfile = NULL ;
  ix[0] = iy[1] = iz[2] = 1.0 ; 
  xs = src = xf = NULL ;
  quad = 0 ;
  wb = 1.0 ;
  level = 2 ;

  while ( (ch = getopt(argc, argv, "1f:i:k:l:n:N:q:t:w:x:")) != EOF ) {
    switch ( ch ) {
    default: g_assert_not_reached() ; break ;
    case 'f': nfld = atoi(optarg) ; break ;
    case 'i': ipfile = g_strdup(optarg) ; break ;
    case 'k': k = atof(optarg) ; break ;
    case 'l': level = atoi(optarg) ; break ;
    case 'n': nsrc = atoi(optarg) ; break ;
    case 'N': N = atoi(optarg) ; break ;
    case 'q': quad = atoi(optarg) ; break ;
    case 't': test = parse_test(optarg) ; break ;
    case 'w': wb = atof(optarg) ; break ;
    case 'x': x = atof(optarg) ; break ;
    }
  }

  timer = g_timer_new() ;

  if ( test == 0 ) {
    legendre_test(N, x) ;

    return 0 ;
  }

  if ( test == 1 ) {
    besselj_test(N, x) ;

    return 0 ;
  }

  if ( test == 2 ) {
    besselh_test(N, x) ;

    return 0 ;
  }

  if ( ipfile != NULL ) {
    read_data(ipfile, xc, x0, &xs, &src, &nsrc, &xstride, &sstride,
	      &xf, &fstride, &nfld, ix, iy, iz) ;
  } else {
    xs = &(defaults[0]) ;
    xs[0] = 0.1 ; xs[1] = 0.2 ; xs[2] = -0.05 ;
    xstride = 3*sizeof(gfloat) ;
    src = &(xs[3]) ; nsrc = 1 ;
    src[0] = 0.5 ; src[1] = -0.9 ;
    xf = &(src[2]) ;
    xf[0] = 1.5 ; xf[1] = 2.5 ; xf[2] = -1.3 ;

    ix[0] = 0.0 ; ix[1] =  1.0/sqrt(2.0) ; ix[2] = 1.0/sqrt(2.0) ;
    iy[0] = 0.0 ; iy[1] = -1.0/sqrt(2.0) ; iy[2] = 1.0/sqrt(2.0) ;
    iz[0] = 1.0 ; iz[1] =  0.0 ;           iz[2] = 0.0 ;
  }    

  if ( test == 3 || test == 4 || test == 5 || test == 6 ||
       test == 19 || test == 21 ) {

    if ( xs == NULL ) {
      src = (gfloat *)g_malloc(nsrc*sstride*sizeof(gfloat)) ;
      xs = (gfloat *)g_malloc(nsrc*xstride*sizeof(gfloat)) ;
      for ( i = 0 ; i < nsrc ; i ++ ) {
	xs[i*xstride+0] = g_random_double_range(-r, r) ;
	xs[i*xstride+1] = g_random_double_range(-r, r) ;
	xs[i*xstride+2] = g_random_double_range(-r, r) ;
	src[i*sstride+0] = g_random_double_range(-1, 1) ;
	src[i*sstride+1] = g_random_double_range(-1, 1) ;
      }
    }

    if ( xf == NULL ) {
      xf = (gfloat *)g_malloc(nfld*xstride*sizeof(gfloat)) ;
      for ( i = 0 ; i < nfld ; i ++ ) {
	xf[i*xstride+0] = g_random_double_range(4*r, 8*r) ;
	xf[i*xstride+1] = g_random_double_range(4*r, 8*r) ;
	xf[i*xstride+2] = g_random_double_range(4*r, 8*r) ;
      }
    }
  }

  if ( test == 3 ) {
    expansion_test(k, N, x0, xs, xstride, src, sstride, nsrc,
		   xf, nfld) ;
    return 0 ;
  }

  if ( test == 4 ) {
    translation_test(k, N, xc, x, xs, xstride, src, sstride, nsrc,
  		     xf, fstride, nfld) ;

    return 0 ;
  }

  if ( test == 5 ) {
    rotation_test(k, N, xc, xs, xstride, src, sstride, nsrc,
		  ix, iy, iz, xf, nfld) ;
    return 0 ;
  }

  if ( test == 6 ) {
    shift_test(k, N, xc, x0, xs, xstride, src, sstride, nsrc,
	       xf, fstride, nfld) ;
    return 0 ;
  }

  if ( test == 7 ) {
    location_test(x0, 1.0) ;

    return 0 ;
  }

  if ( test == 8 ) {
    box_index_test(x0, 1.0, 1) ;

    return 0 ;
  }

  if ( test == 9 ) {
    box_parent_test(x0, 1.0, 2) ;

    return 0 ;
  }

  if ( test == 10 ) {
    box_children_test(x0, 1.0, 1) ;

    return 0 ;
  }

  if ( test == 11 ) {
    tree_test(x0, 1.0, 1000) ;

    return 0 ;
  }

  if ( test == 12 ) {
    child_parent_test(k, N, xc, wb, quad) ;

    return 0 ;
  }

  if ( test == 13 ) {
    shift_local_test(k, N, x0, xc, xs, xstride, src, sstride, nsrc,
		     xf, fstride, nfld) ;
    return 0 ;
  }

  if ( test == 14 ) {
    translation_local_test(k, N, xc, x, xs, xstride, src, sstride, nsrc,
			   xf, fstride, nfld) ;

    return 0 ;
  }

  if ( test == 15 ) {
    neighbour_test(level) ;
    
    return 0 ;
  }

  if ( test == 16 ) {
    ilist4_test(level) ;
    
    return 0 ;
  }

  if ( test == 17 ) {
    interaction_shift_test(wb) ;
    
    return 0 ;
  }

  if ( test == 18 ) {
    parent_child_test(k, N, xc, wb) ;

    return 0 ;
  }

  if ( test == 19 ) {
    xf[0] = 2.0 ;
    expansion_dipole_test(k, N, x0, xs, xstride, src, sstride, nsrc,
			  xf, nfld) ;
    return 0 ;
  }

  if ( test == 20 ) {
    rotations_write(N, ix, iy, iz) ;
    
    return 0 ;
  }

  if ( test == 21 ) {
    xf[0] = 2.0 ; xf[1] = 1.0 ; xf[2] = 3.0 ;
    expansion_normal_test(k, N, x0, xs, xstride, src, sstride, nsrc,
			  xf, nfld) ;
    return 0 ;
  }

  if ( test == 22 ) {
    expansion_gradient_test(k, N, x0, xs, xstride, src, sstride, nsrc,
			    xf, nfld) ;
    return 0 ;
  }

  if ( test == 23 ) {
    local_gradient_test(k, N, xc, x, xs, xstride, src, sstride, nsrc,
			xf, fstride, nfld) ;

    return 0 ;
  }

  
  return 0 ;
}
