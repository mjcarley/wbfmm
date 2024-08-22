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
char *progname ;

gint read_points(char *file,
		 gdouble **xs, gint *xstr,
		 gdouble **q,  gint *qstr, gint *nq,
		 gdouble **n,  gint *nstr, 
		 gdouble **f,  gint *fstr,
		 gint *nsrc) ;

gint read_points(char *file,
		 gdouble **xs, gint *xstr,
		 gdouble **q,  gint *qstr, gint *nq,
		 gdouble **n,  gint *nstr, 
		 gdouble **f,  gint *fstr,
		 gint *nsrc)

{
  FILE *input = stdin ;
  gdouble *s ;
  char code[8] ;
  gint i, j, nqt ;

  if ( file != NULL ) {
    input = fopen(file, "r") ;
    if ( input == NULL ) {
      fprintf(stderr, "%s: cannot open file %s\n", progname, file) ;
      exit(1) ;
    }
  }

  fscanf(input, "%d", nsrc) ;
  fscanf(input, "%d", &nqt) ;
  fscanf(input, "%s", code) ;
  fprintf(stderr, "%s: %d point%c of %d components\n", 
	  progname, *nsrc, (*nsrc > 1 ? 's' : ' '), nqt) ;
  if ( strcmp(code, "M") == 0) {
    *xs = *q = *n = *f = NULL ;

    *nq = nqt ;
    *xstr = 3 + *nq ;
    s = *xs = (gdouble *)g_malloc0((*xstr)*(*nsrc)*sizeof(gdouble)) ;

    for ( i = 0 ; i < *nsrc ; i ++ ) {
      for ( j = 0 ; j < *xstr ; j ++ )
	fscanf(input, "%lg", &(s[(*xstr)*i+j])) ;
    }

    *q = &(s[3]) ; *qstr = *xstr ;
    
    if ( file != NULL ) fclose(input) ;

    return 0 ;
  }

  if ( strcmp(code, "MN") == 0) {
    /*monopole and normals, order:  */
    /* x y z q0 q1 ... qn nx ny nz f0 f1 ... fn  */
    /*   */
    /*   */
    *xs = *q = *n = *f = NULL ;

    *nq = nqt ;
    *xstr = 3 + 3 + 2*(*nq) ;
    s = *xs = (gdouble *)g_malloc0((*xstr)*(*nsrc)*sizeof(gdouble)) ;

    for ( i = 0 ; i < *nsrc ; i ++ ) {
      for ( j = 0 ; j < *xstr ; j ++ )
	fscanf(input, "%lg", &(s[(*xstr)*i+j])) ;
    }

    *q = &(s[3])         ; *qstr = *xstr ;
    *n = &(s[3+(*nq)])   ; *nstr = *xstr ;
    *f = &(s[6+(*nq)])   ; *fstr = *xstr ;
    
    if ( file != NULL ) fclose(input) ;

    return 0 ;
  }

  if ( strcmp(code, "F") == 0) {
    *xstr = 3 ;
    s = *xs = (gdouble *)g_malloc0((*xstr)*(*nsrc)*sizeof(gdouble)) ;

    for ( i = 0 ; i < *nsrc ; i ++ ) {
      for ( j = 0 ; j < *xstr ; j ++ )
	fscanf(input, "%lg", &(s[(*xstr)*i+j])) ;
    }

    if ( file != NULL ) fclose(input) ;

    return 0 ;
  }

  g_assert_not_reached() ;
  
  return 0 ;
}

gint main(gint argc, char **argv)

{
  gdouble k, *xs ;
  gdouble *xf, *f, *q, *normals, *dipoles ;
  gint nsrc, i, j, xstr, fstr, nf, qstr, dstr, nstr, fcstr, nq ;
  char ch, *sfile = NULL, *ffile = NULL ;
  gboolean gradient ;

  k = 1.0 ; nq = 1 ;
  gradient = FALSE ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;
  timer = g_timer_new() ;

  while ( (ch = getopt(argc, argv, "hf:gk:s:")) != EOF ) {
    switch ( ch ) {
    default:
    case 'h':
      fprintf(stderr,
	      "Usage: %s <options> > <output file>\n\n"
	      "Compute field generated by list of sources at specified "
	      "field points using\n"
	      "direct evaluation method\n"
	      "Options:\n\n"
	      "  -f (field point name)\n"
	      "  -g calculate gradient of field\n"
	      "  -k # wavenumber (%lg)\n"
	      "  -s (source file name)\n",
	      progname, k) ;
      return 0 ;
      break ;
    case 'f': ffile = g_strdup(optarg) ; break ;      
    case 'g': gradient = TRUE ; break ;
    case 'k': k = atof(optarg) ; break ;
    case 's': sfile = g_strdup(optarg) ; break ;
    }
  }

  if ( sfile != NULL ) {
    read_points(sfile,
		&xs, &xstr,
		&q,  &qstr, &nq,
		&normals, &nstr,
		&dipoles, &dstr,
		&nsrc) ;
    /* read_points(sfile, &sources, &nsrc, &str) ; */
  } else {
    fprintf(stderr, "%s: source list must be specified (-s)\n",
	    progname) ;
    return 1 ;
  }

  nq /= 2 ;
  fcstr = 2*nq ;
  if ( gradient ) fcstr *= 3 ;
  
  if ( ffile != NULL ) {
    /* read_points(ffile, &xf, &nf, &strf) ; */
    read_points(ffile,
		&xf, &fstr, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &nf) ;
  } else {
    fprintf(stderr, "%s: field point list must be specified (-f)\n",
	    progname) ;
    return 1 ;
  }

  f = (gdouble *)g_malloc0(nf*fcstr*sizeof(gdouble)) ;

  fprintf(stderr, "%s: computing direct field; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  if ( gradient )
    for ( i = 0 ; i < nf ; i ++ ) {
      wbfmm_total_field_grad(k, xs, xstr,
				  q, qstr, normals, nstr, dipoles, dstr,
				  nq, nsrc,
				  &(xf[i*fstr]), &(f[fcstr*i]), 6) ;
    }
  else
    for ( i = 0 ; i < nf ; i ++ ) {
      wbfmm_total_field(k, xs, xstr,
			     q, qstr, normals, nstr, dipoles, dstr, nq,
			     nsrc,
			     &(xf[i*fstr]), &(f[fcstr*i]), nq) ;
    }

  fprintf(stderr, "%s: direct field computed; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  for ( i = 0 ; i < nf ; i ++ ) {
    fprintf(stdout, 
	    "%1.16e %1.16e %1.16e ",
	    xf[i*fstr+0], xf[i*fstr+1], xf[i*fstr+2]) ;
    for ( j = 0 ; j < fcstr/2 ; j ++ ) {
      fprintf(stdout, 
	      " %1.16e %1.16e", f[fcstr*i+2*j+0], f[fcstr*i+2*j+1]) ;
    }
    fprintf(stdout, "\n") ;
  }

  return 0 ;
}