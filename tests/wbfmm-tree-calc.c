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
gchar *progname ;

gint parse_origin(gdouble *x, gchar *str) ;
gint read_points(gchar *file,
		 gdouble **xs, gint *xstr,
		 gdouble **q,  gint *qstr, gint *nq,
		 gdouble **n,  gint *nstr, 
		 gdouble **f,  gint *fstr,
		 gint *nsrc) ;

gint read_points(gchar *file,
		 gdouble **xs, gint *xstr,
		 gdouble **q,  gint *qstr, gint *nq,
		 gdouble **n,  gint *nstr, 
		 gdouble **f,  gint *fstr,
		 gint *nsrc)

{
  FILE *input = stdin ;
  gdouble *s ;
  gchar code[8] ;
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

gint parse_origin(gdouble *x, gchar *str)

{
  sscanf(str, "%lg,%lg,%lg", &(x[0]), &(x[1]), &(x[2])) ;
  
  return 0 ;
}

gint main(gint argc, gchar **argv)

{
  wbfmm_tree_t *tree ;
  wbfmm_shift_operators_t *shifts ;
  gdouble k, D, xtree[3] = {0.0}, xtmax[3], *xs ;
  gdouble del, *x, *work, *xf, *f, tol, *q, *normals, *dipoles ;
  gint nsrc, i, j, fstr, nf, xstr, qstr, nstr, dstr, fcstr, nq ;
  gsize pstr ;
  guint depth, order[48] = {0}, order_s, order_r, order_max, level ;
  guint sizew ;
  guint64 b ;
  gchar ch, *sfile = NULL, *ffile = NULL ;
  gboolean write_sources, fit_box, shift_bw ;
  wbfmm_library_config_t lconfig ;
  wbfmm_field_t field ;
  
  k = 1.0 ; D = 1.0 ; nsrc = 1 ; del = 1e-2 ; tol = 1e-6 ;
  field = WBFMM_FIELD_SCALAR ;
  nq = 1 ;
  depth = 2 ;
  xtree[0] = xtree[1] = xtree[2] = 0.0 ;
  /* order_s = 8 ; order_r = 8 ; */
  order_s = order_r = 0 ;
  order_max = 0 ;
  write_sources = FALSE ;
  fit_box = FALSE ;
  shift_bw = FALSE ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;
  timer = g_timer_new() ;

  while ( (ch = getopt(argc, argv, "hBbcD:d:f:gk:O:R:s:S:t:w")) != EOF ) {
    switch ( ch ) {
    default:
    case 'h':
      fprintf(stderr,
	      "Usage: %s <options> > <output file>\n\n"
	      "Compute field generated by list of sources at specified "
	      "field points using\n"
	      "Wide Band (some day) Fast Multipole Method\n\n"
	      "Options:\n\n"
	      "  -B use backward shift algorithm\n"
	      "  -b fit octree box to sources\n"
	      "  -c write the library configuration to stderr and exit\n"
	      "  -d # depth of octree (%d)\n"
	      "  -D # width of octree (%lg)\n"
	      "  -f (field point name)\n"
	      "  -g calculate gradient of field\n"
	      "  -k # wavenumber (%lg)\n"
	      "  -O #,#,# origin of octree (%lg,%lg,%lg)\n"
	      "  -R # order of regular expansions at leaf level (%u)\n"
	      "  -S # order of singular expansions at leaf level (%u)\n"
	      "  -s (source file name)\n"
	      "  -t # tolerance (%lg)\n"
	      "  -w write source data to stdout\n",
	      progname, depth, D, k, xtree[0], xtree[1], xtree[2],
	      order_r, order_s, tol) ;
      return 0 ;
      break ;
    case 'B': shift_bw = TRUE ; break ;
    case 'b': fit_box = TRUE ; break ;
    case 'c':
      wbfmm_library_config(&lconfig) ;  
      wbfmm_library_config_print(&lconfig, stderr) ;
      return 0 ;
      break ;
    case 'd': depth = atoi(optarg) ; break ;
    case 'D': D = atof(optarg) ; break ;
    case 'f': ffile = g_strdup(optarg) ; break ;      
    case 'g': field = WBFMM_FIELD_GRADIENT ; break ;
    case 'k': k = atof(optarg) ; break ;
    case 'O': parse_origin(xtree, optarg) ; break ;
    case 'R': order_r = atoi(optarg) ; break ;
    case 'S': order_s = atoi(optarg) ; break ;
    case 's': sfile = g_strdup(optarg) ; break ;
    case 't': tol = atof(optarg) ; break ;
    case 'w': write_sources = TRUE ; break ;
    }
  }

  if ( sfile != NULL ) {
    read_points(sfile,
		&xs, &xstr,
		&q,  &qstr, &nq,
		&normals, &nstr,
		&dipoles, &dstr,
		&nsrc) ;
  } else {
    fprintf(stderr, "%s: source list must be specified (-s)\n",
	    progname) ;
    return 1 ;
  }

  nq /= 2 ;
  
  if ( ffile != NULL ) {
    read_points(ffile,
		&xf, &fstr, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &nf) ;
  } else {
    fprintf(stderr, "%s: field point list must be specified (-f)\n",
  	    progname) ;
    return 1 ;
  }

  fcstr = 2*nq ;

  if ( fcstr | WBFMM_FIELD_GRADIENT ) fcstr *= 3 ;
  f = (gdouble *)g_malloc0(nf*fcstr*sizeof(gdouble)) ;
  
  if ( fit_box ) {
    wbfmm_points_origin_width(xs, xstr, nsrc, xtree, xtmax, &D, TRUE) ;
    wbfmm_points_origin_width(xf, fstr, nf, xtree, xtmax, &D, FALSE) ;

    xtree[0] -= del ; xtree[1] -= del ; xtree[2] -= del ;
    D += 2.0*del ;
  }

  pstr = xstr*sizeof(gdouble) ;
  tree = wbfmm_tree_new(xtree, D, 2*nsrc) ;
  wbfmm_tree_source_size(tree) = nq ;
  
  if ( order_s != 0 && order_r != 0 ) {
    order[2*depth+0] = order_s ; 
    order[2*depth+1] = order_r ; 
    order_max = MAX(order_s, order_r) ;
    for ( i = depth-1 ; i > 0 ; i -- ) {
      order[2*i+0] = order[2*(i+1)+0] + 4 ;
      order[2*i+1] = order[2*(i+1)+1] + 4 ;
      order_max = MAX(order_max, order[2*i+0]) ;
      order_max = MAX(order_max, order[2*i+1]) ;
    }
  } else {
    for ( i = 1 ; i <= depth ; i ++ ) {
      order[2*i+0] = order[2*i+1] =
	wbfmm_truncation_number(tree, k, i, tol) ;
      order_max = MAX(order_max, order[2*i+0]) ;
      order_max = MAX(order_max, order[2*i+1]) ;
    }
  }

  sizew = wbfmm_element_number_rotation(2*order_max) ;
  sizew = MAX(sizew, 64*nq*(wbfmm_coefficient_index_nm(order_max+1,0))) ;
  work = (gdouble *)g_malloc0(sizew*sizeof(gdouble)) ;

  x = wbfmm_tree_origin(tree) ;
  fprintf(stderr, "%s: box origin: %lg %lg %lg\n",
	  progname, x[0], x[1], x[2]) ;
  fprintf(stderr, "%s: box width: %lg\n",
	  progname, wbfmm_tree_width(tree)) ;
  for ( i = 1 ; i <= depth ; i ++ ) {
    fprintf(stderr, "%s: expansion order, level %d: singular %u; regular %u\n",
	    progname, i, order[2*i+0], order[2*i+1]) ;
  }
  
  fprintf(stderr, "%s: initializing shift rotation operators; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  wbfmm_shift_angle_table_init() ;

  shifts = wbfmm_shift_operators_new(order_max, shift_bw, work) ;

  fprintf(stderr, "%s: shift rotation operators initialized; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  fprintf(stderr, "%s: initializing shift translation operators; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  for ( i = 1 ; i <= depth ; i ++ ) {
    wbfmm_shift_operators_coaxial_SR_init(shifts, D, i, order[2*i+0], 
					       k, work) ;
  }

  fprintf(stderr, "%s: shift translation operators initialized; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  fprintf(stderr, "%s: initializing upward pass translation operators; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  for ( level = 2 ; level <= depth ; level ++ ) {
    wbfmm_shift_operators_coaxial_SS_init(shifts, D, level, 
					       order[2*(level-1)+0], 
					       k, work) ;
  }

  fprintf(stderr, "%s: upward pass translation operators initialized; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  wbfmm_tree_add_points(tree, (gpointer)xs, nsrc, pstr) ;

  for ( i = 0 ; i < depth ; i ++ ) wbfmm_tree_refine(tree) ;

  if ( write_sources ) {
    wbfmm_tree_write_sources(tree, &(xs[3]), xstr, stderr) ;
    
    return 0 ;
  }

  wbfmm_tree_problem(tree) = WBFMM_PROBLEM_HELMHOLTZ ;
  for ( i = 1 ; i <= depth ; i ++ ) {
    wbfmm_tree_coefficient_init(tree, i, order[2*i+1], order[2*i+0]) ;
  }

  fprintf(stderr, "%s: initializing leaf expansions; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;
  
  wbfmm_tree_leaf_expansions(tree, k,
				  q, qstr, normals, nstr, dipoles, dstr,
				  TRUE, work) ;
  
  fprintf(stderr, "%s: leaf expansions initialized; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  fprintf(stderr, "%s: upward pass; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  for ( level = depth ; level >= 3 ; level -- ) 
    wbfmm_upward_pass(tree, shifts, level, work) ;

  fprintf(stderr, "%s: upward pass completed; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  fprintf(stderr, "%s: downward pass; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  for ( level = 2 ; level <= depth ; level ++ ) {
    wbfmm_downward_pass(tree, shifts, level, work) ;
  }

  fprintf(stderr, "%s: downward pass completed; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  level = depth ;

  fprintf(stderr, "%s: computing fmm field; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  for ( i = 0 ; i < nf ; i ++ ) {
    b = wbfmm_point_box(tree, level, &(xf[i*fstr])) ;
    wbfmm_tree_box_local_field(tree, level, b, k, 
				    &(xf[i*fstr]), &(f[fcstr*i]), fcstr,
				    q, qstr, normals, nstr, dipoles, dstr,
				    TRUE, field, work) ;
  }

  fprintf(stderr, "%s: fmm field computed; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  for ( i = 0 ; i < nf ; i ++ ) {
    fprintf(stdout, 
	    "%1.16e %1.16e %1.16e",
	    xf[i*fstr+0], xf[i*fstr+1], xf[i*fstr+2]) ;
    for ( j = 0 ; j < fcstr ; j ++ ) {
      fprintf(stdout, " %1.16e", f[fcstr*i+j]) ;
    }
    fprintf(stdout, "\n") ;
  }
  
  return 0 ;
}
