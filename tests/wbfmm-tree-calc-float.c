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

gint parse_origin(gfloat *x, gchar *str) ;
gint read_points(gchar *file, gfloat **points, gint *nsrc, gint *str) ;

gint read_points(gchar *file, gfloat **points, gint *nsrc, gint *str)

{
  FILE *input = stdin ;
  gfloat *s ;
  gint i, j ;

  if ( file != NULL ) {
    input = fopen(file, "r") ;
    if ( input == NULL ) {
      fprintf(stderr, "%s: cannot open file %s\n", progname, file) ;
      exit(1) ;
    }
  }

  fscanf(input, "%d", nsrc) ;
  fscanf(input, "%d", str) ;
  fprintf(stderr, "%s: %d point%c\n", 
	  progname, *nsrc, (*nsrc > 1 ? 's' : ' ')) ;
  s = *points = (gfloat *)g_malloc0((*str)*(*nsrc)*sizeof(gfloat)) ;

  for ( i = 0 ; i < *nsrc ; i ++ ) {
    for ( j = 0 ; j < *str ; j ++ ) 
      fscanf(input, "%g", &(s[(*str)*i+j])) ;
  }

  if ( file != NULL ) fclose(input) ;

  return 0 ;
}

gint parse_origin(gfloat *x, gchar *str)

{
  sscanf(str, "%g,%g,%g", &(x[0]), &(x[1]), &(x[2])) ;
  
  return 0 ;
}

gint main(gint argc, gchar **argv)

{
  wbfmm_tree_t *tree ;
  wbfmm_shift_operators_t *shifts ;
  gfloat k, D, xtree[3] = {0.0}, xtmax[3], *sources ;
  gfloat del, *x, *work, *xf, *f, tol, *q, *normals, *dipoles ;
  gint nsrc, i, str, strf, nf, qstr, nstr, dstr ;
  gsize pstr ;
  guint depth, order[48] = {0}, order_s, order_r, order_max, level ;
  guint sizew ;
  guint64 b ;
  gchar ch, *sfile = NULL, *ffile = NULL ;
  gboolean write_sources, fit_box, shift_bw ;
  wbfmm_library_config_t lconfig ;
  
  k = 1.0 ; D = 1.0 ; nsrc = 1 ; del = 1e-2 ; tol = 1e-6 ;
  depth = 2 ; str = 5 ;
  xtree[0] = xtree[1] = xtree[2] = 0.0 ;
  /* order_s = 8 ; order_r = 8 ; */
  order_s = order_r = 0 ;
  order_max = 0 ;
  write_sources = FALSE ;
  fit_box = FALSE ;
  shift_bw = FALSE ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;
  timer = g_timer_new() ;

  while ( (ch = getopt(argc, argv, "hBbcD:d:f:k:O:R:s:S:t:w")) != EOF ) {
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
	      "  -D # width of octree (%g)\n"
	      "  -f (field point name)\n"
	      "  -k # wavenumber (%g)\n"
	      "  -O #,#,# origin of octree (%g,%g,%g)\n"
	      "  -R # order of regular expansions at leaf level (%u)\n"
	      "  -S # order of singular expansions at leaf level (%u)\n"
	      "  -s (source file name)\n"
	      "  -t # tolerance (%g)\n"
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
    read_points(sfile, &sources, &nsrc, &str) ;
    switch ( str ) {
    default:
      fprintf(stderr, "%s: don't know how to interpret stride %d\n",
	      progname, str) ;
      exit(1) ;
      break ;
    case  5: fprintf(stderr, "%s: monopole sources\n", progname) ; break ;
    case  8: fprintf(stderr, "%s: dipole sources\n", progname) ; break ;
    case 10: fprintf(stderr, "%s: mixed sources\n", progname) ; break ;
    }
  } else {
    fprintf(stderr, "%s: source list must be specified (-s)\n",
	    progname) ;
    return 1 ;
  }

  if ( ffile != NULL ) {
    read_points(ffile, &xf, &nf, &strf) ;
  } else {
    fprintf(stderr, "%s: field point list must be specified (-f)\n",
	    progname) ;
    return 1 ;
  }

  f = (gfloat *)g_malloc0(nf*2*sizeof(gfloat)) ;

  if ( fit_box ) {
    wbfmm_points_origin_width_f(sources, str, nsrc, xtree, xtmax, &D, TRUE) ;
    wbfmm_points_origin_width_f(xf, strf, nf, xtree, xtmax, &D, FALSE) ;

    xtree[0] -= del ; xtree[1] -= del ; xtree[2] -= del ;
    D += 2.0*del ;
  }
  
  pstr = str*sizeof(gfloat) ;
  tree = wbfmm_tree_new_f(xtree, D, 2*nsrc) ;
  wbfmm_tree_source_size(tree) = 1 ;
  
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
	wbfmm_truncation_number_f(tree, k, i, tol) ;
      order_max = MAX(order_max, order[2*i+0]) ;
      order_max = MAX(order_max, order[2*i+1]) ;
    }
  }

  sizew = wbfmm_element_number_rotation(2*order_max) ;
  sizew = MAX(sizew, 16*(wbfmm_coefficient_index_nm(order_max+1,0))) ;
  work = (gfloat *)g_malloc0(sizew*sizeof(gfloat)) ;

  x = wbfmm_tree_origin(tree) ;
  fprintf(stderr, "%s: box origin: %g %g %g\n",
	  progname, x[0], x[1], x[2]) ;
  fprintf(stderr, "%s: box width: %g\n",
	  progname, wbfmm_tree_width(tree)) ;
  for ( i = 1 ; i <= depth ; i ++ ) {
    fprintf(stderr, "%s: expansion order, level %d: singular %u; regular %u\n",
	    progname, i, order[2*i+0], order[2*i+1]) ;
  }
  
  fprintf(stderr, "%s: initializing shift rotation operators; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  wbfmm_shift_angle_table_init_f() ;

  shifts = wbfmm_shift_operators_new_f(order_max, shift_bw, work) ;

  fprintf(stderr, "%s: shift rotation operators initialized; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  fprintf(stderr, "%s: initializing shift translation operators; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  for ( i = 1 ; i <= depth ; i ++ ) {
    wbfmm_shift_operators_coaxial_SR_init_f(shifts, D, i, order[2*i+0], 
					       k, work) ;
  }

  fprintf(stderr, "%s: shift translation operators initialized; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  fprintf(stderr, "%s: initializing upward pass translation operators; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  for ( level = 2 ; level <= depth ; level ++ ) {
    wbfmm_shift_operators_coaxial_SS_init_f(shifts, D, level, 
					       order[2*(level-1)+0], 
					       k, work) ;
  }

  fprintf(stderr, "%s: upward pass translation operators initialized; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  wbfmm_tree_add_points_f(tree, (gpointer)sources, nsrc, pstr) ;

  for ( i = 0 ; i < depth ; i ++ ) wbfmm_tree_refine_f(tree) ;

  if ( write_sources ) {
    wbfmm_tree_write_sources_f(tree, &(sources[3]), str, stderr) ;
    
    return 0 ;
  }

  wbfmm_tree_problem(tree) = WBFMM_PROBLEM_HELMHOLTZ ;
  for ( i = 1 ; i <= depth ; i ++ ) {
    wbfmm_tree_coefficient_init_f(tree, i, order[2*i+1], order[2*i+0]) ;
  }

  fprintf(stderr, "%s: initializing leaf expansions; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  qstr = nstr = dstr = 0 ;
  q = normals = dipoles = NULL ;
  switch ( str ) {
  default: g_assert_not_reached() ; break ;
  case 5:
    q = &(sources[3]) ; qstr = str ;
    break ;
  case 8:
    normals = &(sources[3]) ; nstr = str ;
    dipoles = &(sources[6]) ; dstr = str ;
    break ;
  case 10:
    q = &(sources[3]) ; qstr = str ;
    normals = &(sources[5]) ; nstr = str ;
    dipoles = &(sources[8]) ; dstr = str ;
    break ;
  }
  
  wbfmm_tree_leaf_expansions_f(tree, k,
				  q, qstr, normals, nstr, dipoles, dstr,
				  TRUE, work) ;
  
  fprintf(stderr, "%s: leaf expansions initialized; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  fprintf(stderr, "%s: upward pass; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  for ( level = depth ; level >= 3 ; level -- ) 
    wbfmm_upward_pass_f(tree, shifts, level, work) ;

  fprintf(stderr, "%s: upward pass completed; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  fprintf(stderr, "%s: downward pass; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  for ( level = 2 ; level <= depth ; level ++ ) {
    wbfmm_downward_pass_f(tree, shifts, level, work) ;
  }

  fprintf(stderr, "%s: downward pass completed; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  level = depth ;

  fprintf(stderr, "%s: computing fmm field; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  for ( i = 0 ; i < nf ; i ++ ) {
    b = wbfmm_point_box_f(tree, level, &(xf[i*strf])) ;
    wbfmm_tree_box_local_field_f(tree, level, b, k, 
				    &(xf[i*strf]), &(f[2*i]),
				    q, qstr, normals, nstr, dipoles, dstr,
				    TRUE, work) ;
  }

  fprintf(stderr, "%s: fmm field computed; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  for ( i = 0 ; i < nf ; i ++ ) {
    fprintf(stdout, 
	    "%1.16e %1.16e %1.16e %1.16e %1.16e\n",
	    xf[i*strf+0], xf[i*strf+1], xf[i*strf+2], f[2*i+0], f[2*i+1]) ;
  }

  return 0 ;
}
