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

#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include <glib.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

/* #ifdef HAVE_AVX_INSTRUCTIONS */
/* #define WBFMM_USE_AVX */
/* #endif /\*HAVE_AVX_INSTRUCTIONS*\/ */

/*include wbfmm.h here so it picks up the WBFMM_USE_AVX #define */

#include <wbfmm.h>

GTimer *timer ;
char *progname ;

gint parse_origin(gdouble *x, char *str) ;
gint read_points(char *file,
		 gdouble **xs, gint *xstr,
		 gdouble **q,  gint *qstr, gint *nq,
		 gdouble **n,  gint *nstr, 
		 gdouble **f,  gint *fstr,
		 gint *nsrc) ;

static gint print_longer_help(char *progname)

{
  fprintf(stderr, "%s: detailed notes\n\n", progname) ;

  fprintf(stderr,
	  "Input file formats:\n\n"
	  "Input files for source and field points take the same format\n"
	  "consisting of a header followed by point data:\n\n"
	  "[number of points] [number of source components] [type code]\n"
	  "x y z [source data]\n"
	  "x y z [source data]\n"
	  "...\n\n"
	  "The number of source components allows for multiple source\n"
	  "values at each source position, such as in vector problems.\n"
	  "The type code is `M' for monopoles and `F' for field points.\n"
	  "For field data, the number of source components can be zero.\n") ;
  return 0 ;
}

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

  if ( strcmp(code, "N") == 0) {
    *xs = *q = *n = *f = NULL ;

    *nq = nqt ;
    *xstr = 3 + *nq ;
    *nstr = 3 ;
    s = *xs = (gdouble *)g_malloc0((*xstr)*(*nsrc)*sizeof(gdouble)) ;
    *n = (gdouble *)g_malloc0(3*(*nsrc)*sizeof(gdouble)) ;
    *f = &(s[3]) ; *fstr = *xstr ;

    for ( i = 0 ; i < *nsrc ; i ++ ) {
      for ( j = 0 ; j < 3 ; j ++ ) {
	fscanf(input, "%lg", &((*xs)[(*xstr)*i+j])) ;
      }
      for ( j = 0 ; j < 3 ; j ++ ) {
	fscanf(input, "%lg", &((*n)[(*nstr)*i+j])) ;
      }
      for ( j = 0 ; j < (*nq) ; j ++ ) {
	fscanf(input, "%lg", &((*f)[(*fstr)*i+j])) ;
      }
    }
    
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
  
  g_assert_not_reached() ;
  
  return 0 ;
}

gint parse_origin(gdouble *x, char *str)

{
  sscanf(str, "%lg,%lg,%lg", &(x[0]), &(x[1]), &(x[2])) ;
  
  return 0 ;
}

gint main(gint argc, char **argv)

{
  wbfmm_tree_t *tree ;
  wbfmm_target_list_t *targets ;
  wbfmm_shift_operators_t *shifts ;
  gdouble D, xtree[3] = {0.0}, xtmax[3], *xs ;
  gdouble del, *x, *work, *xf, *f, tol, *q, *normals, *dipoles ;
  gint nsrc, nq, i, j, xstr, strf, nf, fstr, qstr, nstr, dstr, fcstr ;
  gint order_inc ;
  gsize pstr, pnstr ;
  guint depth, order[48] = {0}, order_s, order_r, order_max, level ;
  guint sizew, field, source ;
  char ch, *sfile = NULL, *ffile = NULL ;
  gboolean fit_box, shift_bw, target_list, sort_sources ;
  gint nthreads, nproc ;
  
  D = 1.0 ; nsrc = 1 ; del = 1e-2 ; tol = 1e-6 ;
  depth = 2 ;
  xtree[0] = xtree[1] = xtree[2] = 0.0 ;
  nthreads = 0 ;
  /* order_s = 8 ; order_r = 8 ; */
  order_s = order_r = 0 ; order_inc = 0 ;
  order_max = 0 ;
  fit_box = FALSE ;
  shift_bw = FALSE ;
  target_list = FALSE ;
  sort_sources = FALSE ;
  field = WBFMM_FIELD_SCALAR ;
  source = 0 ;
  targets = NULL ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;
  timer = g_timer_new() ;

  nproc = 1 ;
#ifdef _OPENMP
  nproc = g_get_num_processors() ;
#endif
  
  while ( (ch = getopt(argc, argv, "hHBbcD:d:f:gi:lO:pR:s:S:T:t:")) != EOF ) {
    switch ( ch ) {
    default:
    case 'h':
      fprintf(stderr,
	      "Usage: %s <options> > <output file>\n\n"
	      "Compute field generated by list of sources at specified "
	      "field points using\n"
	      "Wide Band (some day) Fast Multipole Method\n\n"
	      "Options:\n\n"
	      "  -h print this message and exit\n"
	      "  -H print some more detailed information and exit\n"
	      "  -B use backward shift algorithm\n"
	      "  -b fit octree box to sources\n"
	      "  -c calculate curl of vector field\n"
	      "  -d # depth of octree (%d)\n"
	      "  -D # width of octree (%lg)\n"
	      "  -f (field point name)\n"
	      "  -g calculate gradient of field\n"
	      "  -i increment of order with level (%d)\n"
	      "  -l use target lists to calculate field at points\n"
	      "  -O #,#,# origin of octree (%lg,%lg,%lg)\n"
	      "  -p sort source points before generating tree\n"
	      "  -R # order of regular expansions at leaf level (%u)\n"
	      "  -S # order of singular expansions at leaf level (%u)\n"
	      "  -s (source file name)\n"
	      "  -T # (number of threads)\n"
	      "  -t # tolerance (%lg)\n",
	      progname, depth, D, order_inc, xtree[0], xtree[1], xtree[2],
	      order_r, order_s, tol) ;
      return 0 ;
      break ;
    case 'H': print_longer_help(progname) ;  return 0 ; break ;
    case 'B': shift_bw = TRUE ; break ;
    case 'b': fit_box = TRUE ; break ;
    case 'c': field = WBFMM_FIELD_CURL ; break ;
    case 'd': depth = atoi(optarg) ; break ;
    case 'D': D = atof(optarg) ; break ;
    case 'f': ffile = g_strdup(optarg) ; break ;      
    case 'g': field = WBFMM_FIELD_GRADIENT ; break ;
    case 'i': order_inc = atoi(optarg) ; break ;
    case 'l': target_list = TRUE ; break ;
    case 'O': parse_origin(xtree, optarg) ; break ;
    case 'p': sort_sources = TRUE ; break ;
    case 'R': order_r = atoi(optarg) ; break ;
    case 'S': order_s = atoi(optarg) ; break ;
    case 's': sfile = g_strdup(optarg) ; break ;
    case 'T': nthreads = atoi(optarg) ; break ;
    case 't': tol = atof(optarg) ; break ;
    }
  }

  fprintf(stderr, "%s: using %d threads (%d processes)\n",
	  progname, nthreads, nproc) ;
  
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

  if ( ffile != NULL ) {
    read_points(ffile,
		&xf, &strf, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &nf) ;
  } else {
    fprintf(stderr, "%s: field point list must be specified (-f)\n",
	    progname) ;
    return 1 ;
  }

  fcstr = nq ;
  if ( field == WBFMM_FIELD_GRADIENT ) fcstr *= 3 ;
  
  f = (gdouble *)g_malloc0(nf*fcstr*sizeof(gdouble)) ;

  /*fitting the bounding box to the field and source points*/
  if ( fit_box ) {
    wbfmm_points_origin_width(xs, xstr, nsrc, xtree, xtmax, &D, TRUE) ;
    wbfmm_points_origin_width(xf, strf, nf, xtree, xtmax, &D, FALSE) ;

    xtree[0] -= del ; xtree[1] -= del ; xtree[2] -= del ;
    D += 2.0*del ;
  }

  /*data strides and tree allocation*/
  pstr  = xstr*sizeof(gdouble) ;
  pnstr = nstr*sizeof(gdouble) ;
  fstr  = strf*sizeof(gdouble) ;
  tree  = wbfmm_tree_new(xtree, D, 2*nsrc) ;

  if ( sort_sources ) 
    wbfmm_tree_sort_points(tree, xs, pstr, nsrc) ;
  
  /*set expansion orders at each level of the tree*/
  if ( order_s != 0 && order_r != 0 ) {
    order[2*depth+0] = order_s ; 
    order[2*depth+1] = order_r ; 
    order_max = MAX(order_s, order_r) ;
    for ( i = depth-1 ; i > 0 ; i -- ) {
      order[2*i+0] = order[2*(i+1)+0] + order_inc ;
      order[2*i+1] = order[2*(i+1)+1] + order_inc ;
      /* order[2*i+1] = order[2*(i+1)+1] ; */
      order_max = MAX(order_max, order[2*i+0]) ;
      order_max = MAX(order_max, order[2*i+1]) ;
    }
  } else {
    for ( i = 1 ; i <= depth ; i ++ ) {
      order[2*i+0] = order[2*i+1] = 8 ;
      order_max = MAX(order_max, order[2*i+0]) ;
      order_max = MAX(order_max, order[2*i+1]) ;
    }
  }

  /*size the workspace: this needs to be tighter to account for
    Laplace case*/
  sizew = wbfmm_element_number_rotation(2*order_max) ;
  sizew = MAX(sizew, (order_max+1)*(order_max+1)*nq*16) ;
  work = (gdouble *)g_malloc0(2*sizew*sizeof(gdouble)) ;

  fprintf(stderr, "%s: %d elements allocated for workspace\n",
	  progname, sizew) ;
  
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

  fprintf(stderr, "%s: initializing coaxial translation coefficients; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;
  wbfmm_laplace_coaxial_translate_init(order_max+2) ;
  fprintf(stderr, "%s: coaxial translation coefficients initialized; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  /*add source points to the tree and refine to allocate sources to
    leaf boxes*/
  wbfmm_tree_add_points(tree, (gpointer)xs, pstr, normals, pnstr,
			     nsrc, sort_sources) ;
  for ( i = 0 ; i < depth ; i ++ ) wbfmm_tree_refine(tree) ;
  
  /*initialize memory for box coefficients at each level*/
  wbfmm_tree_problem(tree) = WBFMM_PROBLEM_LAPLACE ;
  wbfmm_tree_source_size(tree) = nq ;
  for ( i = 1 ; i <= depth ; i ++ ) {
    wbfmm_tree_laplace_coefficient_init(tree, i,
					     order[2*i+1], order[2*i+0]) ;
  }

  if ( target_list ) {
    g_assert(field != WBFMM_FIELD_CURL) ;
    fprintf(stderr, "%s: initializing target point list; %lg\n",
	    progname, g_timer_elapsed(timer, NULL)) ;
    if ( q != NULL ) source |= WBFMM_SOURCE_MONOPOLE ;
    if ( normals != NULL ) source |= WBFMM_SOURCE_DIPOLE ;
    targets = wbfmm_target_list_new(tree, nf) ;
    wbfmm_target_list_add_points(targets, xf, fstr, nf) ;
    wbfmm_target_list_coefficients_init(targets, field) ;
    wbfmm_laplace_target_list_local_coefficients(targets, source, work) ;
    fprintf(stderr, "%s: target point list initialized; %lg\n",
	    progname, g_timer_elapsed(timer, NULL)) ;
  }
  
  fprintf(stderr, "%s: initializing leaf expansions; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;  
  wbfmm_tree_laplace_leaf_expansions(tree,
					   q, qstr,
					   /* normals, nstr, */
					   dipoles, dstr,
					   TRUE, work) ;  
  fprintf(stderr, "%s: leaf expansions initialized; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  fprintf(stderr, "%s: upward pass; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;
  for ( level = depth ; level >= 3 ; level -- ) {
    wbfmm_laplace_upward_pass(tree, shifts, level, work) ;
  }  
  fprintf(stderr, "%s: upward pass completed; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  fprintf(stderr, "%s: downward pass; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;
  for ( level = 2 ; level <= depth ; level ++ ) {
    wbfmm_laplace_downward_pass(tree, shifts, level, work, nthreads) ;
  }
  fprintf(stderr, "%s: downward pass completed; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  fprintf(stderr, "%s: computing fmm field; %lg\n",
  	  progname, g_timer_elapsed(timer, NULL)) ;
  if ( target_list ) {
    wbfmm_target_list_local_field(targets, q, qstr, dipoles, dstr,
				       f, fcstr) ;
  } else {
    switch ( field ) {
    default: g_assert_not_reached() ; break ;
    case WBFMM_FIELD_SCALAR:
      for ( i = 0 ; i < nf ; i ++ ) {
	guint64 box ;
	box = wbfmm_point_box(tree, tree->depth, &(xf[i*strf])) ;
	wbfmm_tree_laplace_box_local_field(tree, tree->depth, box,
						&(xf[i*strf]),
						 &(f[i*fcstr]), 1,
						 q, qstr,
						 dipoles, dstr, TRUE, work) ;
      }
      break ;
    case WBFMM_FIELD_GRADIENT:
      for ( i = 0 ; i < nf ; i ++ ) {
	guint64 box ;
	box = wbfmm_point_box(tree, tree->depth, &(xf[i*strf])) ;
	wbfmm_tree_laplace_box_local_grad(tree, tree->depth, box,
					       &(xf[i*strf]),
					       &(f[i*fcstr]), 3, q, qstr,
					       /* normals, nstr, */
					       dipoles, dstr, TRUE, work) ;
      }
      break ;
    case WBFMM_FIELD_CURL:
      for ( i = 0 ; i < nf ; i ++ ) {
	guint64 box ;
	box = wbfmm_point_box(tree, tree->depth, &(xf[i*strf])) ;
	wbfmm_tree_laplace_box_local_curl(tree, tree->depth, box,
					       &(xf[i*strf]),
					       &(f[i*fcstr]), 3, q, qstr,
					       /* normals, nstr, */
					       dipoles, dstr, TRUE, work) ;
      }
      break ;
    }
  }

  fprintf(stderr, "%s: fmm field computed; %lg\n",
  	  progname, g_timer_elapsed(timer, NULL)) ;

  if ( field != WBFMM_FIELD_CURL ) {
    for ( i = 0 ; i < nf ; i ++ ) {
      fprintf(stdout, 
	      "%1.16e %1.16e %1.16e",
	      xf[i*strf+0], xf[i*strf+1], xf[i*strf+2]) ;
      for ( j = 0 ; j < fcstr ; j ++ ) {
	fprintf(stdout, " %1.16e", f[i*fcstr+j]) ;
      }
      fprintf(stdout, "\n") ;
    }
  } else {
    for ( i = 0 ; i < nf ; i ++ ) {
      fprintf(stdout, 
	      "%1.16e %1.16e %1.16e",
	      xf[i*strf+0], xf[i*strf+1], xf[i*strf+2]) ;
      for ( j = 0 ; j < 3 ; j ++ ) {
	fprintf(stdout, " %1.16e", f[i*fcstr+j]) ;
      }
      fprintf(stdout, "\n") ;
    }
  }

  g_free(ffile) ; g_free(sfile) ; g_free(progname) ;
  g_free(xf) ; g_free(work) ;
  
  return 0 ;
}
