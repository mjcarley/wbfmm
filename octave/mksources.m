function dat = mksources(file, xmin, D, nsrc, nc, code)

  if ( nargin < 5 ) nc = 2 ; endif
  if ( nargin < 6 ) code = "M" ; endif
  
if ( strcmp(code, "M") )
  fid = fopen(file, "w") ;

  fprintf(fid, "%d %d %s\n", nsrc, nc, code) ;

  ## monopoles
  q = randn(nsrc, nc) ;

  dat = [D*rand(nsrc, 3) q] ;
  dat(:,1) += xmin(1) ;
  dat(:,2) += xmin(2) ;
  dat(:,3) += xmin(3) ;

  fmt = "" ;
  for i=1:nc+3
    fmt = [fmt " %f"] ;
  endfor
  fmt = [fmt "\n"] ;

  fprintf(fid, fmt, dat') ;

  fclose(fid) ;

  return ;
endif

if ( strcmp(code, "MN") )
  fid = fopen(file, "w") ;

  fprintf(fid, "%d %d %s\n", nsrc, nc, code) ;

  ## monopoles
  q = randn(nsrc, nc) ;
  ##q *= 0 ;
  ## normals
  n = randn(nsrc, 3) ;
  ## normal dipole strengths
  f = randn(nsrc, nc) ;
  dat = [D*rand(nsrc, 3) q n f] ;
  dat(:,1) += xmin(1) ;
  dat(:,2) += xmin(2) ;
  dat(:,3) += xmin(3) ;

  fmt = "" ;
  for i=1:2*nc+6
    fmt = [fmt " %f"] ;
  endfor
  fmt = [fmt "\n"] ;
  
  fprintf(fid, fmt, dat') ;

  fclose(fid) ;

  return ;
endif

if ( strcmp(code, "N") )
  fid = fopen(file, "w") ;

  fprintf(fid, "%d %d %s\n", nsrc, nc, code) ;

  ## normals
  n = randn(nsrc, 3) ;
  ## normal dipole strengths
  f = randn(nsrc, nc) ;
  dat = [D*rand(nsrc, 3) n f] ;
  dat(:,1) += xmin(1) ;
  dat(:,2) += xmin(2) ;
  dat(:,3) += xmin(3) ;

  fmt = "" ;
  for i=1:2*nc+6
    fmt = [fmt " %f"] ;
  endfor
  fmt = [fmt "\n"] ;
  
  fprintf(fid, fmt, dat') ;

  fclose(fid) ;

  return ;
endif

error(["Unrecognized code " code]) ;
