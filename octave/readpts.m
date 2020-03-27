function [x,q,n,f]=readpts(file)

  fid = fopen(file, "r") ;

  dat = fscanf(fid, "%d", 2) ;
  cstr = strtrim(fscanf(fid, "%[^\n]c")) ;
  nq = dat(2) ;
  ns = dat(1) ;

  if ( strcmp(cstr, "M") )
    dat = fscanf(fid, "%f", ns*(3+nq)) ;
    
    dat = reshape(dat, 3+nq, ns)' ;

    x = dat(:, 1:3) ;
    q = dat(:, 4:end) ;
  endif

  if ( strcmp(cstr, "F") )
    dat = fscanf(fid, "%f", ns*(3+nq)) ;
    
    dat = reshape(dat, 3+nq, ns)' ;

    x = dat(:, 1:3) ;
    q = 0 ;
  endif
    
  fclose(fid) ;
  
