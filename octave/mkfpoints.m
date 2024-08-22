function dat = mkfpoints(file, xmin, D, nf)

  ## MKFPOINTS(FILE, XMIN, D, NF)
  ##
  ## Generate a field point file for testing FMM solvers
  ##
  ## FILE: name of file for input to FMM solver
  ## XMIN: three vector for coordinates of origin of bounding box of
  ##       field points
  ## D:    side length of bounding box
  ## NF:   number of points in file
  ##
  ## The output file can be supplied to the FMM test codes using the -f
  ## command line option

  fid = fopen(file, "w") ;

  fprintf(fid, "%d 0 F\n", nf) ;

  dat = D*rand(nf, 3) ;
  dat(:,1) += xmin(1) ;
  dat(:,2) += xmin(2) ;
  dat(:,3) += xmin(3) ;
  
  fprintf(fid, "%f %f %f\n", dat') ;

  fclose(fid) ;
