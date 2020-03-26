function dat = mksources(file, xmin, D, nsrc, nc)

  if ( nargin < 5 ) nc = 2 ; endif
  
fid = fopen(file, "w") ;

fprintf(fid, "%d %d\n", nsrc, nc+3) ;

dat = [D*rand(nsrc, 3) randn(nsrc,nc)] ;
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