function dat = mkfpoints(file, xmin, D, nsrc)

fid = fopen(file, "w") ;

fprintf(fid, "%d 3\n", nsrc) ;

dat = D*rand(nsrc, 3) ;
dat(:,1) += xmin(1) ;
dat(:,2) += xmin(2) ;
dat(:,3) += xmin(3) ;

fprintf(fid, "%f %f %f\n", dat') ;


fclose(fid) ;
