x = 0.4 ; y = 0.9 ; z = 3.1 ;  dz = 1e-6 ;

ee = [] ;
for n = 0:4
  for m = -n:n

    S = sfunc(n, m, x, y, z) ;
    dS = (sfunc(n, m, x+0.5*dz, y, z) - sfunc(n, m, x-0.5*dz, y, z))/dz ;
    ##dS = (sfunc(n, m, x, y, z+0.5*dz) - sfunc(n, m, x, y, z-0.5*dz))/dz ;
    ##dS = (sfunc(n, m, x, y, z+0.5*dz) - sfunc(n, m, x, y, z-0.5*dz))/dz ;

    d1 = -sfunc(n+1, m+1, x, y, z)*sqrt((2*n+1)/(2*n+3)*(n+m+1)*(n+m+2)) ;
    d2 =  sfunc(n+1, m-1, x, y, z)*sqrt((2*n+1)/(2*n+3)*(n-m+1)*(n-m+2)) ;
    if ( m <= 0 ) d2 = -d2 ; endif
    if ( m < 0) d1 = -d1 ; endif
    ds = 0.5*(d2 + d1) ;
    ##ds = -sfunc(n+1, m, x, y, z)*sqrt((2*n+1)/(2*n+3)*(n+m+1)*(n-m+1)) ;    
    ##ds = -sfunc(n+1, m, x, y, z)*sqrt((2*n+1)/(2*n+3)*(n+m+1)*(n-m+1)) ;
    ee = [ee; abs(dS-ds)] ;
  endfor
endfor
