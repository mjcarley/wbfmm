ph = 0.7 ; mu = -0.4 ; r = 1.9 ;
n = 3 ; m = -1 ;

d = 1e-4 ;

Rnm = rnmfunc(n, m, r, ph, mu) ;

ee = id = rr = [] ;
for n=1:4
  for m=-n:n

    drdr = ...
    (rnmfunc(n, m, r+0.5*d, ph, mu) - rnmfunc(n, m, r-0.5*d, ph, mu))/d ;
    drdmu = ...
    (rnmfunc(n, m, r, ph, mu+0.5*d) - rnmfunc(n, m, r, ph, mu-0.5*d))/d ;
    drdph = ...
    (rnmfunc(n, m, r, ph+0.5*d, mu) - rnmfunc(n, m, r, ph-0.5*d, mu))/d ;

    dxy = ((1-mu^2)*(r*drdr - mu*drdmu) + j*drdph)/r/sqrt(1-mu^2)*exp(j*ph) ;
    dxyb = ((1-mu^2)*(r*drdr - mu*drdmu) - j*drdph)/r/sqrt(1-mu^2)*exp(-j*ph) ;

    if 1
      if ( m >= 0 ) 
	tt = -sqrt((2*n+1)/(2*n-1)*(n-abs(m))*(n-abs(m)-1))*...
	      rnmfunc(n-1, m+1, r, ph, mu) ;
      else
	tt = sqrt((2*n+1)/(2*n-1)*(n+abs(m))*(n+abs(m)-1))*...
	     rnmfunc(n-1, m+1, r, ph, mu) ;
      endif
      ee = [ee; abs(dxy-tt)] ;
      id = [id; n m] ;
      rr = [rr; tt/dxy] ;
    endif

    if 0
      if ( m > 0 ) 
	tt = sqrt((2*n+1)/(2*n-1)*(n+abs(m))*(n+abs(m)-1))*...
	     rnmfunc(n-1, m-1, r, ph, mu) ;
      else
	tt = -sqrt((2*n+1)/(2*n-1)*(n-abs(m))*(n-abs(m)-1))*...
	      rnmfunc(n-1, m-1, r, ph, mu) ;
      endif
      ee = [ee; abs(dxyb-tt)] ;
      id = [id; n m] ;
      rr = [rr; tt/dxyb] ;
    endif
    endfor
endfor

rr = real(rr) ;
