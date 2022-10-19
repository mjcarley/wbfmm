ph = 0.7 ; mu = -0.4 ; r = 1.9 ;
n = 3 ; m = -1 ;

d = 1e-6 ;

Snm = snmfunc(n, m, r, ph, mu) ;

ee = id = rr = [] ;
for n=0:5
  for m=-n:n

    dSdr = ...
    (snmfunc(n, m, r+0.5*d, ph, mu) - snmfunc(n, m, r-0.5*d, ph, mu))/d ;
    dSdmu = ...
    (snmfunc(n, m, r, ph, mu+0.5*d) - snmfunc(n, m, r, ph, mu-0.5*d))/d ;
    dSdph = ...
    (snmfunc(n, m, r, ph+0.5*d, mu) - snmfunc(n, m, r, ph-0.5*d, mu))/d ;

    dxy = ((1-mu^2)*(r*dSdr - mu*dSdmu) + j*dSdph)/r/sqrt(1-mu^2)*exp(j*ph) ;
    dxyb = ((1-mu^2)*(r*dSdr - mu*dSdmu) - j*dSdph)/r/sqrt(1-mu^2)*exp(-j*ph) ;

    #if 0
    tt = -sqrt((2*n+1)/(2*n+3)*(n+(m)+2)*(n+(m)+1))*...
	    snmfunc(n+1, m+1, r, ph, mu) ;
    if ( m < 0 ) tt = -tt ; endif
    ee = [ee; abs(dxy-tt)] ;
    id = [id; n m] ;

    rr = [rr; dxy/tt] ;
    #endif

    if 0
    if ( m > 0 )
      tt = sqrt((2*n+1)/(2*n+3)*(n-(m)+2)*(n-(m)+1))*...
	   snmfunc(n+1, m-1, r, ph, mu) ;
    else
      tt = -sqrt((2*n+1)/(2*n+3)*(n-(m)+2)*(n-(m)+1))*...
	   snmfunc(n+1, m-1, r, ph, mu) ;
    endif
    
    ee = [ee; abs(dxyb-tt)] ;
    id = [id; n m] ;

    rr = [rr; dxyb/tt] ;
    endif
    
  endfor
endfor

rr = real(rr) ;
