function f=rnmfunc(n,m,r,ph,mu)


  if n < abs(m)
    f = 0 ; ##*(0:abs(m)) ;
    return ;
  endif
  
  P = nlegendre(n, mu) ;
  P = P(abs(m)+1) ;

  f = P*r^n*exp(j*m*ph) ;
  
