function P=nlegendre(n, x)


  P = legendre(n, x) ;

  m = (0:n)' ;
  P .*= (-1).^m .*sqrt(gamma(n-abs(m)+1)./gamma(n+abs(m)+1)*(2*n+1)/4/pi);
