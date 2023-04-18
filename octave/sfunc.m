function S=sfunc(n, m, x, y, z)

r = sqrt(x.^2 + y.^2 + z.^2) ;
th = acos(z/r) ;
ph = atan2(y, x) ;

C = cos(th) ;

if ( n < abs(m) )
  S = 0 ;
  return ;
endif

P = nlegendre(n, C) ;
S = exp(j*m*ph)*P(abs(m)+1)/r^(n+1) ;


