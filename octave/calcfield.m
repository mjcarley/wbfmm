function p=calcfield(xs, ns, qs, fs, xf)

  nf = size(xf, 1) ;
  p = zeros(nf, 1) ;

  for i=1:nf
    r = [xf(i,1)-xs(:,1) xf(i,2)-xs(:,2) xf(i,3)-xs(:,3)] ;
    R = sqrt(sum(r.^2, 2)) ;
    p(i) = sum(qs./R) ;

    rn = sum(r.*ns, 2) ;
    p(i) += sum(fs.*rn./R.^3) ;
  endfor

  p /= 4*pi ;
  
