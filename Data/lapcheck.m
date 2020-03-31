grad = 0 ;

[xs,q] = readpts("monopoles.dat") ;

xf = readpts("field.dat") ;

dat = load("monopole-laplace-direct.dat") ;
xf = dat(:,1:3) ;
pd = dat(:,4:end) ;

dat = load("monopole-laplace-fmm.dat") ;
xf = dat(:,1:3) ;
pf = dat(:,4:end) ;

pc = 0*pd ;
if ( ~grad ) 
  for i=1:length(xf)
    r = [xf(i,1) - xs(:,1) xf(i,2) - xs(:,2) xf(i,3) - xs(:,3)] ;
    R = sqrt(sum(r.^2,2)) ;
    for ii=1:size(pd,2)
      pc(i,ii) = sum(q(:,ii)./R)/4/pi ;
    endfor
  endfor
else
  for i=1:length(xf)
    r = [xf(i,1) - xs(:,1) xf(i,2) - xs(:,2) xf(i,3) - xs(:,3)] ;
    R = sqrt(sum(r.^2,2)) ;
    kR = k*R ;
    nR = r./[kR kR kR] ;
    E = [exp(j*kR) exp(j*kR) exp(j*kR)]./R.^2.*(j*kR-1).*nR/4/pi/j ;
    tt = [] ;
    for ii=1:nq
      tt = [tt sum([q(:,ii) q(:,ii) q(:,ii)].*E)] ;
    endfor
    pc(i,:) = tt ;
  endfor
endif
