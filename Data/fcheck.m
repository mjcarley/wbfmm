k = 1.5 ;

[xs,qq] = readpts("monopoles.dat") ;
q = [] ; nq = size(qq,2) ;
for i=1:2:nq
  q = [q qq(:,i:i+1)*[1; j]] ;
endfor
xf = readpts("field.dat") ;

dat = load("monopole-direct.dat") ;
xf = dat(:,1:3) ;
pd = [] ;
for i=4:2:size(dat,2)
  pd = [pd dat(:,i:i+1)*[1; j]] ;
endfor

dat = load("monopole-fmm.dat") ;
xf = dat(:,1:3) ;
pf = [] ;
for i=4:2:size(dat,2)
  pf = [pf dat(:,i:i+1)*[1; j]] ;
endfor

pc = 0*pd ;
for i=1:length(xf)
  r = [xf(i,1) - xs(:,1) xf(i,2) - xs(:,2) xf(i,3) - xs(:,3)] ;
  kR = k*sqrt(sum(r.^2,2)) ;
  for ii=1:size(pd,2)
    pc(i,ii) = sum(q(:,ii).*exp(j*kR)./kR)/4/pi/j ;
  endfor
endfor
