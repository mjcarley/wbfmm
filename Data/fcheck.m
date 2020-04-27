k = 1.5 ; grad = 0 ;
sfile = "mixed" ;

[xs,qq,n,ff] = readpts([sfile ".dat"]) ;
q = [] ; nq = size(qq,2) ;
for i=1:2:nq
  q = [q qq(:,i:i+1)*[1; j]] ;
endfor
if ( ff ~= 0 )
  f = [] ;
  for i=1:2:nq
    f = [f ff(:,i:i+1)*[1; j]] ;
  endfor
else
  n = zeros(size(xs)) ;
  f = zeros(size(q)) ;
endif

nq /= 2 ;

xf = readpts("field.dat") ;

dat = load([sfile "-direct.dat"]) ;
xf = dat(:,1:3) ;
pd = [] ;
for i=4:2:size(dat,2)
  pd = [pd dat(:,i:i+1)*[1; j]] ;
endfor

if 0
dat = load("monopole-fmm.dat") ;
xf = dat(:,1:3) ;
pf = [] ;
for i=4:2:size(dat,2)
  pf = [pf dat(:,i:i+1)*[1; j]] ;
endfor
endif

dat = load([sfile "-avx.dat"]) ;
xf = dat(:,1:3) ;
pa = [] ;
for i=4:2:size(dat,2)
  pa = [pa dat(:,i:i+1)*[1; j]] ;
endfor

if 0
pc = 0*pd ;
if ( ~grad ) 
  for i=1:length(xf)
    r = [xf(i,1) - xs(:,1) xf(i,2) - xs(:,2) xf(i,3) - xs(:,3)] ;
    R = sqrt(sum(r.^2,2)) ;
    nR = -sum(r.*n,2)./R ;
    kR = k*R ;
    for ii=1:size(pd,2)
      pc(i,ii)  = sum(q(:,ii).*exp(j*kR)./R)/4/pi/j/k ;
      pc(i,ii) -= sum(f(:,ii).*exp(j*kR)./R.^2.*(j*k*R-1).*nR)/4/pi/j/k ;
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
endif
