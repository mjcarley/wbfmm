grad = 0 ;
sfile = "monopoles" ;

[xs,q,n,f] = readpts([sfile ".dat"]) ;

nq = size(q, 2) ;

xf = readpts("field.dat") ;

dat = load([sfile "-laplace-fmm.dat"]) ;
xf = dat(:,1:3) ;
pf = dat(:,4:end) ;

dat = load("monopole-laplace-avx.dat") ;
xf = dat(:,1:3) ;
pa = dat(:,4:end) ;

dat = load([sfile "-laplace-direct.dat"]) ;
xf = dat(:,1:3) ;
pd = dat(:,4:end) ;
