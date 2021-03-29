##[xs, qs, ns, fs] = readpts("monopoles.dat") ;
[xs, qs, ns, fs] = readpts("dipoles.dat") ;
xf = readpts("field.dat") ;

p = calcfield(xs, ns, qs, fs, xf) ;

pc = load("test.dat") ;
pc = pc(:,end) ;

