N = 5 ;

id = [] ;
i = 0 ;
#for n=0:N
#  for m=-n:n
#    id = [id; n m i n*(n+1)+m] ;
#    i ++ ;
#  endfor
#endfor

for n=0:N
  for m=0:n
    id = [id; n m i n*(n+1)/2+m] ;
    i ++ ;
  endfor
endfor
