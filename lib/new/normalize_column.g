Read("normalize_pair.g");

normalize_column := function(col)
  local a,b,j,n,nn,u,A,U,R;

  R := HomalgRing(col);
  
  n := NrRows(col);  

  U := HomalgIdentityMatrix(n, R);
  
  nn := NonZeroRows(col);

  for j in nn do 
    
    a := MatElm(col,1,1);
    b := MatElm(col,j,1);

    u := normalize_pair_inflated(n,a,b,1,j,R);

    A := u * A;
    U := u * U;
    
  od;
  

end;
