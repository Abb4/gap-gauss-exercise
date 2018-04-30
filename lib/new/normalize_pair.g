LoadPackage("RingsForHomalg");

ZZ := HomalgRingOfIntegers();
QQ := HomalgFieldOfRationals();

Qx := QQ * "x";

AssignGeneratorVariables(Qx);

mulmat := function(a, R) 
  return HomalgMatrix([a^-1, 0, 0, 1], 2, 2, R);
end;

addmat := function(a, R)
  return HomalgMatrix([1, 0, -a, 1], 2, 2, R);
end;

vermat := function( m, i, j, R)
  local M;
  M := HomalgInitialIdentityMatrix( m, R);

  SetMatElm( M, i, i, "0");
  SetMatElm( M, i, j, "1");
  SetMatElm( M, j, j, "0");
  SetMatElm( M, j, i, "1");
  
  return M;
end;

#--------------------------------------------------------------------

normalize_pair_Field := function(a, b, R)
  local U, P; 

  U := HomalgIdentityMatrix(2, R);

  if IsZero(a) and IsZero(b) then
     return U;
  elif IsZero(b) then
     return mulmat(a, R);
  elif IsZero(a) then
     return mulmat(b, R) * P;
  fi;

  return addmat(b, R) * mulmat(a, R);
end;

#----------------------------------------------------------------------------

normalize_pair_Ring := function(a, b, R)
  local q, r, P, U, res;
  
  a := One(R)*a;
  b := One(R)*b;
  
  U := HomalgIdentityMatrix(2, R);

  while not IsZero(b) do
    q := EuclideanQuotient(R, a, b);
    
    U := HomalgMatrix( [ 0, 1, 1, -q ], 2, 2, R) * U;

    r := a - (q*b);
    a := b;
    b := r;
  od;

  return U;
end;

#----------------------------------------------------------------------------

normalize_pair := function(a, b, R)
  if IsFieldForHomalg(R) then
    return normalize_pair_Field(a, b, R);
  else
    return normalize_pair_Ring(a, b, R);
  fi;
end;

#----------------------------------------------------------------------------

normalize_pair_inflated := function(a, b, i, j, n, R)
  local U,u;

  U := HomalgInitialIdentityMatrix( n, R );

  u := normalize_pair( a, b, R );

  SetMatElm( U, i, i, MatElm( u, 1, 1 ) );
  SetMatElm( U, i, j, MatElm( u, 1, 2 ) );
  SetMatElm( U, j, i, MatElm( u, 2, 1 ) );
  SetMatElm( U, j, j, MatElm( u, 2, 2 ) );

  return U;
end;

#----------------------------------------------------------------------------Testing

test_pair := function(a,b,R) 
  return normalize_pair(a,b,R)*HomalgMatrix([a,b],2,1,R);
end;

