LoadPackage( "RingsForHomalg" );

ZZ := HomalgRingOfIntegers( );
QQ := HomalgFieldOfRationals( );

mulmat := function(a, k) 
  return HomalgMatrix([a^-1, 0, 0, 1], 2, 2, k);
end;

addmat := function(a, k)
  return HomalgMatrix([1, 0, -a, 1], 2, 2, k);
end;

addmatn := function(a, i, j, dim, K)
  local M;
  
  M := HomalgInitialIdentityMatrix( dim, K);
  
  SetMatElm(M, j, i, a);
  
  return M;
end;

vermat := function( dim, i, j, K)
  local M;
  M := HomalgInitialIdentityMatrix( dim, K);

  SetMatElm( M, i, i, 0);
  SetMatElm( M, i, j, 1);
  SetMatElm( M, j, j, 0);
  SetMatElm( M, j, i, 1);
  
  return M;
end;

#--------------------------------------------------------------------

normalize_pair_Q := function(mat)
  
  local a, b, U, P, R;
  a := MatElm(mat, 1, 1);
  b := MatElm(mat, 2, 1);

  U := HomalgIdentityMatrix(2, QQ);
  P := HomalgMatrix( [ 0, 1, 1, 0 ], 2, 2, QQ);

  R := HomalgRing( mat );

  if a = 0 and b = 0 then
     return U;
  elif b = 0 then
     return mulmat( a, R);
  elif a = 0 then
     return mulmat( b, R) * P;
  fi;

  return addmat(b, R) * mulmat(a, R);
end;

#-----------------------------------------------------------------------

normalize_pair_Z := function(m)
  
  local a, b, q, r, P, U, res;
  a := MatElm(m, 1, 1);
  b := MatElm(m, 2, 1);

  U := HomalgIdentityMatrix(2, ZZ);
  P := HomalgMatrix( [ 0, 1, 1, 0 ], 2, 2, ZZ);

  if a = 0 and b = 0 then
     return U;
  elif b = 0 then
     return U;
  elif a = 0 then
     return P;
  fi;

  q := QuoInt(a,b);
  r := a - (b*q);

  res := HomalgIdentityMatrix(2, ZZ);

  if r = 0 then
    return HomalgMatrix( [ 0, 1, 1, -q ], 2, 2, ZZ) * res;
  fi;

  while not r = 0 do
    q := QuoInt(a,b);
    r := a - (b*q);
    a := b;
    b := r; 

    res := HomalgMatrix( [ 0, 1, 1, -q ], 2, 2, ZZ) * res;
  od;

  return res;
end;

#----------------------------------------------------------------------------

normalize_pair := function(m)
  local R;
  R := HomalgRing( m );

  if HasIsIntegersForHomalg(R) and IsIntegersForHomalg(R) then
    return normalize_pair_Z(m);
  else
    return normalize_pair_Q(m);
  fi;
end;

#----------------------------------------------------------------------------

HomalgRandomColumn := function(r, K)
  return HomalgMatrix( RandomMat(r,1),r,1,K);
end;

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------

normalize_column_Q := function(A)

  local U, i, j, mat, mati, swapped, swappedJ, K;
  
  swapped := false;
  swappedJ := 1;

  K := HomalgRing( A );

  U := HomalgInitialIdentityMatrix(NrRows(A), K);
  
  if NrRows(A) = 0 then
     return HomalgInitialIdentityMatrix(NrRows(A), K);
  elif NrRows(A) = 1 and MatElm(A, 1, 1) = 0 then
     return U;
  elif NrRows(A) = 1 then
     return HomalgMatrix([MatElm(A, 1, 1)^-1], 1, 1, K);
  fi;

  if MatElm(A, 1, 1) = 0 then
     for j in [2..NrRows(A)] do
        if not MatElm(A, j, 1) = 0 then
          A := vermat(NrRows(A), 1, j, K) * A;
          swappedJ := j;
          swapped := true;
          break;
        fi;
     od;
     if swapped = false then return U;
     fi;
  fi;

  for i in [1..NrRows(A)-1] do
    mat := HomalgMatrix( [ MatElm(A, 1, 1), MatElm(A, i + 1, 1) ], 2, 1, K);

    mati := normalize_pair(mat);

    SetMatElm(U, 1, 1, (MatElm(mat, 1, 1)^-1));

    SetMatElm(U, i+1, 1, MatElm(mati, 2, 1));
  od;  
  
  if swapped = true then return U * vermat(NrRows(A), 1, swappedJ, K);
  else return U;
  fi;
end;


#------------------------------------------------------------------------------

normalize_column_Z := function (A)
  
  local a, i, j, k, r, mat, mati, mata, matai, U, K, H;

  K := ZZ;

  U := HomalgInitialIdentityMatrix(NrRows(A), K);
  
  if NrRows(A) = 0 then return U;
  fi;

  H := HomalgIdentityMatrix(NrRows(A), K);

  for i in [1..NrRows(A)] do
    for j in [i..NrRows(A)] do
      if not MatElm(A, j, 1) = 0 then
        H := vermat(NrRows(A), i, j, K) * H;
        A := vermat(NrRows(A), i, j, K) * A;
        break;
      fi;
    od;
  od;

  r := MatElm(A, 1, 1);

  for i in [1..NrRows(A)-1] do

    mat := HomalgMatrix( [ MatElm(A, i, 1), MatElm(A, i+1, 1) ], 2, 1, K);
    mati := normalize_pair_Z(mat);

    mata := HomalgMatrix( [ r, MatElm(A, i+1, 1) ], 2, 1, K);

    matai := normalize_pair_Z(mata);
    
    if not r = 1 then
      r := MatElm(mati * mat, 1, 1);

      a := MatElm(matai, 1, 1);

      SetMatElm(U, 1, i+1, MatElm(matai, 1, 2));

      if i = 1 then 
        SetMatElm(U, 1, 1, MatElm(matai, 1, 1));
      fi;

      if i > 1 then
        for j in [1..i] do
          SetMatElm(U, 1, j, MatElm(U, 1, j) * a);
        od;
      fi;
    fi;
    
    SetMatElm(U, i+1, i, MatElm(mati, 2, 1));
    SetMatElm(U, i+1, i+1, MatElm(mati, 2, 2));
  od;

  return U * H;
end;

#------------------------------------------------------------------------------


normalize_column := function (A)
  local R;
  R := HomalgRing( A );

  if HasIsIntegersForHomalg(R) and IsIntegersForHomalg(R) then
    return normalize_column_Z(A);
  else
    return normalize_column_Q(A);
  fi;
end;

#------------------------------------------------------------------------------

normalize_matrix := function (A)
   local i, j, B, c, s, si, Bi, Be, K, U, NZC; 

   K := HomalgRing( A );

   U := HomalgInitialIdentityMatrix(NrRows(A), K);

   if NrRows(A) = 0 then 
      return HomalgInitialIdentityMatrix(NrRows(A), K);
   fi;

   NZC := NonZeroColumns(A);
   i := 1;

   while not NZC = [] do
     
     Bi := normalize_column(CertainColumns(A, [NZC[1]]));
     #Display("Bi");
     #Display(Bi);
     #Be := matrix_expand(Bi, NrRows(A), i, i);
     Be := DiagMat([HomalgIdentityMatrix(i - 1, K), Bi]);

     #Display("Be");
     #Display(Be);

     A := Bi * A;
     U := Be * U;
     
     A := CertainRows(A,[2..NrRows(A)]);
   
     NZC := NonZeroColumns(A);
     
     i := i + 1;
   od;

   return U;
end;

#--------------------------------------------------------------------------------------------

sindex := function(M)
   local res, i,j;
   res := [];
   i := 1;

   for j in [1..NrColumns(M)] do
      if i < NrRows(M) + 1 then
        if not MatElm(M, i, j) = 0 then
          Add(res, j);
          i := i+1;
        fi;
      fi;
   od;

   return res;
end;

#--------------------------------------------------------------------------------------------

strictly_normalize_matrix := function (A)
   local i, j, k, si, U, u, a, b, R;

   R := HomalgRing( A );
   
   U := normalize_matrix(A);

   A := U * A;

   si := sindex(A);
   
   i := 1;
   
   for j in si do
     a := MatElm(A,i,j);
     for k in [1..i-1] do
       b := MatElm(A,k,j);
       
       if HasIsIntegersForHomalg(R) and IsIntegersForHomalg(R) then
         u := addmatn(-QuoInt(b,a), i, k, NrRows(A), R);
       else
         u := addmatn(-b, i, k, NrRows(A), R);
       fi;

       U := u * U;
       A := u * A;
     od;
     i := i + 1;
   od;

   return U;
end;

#--------------------------------------------------------------------------------------------

decide_zero_rows := function(B,A)
  local r, rs, c, C, K, E, L, R, X;

  K := HomalgRing( A );

  r  := NrRows(A);
  rs := NrRows(B);
  c  := NrColumns(A);

  Display(["r: ",r,"rs: ", rs, "c: ", c]);

  E := HomalgIdentityMatrix(r, K);

  
  C := UnionOfRows(B,A);
  
  L := UnionOfRows(HomalgIdentityMatrix(rs, K), HomalgZeroMatrix(r,rs,K));
  R := UnionOfRows(HomalgZeroMatrix(rs,r,K),HomalgIdentityMatrix(r, K));
  
  C := UnionOfColumns(L,C);
  C := UnionOfColumns(C,R);
  
  Display("C:");
  Display(C);

  X := strictly_normalize_matrix(C);
  
  Display("Xnorm:");
  Display(X);  

  X := CertainColumns(X, [(NrColumns(X) - r + 1)..NrColumns(X)]);
  
  Display("Xno columns:");
  Display(X); 

  X := CertainRows(X, [1..rs]);
  
  X := -X;  

  Display("Xno rows:");
  Display(X); 
  
  if X*A = B then
     return true;
  else
     return false;
  fi;
end;

#--------------------------------------------------------------------------------------------

decide_zero_rows_effectively := function(B,A)
  local r, rs, c, C, K, E, L, R, X, BS;

  K := HomalgRing( A );

  r  := NrRows(A);
  rs := NrRows(B);
  c  := NrColumns(A);

  Display(["r: ",r,"rs: ", rs, "c: ", c]);

  E := HomalgIdentityMatrix(r, K);

  
  C := UnionOfRows(B,A);
  
  L := UnionOfRows(HomalgIdentityMatrix(rs, K), HomalgZeroMatrix(r,rs,K));
  R := UnionOfRows(HomalgZeroMatrix(rs,r,K),HomalgIdentityMatrix(r, K));
  
  C := UnionOfColumns(L,C);
  C := UnionOfColumns(C,R);
  
  Display("C:");
  Display(C);

  X := strictly_normalize_matrix(C);
  
  Display("Xnorm:");
  Display(X);  

  X := CertainColumns(X, [(NrColumns(X) - r + 1)..NrColumns(X)]);
  
  Display("Xno columns:");
  Display(X); 

  X := CertainRows(X, [1..rs]);
  
  X := -X;  

  Display("Xno rows:");
  Display(X); 
  
  BS := B - X*A;

  Display("BS");
  Display(BS);
  
  if X*A = B then
     return BS;
  else
     return fail;
  fi;
end;

#--------------------------------------------------------------------------------------------

syzygies_of_rows := function(A)
  local r, rs, c, C, K, L, R, X, B, BS;

  K := HomalgRing( A );

  r  := NrRows(A);
  c  := NrColumns(A);

  B := HomalgZeroMatrix(1, c, K);

  rs := NrRows(B);

  Display(["r: ",r ,"c: ", c]);
  
  C := UnionOfRows(B,A);
  
  L := UnionOfRows(HomalgIdentityMatrix(rs, K), HomalgZeroMatrix(r,rs,K));
  R := UnionOfRows(HomalgZeroMatrix(rs,r,K),HomalgIdentityMatrix(r, K));
  
  C := UnionOfColumns(L,C);
  C := UnionOfColumns(C,R);
  
  Display("C:");
  Display(C);

  X := strictly_normalize_matrix(C);
  
  Display("Xnorm:");
  Display(X);  

  X := CertainColumns(X, [(NrColumns(X) - r + 1)..NrColumns(X)]);
  
  Display("Xno columns:");
  Display(X); 

  X := CertainRows(X, [1..rs]);
  
  X := -X;  

  Display("Xno rows:");
  Display(X); 
  
  BS := B - X*A;

  Display("BS");
  Display(BS);
  
  if X*A = B then
     return BS;
  else
     return fail;
  fi;
end;




#--------------------------------------------------------------------------------------------
HomalgRandomMatrix := function(r,c,R)
  return HomalgMatrix(RandomMat(r,c),r,c,R);
end;

#r := 1;; c:= 4;; R := QQ;; mat := HomalgRandomMatrix(r,c,R);; Display(mat); U := normalize_matrix(mat);; Display(U); Display(U*mat);



