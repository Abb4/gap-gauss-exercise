LoadPackage( "RingsForHomalg" );;

ZZ := HomalgRingOfIntegers( );
QQ := HomalgFieldOfRationals( );

mulmat := function(a, k) 
  return HomalgMatrix([a^-1, 0, 0, 1], 2, 2, k);
end;

addmat := function(a, k)
  return HomalgMatrix([1, 0, -a, 1], 2, 2, k);
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
strictly_normalize_matrix := function (A)
  local NZC, nr, nc, SF, SSF;

  nr := NrRows(A);
  nc := NrColumns(A);

  SF := normalize_matrix(A);

  if nr = 1 or nc = 1 then
    return SF;
  fi;

  A := SF * A;

  if nr < nc then
    NZC := NonZeroColumns(A);

    A := CertainColumns(A, NZC{[1..Minimum(Length(NZC),nr)]});
  fi;
  
  A := Involution(A);
  SSF := normalize_matrix(A);

  if NrRows(SF) > NrRows(SSF) then
    SSF := DiagMat([SSF,HomalgIdentityMatrix(NrRows(SF) - NrRows(SSF), HomalgRing(A))]);
  fi;

  return Involution(SSF) * SF;
end;
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
  local Bi, PartialMatI, U, NZC, R; 

  R := HomalgRing(A);
  NZC := NonZeroColumns(A);
  
  if NZC = [] then
    return HomalgIdentityMatrix(NrRows(A), R);
  else
    Bi := normalize_column(CertainColumns(A, [NZC[1]]));
    
    # new matrix
    A := Bi * A;
    # cut away top row
    A := CertainRows(A,[2..NrRows(A)]);
    
    # normalize patrial matrix
    PartialMatI := normalize_matrix(CertainColumns(A,NZC));
    # expand to original size
    PartialMatI := DiagMat([HomalgIdentityMatrix(1, R), PartialMatI]);

    return PartialMatI * Bi;
  fi;
end;

#--------------------------------------------------------------------------------------------

strictly_normalize_matrix_qq := function (A)
  local NZC, nr, nc, SF, SSF;

  nr := NrRows(A);
  nc := NrColumns(A);

  SF := normalize_matrix(A);

  if nr = 1 or nc = 1 then
    return SF;
  fi;

  A := SF * A;

  if nr < nc then
    NZC := NonZeroColumns(A);

    A := CertainColumns(A, NZC{[1..Minimum(Length(NZC),nr)]});
  fi;
  
  A := Involution(A);
  SSF := normalize_matrix(A);

  if NrRows(SF) > NrRows(SSF) then
    SSF := DiagMat([SSF,HomalgIdentityMatrix(NrRows(SF) - NrRows(SSF), HomalgRing(A))]);
  fi;

  return Involution(SSF) * SF;
end;

#--------------------------------------------------------------------------------------------

strictly_normalize_matrix_zz := function (A)
  local NZC, nr, nc, SF, SSF;

  nr := NrRows(A);
  nc := NrColumns(A);

  SF := normalize_matrix(A);

  A := SF * A;

  NZC := NonZeroColumns(A);

  return fail;
end;




HomalgRandomMatrix := function(r,c,R)
  return HomalgMatrix(RandomMat(r,c),r,c,R);
end;

#r := 1;; c:= 4;; R := QQ;; mat := HomalgRandomMatrix(r,c,R);; Display(mat); U := normalize_matrix(mat);; Display(U); Display(U*mat);
