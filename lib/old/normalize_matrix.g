LoadPackage( "RingsForHomalg" );

ZZ := HomalgRingOfIntegers( );
QQ := HomalgFieldOfRationals( );

mulmat := function(a, k) 
  return HomalgMatrix([a^-1, 0, 0, 1], 2, 2, k);
end;

addmat := function(a, k)
  return HomalgMatrix([1, 0, -a, 1], 2, 2, k);
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

vermat := function( m, i, j, K)
  local M;
  M := HomalgInitialIdentityMatrix( m, K);

  SetMatElm( M, i, i, "0");
  SetMatElm( M, i, j, "1");
  SetMatElm( M, j, j, "0");
  SetMatElm( M, j, i, "1");
  
  return M;
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

#r := 1;; R := QQ;; col:= HomalgRandomColumn( r, R);; Display(col); U := normalize_column_Z(col);; Display(U); Display(U * col);


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

matrix_expand := function(A, n, i, j)
# gegeben eine Matrix A
# setze A in neue Einheitsmatrix von Groesse n mit rein, wobei I_ij = a_11

  local k, l, U, K;
  
  K := HomalgRing( A );
  # schreibschutz von A umgehen
  U := HomalgInitialIdentityMatrix(n, K);

  for k in [1..NrRows(A)] do
    for l in [1..NrColumns(A)] do
      SetMatElm(U, i + k - 1, j + l - 1, MatElm(A, k, l));
    od;
  od;

  return U;
end;

#A  := HomalgMatrix( [ 1, 2, 3, 4, 5, 6, 7, 8, 9], 3, 3, QQ);

#Display(matrix_expand(A,4,2,2));

#------------------------------------------------------------------------------

normalize_matrix2 := function (A)
   local i, j, b, c, s, si, K, U; 

   K := HomalgRing( A );

   U := HomalgInitialIdentityMatrix(NrRows(A), K);

   if NrRows(A) = 0 then 
      return U;
   fi;
   
   # todo skip some rows columns here
   i := 1;
   j := 1;

   s := CertainColumns(A, [i]);
   si := normalize_column(s);
   
   U := U * si;

   b := CertainRows(CertainColumns(A, [(j + 1)..NrRows(A)]), [(i + 1)..NrRows(A)]);

   #return matrix_union(U, normalize_matrix(b), i + 1, j + 1);

   return U;

   #for i in [1..NrRows(A)] do
   #  b := CertainRows(CertainColumns(A, [i..NrRows(A)]), [i..NrRows(A)]);

   #  Display(b);
   #  s := CertainColumns(A, [i]);

   #  si := normalize_column(s);

     #A := si * A;
   #od;

   #Display(si);

   #return A;

   #for i in [1..NrRows(A)] do

    # B := CertainRows(CertainColumns(A, [i..NrRows(A)]), [i..NrRows(A)]);
     
    # Bi := normalize_column(CertainColumns(B, [1]));
     
    # Be := matrix_expand(Bi, NrRows(A), i, i);
     
    # A := Be * A;
    # U := Be * U;
   #od;

   return U;

end;



normalize_matrix := function (A)
   local i, Bi, Be, U, NZC, R; 

   R := HomalgRing( A );

   U := HomalgInitialIdentityMatrix(NrRows(A), R);

   if NrRows(A) = 0 or NrColumns(A) = 0 then 
    return HomalgIdentityMatrix(NrRows(A), HomalgRing(A));
  fi;

   NZC := NonZeroColumns(A);
   i := 1;

   while not NZC = [] do
     
     Bi := normalize_column(CertainColumns(A, [NZC[1]]));
     #Display("Bi");
     #Display(Bi);
     #Be := matrix_expand(Bi, NrRows(A), i, i);
     Be := DiagMat([HomalgIdentityMatrix(i - 1, R), Bi]);

     #Display("Be");
     #Display(Be);

     A := Bi * A;
     U := Be * U;
     
     A := CertainRows(A,[2..NrRows(A)]);
   
     NZC := NonZeroColumns(A);
     
     i := i + 1;
   od;
   
   MakeImmutable(U);
   return U;
end;

# TODO: Rem this
strictly_normalize_matrix := function (A)
  local i, j, B, c, s, si, Bi, Be, NZC, nr, nc, SF, SSF, Aorig, R; 

  Aorig := A;

  #if NrRows(A) = 0 or NrColumns(A) = 0 then 
  #  return HomalgIdentityMatrix(NrRows(A), HomalgRing(A));
  #fi;

  #Display("A");
  #Display(A);

  R := HomalgRing(A);

  nr := NrRows(A);
  nc := NrColumns(A);

  SF := normalize_matrix(A);

  if nr = 1 or nc = 1 then
    return SF;
  fi;
  
  #Display("SF");
  #Display(SF);

  A := SF * A;

  #Display("SF * A");
  #Display(A);

  if nr < nc then
    NZC := NonZeroColumns(A);
    #Display("NZC");
    #Display(NZC);
    #Display("A NZC Squared");
    #Display(CertainColumns(A, NZC{[1..Minimum(Length(NZC),nr)]}));
    A := CertainColumns(A, NZC{[1..Minimum(Length(NZC),nr)]});
    #Display("Areduced");
    #Display(A);
  fi;
  
  A := Involution(A);
  SSF := normalize_matrix(A);
  #Display("Ainvoluted");
  #Display(A);
  #Display("SSF(Ainvoluted)");
  #Display(SSF);

  #Display("SSF(Ainvolued) * A(involuted)");
  #Display(SSF * A);

  #Display("SSF(Ainvolued)Involuted 2 ");
  #Display(SSF);

  #Display("SSF(Ainvolued)Involuted * SF");

  if NrRows(SF) > NrRows(SSF) then
    SSF := DiagMat([SSF,HomalgIdentityMatrix(NrRows(SF) - NrRows(SSF), R)]);
  fi;

  SSF := Involution(SSF) * SF;
  #Display(SSF);

  #Display("Aorig");
  #Display(Aorig);

  #Display("SF * SSF(Ainvolued)Involuted * Aorig");
  #Display(SSF * Aorig);

  return SSF;
end;


HomalgRandomMatrix := function(r,c,R)
  return HomalgMatrix(RandomMat(r,c),r,c,R);
end;

#r := 1;; c:= 4;; R := QQ;; mat := HomalgRandomMatrix(r,c,R);; Display(mat); U := normalize_matrix(mat);; Display(U); Display(U*mat);