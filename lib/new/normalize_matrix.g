Read("normalize_column.g");

normalize_matrix := function (mat)
   local i, j, B, c, s, si, Bi, Be, R, U, NZC; 

   R := HomalgRing(mat);

   U := HomalgInitialIdentityMatrix(NrRows(mat), R);

   if NrRows(mat) = 0 then 
      return HomalgInitialIdentityMatrix(NrRows(mat), R);
   fi;

   NZC := NonZeroColumns(mat);
   i := 1;

   while not NZC = [] do
     
     Bi := normalize_column(CertainColumns(mat, [NZC[1]]));

     Be := DiagMat([HomalgIdentityMatrix(i - 1, R), Bi]);

     mat := Bi * mat;
     U := Be * U;
     
     mat := CertainRows(mat,[2..NrRows(mat)]);
   
     NZC := NonZeroColumns(mat);
     
     i := i + 1;
   od;

   return U;
end;
