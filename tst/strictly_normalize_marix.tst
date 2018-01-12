gap> START_TEST( "tests for gauss matrix normalisation (reduced row echelon form)" );

#setup
gap> ZZ := HomalgRingOfIntegers( );
Z
gap> QQ := HomalgFieldOfRationals( );
Q

# actual testing
# R = Q (Field)
# Q 1x1 matrix
gap> A := HomalgMatrix([-1], 1, 1, QQ);
<A 1 x 1 matrix over an internal ring>
gap> strictly_normalize_matrix(A);
<An unevaluated 1 x 1 matrix over an internal ring>
gap> Display(last*A);
[ [  1 ] ]

# Q 1x3 matrix
gap> A := HomalgMatrix([-1, 0, 5], 1, 3, QQ);
<A 1 x 3 matrix over an internal ring>
gap> strictly_normalize_matrix(A);
<An unevaluated 3 x 3 matrix over an internal ring>
gap> Display(last*A);
[ [  1 ] ]

# Q 3x1 matrix
gap> A := HomalgMatrix([2, 3, 5], 3, 1, QQ);
<A 3 x 1 matrix over an internal ring>
gap> normalize_matrix(A);
<An unevaluated 3 x 3 matrix over an internal ring>
gap> Display(last*A);
[ [  1 ],
  [  0 ],
  [  0 ] ]

# Q 2x2 matrix
gap> A := HomalgMatrix([[-3, -6], [-6, -13]], 2, 2, QQ);
<A 2 x 2 matrix over an internal ring>
gap> strictly_normalize_matrix(A);
<An unevaluated 2 x 2 matrix over an internal ring>
gap> Display(last*A);
[ [  1,  0 ],
  [  0,  1 ] ]

# Q 3x3 matrix
gap> A := HomalgMatrix([[1, 2 , 3], [4 , 5 , 6], [7 , 8 , 8 ]], 3, 3, QQ);
<A 3 x 3 matrix over an internal ring>
gap> strictly_normalize_matrix(A);
<An unevaluated 3 x 3 matrix over an internal ring>
gap> Display(last*A);
[ [  1,  0,  0 ],
  [  0,  1,  0 ],
  [  0,  0,  1 ] ]

# Q 0x0 matrix
gap> A := HomalgMatrix([], 0, 0, QQ);
<An unevaluated 0 x 0 zero matrix over an internal ring>
gap> strictly_normalize_matrix(A);
<An unevaluated 0 x 0 zero matrix over an internal ring>
gap> Display(last*A);
(an empty 0 x 0 matrix)

# Q 1x0 matrix
gap> A := HomalgMatrix([], 1, 0, QQ);
<An unevaluated 1 x 0 zero matrix over an internal ring>
gap> strictly_normalize_matrix(A);
<An unevaluated 1 x 1 identity matrix over an internal ring>
gap> Display(last*A);
(an empty 1 x 0 matrix)

# Q 0x1 matrix
gap> A := HomalgMatrix([], 0, 1, QQ);
<An unevaluated 0 x 1 zero matrix over an internal ring>
gap> strictly_normalize_matrix(A);
<An unevaluated 0 x 0 zero matrix over an internal ring>
gap> Display(last*A);
(an empty 0 x 1 matrix)

#
gap> STOP_TEST( "strictly_normalize_marix.tst", 10000 );
