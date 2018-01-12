gap> START_TEST( "tests for gauss matrix normalisation (row echelon form)" );

#setup
gap> ZZ := HomalgRingOfIntegers( );
Z
gap> QQ := HomalgFieldOfRationals( );
Q
gap> A1 := HomalgMatrix( [-1], 1, 1, QQ);
<A 1 x 1 matrix over an internal ring>
gap> A2 := HomalgMatrix( [-3, -6, -6, -13], 2, 2, QQ);
<A 2 x 2 matrix over an internal ring>
gap> A3 := HomalgMatrix( [ [1, 2 , 3] , [4 , 5 , 6], [7 , 8 , 8 ]], 3, 3, QQ);
<A 3 x 3 matrix over an internal ring>

# actual testing
# Q 1x1 matrix
gap> normalize_matrix(A1);
<An unevaluated 1 x 1 matrix over an internal ring>
gap> Display(last*A1);
[ [  1 ] ]

# Q 2x2 matrix
gap> normalize_matrix(A2);
<An unevaluated 2 x 2 matrix over an internal ring>
gap> Display(last*A2);
[ [  1,  2 ],
  [  0,  1 ] ]

# Q 3x3 matrix
gap> normalize_matrix(A3);
<An unevaluated 3 x 3 matrix over an internal ring>
gap> Display(last*A3);
[ [  1,  2,  3 ],
  [  0,  1,  2 ],
  [  0,  0,  1 ] ]

# Q 2x2 with swap
gap> A2 := HomalgMatrix( [0, 2, 2, 3], 2, 2, QQ);
<A 2 x 2 matrix over an internal ring>
gap> normalize_matrix(A2);
<An unevaluated 2 x 2 matrix over an internal ring>
gap> Display(last*A2);
[ [    1,  3/2 ],
  [    0,    1 ] ]

#
gap> STOP_TEST( "strictly_normalize_marix.tst", 10000 );
