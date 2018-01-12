Display("Reading library..");
if not IsExistingFile( "lib/strictly_normalize_matrix.g" ) then
  Display(LastSystemError().message);
else
  Display("Attempting to read strictly_normalize_matrix.g");
  Read("lib/strictly_normalize_matrix.g");
fi;

Display("Running tests..");
TestDirectory("tst");


