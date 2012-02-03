#!/usr/bin/perl

# run code in directorites TEST1, TEST2, ... TESTN

# enter full path for code directory
my $CODE_DIR = '/scratch2/maerzkek/MCCCS-2010-SERIAL';
# enter number of TEST directories
my $NUM_TESTS = 9;

# change this if want to run script from other directory
chomp(my $WORK_DIR = `pwd`);
my @TEST_DIR;
my @EX;
my $diff_file = DIFF;
my @diff_lines;

for ($i=1;$i<=$NUM_TESTS;$i++){
    $TEST_DIR[$i] = $WORK_DIR ."/TEST" . "$i";
#    print "$TEST_DIR[$i] \n";
}

# assign proper executable to each TEST_DIR

# TEST1 = NpT
$EX[1] = $CODE_DIR . "/NpTewald";
# TEST2 = NVT field
$EX[2] = $CODE_DIR . "/NVTfield";
# TEST3 = NpT tab
$EX[3] = $CODE_DIR . "/NpTtab";
# TEST4 = NVT-Gibbs ewald (1,2 butanediol)
$EX[4] = $CODE_DIR . "/NVT-Gibbs";
# TEST5 = NpT garofalini
$EX[5] = $CODE_DIR . "/NpTgaro";
# TEST6 = NVT-Gibbs ewald (TIP4P)
$EX[6] = $CODE_DIR . "/NVT-Gibbs";
# TEST7 = RPLC (RPLC)
$EX[7] = $CODE_DIR . "/RPLC";
# TEST8 = NVT-Gibbs no ewald (n-octane)
$EX[8] = $CODE_DIR . "/NVT-Gibbs-noq";
# TEST9 = NpT solid TATB
$EX[9] = $CODE_DIR . "/NpTsolid";


# remove existing DIFF file, if it exists
unlink "$diff_file";

# open DIFF with append access (creates file, then appends)
open DIFF, ">>", "$diff_file";

# loop over test directories and run executable
# then diff previous results
for ($i=1; $i<=$NUM_TESTS; $i++){
    chdir $TEST_DIR[$i];
    print "now running TEST $i\n";
    system "$EX[$i]";
#   move DIFF to current TEST_DIR
    rename "$WORK_DIR"."/$diff_file", "$TEST_DIR[$i]"."/$diff_file";
    print DIFF "TEST $i \n\n";
    @diff_lines = `diff run1a.dat prev.out`;
    print DIFF @diff_lines;
#   print space and line to separate different TEST directories
    print DIFF "\n\n"; print DIFF ("-") x 95; print DIFF "\n\n";
#   move DIFF back to WORK_DIR
    rename "$TEST_DIR[$i]"."/$diff_file", "$WORK_DIR"."/$diff_file";
}



