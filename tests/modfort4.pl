#!/usr/bin/perl

# remember to include all fort.4 files in the arguments!
# i.e., ./modfort4.pl fort.4

#back-up fort.4 files using ~
$^I = "~";

my %todo;
my $temp = undef;
my @lines;
my @inputfiles;
my $i = 0;

while (<>){

    $i = $.+2;
    $lines[$i] = $_;

# check previous 2 lines to avoid unintentional matches
    if (($lines[$i-2] =~ /mseed/)&& ($lines[$i-1] =~ /\d+/)){
	if ($lines[$i] =~ /io unit/){
##           do nothing
	} else {
	    print "io unit\n";
	    print "2\n";
	}
    }

    if (($lines[$i-2] =~ /L_tor_table/)&& ($lines[$i-1] =~ /false|true/)){
	if ($lines[$i] =~ /L_vib_table/){
##           do nothing
	} else {
	    print "L_vib_table L_bend_table L_vdW_table L_elect_table\n";
	    print ".false.     .false.      .false.    .false.\n";
	}
    }

# change current line (Elect_field(ibox=1,nbox))
    if ($lines[$i] =~/temp\s+fqtemp/){
	chomp;
	$_ = $_ . "(ibox=1,nbox)\n";
    }
    if (($lines[$i-1] =~ /temp\s+fqtemp/)&& ($lines[$i] =~ /\d+/)){
	if ($lines[$i] =~ /L_vib_table/){
##           do nothing
	} else {
	    chomp;
	    $_ = $_ . " 0.0d0\n";
	}
    }

# change current line (Elect_field(ibox=1,nbox))
    if ($lines[$i] =~/temp\s+fqtemp/){
	chomp;
	$_ = $_ . "(ibox=1,nbox)\n";
    }
    if (($lines[$i-1] =~ /temp\s+fqtemp/)&& ($lines[$i] =~ /\d+/)){
	chomp;
	$_ = $_ . " 0.0d0\n";
    }

# change current line (iheatcapacity)
    if ($lines[$i] =~/imv\s+iratio/){
	chomp;
	$_ = $_ . "iheatcapacity\n";
    }
    if (($lines[$i-1] =~ /imv\s+iratio/)&& ($lines[$i] =~ /\d+/)){
	chomp;
	$_ = $_ . "   1000000\n";
    }

# match anything except newline (remove blank line)
    if(/./){
	print $_;
    }
}
