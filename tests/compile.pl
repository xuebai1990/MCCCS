#!/usr/bin/perl

# modify control.inc, external.inc, and Makefile

#use Shell qw(make);

my $CODE_DIR = '/scratch2/maerzkek/MCCCS-2010-PARALLEL';
chomp(my $WORK_DIR = `pwd`);

chdir $CODE_DIR;

if (! open M, "<", "Makefile"){
    die "cannot open Makefile: $!";
}
if (! open M2, ">", "Makefile2"){
    die "cannot open Makefile2: $!";
}
open C, "<", "control.inc";
open C2, ">", "control.inc2";
open E, "<", "external.inc";
open E2, ">", "external.inc2";

my @splits;

system 'make clean';

system 'cp -f control.inc.TEST control.inc';

##########################################################
#NpT w/Ewald 

while (<M>){
#    if (/PROGRAM\s+=\s+\w+/){  # more complicated match pattern
    if ((/PROGRAM/)&&(/=/)){
# use split to change the last 'word' regardless of its initial value
	@splits = split /\s+/, $_;
	$_ = "$splits[0]" . "       $splits[1]". " NpTewald\n";
#	print $_;
    }
    print M2 $_;
}

rename "Makefile2", "Makefile";

while (<C>){
    if ((/parameter/)&&(/lnpt/)){
	s/false/true/;
    }
    if ((/parameter/)&&(/lgibbs/)){
	s/true/false/;
    }
    if ((/parameter/)&&(/lgrand/)){
	s/true/false/;
    }
    if ((/parameter/)&&(/lcutcm/)){
        s/false/true/;
    }
    if ((/parameter/)&&(/lewald/)){
	s/false/true/;
    }
    if ((/parameter/)&&(/ltailc/)){
	s/false/true/;
    }
    if ((/parameter/)&&(/lshift/)){
	s/true/false/;
    }
    if ((/parameter/)&&(/llj/)){
	s/false/true/;
    }
    if ((/parameter/)&&(/lgaro/)){
	s/true/false/;
    }
    if ((/parameter/)&&(/lionic/)){
	s/true/false/;
    }
    print C2 $_;
}

rename "control.inc2", "control.inc";

while (<E>){
    if ((/parameter/)&&(/lelect_field/)){
	s/true/false/;
    }

    print E2 $_;
}

rename "external.inc2", "external.inc";

# close files
close M;
close M2;
close C;
close C2;
close E;
close E2;

# make the executable
system 'make';

##########################################################
#NpT w/o Ewald (for tab potentials) 

# re-open files for further editing
open M, "<", "Makefile";
open M2, ">", "Makefile2";
open C, "<", "control.inc";
open C2, ">", "control.inc2";
open E, "<", "external.inc";
open E2, ">", "external.inc2";

system 'make clean';

while (<M>){
#    if (/PROGRAM\s+=\s+\w+/){  # more complicated match pattern
    if ((/PROGRAM/)&&(/=/)){
# use split to change the last 'word' regardless of its initial value
        @splits = split /\s+/, $_;
        $_ = "$splits[0]" . "       $splits[1]". " NpTtab\n";
#       print $_;
    }
    print M2 $_;
}

rename "Makefile2", "Makefile";

while (<C>){
    if ((/parameter/)&&(/lnpt/)){
        s/false/true/;
    }
    if ((/parameter/)&&(/lgibbs/)){
        s/true/false/;
    }
    if ((/parameter/)&&(/lgrand/)){
        s/true/false/;
    }
    if ((/parameter/)&&(/lewald/)){
        s/true/false/;
    }
    if ((/parameter/)&&(/lcutcm/)){
        s/false/true/;
    }
    if ((/parameter/)&&(/ltailc/)){
        s/false/true/;
    }
    if ((/parameter/)&&(/lshift/)){
        s/true/false/;
    }
    if ((/parameter/)&&(/llj/)){
        s/false/true/;
    }
    if ((/parameter/)&&(/lgaro/)){
        s/true/false/;
    }
    if ((/parameter/)&&(/lionic/)){
        s/true/false/;
    }
    print C2 $_;
}

rename "control.inc2", "control.inc";

while (<E>){
    if ((/parameter/)&&(/lelect_field/)){
        s/true/false/;
    }

    print E2 $_;
}

rename "external.inc2", "external.inc";

# close files
close M;
close M2;
close C;
close C2;
close E;
close E2;

# make the executable
system 'make';

##########################################################
#NpT for solid 

# re-open files for further editing
open M, "<", "Makefile";
open M2, ">", "Makefile2";
open C, "<", "control.inc";
open C2, ">", "control.inc2";
open E, "<", "external.inc";
open E2, ">", "external.inc2";

system 'make clean';

while (<M>){
#    if (/PROGRAM\s+=\s+\w+/){  # more complicated match pattern
    if ((/PROGRAM/)&&(/=/)){
# use split to change the last 'word' regardless of its initial value
        @splits = split /\s+/, $_;
        $_ = "$splits[0]" . "       $splits[1]". " NpTsolid\n";
#       print $_;
    }
    print M2 $_;
}

rename "Makefile2", "Makefile";

while (<C>){
    if ((/parameter/)&&(/lnpt/)){
        s/false/true/;
    }
    if ((/parameter/)&&(/lgibbs/)){
        s/true/false/;
    }
    if ((/parameter/)&&(/lgrand/)){
        s/true/false/;
    }
    if ((/parameter/)&&(/lcutcm/)){
        s/true/false/;
    }
    if ((/parameter/)&&(/lewald/)){
        s/false/true/;
    }
    if ((/parameter/)&&(/ltailc/)){
        s/false/true/;
    }
    if ((/parameter/)&&(/lshift/)){
        s/true/false/;
    }
    if ((/parameter/)&&(/llj/)){
        s/false/true/;
    }
    if ((/parameter/)&&(/lgaro/)){
        s/true/false/;
    }
    if ((/parameter/)&&(/lionic/)){
        s/true/false/;
    }
    print C2 $_;
}

rename "control.inc2", "control.inc";

while (<E>){
    if ((/parameter/)&&(/lelect_field/)){
        s/true/false/;
    }

    print E2 $_;
}

rename "external.inc2", "external.inc";

# close files
close M;
close M2;
close C;
close C2;
close E;
close E2;

# make the executable
system 'make';

##########################################################
#NVT field 

# re-open files for further editing
open M, "<", "Makefile";
open M2, ">", "Makefile2";
open C, "<", "control.inc";
open C2, ">", "control.inc2";
open E, "<", "external.inc";
open E2, ">", "external.inc2";

system 'make clean';

while (<M>){
    if (/PROGRAM\s+=\s+\w+/){
# use split to change the last 'word' regardless of its initial value
	@splits = split /\s+/, $_;
	$_ = "$splits[0]" . "       $splits[1]". " NVTfield\n";
#	print $_;
    }
    print $_;
    print M2 $_;
}

rename "Makefile2", "Makefile";

while (<C>){
    if ((/parameter/)&&(/lnpt/)){
	s/true/false/;
    }
    if ((/parameter/)&&(/lgibbs/)){
	s/true/false/;
    }
    if ((/parameter/)&&(/lgrand/)){
	s/true/false/;
    }
    if ((/parameter/)&&(/lcutcm/)){
        s/false/true/;
    }
    if ((/parameter/)&&(/lewald/)){
	s/false/true/;
    }
    if ((/parameter/)&&(/ltailc/)){
	s/false/true/;
    }
    if ((/parameter/)&&(/lshift/)){
	s/true/false/;
    }
    if ((/parameter/)&&(/llj/)){
	s/false/true/;
    }
    if ((/parameter/)&&(/lgaro/)){
	s/true/false/;
    }
    if ((/parameter/)&&(/lionic/)){
	s/true/false/;
    }
    print C2 $_;
}

rename "control.inc2", "control.inc";

while (<E>){
    if ((/parameter/)&&(/lelect_field/)){
	s/false/true/;
    }

    print E2 $_;
}

rename "external.inc2", "external.inc";

# make the executable
system 'make';

close M;
close M2;
close C;
close C2;
close E;
close E2;

##########################################################
#NVT-Gibbs w/Ewald

# re-open files for further editing
open M, "<", "Makefile";
open M2, ">", "Makefile2";
open C, "<", "control.inc";
open C2, ">", "control.inc2";
open E, "<", "external.inc";
open E2, ">", "external.inc2";

system 'make clean';

while (<M>){
    if (/PROGRAM\s+=\s+\w+/){
# use split to change the last 'word' regardless of its initial value
	@splits = split /\s+/, $_;
	$_ = "$splits[0]" . "       $splits[1]". " NVT-Gibbs\n";
#	print $_;
    }
    print $_;
    print M2 $_;
}

rename "Makefile2", "Makefile";

while (<C>){
    if ((/parameter/)&&(/lnpt/)){
	s/true/false/;
    }
    if ((/parameter/)&&(/lgibbs/)){
	s/false/true/;
    }
    if ((/parameter/)&&(/lgrand/)){
	s/true/false/;
    }
    if ((/parameter/)&&(/lcutcm/)){
        s/false/true/;
    }
    if ((/parameter/)&&(/lewald/)){
	s/false/true/;
    }
    if ((/parameter/)&&(/ltailc/)){
	s/false/true/;
    }
    if ((/parameter/)&&(/lshift/)){
	s/true/false/;
    }
    if ((/parameter/)&&(/llj/)){
	s/false/true/;
    }
    if ((/parameter/)&&(/lgaro/)){
	s/true/false/;
    }
    if ((/parameter/)&&(/lionic/)){
	s/true/false/;
    }
    print C2 $_;
}

rename "control.inc2", "control.inc";

while (<E>){
    if ((/parameter/)&&(/lelect_field/)){
	s/true/false/;
    }

    print E2 $_;
}

rename "external.inc2", "external.inc";

# make the executable
system 'make';

close M;
close M2;
close C;
close C2;
close E;
close E2;

##########################################################
#NVT-Gibbs w/o Ewald

# re-open files for further editing
open M, "<", "Makefile";
open M2, ">", "Makefile2";
open C, "<", "control.inc";
open C2, ">", "control.inc2";
open E, "<", "external.inc";
open E2, ">", "external.inc2";

system 'make clean';

while (<M>){
    if (/PROGRAM\s+=\s+\w+/){
# use split to change the last 'word' regardless of its initial value
	@splits = split /\s+/, $_;
	$_ = "$splits[0]" . "       $splits[1]". " NVT-Gibbs-noq\n";
#	print $_;
    }
    print $_;
    print M2 $_;
}

rename "Makefile2", "Makefile";

while (<C>){
    if ((/parameter/)&&(/lnpt/)){
	s/true/false/;
    }
    if ((/parameter/)&&(/lgibbs/)){
	s/false/true/;
    }
    if ((/parameter/)&&(/lgrand/)){
	s/true/false/;
    }
    if ((/parameter/)&&(/lcutcm/)){
        s/false/true/;
    }
    if ((/parameter/)&&(/lewald/)){
	s/true/false/;
    }
    if ((/parameter/)&&(/ltailc/)){
	s/false/true/;
    }
    if ((/parameter/)&&(/lshift/)){
	s/true/false/;
    }
    if ((/parameter/)&&(/llj/)){
	s/false/true/;
    }
    if ((/parameter/)&&(/lgaro/)){
	s/true/false/;
    }
    if ((/parameter/)&&(/lionic/)){
	s/true/false/;
    }
    print C2 $_;
}

rename "control.inc2", "control.inc";

while (<E>){
    if ((/parameter/)&&(/lelect_field/)){
	s/true/false/;
    }

    print E2 $_;
}

rename "external.inc2", "external.inc";

# make the executable
system 'make';

close M;
close M2;
close C;
close C2;
close E;
close E2;


##########################################################
# RPLC

# re-open files for further editing
open M, "<", "Makefile";
open M2, ">", "Makefile2";
open C, "<", "control.inc";
open C2, ">", "control.inc2";
open E, "<", "external.inc";
open E2, ">", "external.inc2";

system 'make clean';

while (<M>){
    if (/PROGRAM\s+=\s+\w+/){
# use split to change the last 'word' regardless of its initial value
	@splits = split /\s+/, $_;
	$_ = "$splits[0]" . "       $splits[1]". " RPLC\n";
#	print $_;
    }
    print $_;
    print M2 $_;
}

rename "Makefile2", "Makefile";

while (<C>){
    if ((/parameter/)&&(/lnpt/)){
	s/false/true/;
    }
    if ((/parameter/)&&(/lgibbs/)){
	s/false/true/;
    }
    if ((/parameter/)&&(/lgrand/)){
	s/true/false/;
    }
    if ((/parameter/)&&(/lcutcm/)){
        s/false/true/;
    }
    if ((/parameter/)&&(/lewald/)){
	s/false/true/;
    }
    if ((/parameter/)&&(/ltailc/)){
	s/true/false/;
    }
    if ((/parameter/)&&(/lshift/)){
	s/true/false/;
    }
    if ((/parameter/)&&(/llj/)){
	s/false/true/;
    }
    if ((/parameter/)&&(/lgaro/)){
	s/true/false/;
    }
    if ((/parameter/)&&(/lionic/)){
	s/false/true/;
    }
    print C2 $_;
}

rename "control.inc2", "control.inc";

while (<E>){
    if ((/parameter/)&&(/lelect_field/)){
	s/true/false/;
    }

    print E2 $_;
}

rename "external.inc2", "external.inc";

# make the executable
system 'make';

close M;
close M2;
close C;
close C2;
close E;
close E2;

##########################################################
#NpT w/Garofalini

# re-open files for further editing
open M, "<", "Makefile";
open M2, ">", "Makefile2";
open C, "<", "control.inc";
open C2, ">", "control.inc2";
open E, "<", "external.inc";
open E2, ">", "external.inc2";

system 'make clean';

while (<M>){
    if (/PROGRAM\s+=\s+\w+/){
# use split to change the last 'word' regardless of its initial value
	@splits = split /\s+/, $_;
	$_ = "$splits[0]" . "       $splits[1]". " NpTgaro\n";
#	print $_;
    }
    print $_;
    print M2 $_;
}

rename "Makefile2", "Makefile";

while (<C>){
    if ((/parameter/)&&(/lnpt/)){
	s/false/true/;
    }
    if ((/parameter/)&&(/lgibbs/)){
	s/true/false/;
    }
    if ((/parameter/)&&(/lgrand/)){
	s/true/false/;
    }
    if ((/parameter/)&&(/lcutcm/)){
        s/false/true/;
    }
    if ((/parameter/)&&(/lewald/)){
	s/true/false/;
    }
    if ((/parameter/)&&(/ltailc/)){
	s/true/false/;
    }
    if ((/parameter/)&&(/lshift/)){
	s/false/true/;
    }
    if ((/parameter/)&&(/llj/)){
	s/true/false/;
    }
    if ((/parameter/)&&(/lgaro/)){
	s/false/true/;
    }
    if ((/parameter/)&&(/lionic/)){
	s/false/true/;
    }
    print C2 $_;
}

rename "control.inc2", "control.inc";

while (<E>){
    if ((/parameter/)&&(/lelect_field/)){
	s/true/false/;
    }

    print E2 $_;
}

rename "external.inc2", "external.inc";

# make the executable
system 'make';

close M;
close M2;
close C;
close C2;
close E;
close E2;

# re-open Makefile to rename executable
open M, "<", "Makefile";
open M2, ">", "Makefile2";

system 'make clean';

while (<M>){
    if (/PROGRAM\s+=\s+\w+/){
# use split to change the last 'word' regardless of its initial value
	@splits = split /\s+/, $_;
	$_ = "$splits[0]" . "       $splits[1]". " topmon4\n";
#	print $_;
    }
    print $_;
    print M2 $_;
}

rename "Makefile2", "Makefile";
close M;
close M2;


system 'cp -f control.inc.TEST control.inc';
