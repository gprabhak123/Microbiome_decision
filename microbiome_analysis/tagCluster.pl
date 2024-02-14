#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my $accFile = $ARGV[0];
open (my $fh, "<", $accFile) || die $!;
chomp(my @fams = <$fh>);
close $fh;

#print Data::Dumper->Dump([\@fams],[qw(*fams)]);
my $treeFile = $ARGV[1];

my $outTree  = $treeFile;
#$outTree =~ s/\.tree/tagged\.tree/;


open (my $fh2, "<", $treeFile) || die $!;
my $tree = <$fh2>;
close $fh2;


foreach my $fam (@fams) {
  my $fam2 = $fam;
  $fam2 =~ s/\./\\\./g;
  #print "$fam2\n";
  $tree =~ s/${fam2}:/CL3_${fam}:/g;
}

open (my $ofh, ">", $outTree) || die $!;
print $ofh $tree;
close $ofh;

