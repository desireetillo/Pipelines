#!/bin/perl
use strict;
use File::Basename;

my $ins = $ARGV[0];
my @files = split("\,", $ins);
my @header=qw(Sample NReads NReadsTrim MappedReads(All) mm10_MappedReads(All) mm10_Reads(>120) spikein_Reads NRF PBC1 PBC2);

print join("\t", @header), "\n";

foreach my $file(@files){
    chomp $file;
    my $sampName=basename($file);
    $sampName=~ s/\.qcmetrics//;
    my $out =`cat $file | ./Scripts/perl_lib/transpose.pl -q | tail -n 1`;
    chomp $out;
    print "$sampName\t$out\n";
}

