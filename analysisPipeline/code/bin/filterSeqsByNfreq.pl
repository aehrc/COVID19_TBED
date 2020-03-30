#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;

my $input;
my $freq;
my $output;
my $nlength;

GetOptions('input=s'   => \$input,
           'freq=s'    => \$freq,
           'nlength=s' => \$nlength,
           'output=s'  => \$output);
           
open (INPUT, $input) || die "Could not open input: $input\n";
open (OUTPUT, ">$output") || die "Could not open output: $output\n";

my @in   = <INPUT>;
my $in   = join ('', @in);
my @seqs = split(">", $in);
shift @seqs;

foreach my $entry (@seqs) {
    my @lines  = split("\n", $entry);
    my $header = shift(@lines);
    my $seq    = join('', @lines);

    my $length = length($seq);
    my $counts = () = $seq =~ /[Nn-]/g;
    my $f      = $counts/$length*100;

    my $flag = 0;
    if ($seq =~ /N{$nlength}/i) { 
        $flag = 1; 
    }

    if ($f < $freq and $flag == 0) { 
        print OUTPUT ">${header}\n${seq}\n"; 
    }
}
close OUTPUT;
close INPUT;

exit;
