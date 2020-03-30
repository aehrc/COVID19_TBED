#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;

my $input;
my $freqs;
my $alignment;
my $alignout;
my $sequences;

GetOptions('input=s'     => \$input,
           'freqs=s'     => \$freqs,
           'alignout=s'  => \$alignout,
           'sequences=s' => \$sequences);

open (FREQ, $freqs) || die "Could not open file: $freqs\n";
my @lines = <FREQ>;
my $i    = 0;
my $pos  = 0;
my $gap  = 100;
do {
    my $line = $lines[$i];
    chomp($line);
    ($pos, $gap) = split("\t",$line);
    $i++;
} until($gap<0.05);

my $upPos = $pos;

$i   = scalar(@lines) - 1;
$pos = 0;
$gap = 100;
do {
    my $line = $lines[$i];
    chomp($line);
    ($pos, $gap) = split("\t",$line);
    $i--;
} until($gap<0.05);

close FREQ;

my $downPos = scalar(@lines) - $pos;

open (INPUT, $input)            || die "Could not open input: $input\n";
open (ALIGNOUT, ">$alignout")   || die "Could not open output: $alignout\n";
open (SEQUENCES, ">$sequences") || die "Could not open output: $sequences\n";

my @in   = <INPUT>;
my $in   = join('',@in);
my @seqs = split(">", $in);
shift @seqs;

foreach my $entry (@seqs) {
    my @lines    = split("\n", $entry);
    my ($header) = shift(@lines);
    my $seq      = join('',@lines);

    $seq = substr($seq, $upPos);
    $seq = substr($seq, 0, -${downPos});

    print ALIGNOUT ">${header}\n${seq}\n";

    ($seq) =~ s/-//g;

    print SEQUENCES ">${header}\n${seq}\n";
}
close INPUT;
close ALIGNOUT;
close SEQUENCES;

exit;