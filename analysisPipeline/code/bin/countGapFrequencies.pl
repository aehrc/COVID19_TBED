#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;

my $input;
my $output;

GetOptions('input=s'  => \$input,
           'output=s' => \$output);

open (INPUT, $input) || die "Could not open file: $input\n";
open (OUTPUT, ">$output") || die "Could not open file: $output\n";

my $no = 0;

my @fileLines = <INPUT>;
my $fileData  = join('', @fileLines);

my @seqs = split(">", $fileData);
shift @seqs;

my @counts     = ();
$counts[29947] = 0;

foreach my $entry (@seqs) {
    my @lines    = split("\n", $entry);
    my ($header) = shift(@lines);
    my $seq      = join('', @lines);

    while ($seq =~ m/-/g) {
        my $position = pos($seq);
        if ($counts[$position-1]) { $counts[$position-1] = $counts[$position-1] + 1; }
        else                      { $counts[$position-1] = 1; }
    }
}
close INPUT;

for (my $i = 0; $i<@counts; $i++) {
    my $freq = '';
    if ($counts[$i]) { $freq = $counts[$i]/scalar@seqs; }
    else             { $freq = 0; }
    print OUTPUT "$i\t$freq\n";
}
close OUTPUT;

exit;
