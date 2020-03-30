#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;

my $input;
my $output;
my $record;

GetOptions('input=s'  => \$input,
           'output=s' => \$output,
           'record=s' => \$record);
           
my %seq2id;

open (INPUT, $input) || die "Could not open file: $input\n";
my $header = '';
my $seq    = '';
while (my $line = <INPUT>) {
    chomp $line;
    if ($line =~ /^>/) { $header = $line; }
    else               { $seq    = $line; }

    if ($header and $seq) {
        ($header) =~ s/>//g;
        if (exists $seq2id{$seq}) { $seq2id{$seq} = join(', ', ($seq2id{$seq}, $header)); }
        else                      { $seq2id{$seq} = $header; }
        $header = '';
        $seq    = '';
    }
}
close INPUT;

open (OUTPUT, ">$output") || die "Could not open file: $output\n";
open (RECORD, ">$record") || die "Could not open file: $record\n";

my $combo=1;
my @keys = keys %seq2id;
foreach my $key (@keys) {
    my @entries = ('0');
    if ($seq2id{$key} =~ /,/) { @entries = split(",",$seq2id{$key}); }
    if (scalar @entries > 1) {
        my $newHead = join('',('CombinedSequences_',$combo));
        print OUTPUT ">${newHead}\n$key\n";
        my $count = scalar(@entries);
        print RECORD "$newHead\t$count\t$seq2id{$key}\n";
        $combo++;
    }
    else { print OUTPUT ">$seq2id{$key}\n$key\n"; }
}
close OUTPUT;
close RECORD;

exit;
