#!/usr/local/bin/perl

use strict;
use warnings;
use Carp;

my $usage = "usage: $0 min_chain_length < file.aligncoords\n\n";

my $min_chain_len = $ARGV[0] or die $usage;

my $chain_text = "";

my %gene_models;

while (<STDIN>) {
    unless (/\w/) { next; }
    
    if (/^\#/) {
        if ($chain_text) {
            &process_chain_text($chain_text);
        }
        $chain_text = $_;
    }
    else {
        $chain_text .= $_;
    }
}
if ($chain_text) {
    &process_chain_text($chain_text); # get last one.
}

foreach my $gene_model (keys %gene_models) {
    print "$gene_model\n";
}


exit(0);


####
sub process_chain_text {
    my $text = shift;

    my @entries = split (/\n/, $text);
    
    my $header = shift @entries;
    $header =~ /num aligned pairs: (\d+)/ or confess "Error, cannot determine num aligned pairs from header\n$header\n\n$text\n\n!!";

    my $num_pairs = $1;

    if ($num_pairs < $min_chain_len) { 
        return;
    }

    foreach my $line (@entries) {
        my @x = split (/\t/, $line);
        my ($chr_A, $gene_A, $chr_B, $gene_B) = ($x[0], $x[1], $x[4], $x[5]);

        $gene_models{"$chr_A; $gene_A"} = 1;
        $gene_models{"$chr_B; $gene_B"} = 1;
        
    }


}
