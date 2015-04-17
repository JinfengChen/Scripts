#!/usr/local/bin/perl

use strict;
use warnings;
use Carp;

my $usage = "usage: $0 min_chain_length < file.aligncoords\n\n";

my $min_chain_len = $ARGV[0] or die $usage;

my $chain_text = "";

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
    
    print $text . "\n";
    

}
