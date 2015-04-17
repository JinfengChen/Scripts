#!/usr/local/bin/perl -w

use strict;

use Bio::SearchIO;

my $in = new Bio::SearchIO(-format => 'blast', -file => $ARGV[0]);
while (my $result = $in->next_result ) {
    while (my $hit = $result->next_hit ) {
        while (my $hsp = $hit->next_hsp ) {
            my $pcid = $hsp->percent_identity;
            my $len = $hsp->length('total');
            my $mismatches = $len - $hsp->num_identical -
$hsp->gaps('total');
            #Find runs of '-' in query, hit strings
            my @gap_list = (($hsp->query_string . $hsp->hit_string) =~
/-+/g);
            my $ngaps = @gap_list;
            my @fields = (
                $result->query_name, $hit->name, sprintf("%.2f",$pcid),
                $len, $mismatches, $ngaps,
                $hsp->start('query'), $hsp->end('query'),
                $hsp->start('hit'), $hsp->end('hit'),
                $hsp->evalue, $hsp->bits
            );
            print join("\t", @fields), "\n";
        }
    }
}
