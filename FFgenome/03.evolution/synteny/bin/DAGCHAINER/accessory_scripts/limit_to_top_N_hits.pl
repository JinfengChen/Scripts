#!/usr/local/bin/perl

use strict;

my $usage = "usage: $0 num < dataFile \n\n";
my $num = $ARGV[0] or die $usage;

my @all_data;

while (my $line = <STDIN>) {
    chomp $line;
    my @x = split (/\t/, $line);
    my ($acc, $evalue) = ($x[1], $x[8]);
    push (@all_data, [$acc, $evalue, $line]);
}

## Sort by accession and evalue
@all_data = sort {
    $a->[0] cmp $b->[0]
	||
	$b->[1] <=> $a->[1] 
    } @all_data;

my $curr_acc = "";
my $count = 0;
while (my $entry_aref = shift @all_data) {
    my ($acc, $evalue, $line) = @$entry_aref;
    if ($acc eq $curr_acc) {
	$count++;
	if ($count > $num) {
	    next;
	} else {
	    print "$line\n";
	}
    } else {
	$curr_acc = $acc;
	$count = 1;
	print "$line\n";
    }
}

exit(0);

