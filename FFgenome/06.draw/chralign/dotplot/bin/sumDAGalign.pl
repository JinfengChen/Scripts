#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"help");


my $help=<<USAGE;


USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

 
