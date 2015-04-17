#!/usr/bin/perl
use strict;

die "\nFunction: make simple name for swissprot protein ID\n\nUsage: perl simplify_swissprot.pl uniprot_sprot.fasta > uniprot_sprot.fasta.simple\n\n" if(@ARGV==0);

my $uniprot_file = shift;

open(IN, $uniprot_file) || die ("can not open $uniprot_file\n");
$/=">"; <IN>; $/="\n";
while (<IN>) {
	chomp;
	my $head = $1 if(/^(\S+)/);
		
	$/=">";
	my $seq = <IN>;
	chomp $seq;
	$/="\n";

	print ">".$head."\n".$seq;
}
close(IN);
