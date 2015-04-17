#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"go:s","species:s","help");


my $help=<<USAGE;
perl $0 --go ../input/tigr.all.final.iprscan.gene.GO--species Os > single_gene_Os.txt

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

open(IN,"$opt{go}");
while($line=<IN>){
	@seg=split(/\s+/,$line);
	$gene_name=shift @seg;
        $gene_name=$1 if ($gene_name=~/LOC_(.*)$/);
	shift @seg;
	$go=join(" ",@seg);
	$hash_go{$gene_name}=$go;
}


open(IN,"all_orthomcl.out");
while($line=<IN>){
	@gene=split(/\s+/,(split(/\:\s+/,$line))[1]);
	foreach $gene(@gene){
		$gene=~/(.*?)\((\w+)\)/;
		($gene_name,$spec)=($1,$2);
		#$gene_name=~/(.*)_\w+/;
		#$gene_name=$1;
                #print "$gene_name\n";
		if($spec eq "$opt{species}"){
			$hash_gene{$gene_name}++;
		}
	}
}

open(IN,"all.gg");
while($line=<IN>){
	$spec=(split(/\:\s+/,$line))[0];
	if($spec eq "$opt{species}"){
		last;
	}
}
@gene=split(/\s+/,(split(/\:\s+/,$line))[1]);
foreach $gene(@gene){
	if(!$hash_gene{$gene}){
		$gene=~/(.*?)_\w+/;
		$gene=$1;
		$go=$hash_go{$gene};
		print "$gene\t$go\n";
	}
}

