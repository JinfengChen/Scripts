#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"go:s","species:s","help");


my $help=<<USAGE;
perl $0 --go ../input/tigr.all.final.iprscan.gene.GO--species Os > single_family_Os.txt

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

open(IN,"$opt{go}");
while($line=<IN>){
        @seg=split(/\t|\n/,$line);
        $gene_name=shift @seg;
        shift @seg;
        $go=join("\t",@seg);
        $gene_name=$1 if ($gene_name=~/LOC\_(.*?)$/);
        #print "$gene_name\t$go\n";
        $hash_go{$gene_name}=$go;
}

open(IN,"all_orthomcl.out");
while($line=<IN>){
	$go=$genes="";
	%hash=();
	$line=~/ORTHOMCL(\d+)\(\d+ genes,(\d+) taxa\):/;
	($fam_num,$taxa)=($1,$2);
	if($taxa == 1 && $line=~/$opt{species}/){
	        @gene=split(/\s+/,(split(/\:\s+/,$line))[1]);
	        foreach $gene(@gene){
	                $gene=~/(.*?)\((\w+)\)/;
        	        ($gene_name,$spec)=($1,$2);
	                $gene_name=~/(\S+)\_\w+/;
                	$gene_name=$1;
			$genes.="$gene_name ";
	                if(defined $hash_go{$gene_name}){
	                        @seg=split(/\t/,$hash_go{$gene_name});
				foreach $seg(@seg){
					$seg=~/(.*?);/;
					$hash{$1}=$seg;
				}
	                }
	        }
		$go=join(" ",values %hash);
		print "$fam_num\t$genes\t$go\n";			
		
	}
}
