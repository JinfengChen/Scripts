open(IN,"Tobacco.filled.filter_2k.glean.pep.iprscan.gene.GO");
while($line=<IN>){
        @seg=split(/\t|\n/,$line);
        $gene_name=shift @seg;
        shift @seg;
        $go=join("\t",@seg);
        $hash_go{$gene_name}=$go;
}

open(IN,"all_orthomcl.out");
while($line=<IN>){
	$go=$genes="";
	%hash=();
	$line=~/ORTHOMCL(\d+)\(\d+ genes,(\d+) taxa\):/;
	($fam_num,$taxa)=($1,$2);
	if($taxa == 1 && $line=~/NICTA/){
	        @gene=split(/\s+/,(split(/\:\s+/,$line))[1]);
	        foreach $gene(@gene){
	                $gene=~/([^\(]+)\(([A-Z]{5})\)/;
        	        ($gene_name,$spec)=($1,$2);
	                $gene_name=~/(\S+)_[A-Z]{5}/;
                	$gene_name=$1;
			$genes.="$gene_name ";
	                if(defined $hash_go{$gene_name}){
	                        @seg=split(/\t/,$hash_go{$gene_name});
				foreach $seg(@seg){
					$seg=~/GO:(\d+);/;
					$hash{$1}=$seg;
				}
	                }
	        }
		$go=join(" ",values %hash);
		print "$fam_num\t$genes\t$go\n";			
		
	}
}
