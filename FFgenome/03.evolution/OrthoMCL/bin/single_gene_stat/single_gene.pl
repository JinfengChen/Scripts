open(IN,"Tobacco.filled.filter_2k.glean.pep.iprscan.gene.GO");
while($line=<IN>){
	@seg=split(/\s+/,$line);
	$gene_name=shift @seg;
	shift @seg;
	$go=join(" ",@seg);
	$hash_go{$gene_name}=$go;
}


open(IN,"all_orthomcl.out");
while($line=<IN>){
	@gene=split(/\s+/,(split(/\:\s+/,$line))[1]);
	foreach $gene(@gene){
		$gene=~/([^\(]+)\(([A-Z]{5})\)/;
		($gene_name,$spec)=($1,$2);
		#$gene_name=~/(\S+)_[A-Z]{5}/;
		#$gene_name=$1;
		if($spec eq "NICTA"){
			$hash_gene{$gene_name}++;
		}
	}
}

open(IN,"all.gg");
while($line=<IN>){
	$spec=(split(/\:\s+/,$line))[0];
	if($spec eq "NICTA"){
		last;
	}
}
@gene=split(/\s+/,(split(/\:\s+/,$line))[1]);
foreach $gene(@gene){
	if(!$hash_gene{$gene}){
		$gene=~/(\S+)_[A-Z]{5}/;
		$gene=$1;
		$go=$hash_go{$gene};
		print "$gene\t$go\n";
	}
}

