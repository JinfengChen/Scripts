opendir DIR, "./" or die "can not open my dir";
foreach my $file (readdir DIR){
    if ($file=~/(.*)\.blast/){
	   system "perl /home/jfchen/FFproject/FFgenome/02.annotation/09.validiation_geneset/collinearity_validiation/bin/format_blastn.pl $1";
	   system "perl /home/jfchen/FFproject/FFgenome/02.annotation/09.validiation_geneset/collinearity_validiation/bin/toACT.pl $1";
       print "$file\n";
       push (@seq,$1);
    }
}
close DIR;










