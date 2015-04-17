opendir DIR, "./" or die "can not open my dir";
foreach my $file (readdir DIR){
    if ($file=~/(.*)\.blast/){
	   system "perl format_blastn.pl $1";
	   system "perl toACT.pl $1";
       print "$file\n";
       push (@seq,$1);
    }
}
close DIR;










