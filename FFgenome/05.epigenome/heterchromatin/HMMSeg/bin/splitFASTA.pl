#!/usr/bin/perl
## FASTA file will be split into files for each chromosome in a chr directory.

###split rice fasta file of gene anontation and TE annotation into chromosome files
my $fasta=$ARGV[0];
dealfasta($fasta);


sub dealfasta {
my ($fasta)=@_;
chomp $fasta;
my $dir=$fasta.".chr";
`mkdir $dir` unless (-d $dir);
$/=">";
open IN, "$fasta" or die "$!";
while(<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $head=shift @unit;
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    open OUT, ">$dir/$head";
         print "$head\n"; 
         print OUT ">$head\n$seq\n"; 
    close OUT;
}
close IN;
$/="\n";
}
