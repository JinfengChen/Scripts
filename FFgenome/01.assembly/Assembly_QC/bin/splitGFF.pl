#!/usr/bin/perl
## Change the file name in this script and run it in bin. 
## GFF file or BED file will be split into files for each chromosome in a chr directory.

###split rice gff file of gene anontation and TE annotation into chromosome files
my $gff=$ARGV[0];
dealgff($gff);

######################################
sub dealgff {
my ($gff)=@_;
chomp $gff;
my $dir=$gff.".chr";
`mkdir $dir` unless (-d $dir);
open IN, "$gff" or die "$!";
while(<IN>){
    next if ($_ eq "");
    next if ($_ =~/^##/);
    my @unit=split("\t", $_);
    my $line=join("\t",@unit);
    $chr=$unit[0];
    open OUT, ">>$dir/$chr.gff"; 
         print OUT "$line"; 
    close OUT;
}
close IN;
}
