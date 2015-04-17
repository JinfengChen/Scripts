#!/usr/bin/perl
## Change the file name in this script and run it in bin. 
## GFF file or BED file will be split into files for each chromosome in a chr directory.

###split rice gff file of gene anontation and TE annotation into chromosome files
my $align=$ARGV[0];
splitDAG($align);

######################################
sub splitDAG
{
my ($align)=@_;
chomp $align;
my $dir=$align.".chr";
`mkdir $dir` unless (-d $dir);
my %hash;
my $block;
my $qry;
my $ref;
my $chr;
open IN, "$align" or die "$!";
while(<IN>){
    chomp $_;
    my $line=$_;
    if ($_=~/\#\# alignment\s+(\w+)\s+vs\.\s+(\w+)\s+Alignment/) {
        $qry=$1;
        $ref=$2;
        if ($ref=~/(\d+)/){
           if (length $1 == 2){
             $chr="chr".$1;
           }else{
             $chr="chr0".$1;
           }
        }
        open OUT, ">>$dir/$chr.align";
        print OUT "$line\n";
        close OUT;    
        #print "$block\t$qry\t$ref\n";
    }elsif($_=~/$ref/){
        open OUT, ">>$dir/$chr.align";
        print OUT "$line\n";
        close OUT;
    }
}
close IN;
return \%hash;
}


