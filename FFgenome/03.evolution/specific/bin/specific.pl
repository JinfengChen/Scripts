#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"query:s","target:s","help");


my $help=<<USAGE;
Find out lineage specific genes for one species from blastp result.
The input names of query and targe are used to find fasta gene seq and blasttable in ../input/gene ../input/blasttable directory.
The only cutoff used here is blastp e-value=1-e10. These gene with no blastp hit will be output to query.specific.fa, So we can use these protein to blastp with nr for further filter.

Run: perl specific.pl -query OB -target OS,SB,BR
-query : the species which you want to find lineage specific genes.
-target: the sepcies which you want to compared with the query
the blasttable file will be OB_BR.blasttable, OB_SB.blasttable, OB_OS.blasttable.

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}
my $queryfa="../input/gene/$opt{query}";
my $refquery=getfastaseq($queryfa);
my @target=split(",",$opt{target});


####delete gene key from hash if blastp hit found in one genome. 
foreach(@target){
    my $blasttable="../input/blasttable/$opt{query}_$_.blasttable";
    print "Deal with $blasttable\n";
    open IN, "$blasttable" or die "$!";
    while (<IN>){
       my @unit=split("\t",$_);
       if (exists $refquery->{$unit[0]}){
           delete $refquery->{$unit[0]};
       }
    }
    close IN;
}

my $counter;
####output the left protein sequence#####################
open OUT, ">$opt{query}.specific.fa" or die "$!";
foreach(keys %$refquery){
    $counter++;
    print OUT ">$_\n$refquery->{$_}\n";
}
close OUT;
print "$counter\n";
########################sub function##################33
sub getfastaseq
{
$/=">";
my %hash;
my ($file)=@_;
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my @temp1=split(" ",$temp);
    my $head=$temp1[0];
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    #print "$head\n";
    $hash{$head}=$seq;
}
$/="\n";
return \%hash;
}
 
