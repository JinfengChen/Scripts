#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"table:s","fasta:s","help");


my $help=<<USAGE;
perl $0 --table --fasta
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my $ref=readtable($opt{table});
my ($refseq,$refseq2)=getfastaseq($opt{fasta});

my $total;
open OUT2, ">Rice.RLK-LRR.type" or die "$!";
foreach(sort keys %$ref){
    my $n=@{$ref->{$_}};
    $total+=$n;
    my $type=$_;
    print "$type\t$n\n";
    open OUT, ">Rice.$type.fa" or die "$!";
    foreach my $g (@{$ref->{$_}}){
       my $locus= $1 if ($g=~/(.*)\.\d+/);
       if (exists $refseq->{$g}){
          print OUT ">$g\n$refseq->{$g}\n";
          print OUT2 "$g\t$type\n";
       }elsif(exists $refseq2->{$locus}){
          print OUT ">$refseq2->{$locus}->[0]\n$refseq2->{$locus}->[1]\n";
          print OUT2 "$refseq2->{$locus}->[0]\t$type\n";
       }else{
          print "$g\tNot found in fasta\n";
       }
    }
    close OUT;
}
close OUT2;
print "Total\t$total\n";

##############################
sub readtable
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    #print "$unit[0]\n";
    next unless ($unit[0]=~/Orysa4\|(.*)/);
    #print "$unit[0]";
    $gene="LOC"."_".$1;
    if ($unit[1]=~/LRR/){
        #print "$unit[1]\n";
        push (@{$hash{$unit[1]}},$gene);
    } 
}
close IN;
return \%hash;
}

sub getfastaseq
{
$/=">";
my %hash;
my %hash2;
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
    my $locus= $1 if ($head=~/(.*)\.\d+/);
    $hash2{$locus}=[$head,$seq];
    $hash{$head}=$seq;
}
$/="\n";
return (\%hash,\%hash2);
}

