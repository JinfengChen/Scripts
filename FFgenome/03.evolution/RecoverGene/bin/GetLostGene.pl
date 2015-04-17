#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"domain:s","fasta:s","bed:s","pfam:s","project:s","help");


my $help=<<USAGE;
perl $0 --bed --fasta --gff
--domain: file of pfam domain match on scaffold
--fasta: fasta file of scaffold
--bed: bed format of gene annotation, contains only mRNA feature here
--pfam: the pfam domain which we need to find
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my $fgenesh="/home/biosoftware/fgenesh/fgenesh";
my $convert="/home/jfchen/FFproject/tools/bin/predict_convert.pl";
my $refseq=getfastaseq($opt{fasta});

my $domainbed="$opt{domain}.bed";
& domain2bed($opt{domain},"$opt{domain}.bed");

=pod
my $BEDtools="/home/jfchen/FFproject/tools/BEDTools/bin";
system("$BEDtools/windowBed -a $domainbed -b $opt{bed} -w 500 > $opt{pfam}.windowsBED");
my $refBED=windowsBED("$opt{pfam}.windowsBED");
foreach(keys %$refBED){
   print "$_\t$refBED->{$_}\n";
   
}
=cut

& getflankingseq($domainbed,$refseq);


##############################
sub getflankingseq
{
my ($bed,$refseq)=@_;
open IN, "$bed" or die "$!";
while(<IN>){
   chomp $_;
   next if ($_=~/^$/);
   my @unit=split("\t",$_);
   my $start=$unit[1]-2500 > 0 ? $unit[1]-2500 : 0 ;
   my $len=$unit[2]-$unit[1]+1+5000;
   my $subseq=substr($refseq->{$unit[0]},$start,$len);
   my $line=">$unit[3]\n$subseq\n";
   writefile("temp.fa",$line);
   `$fgenesh/fgenesh $fgenesh/Monocots temp.fa > temp.fgenesh`;
   `perl $convert --predict fgenesh temp.fgenesh temp.fa`;
   `cat temp.fgenesh.gff >> all.fgenesh.gff`;
   `cat temp.fgenesh >> all.fgenesh`; 
   `cat temp.fgenesh.cds >> all.fgenesh.cds`;
   `cat temp.fgenesh.pep >> all.fgenesh.pep`;
   `rm temp.*`;
}
close IN;
}

sub writefile
{
my ($file,$line)=@_;
open OUT, ">$file" or die "$!";
    print OUT "$line";
close OUT;
}


sub domain2bed
{
my ($file,$bed)=@_;
my $count;
open OUT, ">$bed" or die "$!";
open IN, "$file" or die "$!";
while(<IN>){
   chomp $_;
   next if ($_=~/^$/);
   my $line=$_;
   if ($line=~/$opt{pfam}/){
     $count++;
     my @unit=split("\t",$_);
     my $fbox=$unit[2].".$count";
     print OUT "$unit[0]\t$unit[3]\t$unit[4]\t$fbox\t$unit[5]\t$unit[1]\n";

   }
}
close IN;
close OUT;
}

sub windowsBED
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    my $len=$unit[10]-$unit[9]+1;
    if (exists $hash{$unit[3]}){
        $hash{$unit[3]}+=$len;
    }else{
        $hash{$unit[3]} =$len;
    }
}
close IN;
return \%hash;
}


sub revcom
{
my ($seq)=@_;
my $rev=reverse $seq;
$rev=~tr/ATGCatgc/TACGtacg/;
return $rev;
}


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



