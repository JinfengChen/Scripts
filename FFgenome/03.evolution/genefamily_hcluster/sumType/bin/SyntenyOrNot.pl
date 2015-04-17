#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"syntid:s","orthid:s","outdir:s","project:s","help");


my $help=<<USAGE;
split ortholog gene into synteny or non-synteny
perl $0 --syntid --orthid --outdir
--syntid: gene id file preduced by synteny_pipeline
--orthid: gene id file of ortholog gene in this species.
--outdir: output dir for synteny and non-synteny gene
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

`mkdir $opt{outdir}`;

my $refsyn=readtable($opt{syntid});
& ort($refsyn,$opt{orthid},$opt{outdir},$opt{project});
###################
sub ort
{
my ($ref,$orth,$out,$name)=@_;
my ($n1,$n2);
my $syn="$out/$name.syn.ortholog.table";
my $nonsyn="$out/$name.nonsyn.ortholog.table";
open OUT1, ">$syn" or die "$!";
open OUT2, ">$nonsyn" or die "$!";
open IN, "$orth" or die "$!";
while(<IN>){
   chomp $_;
   next if ($_=~/^$/);
   my @unit=split("\t",$_);
   my $gene= $1 if ($unit[0]=~/(.*)\_(\w+)$/);
   #print "$gene\n";
   if (exists $ref->{$gene}){
      $n1++;
      print OUT1 "$_\n";
   }else{
      $n2++;
      print OUT2 "$_\n";
   }
}
close IN;
close OUT1;
close OUT2;
print "$name\tSynteny\tNonSynteny\n";
print "$name\t$n1\t$n2\n";
}

sub readtable
{
my ($file,$col)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash{$unit[$col]}=1;
}
close IN;
return \%hash;
}


