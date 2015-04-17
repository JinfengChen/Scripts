#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"gff:s","domain:s","iprscan:s","help");


my $help=<<USAGE;
This script deals with super_scaffold.domain. Find out domains which are not in gene set.
Run: perl FindDomainGene.pl -gff FF.gff -domain super_scaffold.domain -iprscan FF.iprscan > log 2> log2 &
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my $refgff =parseGFF($opt{gff});
our $refpfam=pfam($opt{iprscan});

open IN, "$opt{domain}" or die "$!";
while(<IN>){
    chomp $_;
    my $domain=$_;
    my $geneinf;
    my @unit=split("\t",$_);
    my $pos =[$unit[3],$unit[4]];
    if (exists $refgff->{$unit[0]}){
       my $refarray=$refgff->{$unit[0]};
       $geneinf=findoverlap($pos,$refarray);
    }else{
       print "$domain\tScaffold not in GFF file\n";
    }
    print "$domain\t$geneinf\n";

}
close IN;
#############################################################
sub pfam
{
my ($iprscan)=@_;
my %hash;
open IN, "$iprscan" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    if ($unit[3] eq "HMMPfam"){ 
       if (exists $hash{$unit[0]}){
          $hash{$unit[0]}.="\t$unit[4]\t$unit[5]";
       }else{
          $hash{$unit[0]}="$unit[4]\t$unit[5]";
       }
    }
}
close IN;
return \%hash;
}



sub findoverlap
{
my ($pair,$array)=@_;
my $gene="No Hit";
#print "$pair->[0]\t$pair->[1]\n";
foreach(sort {$a->[3] <=> $b->[3]} @$array){
   #print "$_->[3]\t$_->[4]\n";
   if ($_->[3] >= $pair->[1]){
       last;
   }elsif($_->[4] <= $pair->[0]){
       next;
   }else{
       #$gene=join("\t",@$_);
       my $id;
       my $pfam="NA";
       if ($_->[8]=~/ID=(.*);/ or $_->[8]=~/ID=(.*)/){
          $id=$1;
          if (exists $refpfam->{$id}){
             $pfam=$refpfam->{$id};
          }
       }
       $gene="$id\t$_->[3]\t$_->[4]\t$_->[6]\t$pfam";
       return $gene;
   }
}
return $gene;
}


sub parseGFF
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $seq;   ##Scaffold
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/mRNA/){
        $seq=$unit[0];
        if ($unit[8]=~/ID=(.*);/ or $unit[8] =~/ID=(.*)/){
            $id=$1;
        }
        if (exists $hash{$seq}){
           my $temp=$hash{$seq};
           push (@$temp,[@unit]);
           $hash{$seq}=$temp;
        }else{
           my @temp;
           push (@temp,[@unit]);
           $hash{$seq}=\@temp;
        }
    }

}
close IN;
return \%hash;
}

