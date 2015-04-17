#!/usr/bin/perl
use Getopt::Long;
use strict;
use warnings;

my %opt;

GetOptions (\%opt,"gene:s","blast:s","project:s","help");


my $help=<<USAGE;
perl $0 --gene gene.gff --project opposite

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my $blast=getblast($opt{blast});
my $exon =parseGFF($opt{gene});
& gff2table($opt{gene},"mRNA","$opt{project}.mRNA");
`perl findOverlap.pl $opt{project}.mRNA $opt{project}.mRNA > $opt{project}.mRNA.overlap`;
my $kill=findopposite("$opt{project}.mRNA.overlap",$blast,$exon);
my $n=keys %$kill;
print "opposite gene number: $n\n";
foreach (keys %$kill){
   my $en=@{$exon->{$_}}-1;
   print "$_\t$en\n";
}
cleanGFF($opt{gene},$kill);


##################
sub findopposite
{
my ($file,$blast,$exon)=@_;
my $count=0;
my %gene2cluster;
my %cluster2gene;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    my @unit=split("\t",$_);
    next if ($unit[3] < 2);
    my @gene;
    for(my $i=4;$i<@unit;$i++){
       my @over=split(",",$unit[$i]);
       #print "$over[0]\n";
       push (@gene,$over[0]);
    }
    if (exists $gene2cluster{$unit[0]}){
       foreach (@gene){
           unless (exists $gene2cluster{$_}){
               $gene2cluster{$_}=$gene2cluster{$unit[0]};
               push (@{$cluster2gene{$gene2cluster{$unit[0]}}},$_);
           }
       }
    }else{
       $count++;
       foreach (@gene){
           $gene2cluster{$_}=$count;
           push (@{$cluster2gene{$count}},$_);
       }
    }
}
close IN;

my %kill;
foreach my $c (sort {$a <=> $b} keys %cluster2gene){
   my $gene=join("\t",@{$cluster2gene{$c}});
   #print "$c\t$gene\n";
   foreach my $g (@{$cluster2gene{$c}}){
     if (!exists $blast->{$g} and @{$exon->{$g}} < 4){
        $kill{$g}=1;
     }
   }
}
return \%kill;
}
###################
sub gff2table
{
my ($gff,$feature,$table)=@_;
open IN, "$gff" or die "$!";
open OUT, ">$table.unsort" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    if ($unit[2] eq $feature and ($unit[8] =~/ID=(.*?);/ or $unit[8] =~/Parent=(.*?);/)){
        print OUT "$unit[0]\t$1\t$unit[3]\t$unit[4]\n";
    }
}
close IN;
close OUT;
system ("msort -k 1,n3 $table.unsort > $table");
system ("rm $table.unsort");
}


sub cleanGFF
{
my ($gff,$kill)=@_;
my %hash;  ##hash to store every record by key of Seq_id
open IN, "$gff" or die "$!";
open OUT, ">clean.gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my $line=$_;
    my @unit=split("\t",$_);
    if($unit[2]=~/gene/){
        if ($unit[8]=~/ID=(.*?);/){
            my $id=$1.".1";
            print OUT "$line\n" unless exists ($kill->{$id});
        }
    }elsif($unit[2]=~/mRNA/){
        if ($unit[8]=~/ID=(.*?);/){
            print OUT "$line\n" unless exists ($kill->{$1});
        }
    }elsif($unit[2]=~/CDS/){
        if ($unit[8]=~/Parent=(.*?);/){
            print OUT "$line\n" unless exists ($kill->{$1});
        }
    }
}
close IN;
close OUT;
return \%hash;
}

sub parseGFF
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $id;
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/mRNA/){
        if ($unit[8]=~/ID=(.*?);/ or $unit[8] =~/ID=(.*)/){
            $id=$1;
        }
        push (@{$hash{$id}}, $_);
    }elsif($unit[8] =~/Parent=$id/){
        push (@{$hash{$id}}, $_);
    }

}
close IN;
return \%hash;
}


sub getblast
{
my ($file)=@_;
my %hash;
open IN , "$file" or die "$!";
while(<IN>){
    chomp $_;
    my @unit=split("\t",$_);
    next if exists $hash{$unit[0]};
    next if ($unit[0] eq "$unit[1]");
    next unless ($unit[0]=~/(.*?)\_OBRACH/);
    my $gene=$1;
    $hash{$gene}=$unit[1];
}
close IN; 
return \%hash;
} 
