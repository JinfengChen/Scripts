#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"cluster:s","gff:s","help");


my $help=<<USAGE;
perl $0 --cluster --gff 

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my $refgff=parseGFFpos($opt{gff});
analyze($opt{cluster},$refgff);


####################
sub analyze
{
my ($file,$gff)=@_;
$/=">";
my $title=$1 if ($file=~/(.*)\.cluster/);
open OUT, ">$title.cluster.position" or die "$!";
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if (length $_ < 2);
    my @line=split("\n",$_);
    my $head=shift @line;
    my ($id,$number);
    if ($head=~/ID\=(\d+)\D+Number\=(\d+)/){
       $id=$1;
       $number=$2;
    }
    my %type;
    my $chr;
    my @pos;
    foreach my $clone (@line){
       my @unit=split("\t",$clone);
       $type{$unit[1]}++;
       $chr=$unit[2];
       push @pos, ($unit[3],$unit[4],$unit[6],$unit[7]);
    }
    my @t=sort {$type{$a} <=> $type{$b}} keys %type;
    my $class;
    if ($t[$#t]=~/SAME/){
       $class="inversion";
    }elsif($t[$#t]=~/LONG/){
       $class="contraction";
    }elsif($t[$#t]=~/SHORT/){
       $class="expansion";
    }
    @pos=sort {$a <=> $b} @pos;
    my $len=$pos[$#pos]-$pos[0];
    my @gene;
    my $refgff=$gff->{$chr};
    foreach my $g (keys %$refgff){
       if ($refgff->{$g}->[0] >= $pos[0] and $refgff->{$g}->[1] <= $pos[$#pos] ){
          push @gene,$g;
       }
    }
    @gene =sort @gene;
    my $n=@gene;
    my $l=join(",",@gene);
    print OUT "$id\t$number\t$chr\t$class\t$pos[0]\t$pos[$#pos]\t$len\t$gene[0]\t$gene[$#gene]\t$n\n";
}
close IN;
close OUT;
$/="\n";
}




sub parseGFFpos
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $id;    ##ID for element
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/mRNA/){
        if ($unit[8]=~/ID=(.*?);/ or $unit[8] =~/ID=(.*)/){
            $id=$1;
            $id=$1 if ($id=~/LOC_(.*)/);
        }
        $hash{$unit[0]}{$id}=[$unit[3],$unit[4]]; 
    }
}
close IN;
return \%hash;
}

