#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);


GetOptions (\%opt,"gff:s","help");


my $help=<<USAGE;
perl /rhome/cjinfeng/software/bin/gff2bed.pl --gff MSU_r7.fa.RepeatMasker.out.gff > MSU_r7.fa.RepeatMasker.out.bed

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my @chr=("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12");
my $refgff=parseGFF($opt{gff});
for(my $i=0;$i<@chr;$i++){
   my @te=@{$refgff->{$chr[$i]}};
   my @sort=sort{$a->[0] <=> $b->[0]} @te;
   for(my $j=0;$j<@sort;$j++){
      #print "$chr[$i]\t$sort[$j]->[0]\t$sort[$j]->[1]\t$sort[$j]->[2]\t$sort[$j]->[3]\t$sort[$j]->[4]\n";
      print "$sort[$j]->[5]\t$chr[$i]:$sort[$j]->[0]..$sort[$j]->[1]\n";
   }
}

sub readtable
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash{$unit[0]}=1;
}
close IN;
return \%hash;
}


sub parseGFF
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $seq;
my $id;    ##ID for element
my $record;##all line for this record
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/Trans/){
        $seq=$unit[0];
        if ($unit[8]=~/ID=(.*?);/ or $unit[8] =~/ID=(.*)/){
            $id=$1;
            my $name = $1 if ($unit[8]=~/.*;Target=(\S+)\s+/);
            push @{$hash{$unit[0]}},[$unit[3],$unit[4],$id,".",$unit[6],$name];
            #print "$unit[0]\t$unit[3]\t$unit[4]\t$id\t.\t$unit[6]\n";
        }
    }

}
close IN;

return \%hash;
}
 
