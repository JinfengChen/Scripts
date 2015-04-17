#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"table:s","gff1:s","gff2:s","help");


my $help=<<USAGE;

Run: perl connect.pl -table ob_os.blast -gff1 ff.gene.gff -gff2 rice.gene.gff
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my $refgff1=parseGFF($opt{gff1});
my $refgff2=parseGFF($opt{gff2});


my (@start1,@start2,@end1,@end2);
my ($chr1,$chr2);
open IN, "$opt{table}" or die "$!";
open OUT1, ">all.gene.position" or die "$!";
open OUT2, ">connect.txt" or die "$!";
while(<IN>){
   chomp $_;
   next if ($_ eq "");
   my @unit=split("\t",$_);
   #print "$unit[0]\t$unit[1]\n";
   if (exists $refgff1->{$unit[0]} and exists $refgff2->{$unit[1]}){
      print "$unit[0]\t$unit[1]\n";
      my $tempref1=$refgff1->{$unit[0]};
      my $tempref2=$refgff2->{$unit[1]};
      $chr1=$tempref1->[0][0];
      $chr2=$tempref2->[0][0];
      push (@start1,$tempref1->[0][1]);
      push (@end1,$tempref1->[0][3]);
      push (@start2,$tempref2->[0][1]);
      push (@end2,$tempref2->[0][3]);
      #print OUT2 "$tempref1->[0][0]\t$tempref1->[0][1]\t$tempref1->[0][3]\t$tempref2->[0][0]\t$tempref2->[0][1]\t$tempref2->[0][3]\n";
      my $temp="$tempref1->[0][0]"."__"."$tempref2->[0][0]"."_".$unit[0]."_".$unit[1];
      print OUT1 "$temp\t$tempref1->[0][1]\t$tempref1->[0][2]\t$tempref1->[0][3]\t$tempref2->[0][1]\t$tempref2->[0][2]\t$tempref2->[0][3]\n";
   }
}
@start1 = sort {$a <=> $b} @start1;
@start2 = sort {$a <=> $b} @start2;
@end1 = sort {$a <=> $b} @end1;
@end2 = sort {$a <=> $b} @end2;
my $start1=shift @start1;
my $start2=shift @start2;
my $end1=pop @end1;
my $end2=pop @end2;
print OUT2 "$chr1\t$start1\t$end1\t$chr2\t$start2\t$end2\n";

close OUT1;
close OUT2;
close IN;
###########################################

sub parseGFF
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $seq;
my $id;    ##ID for element
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
        my $mid=$unit[3]+($unit[4]-$unit[3])/2;
        my @record;
        push (@record,[$seq,$unit[3],$mid,$unit[4]]);
        $hash{$id}=\@record;
    }

}
close IN;
return \%hash;
}

