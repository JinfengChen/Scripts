#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"distance:s","typedir:s","help");


my $help=<<USAGE;
perl $0 --distance --typedir
--typedir: dir of cluster types, shared.table/unique.table/noncluster.table.
--distance: ortholog distance file, OBRACH_TIGR6.distance.txt
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my ($ref1,$ref2,$com)=parsedistance($opt{distance});
`mkdir $com`;
print "Species\tOrtholog\tNonOrtholog\n";
& ort($ref1,$opt{typedir},$com);
& ort($ref2,$opt{typedir},$com);

###################
sub ort
{
my ($ref,$dir,$com)=@_;
my ($n1,$n2);
my @spe=values %$ref;
my $shared="$dir/$spe[0].shared.table";
my $ort="$com/$spe[0].ortholog.table";
my $nonort="$com/$spe[0].nonortholog.table";
open OUT1, ">$ort" or die "$!";
open OUT2, ">$nonort" or die "$!";
open IN, "$shared" or die "$!";
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
print "$spe[0]\t$n1\t$n2\n";
}


sub parsedistance
{
my ($file)=@_;
my ($spe1,$spe2,%hash1,%hash2,$com);
if ($file=~/((\w+)\_(\w+))\.distance\.txt/){
   $spe1=$2;
   $spe2=$3;
   $com=$1;
   #print "$1\t$2\t$3\n";
}
#$com="$spe1\_$spe2";
open IN, "$file" or die "$!";
<IN>;
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash1{$unit[0]}=$spe1;
    $hash2{$unit[1]}=$spe2;
}
close IN;
return (\%hash1,\%hash2,$com)
}
