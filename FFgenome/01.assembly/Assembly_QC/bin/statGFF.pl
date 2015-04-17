#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"fasta:s","genegff:s","tegff:s","help");


my $help=<<USAGE;
perl $0 --fasta --genegff --tegff

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my $length=getfastalen($opt{fasta});

my @gene=glob("$opt{genegff}/*.gff");
my @te  =glob("$opt{tegff}/*.gff");


#####count gene
my %gene;
foreach my $gff1 (@gene){
    my $bac=$1 if $gff1=~/(\w+).gff$/;
    my $genenumber=`grep mRNA -c $gff1`;
    chomp $genenumber;
    $gene{$bac}=$genenumber > 0 ? $genenumber : 0;
}
#####repeat 
my %te;
foreach my $gff2 (@te){
    my $bac=$1 if $gff2=~/(\w+).gff$/;
    `perl /home/jfchen/FFproject/tools/bin/stat_TE.pl --rank subtype -gff $gff2 > temp.te`;
    my $refhash=statTE("temp.te");
    $te{$bac}=[$refhash->{"DNA"},$refhash->{"MULE"},$refhash->{"RNA"},$refhash->{"LTR"},$refhash->{"TE"}];
}

foreach my $bac (sort keys %$length){
    #print "$bac\t$te{$bac}[0]\n";
    my $genedensity=$gene{$bac} > 0 ? $length->{$bac}/$gene{$bac} : "NA";
    my $dna=$te{$bac}[0]/$length->{$bac};
    my $mule=$te{$bac}[1]/$length->{$bac};
    my $rna=$te{$bac}[2]/$length->{$bac};
    my $ltr=$te{$bac}[3]/$length->{$bac};
    my $te =$te{$bac}[4]/$length->{$bac};
    print "$bac\t$length->{$bac}\t$gene{$bac}\t$genedensity\t$te\t$dna\t$mule\t$rna\t$ltr\n";
}

##########################
sub getfastalen
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
    $hash{$head}=length $seq;
}
close IN;
$/="\n";
return \%hash;
}


sub statTE
{
my ($te)=@_;
my %hash;
my ($dna,$rna,$mule,$ltr);
open IN ,"$te" or die "$!";
while(<IN>){
   chomp $_;
   my @unit=split(" ",$_);
   $hash{"TE"}+=$unit[1];
   if ($unit[0]=~/DNA/){
      $hash{"DNA"}+=$unit[1];
   }elsif($unit[0]=~/LTR/ or $unit[0]=~/SINE/ or $unit[0]=~/LINE/){
      $hash{"RNA"}+=$unit[1];
   }
   if ($unit[0]=~/MULE/){
      $hash{"MULE"}+=$unit[1];
   }
   if ($unit[0]=~/LTR/){
      $hash{"LTR"}+=$unit[1];
   }
}
close IN;
return \%hash;
}
