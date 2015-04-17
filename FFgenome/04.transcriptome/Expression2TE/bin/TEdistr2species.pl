#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"bed1:s","gff1:s","bed2:s","gff2:s","ortholog:s","windows:s","project:s","help");


my $help=<<USAGE;
Distribution of gene number along different levels of TE% in two species. 
-bed1: BED format gene annotation file, generate by GFF2BED.pl from GFF3 format gene annotation.
-gff1: GFF format TE annotation file, should be contain only TE class that interested to calculate correlation.
-bed2: BED format gene annotation file, generate by GFF2BED.pl from GFF3 format gene annotation.
-gff2: GFF format TE annotation file, should be contain only TE class that interested to calculate correlation.
-ortholog: ortholog pair between 2 species.
-windows: windows size upstream/downstream of gene, in which we calculate the fraction of TE.
-project: project name.

Run: perl $0 -ortholog OB-OS.orth -bed1 FF.mRNA.bed -gff1 ../input/chr01/FF.COPIA.gff -bed2 Rice.mRNA.bed -gff2 ../input/chr01/Rice.COPIA.gff -windows 50000 -project A01_Gene_COPIA_2
 
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}
my $BEDtools="/home/jfchen/FFproject/tools/BEDTools/bin";
system("$BEDtools/windowBed -a $opt{bed1} -b $opt{gff1} -w $opt{windows} > $opt{project}.1.windowsBED");
system("$BEDtools/windowBed -a $opt{bed2} -b $opt{gff2} -w $opt{windows} > $opt{project}.2.windowsBED");
my $refBED1=windowsBED("$opt{project}.1.windowsBED");
my $refBED2=windowsBED("$opt{project}.2.windowsBED");
`rm *.windowsBED`;
my $reforth=orth($opt{ortholog});

open OUT, ">$opt{project}.4r" or die "$!";
foreach(keys %$reforth){
   #print "$_\t$reforth->{$_}\n";
   if (exists $refBED1->{$_} and exists $refBED2->{$reforth->{$_}}){  
      my $freq1=$refBED1->{$_}/(2*$opt{windows});
      my $freq2=$refBED2->{$reforth->{$_}}/(2*$opt{windows});
      $freq1=$freq1 > 1? 1 : $freq1;
      $freq2=$freq2 > 1? 1 : $freq2;   
      print OUT "$_\t$freq1\t$reforth->{$_}\t$freq2\n";
   }
}
close OUT;



open OUT, ">$opt{project}.r" or die "$!"; 
print OUT <<"END.";
read.table("$opt{project}.4r") -> x
pdf("$opt{project}.pdf")
hist(x[,2],breaks=20,col=2,xlim=c(0,1),xlab="TE%",ylab="Frequency",main="FF")
hist(x[,4],breaks=20,col=2,xlim=c(0,1),xlab="TE%",ylab="Frequency",main="Rice")
plot(x[,2]~x[,4],xlab="TE% Rice",ylab="TE% FF");
dev.off()
END.
close OUT;


system ("cat $opt{project}.r | R --vanilla --slave");
################################################
sub windowsBED
{
#### sum len of feature that fill in windows of gene
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


sub orth
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
     chomp $_;
     next if ($_ eq "");
     my @unit=split("\t",$_);
     unless (exists $hash{$unit[0]}){
        $hash{$unit[0]}=$unit[1];
     }
}
close IN;
return \%hash;
}


sub rpkmexpr
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split(" ",$_);
    if ($unit[3]=~/\-/){
       $hash{$unit[0]}="NA";
    }else{
       $hash{$unit[0]}=$unit[3];
    }
}
close IN;
return \%hash;
}




 
