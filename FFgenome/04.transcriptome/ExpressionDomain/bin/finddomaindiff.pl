#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"refbed:s","refrpkm:s","refmeth:s","qryrpkm:s","qrymeth:s","orth:s","ks:s","help");


my $help=<<USAGE;
perl finddomaindiff.pl --refbed ../input/rice.mRNA.bed --refrpkm ../input/normalization/normalization.2.4r --refmeth ../input/rice.gene.CG.level.part --qryrpkm ../input/normalization/normalization.1.4r --qrymeth ../input/FF.gene.CG.level.part --orth ../input/OB-OS.orth > log &
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 
my $refbed=getbed($opt{refbed});
my $refrpkm=getrpkm($opt{refrpkm});
my $qryrpkm=getrpkm($opt{qryrpkm});
my $refmeth=getmeth($opt{refmeth});
my $qrymeth=getmeth($opt{qrymeth});
my $ks=getks($opt{ks});
my $ort=getorth($opt{orth});

foreach my $chr (sort keys %$refbed){
   foreach my $g (sort {$refbed->{$chr}->{$a} <=> $refbed->{$chr}->{$b}} keys %{$refbed->{$chr}}){
       my $ortg=$ort->{$g} ? $ort->{$g} : "NA";
       next if ($ortg eq "NA");
       my $refr=$refrpkm->{$g} ? $refrpkm->{$g} : "NA";
       my $qryr=$qryrpkm->{$ortg} ? $qryrpkm->{$ortg} : "NA";
       my $refm=$refmeth->{$g} ? $refmeth->{$g} : "NA"; 
       my $qrym=$qrymeth->{$ortg} ? $qrymeth->{$ortg} : "NA";
       my $refks=$ks->{$ortg};
       next if ($refr eq "NA" or $qryr eq "NA");
       my $ratio=$qryr-$refr;      
       print "$chr\t$refbed->{$chr}->{$g}\t$g\t$refr\t$refm\t$ortg\t$qryr\t$qrym\t$ratio\t$refks\n";
   }
}



#chr01   7355    7553    OBR_GLEAN_10003592      1       -
sub getbed
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash{$unit[0]}{$unit[3]}=$unit[1];
}
close IN;
return \%hash;
}

#OBR_GLEAN_10000026      0.0374312       0.21342 0.175388
sub getks
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash{$unit[0]}=$unit[2];
}
close IN;
return \%hash;
}


#OBR_GLEAN_10017382      1.815706324284  1.01114267114563
sub getrpkm
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash{$unit[0]}=$unit[2];
}
close IN;
return \%hash;
}

sub getorth
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash{$unit[1]}=$unit[0];
}
close IN;
return \%hash;
}


#Os01t0953801-00 233     27      108     0       164     115
sub getmeth
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/ or $_=~/^ID/);
    my @unit=split("\t",$_);
    my $g=shift @unit;
    $hash{$g}=join(":",@unit);;
}
close IN;
return \%hash;
}

