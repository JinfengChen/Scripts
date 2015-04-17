#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"gff:s","help");


my $help=<<USAGE;
Split TE gff file into gff of different families.
perl $0 -gff ../input/OBa.TE.gff
--gff: the input gff is a file of gff for genome
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my $type=$opt{gff}.".type";
`mkdir $type` unless (-e $type); 

my $prefix=$1 if ($opt{gff}=~/.*\/(.+?)\.gff/);
$prefix=$type."/".$prefix;
print "$prefix\n";

& splitTE($opt{gff},$prefix);

########
sub splitTE
{
my ($gff,$prefix)=@_;
open IN, "$gff" or die "$!";
open DNA, ">$prefix.DNA.gff" or die "$!";
open RT, ">$prefix.RT.gff" or die "$!";
open CACTA, ">$prefix.CACTA.gff" or die "$!";
open HELI, ">$prefix.HELITRON.gff" or die "$!";
open HAT, ">$prefix.HAT.gff" or die "$!";
open MITE, ">$prefix.MITE.gff" or die "$!";
open MUDR, ">$prefix.MUDR.gff" or die "$!";
open GYPSY, ">$prefix.GYPSY.gff" or die "$!";
open COPIA, ">$prefix.COPIA.gff" or die "$!";
open SINE, ">$prefix.SINE.gff" or die "$!";
open LINE, ">$prefix.LINE.gff" or die "$!";
while(<IN>){
    next if ($_ eq "");
    next if ($_ =~/^#/);
    my @temp=split("\t",$_);
    my $line=join("\t",@temp);
    my $TE_class= $1 if ($temp[8] =~ /Class=([^;]+);*/);
    #print "$TE_class\n";
    if ($TE_class =~/^DNA/){
       print DNA "$line";
    }elsif($TE_class=~/^LTR/){
       print RT "$line";
    }
    if ($TE_class=~/En-Spm/i){
       print CACTA "$line";
    }elsif($TE_class=~/hAT/i){
       print HAT "$line";
    }elsif($TE_class=~/Helitron/i){
       print HELI "$line";
    }elsif($TE_class=~/MuDR/i){
       print MUDR "$line";
    }elsif($TE_class=~/Stowaway/i or $TE_class=~/Tourist/i or $TE_class=~/MITE/i){
       print MITE "$line";
    }elsif($TE_class=~/Gypsy/i){
       print GYPSY "$line";
    }elsif($TE_class=~/Copia/i){
       print COPIA "$line";
    }elsif($TE_class=~/SINE/i){
       print SINE "$line";
    }elsif($TE_class=~/LINE/i){
       print LINE "$line";
    }
}
close SINE;
close LINE;
close GYPSY;
close COPIS;
close MUDR;
close MITE;
close HAT;
close HELI;
close CACTA;
close DNA;
close RT;
close IN;
}
 
