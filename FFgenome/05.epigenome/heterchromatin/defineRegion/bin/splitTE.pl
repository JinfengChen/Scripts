#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"gff:s","help");


my $help=<<USAGE;
Split TE gff file into gff of different families.
perl $0 -gff ../input/OBa.TE.gff
--gff: the input gff should be directory, which contein gff of each chromosomes, named chr01 ......
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my @gff=glob("$opt{gff}/chr*");
my $type=$opt{gff}.".type";
`mkdir $type` unless (-e $type); 

foreach my $file (@gff){
    my $chr;
    if ($file=~/\/(chr\d+)$/){
       $chr=$1;
    }else{
       next;
    }
    my $prefix=$type."/"."$chr";
    splitTE($file,$prefix);
}

########
sub splitTE
{
my ($gff,$prefix)=@_;
open IN, "$gff" or die "$!";
open DNA, ">$prefix.DNA.gff" or die "$!";
open RT, ">$prefix.RT.gff" or die "$!";
open CACTA, ">$prefix.CACTA.gff" or die "$!";
open MITE, ">$prefix.MITE.gff" or die "$!";
open MUDR, ">$prefix.MUDR.gff" or die "$!";
open GYPSY, ">$prefix.GYPSY.gff" or die "$!";
open COPIA, ">$prefix.COPIA.gff" or die "$!";
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
    }elsif($TE_class=~/MuDR/i){
       print MUDR "$line";
    }elsif($TE_class=~/Stowaway/i or $TE_class=~/Tourist/i or $TE_class=~/MITE/i){
       print MITE "$line";
    }elsif($TE_class=~/Gypsy/i){
       print GYPSY "$line";
    }elsif($TE_class=~/Copia/i){
       print COPIA "$line";
    }
}
close GYPSY;
close COPIS;
close MUDR;
close MITE;
close CACTA;
close DNA;
close RT;
close IN;
}
 
