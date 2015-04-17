#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"gff:s","project:s","help");


my $help=<<USAGE;

Convert repeat gff format file to bed format.
Example:
Scaffold000021  4057231 4057449 repeat    TE00001      +

Run: perl $0 -gff OBa.repeat.gff -project FF_repeat_bed
-gff: GFF file 
-project: prefix for output file

USAGE

if (keys %opt < 1){
   print "$help\n";
   exit();
}

`mkdir $opt{project}` unless (-f $opt{project});
 
my $chr=splitgff($opt{gff},"./$opt{project}");
foreach(keys %$chr){
    print "$_\n";
    my $file="./$opt{project}/$_";
    my $pre="./$opt{project}/$_";
    & splitTE($file,$pre);
}

##########################################
sub splitgff {
my ($gff,$dir)=@_;
my %chr;
chomp $gff;
`mkdir $dir` unless (-d $dir);
open IN, "$gff" or die "$!";
while(<IN>){
    next if ($_ eq "");
    next if ($_ =~/^##/);
    my @unit=split("\t", $_);
    $chrgff=$unit[0];
    unless (exists $chr{$unit[0]}){
       $chr{$unit[0]}=1;
    }
    my $line=join("\t",@unit);
    open OUT, ">>$dir/$chrgff"; 
         print OUT "$line"; 
    close OUT;
}
close IN;
return \%chr;
}


######
sub splitTE
{
my ($gff,$prefix)=@_;
open IN, "$gff" or die "$!";
open ALL, ">$prefix.REPEAT.bed" or die "$!";
open DNA, ">$prefix.DNA.bed" or die "$!";
open RT, ">$prefix.RT.bed" or die "$!";
open CACTA, ">$prefix.CACTA.bed" or die "$!";
open MITE, ">$prefix.MITE.bed" or die "$!";
open MUDR, ">$prefix.MUDR.bed" or die "$!";
open GYPSY, ">$prefix.GYPSY.bed" or die "$!";
open COPIA, ">$prefix.COPIA.bed" or die "$!";
while(<IN>){
    next if ($_ eq "");
    next if ($_ =~/^#/);
    my @temp=split("\t",$_);
    my $TE_class= $1 if ($temp[8] =~ /Class=([^;]+);*/);
    my $TE_ID=$1 if ($temp[8]=~ /ID=(\w+);/);
    my $line="$temp[0]\t$temp[3]\t$temp[4]\trepeat\t$TE_ID\t$temp[6]\n";
    #print "$TE_class\n";
    print ALL "$line";
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
close ALL;
close IN;
}


