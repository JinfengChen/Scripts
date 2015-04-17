#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);


GetOptions (\%opt,"vcf:s", "gff:s","help");


my $help=<<USAGE;
perl orderchr.pl --vcf ALL.gatk.snp.raw.vcf --gff Saccharomyces.fa.out.gff
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my @chrs = readvcf($opt{vcf});

#chrI    RepeatMasker    Transposon      2       62      250     +       .       ID=TE0001
my %hash;
open IN, "$opt{gff}" or die "$!";
while(<IN>){
    chomp $_;
    next unless ($_=~/^chr/);
    my @unit=split("\t", $_);
    my $id = $1 if ($unit[8]=~/ID=(\w+?);/);
    #print "$unit[0]\t$unit[3]\t$unit[4]\t$id\t$unit[7]\t$unit[6]\n";
    push (@{$hash{$unit[0]}}, "$unit[0]\t$unit[3]\t$unit[4]\t$id\t$unit[7]\t$unit[6]");
}
close IN;

my $prefix = basename($opt{gff}, '.gff');
open OUT, ">$prefix.bed" or die "$!";
for(my $i=0; $i<@chrs; $i++){
    #print "$chrs[$i]\t$prefix\n";
    if (exists $hash{$chrs[$i]}){
       my $lines = join("\n", @{$hash{$chrs[$i]}});
       print OUT "$lines\n";
    }
}
close OUT;



###
##contig=<ID=chrI,length=230218>
##contig=<ID=chrII,length=813184>
sub readvcf
{
my ($file)=@_;
my @chrs;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next unless ($_=~/^#/);
    my $line = $_;
    if ($line=~/^##contig=<ID=(chr.*?),length.*/){
        push (@chrs, $1);
    }
}
close IN;
return @chrs;
}
 
