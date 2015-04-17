#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"help");


my $help=<<USAGE;
Run the program in the directory of fasta cds sequence (*.fa). It will tranlate the cds into aa and align aa using muscle. Transform the aa align into cds align. Run Kaks_calculator.

Run: perl pipe_kaks.pl ./ > KaKs.summary
USAGE


if ($opt{help} or @ARGV < 1){
    print "$help\n";
    exit();
} 

my $dir=$ARGV[0];

my $cds2aa="/home/jfchen/FFproject/FFgenome/03.evolution/adaptive_selection_analysis/bin/cds2aa.pl";
my $aa2cds="/home/jfchen/FFproject/FFgenome/03.evolution/adaptive_selection_analysis/bin/pepMfa_to_cdsMfa.pl";
my $muscle="muscle";
my $calculator="/home/jfchen/FFproject/FFgenome/03.evolution/kaks_pairwise/bin/script_KaKs_calculator/calculate_kaks.pl";
my $calidentity="/home/jfchen/FFproject/FFgenome/03.evolution/kaks_pairwise/bin/script_KaKs_calculator/calculate_cds_aa_identity.pl";
my $sumkaks="/home/jfchen/FFproject/FFgenome/03.evolution/kaks_pairwise/bin/script_KaKs_calculator/sumKaKs.pl";
my @fasta=glob("$dir/*.fa");
print "Gene1	Gene2	Ka	Ks	Ka/Ks	Myr	CDS identity	Protein identity\n";
foreach(@fasta){
   if ($_=~/(.*)\.fa/){
      my $cdsfa=$_;
      my $head=$1;
      my $pepfa=$head."aa".".fa";
      my $pepmuscle=$pepfa.".muscle";
      my $cdsalign=$head.".fas";
      my $axtalign=$cdsalign.".axt";
      my $kaks=$axtalign.".kaks";
      my $identity=$axtalign.".identity";
      #print "$cdsfa\n$pepfa\n$pepmuscle\n$cdsalign\n$axtalign\n$kaks\n$identity\n";
      system ("perl $cds2aa $cdsfa > $pepfa");
      system ("$muscle -in $pepfa -out $pepmuscle -stable > muscle.log 2> muscle.log2");
      system ("perl $aa2cds $pepmuscle $cdsfa > $cdsalign");
      system ("perl $calculator --outdir $dir $cdsalign");
      system ("perl $calidentity $axtalign > $identity");
      system ("perl $sumkaks $kaks $identity");
      #system ("rm $cdsfa $pepfa $pepmuscle $cdsalign $axtalign $kaks $identity");
   }
}

