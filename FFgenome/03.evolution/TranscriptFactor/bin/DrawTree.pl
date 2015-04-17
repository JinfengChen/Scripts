#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"dna:s","protein:s","align","phyml","nj","help");


my $help=<<USAGE;
Construct tree using Phyml and draw a simply figure using tree_plot.pl
perl $0 --protien protein.fa --align --nj > log 2> log2 &
--protein: protein sequence fasta
--align: do alignment in this script
--phyml: draw ml draw, slow
--nj: draw nj tree, fast
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 
my $bin="/home/jfchen/FFproject/FFgenome/03.evolution/Family_analysis/bin";
my $title=$1 if ($opt{protein}=~/(.*)\.fa/ or $opt{protein}=~/(.*)\.phy/);

if ($opt{dna}){
    `phyml --quiet -i $opt{dna} -d nt > log 2> log2`;  
    `perl /home/jfchen/software/tree/tree_plot.pl --png $opt{dna}_phyml_tree.txt`;
}
if ($opt{protein}){
   if ($opt{align}){
      `muscle -maxiters 2 -in $opt{protein} -out $opt{protein}.muscle > $title.log 2> $title.log2` if $opt{align};
      `perl $bin/fasta2phylip.pl $opt{protein}.muscle` if $opt{align};
      `phyml --quiet -i $opt{protein}.muscle.phy -d aa > $title.log 2> $title.log2` if $opt{phyml};
      `perl /home/jfchen/software/tree/tree_plot.pl --png $opt{protein}.muscle.phy_phyml_tree.txt` if $opt{phyml};
      `treebest nj $opt{protein}.muscle > $opt{protein}.muscle.nj.tree` if $opt{nj};
      `perl /home/jfchen/software/tree/tree_plot.pl --png $opt{protein}.muscle.nj.tree` if $opt{nj};
   }else{
      `phyml --quiet -i $opt{protein} -d aa > $title.log 2> $title.log2` if $opt{phyml};
      `perl /home/jfchen/software/tree/tree_plot.pl --png $opt{protein}_phyml_tree.txt` if $opt{phyml};
      `treebest nj $opt{protein} > $opt{protein}.nj.tree 2> $title.log2` if $opt{nj};
      `perl /home/jfchen/software/tree/tree_plot.pl --png $opt{protein}.nj.tree` if $opt{nj};
   }
}

