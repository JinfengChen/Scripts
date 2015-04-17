#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"chromatin:s","meth:s","density:s","project:s","help");


my $help=<<USAGE;
Prepare information of gene length/exon, intron length/intergenic length/gene expression/methylation/flanking TE of gene for euchromatin and heterchromatin for one species. 
perl $0 --chromatin OBa_manual --meth OBa.gene.CG.level.part --density OBa_50000.gene_rpkm_TE.density > log 2> log2 &
--chromatin: project name of chromatin which used by splitChromatin.pl to generate results, (OBa_manual)_chromatin
--meth: methylation information,
   ID      UpC     UpmC    BodyC   BodymC  DownC   DownmC
   OBR_GLEAN_10004303      216     43      76      0       176     153
--density: TE proportion of gene flanking sequence and normalized gene expression level
   geneid	                rpkm    te%
   OBR_GLEAN_10015563      NA      0.39674
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}
my $chromatin=$opt{chromatin}."_chromatin";
my $eudir="$chromatin/euchromatin";
my $hedir="$chromatin/heterochromatin";
`mkdir $eudir` unless (-f $eudir);
`mkdir $hedir` unless (-f $hedir);
my $eugff="$chromatin/$opt{chromatin}.gene.euchromatin.gff";
my $hegff="$chromatin/$opt{chromatin}.gene.heterochromatin.gff";
#######gene structure infor
& genestructure($eugff,$eudir);
& genestructure($hegff,$hedir);

my ($rpkm,$te)=densityinf($opt{density});
my ($up,$body,$down)=methinf($opt{meth});
#######gene expression infor
writeinf($eugff,$rpkm,$eudir,"expression");
writeinf($hegff,$rpkm,$hedir,"expression");
#######te proportion infor
writeinf($eugff,$te,$eudir,"TEdensity");
writeinf($hegff,$te,$hedir,"TEdensity");
#######methylation level
##upstream
writeinf($eugff,$up,$eudir,"Upstream");
writeinf($hegff,$up,$hedir,"Upstream");
##body
writeinf($eugff,$body,$eudir,"Body");
writeinf($hegff,$body,$hedir,"Body");
##downstream
writeinf($eugff,$down,$eudir,"Downstream");
writeinf($hegff,$down,$hedir,"Downstream");

#######################
sub writeinf
{
my ($gff,$hash,$dir,$title)=@_;
my $refid=parseGFFid($gff);
open OUT, ">$dir/$title\.txt" or die "$!";
foreach(keys %$refid){
    #print "$_\t$hash->{$_}\n";
    if (exists $hash->{$_} and $hash->{$_} ne "NA"){
       print OUT "$_\t$hash->{$_}\n";
    }
}
close OUT;
}


sub genestructure
{
my ($gff,$dir)=@_;
#print "$gff\t$dir\n";
my $gffstat="./gff_statistic.pl";
`perl $gffstat $gff $dir`;

}

sub methinf
{
my ($file)=@_;
my %up;
my %body;
my %down;
open IN, "$file" or die "$!";
   <IN>;
   while(<IN>){
      chomp $_;
      my @unit=split("\t");
      $up{$unit[0]}= $unit[1] > 0 ? $unit[2]/$unit[1] : "NA";
      $body{$unit[0]}=$unit[3] > 0 ? $unit[4]/$unit[3] : "NA";
      $down{$unit[0]}=$unit[5] > 0 ? $unit[6]/$unit[5] : "NA";
   }
close IN;
return (\%up,\%body,\%down);
}


sub densityinf
{
my ($file)=@_;
my %rpkm;
my %te;
open IN, "$file" or die "$!";
   while(<IN>){
      chomp $_;
      my @unit=split("\t");
      $rpkm{$unit[0]}=$unit[1];
      $te{$unit[0]}=$unit[2];
   }
close IN;
return (\%rpkm,\%te);
}


sub parseGFFid
{
my ($gff)=@_;
#print "$gff\n";
my %hash;  ##hash to store every record by key of Seq_id
my $id;    ##ID for element
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/mRNA/){
        if ($unit[8]=~/ID=(.*);/ or $unit[8] =~/ID=(.*)/){
            $id=$1;
        }
        $hash{$id}=1;
    }

}
close IN;
return \%hash;
}

