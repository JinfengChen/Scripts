#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"iprscan:s","project:s","help");


my $help=<<USAGE;
Read in all *.iprscan file in input.
Generate a txt file
Gene id for each Pfam will be write into ../input/Pfam/.  

Example output file:
Pfam    OB      OS      Total   Annotation
PF00004 103     106     209     ATPase, AAA-type, core
PF00005 119     110     229     ABC transporter-like
PF00006 8       6       14      ATPase, F1/V1/A1 complex, alpha/beta subunit, nucleotide-binding domain
PF00008 0       1       1       EGF-like


Run: perl $0 --iprscan TIGR6 --project ITGR6
--iprscan: dir of iprscan to parse
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my %pfam2anno;
my %pfam2spec;
my @iprscan=glob("$opt{iprscan}/*.iprscan");
foreach(@iprscan){
   my $species;
   if ($_=~/$opt{iprscan}\/(.*)\.iprscan/){
      $species=$1;
      open IN, "$_" or die "$!";
         while(<IN>){
            chomp $_;
            my @unit=split("\t",$_);
            next unless ($unit[3] eq "HMMPfam");
            unless (exists $pfam2anno{$unit[4]}){
                $pfam2anno{$unit[4]}=$unit[12];
            }  
            if (exists $pfam2spec{$species}){
                my $refgene=$pfam2spec{$species};
                if (exists $refgene->{$unit[4]}){
                    my $refhash=$refgene->{$unit[4]};
                    $refhash->{$unit[0]}=1; ### make sure to kick out repeat annotation for one gene
                }else{
                    my %temp;
                    $temp{$unit[0]}=1;
                    $refgene->{$unit[4]}=\%temp;
                }
            }else{
                my %pfam2gene;
                my %temp;
                $temp{$unit[0]}=1;
                $pfam2gene{$unit[4]}=\%temp;
                $pfam2spec{$species}=\%pfam2gene;
            }
         }
      close IN;
   }
}

my @spec=sort keys %pfam2spec;
my $spec=join("\t",@spec);

open OUT, ">$opt{project}.Pfam.summary" or die "$!";
print OUT "Pfam\t$spec\tTotal\tAnnotation\n";
foreach(sort keys %pfam2anno){
     my $pfam=$_;
     my $total;
     print OUT "$pfam\t"; 
     foreach(sort keys %pfam2spec){
         my $refgene =$pfam2spec{$_};
         my $refhash=$refgene->{$pfam};
         my $number=keys %$refhash;
=pod
         if ($number >= 1){
            my $outdir="../input/Pfam";
            mkdir $outdir unless (-e $outdir);
            my $pfamfile=$outdir."/".$pfam."_".$_."_".$number;
            open PFAM, ">$pfamfile" or die "$!";
                 foreach(sort keys %$refhash){
                    print PFAM "$_\n";
                 }
            close PFAM;
         }
=cut
         $total+=$number;
         print OUT "$number\t";
     }
     print OUT "$total\t$pfam2anno{$pfam}\n";
}
close OUT;


