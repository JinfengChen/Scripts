#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"chrlen:s","tegff:s","genegff:s","segment:s","project:s","help");


my $help=<<USAGE;
Split Gene/TE gff into heterchromatin and euchromatin.
perl $0 --chrlen ../input/OBa.chrlen --tegff ../input/OBa.all.fa.RepeatMasker.out.gff.chr --genegff ../input/OBa.all.gff --segment ../input/OBa_methylation_wave --project OBa > log 2> log2 &
--tegff: chromosome file of te gff 
--genegff: chromsome file of gene gff
--segment: chromsome file of segment file for heterchromatin
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my $chromatin=$opt{project}."_chromatin";
`mkdir $chromatin` unless (-f $chromatin);

####summary statistics for every chromosome
my $refchr=chrlen("$opt{chrlen}");
my %helen; ## heterchromatin length for each chromosome
my %gene; ## accumulated gene number and length for each chromosome
my %hegene; ## accumulated heterchromatin gene number and length for each chromosome
my %te; ## accumulated TE number for each chromosome
my %hete; ## accumulated heterchromatin TE number for each chromosome


my @file=glob("$opt{segment}/*.wig.segment.txt");
open GHE, ">$chromatin/$opt{project}.gene.heterochromatin.gff" or die "$!";
open GEU, ">$chromatin/$opt{project}.gene.euchromatin.gff" or die "$!";
open THE, ">$chromatin/$opt{project}.te.heterochromatin.gff" or die "$!";
open TEU, ">$chromatin/$opt{project}.te.euchromatin.gff" or die "$!";
foreach my $file (@file){
    if ($file=~/(chr\d+).*.wig.segment.txt/){
        my $chr=$1;
        my $eunum;
        my $henum;
        my $refgenegff=parseGFFgene("$opt{genegff}/$chr",\%gene);
        my $reftegff=parseGFFte("$opt{tegff}/$chr",\%te);
        open IN, "$file" or die "$!";
        while(<IN>){
           chomp $_;
           next if ($_=~/^$/);
           my @unit=split("\t",$_);
           if ($unit[3] == 1){
               $henum++;
               my $helen=$unit[2]-$unit[1]+1;
               $helen{$chr}+=$helen;
               my $ref=$chr."H".$henum;
               foreach(sort {$a <=> $b} keys %$refgenegff){
                   if ($_ >= $unit[1] and $_ <= $unit[2]){
                        my $new=changerefgene($refgenegff->{$_},$ref,\%hegene);
                        print GHE "$new\n";
                        delete $refgenegff->{$_};
                   }else{
                        last;
                   }
               }
               foreach(sort {$a <=> $b} keys %$reftegff){
                   if ($_ >= $unit[1] and $_ <= $unit[2]){
                        my $new=changerefte($reftegff->{$_},$ref,\%hete);
                        print THE "$new\n";
                        delete $reftegff->{$_};
                   }else{
                        last;
                   }
               }
           }else{
               $eunum++;
               my $ref=$chr."E".$eunum;
               foreach(sort {$a <=> $b} keys %$refgenegff){
                   if ($_ >= $unit[1] and $_ <= $unit[2]){
                        my $new=changerefgene($refgenegff->{$_},$ref);
                        print GEU "$new\n";
                        delete $refgenegff->{$_};
                   }else{
                        last;
                   }
               }
               foreach(sort {$a <=> $b} keys %$reftegff){
                   if ($_ >= $unit[1] and $_ <= $unit[2]){
                        my $new=changerefte($reftegff->{$_},$ref);
                        print TEU "$new\n";
                        delete $reftegff->{$_};
                   }else{
                        last;
                   }
               }
           }
        }
        close IN;
    }
}
close GHE;
close GEU;
close THE;
close TEU;
#########################################
open SUM, ">$chromatin/$opt{project}.chromatin.summary" or die "$!";
print SUM "Chr\tLength\tHeterLength\t#Gene\tCodingLength\t#HetroGene\tHeterCodingLength\t#TE\tTELength\t#HeterTE\tHeterTELength\n";
foreach(sort {$a cmp $b} keys %$refchr){
    my $chr=$_;
    print SUM "$chr\t$refchr->{$chr}\t$helen{$chr}\t$gene{$chr}->[0]\t$gene{$chr}->[1]\t$hegene{$chr}->[0]\t$hegene{$chr}->[1]\t$te{$chr}->[0]\t$te{$chr}->[1]\t$hete{$chr}->[0]\t$hete{$chr}->[1]\n";
}
close SUM;

#########################################

sub changerefgene
{
my ($gffline,$ref,$refhash)=@_;
my @unit=split("\n",$gffline);
my @new;
foreach my $line (@unit){
   if ($line =~/^chr/){
       my @array=split("\t",$line);
       if ($array[2] =~/mRNA/i){
          $refhash->{$array[0]}->[0]++;
       }else{
          my $len=$array[4]-$array[3]+1;
          $refhash->{$array[0]}->[1]+=$len;   
       }
       $array[0]=$ref;
       my $newline=join("\t",@array);
       push (@new,$newline);
   } 
}
my $newgff=join("\n",@new);
return $newgff;
}

sub changerefte
{
my ($gffline,$ref,$refhash)=@_;
my @unit=split("\n",$gffline);
my @new;
foreach my $line (@unit){
   if ($line =~/^chr/){
       my @array=split("\t",$line);
       $refhash->{$array[0]}->[0]++;
       my $len=$array[4]-$array[3]+1;
       $refhash->{$array[0]}->[1]+=$len;
       $array[0]=$ref;
       my $newline=join("\t",@array);
       push (@new,$newline);
   }
}
my $newgff=join("\n",@new);
return $newgff;
}



sub parseGFFgene
{
my ($gff,$gene)=@_;
#print "$gff\n";
my %hash;  ##hash to store every record by key of Seq_id
my $seq;
my $start;
my $id;    ##ID for element
my $record;##all line for this record
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/mRNA/){
        $seq=$unit[0];
        $start=$unit[3];
        $gene->{$seq}->[0]++;
        if ($unit[8]=~/ID=(.*);/ or $unit[8] =~/ID=(.*)/){
            $id=$1;
        }
        $record="$_\n";
        $hash{$start}=$record;
    }elsif($unit[0] eq $seq and $unit[8] =~ /Parent=$id/){
        $hash{$start}.="$_\n";
        my $len=$unit[4]-$unit[3]+1;
        $gene->{$seq}->[1]+=$len;
    }

}
close IN;
return \%hash;
}
 

sub parseGFFte
{
my ($gff,$te)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $seq;
my $start;
my $id;    ##ID for element
my $record;##all line for this record
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    $seq=$unit[0];
    $start=$unit[3];
    $te->{$seq}->[0]++;
    my $len=$unit[4]-$unit[3]+1;
    $te->{$seq}->[1]+=$len;
    $hash{$start}="$_\n";
}
close IN;
return \%hash;
}


### get chr length hash
sub chrlen
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
   chomp $_;
   my @unit=split("\t",$_);
   $hash{$unit[0]}=$unit[1];
}
close IN;
return \%hash;
}


