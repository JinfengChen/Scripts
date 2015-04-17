#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"chrlen:s","tegff:s","genegff:s","segment:s","project:s","help");


my $help=<<USAGE;
perl $0 --chrlen ../input/OBa.chrlen --tegff ../input/OBa.all.fa.RepeatMasker.out.gff.chr --genegff ../input/OBa.all.gff --segment ../input/OBa_methylation_wave --project OBa > log 2> log2 &
--tegff: chromosome file of te gff 
--genegff: chromsome file of gene gff
--segment: chromsome file of segment file for heterchromatin
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my @file=glob("$opt{segment}/*.wig.segment.txt");
open GHE, ">$opt{project}.gene.heterochromatin.gff" or die "$!";
open GEU, ">$opt{project}.gene.euchromatin.gff" or die "$!";
open THE, ">$opt{project}.te.heterochromatin.gff" or die "$!";
open TEU, ">$opt{project}.te.euchromatin.gff" or die "$!";
foreach my $file (@file){
    if ($file=~/(chr\d+).*.wig.segment.txt/){
        my $chr=$1;
        my $eunum;
        my $henum;
        my $refgenegff=parseGFFgene("$opt{genegff}/$chr");
        my $reftegff=parseGFFte("$opt{tegff}/$chr");
        open IN, "$file" or die "$!";
        while(<IN>){
           chomp $_;
           next if ($_=~/^$/);
           my @unit=split("\t",$_);
           if ($unit[3] == 1){
               $henum++;
               my $ref=$chr."H".$henum;
               foreach(sort {$a <=> $b} keys %$refgenegff){
                   if ($_ >= $unit[1] and $_ <= $unit[2]){
                        my $new=changeref($refgenegff->{$_},$ref);
                        print GHE "$new\n";
                        delete $refgenegff->{$_};
                   }else{
                        last;
                   }
               }
               foreach(sort {$a <=> $b} keys %$reftegff){
                   if ($_ >= $unit[1] and $_ <= $unit[2]){
                        my $new=changeref($reftegff->{$_},$ref);
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
                        my $new=changeref($refgenegff->{$_},$ref);
                        print GEU "$new\n";
                        delete $refgenegff->{$_};
                   }else{
                        last;
                   }
               }
               foreach(sort {$a <=> $b} keys %$reftegff){
                   if ($_ >= $unit[1] and $_ <= $unit[2]){
                        my $new=changeref($reftegff->{$_},$ref);
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


sub changeref
{
my ($gffline,$ref)=@_;
my @unit=split("\n",$gffline);
my @new;
foreach my $line (@unit){
   if ($line =~/^chr/){
       my @array=split("\t",$line);
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
my ($gff)=@_;
print "$gff\n";
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
        if ($unit[8]=~/ID=(.*);/ or $unit[8] =~/ID=(.*)/){
            $id=$1;
        }
        $record="$_\n";
        $hash{$start}=$record;
    }elsif($unit[0] eq $seq and $unit[8] =~ /Parent=$id/){
        $hash{$start}.="$_\n";
    }

}
close IN;
return \%hash;
}
 

sub parseGFFte
{
my ($gff)=@_;
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


