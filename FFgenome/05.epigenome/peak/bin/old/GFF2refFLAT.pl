#!/usr/bin/perl
use Getopt::Long;

my %opt;
GetOptions(\%opt,"gff:s","refFLAT:s","help");

if ($opt{help} or keys %opt < 1){
   print "Usage: perl $0 -g all.gff -r refFlat.txt\n";
   exit ();
}

my @record;
my $refgff=parseGFF($opt{gff});
open OUT, ">$opt{refFLAT}" or die "$!";
foreach(keys %$refgff){
    print "$_\n";
    my $chr;
    my $strand;
    my $function;
    my $gene; ### name of locus
    my @unit=split("\n",$refgff->{$_});
    my $name;### name of trancript
    my $tss; ### transcription start site
    my $tes; ### transcription end site
    my $css; ### coding start site
    my $ces; ### coding end site
    my @exons; ### exon start array
    my @exone; ### exon end array
    my $exonnumber;
    foreach(@unit){
       if ($_=~/^##/){
          next;
       }else{
          my @word=split("\t",$_);
          if ($word[2] eq "mRNA"){
             $chr=$word[0];
             $strand=$word[6];
             if($word[8]=~/ID=(.*);/){
                 $gene=$1;
                 $function="NA";
             } 
             $name=$gene;
             $tss =$word[3];
             $tes =$word[4];
          }elsif($word[2] eq "CDS"){
             $exonnumber++;
             push (@exons,$word[3]);
             push (@exone,$word[4]);
          }
       }
    }
    @exons=sort {$a <=> $b} @exons;
    @exone=sort {$a <=> $b} @exone;
    $css=$exons[0];
    $ces=$exone[$#exone];
    my $exonsline=join(",",@exons);
    my $exoneline=join(",",@exone);   
    my $temp="$name:$function";
    print OUT "$gene\t$name:$function\t$chr\t$strand\t$tss\t$tes\t$css\t$ces\t$exonnumber\t$exonsline\t$exoneline\n";
}
close OUT;

########################################################

sub parseGFF
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $seq;
my $id;    ##ID for element
my $record;##all line for this record
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/mRNA/){
        $seq=$unit[0];
        if ($unit[8]=~/ID=(.*);/ or $unit[8] =~/ID=(.*)/){
            $id=$1;
        }
        $record="$_\n";
        $hash{$id}=$record;
    }elsif($unit[0] eq $seq and $unit[8] =~ /Parent=$id/){
        $hash{$id}.="$_\n";
    }

}
close IN;

return \%hash;
}

