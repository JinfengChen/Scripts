#!/usr/bin/perl
use Getopt::Long;
use Data::Dumper;

GetOptions (\%opt,"gff:s","help");


my $help=<<USAGE;
Format Gramene GFF file into BGI GFF, which was used more frequently.
For each gene, we select the longest mRNA(CDS) to represent this gene.
psuedogene and untranslated gene were parsed out.
perl $0 -gff 
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

changeGFFgramene($opt{gff});

my %hash;
my ($refgff,$refpos)=parseGFF("OBa.gramene.raw.gff");
open OUT , ">OBa.gramene.unsort.gff" or die "$!";
foreach(keys %$refgff){
   unless ($hash{$refpos->{$_}}){
      print OUT "$refgff->{$_}";
      $hash{$refpos->{$_}}=1;
   }
}
close;
`msort -k 1,n4 OBa.gramene.unsort.gff > OBa.gramene.gff`;

#####################
sub changeGFFgramene
{
my ($gff)=@_;
my $prefix="OBa_g";
####read gramene gff and store inf as gene->mrna->CDS
my %gene;
my %genetype;
my %evidence;
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/ or $_=~/Component/);
    my $record=$_;
    my @unit=split("\t",$_);
    if ($unit[2]=~/gene/){
        if ($unit[8]=~/ID=(.*?);Description=(.*)$/){
           $genetype{$1}=$2;
        }
    }elsif($unit[2]=~/mRNA/ or $unit[2]=~/pseudoRNA/){
        if ($unit[8]=~/ID=(.*?);Parent=(.*?);Evidence=(.*)$/){
            $evidence{$1}=[$2,$3];
        }
    }elsif($unit[2] eq "CDS"){
        if ($unit[8]=~/ID=(.*?);Parent=(.*?)$/){
            push (@{$gene{$evidence{$2}->[0]}->{$2}},$record);
        }
    }

}
close IN;


####format GFF into BGI GFF and select longest mRNA/CDS
my %count;
my %final;
open OUT, ">OBa.gramene.alltranscript.gff" or die "$!";
foreach my $g (keys %gene){
    $count{$genetype{$g}}++;
    foreach my $t (keys %{$gene{$g}}){
       my $len;
       my @mrna;
       my $ref;
       my $strand;
       foreach my $e (@{$gene{$g}->{$t}}){
          my @unit=split("\t",$e);
          $ref=$unit[0];
          $strand=$unit[6];
          $unit[1]="Gramene";
          $unit[8]="Parent=$prefix"."$t;";
          push (@mrna,[@unit]);
       } 
       my @sortmrna=sort {$a->[3] <=> $b->[3]} @mrna;
       my $start=$sortmrna[0][3];
       my $end  =$sortmrna[$#sortmrna][4];
       $len=$end-$start+1;
       my $anno ="ID=$prefix"."$t;";
       unshift (@sortmrna,[$ref,Gramene,mRNA,$start,$end,".",$strand,".",$anno]);
       #######write transcript to file
       foreach my $templine (@sortmrna){
          my $line=join("\t",@$templine);
          print OUT "$line\n";
       }
       #######
       unless (exists $final{$g}) {
           $final{$g}->[1]=\@sortmrna; 
	   $final{$g}->[0]=$len;
       }elsif($len > $final{$g}->[0]){
	   $final{$g}->[1]=\@sortmrna; 
	   $final{$g}->[0]=$len;
       }
    }
}
close OUT;
####write into file
open OUT, ">OBa.gramene.summary" or die "$!";
foreach(keys %count){
    print OUT "$_\t$count{$_}\n";
}
close OUT;

open OUT, ">OBa.gramene.raw.gff" or die "$!";
foreach my $g (keys %final) {
	foreach my $t (@{$final{$g}->[1]}) {
	    my $e=join("\t",@$t);
            print OUT "$e\n";
        }
}
close OUT;
}



sub parseGFF
{
my ($gff)=@_;
my %pos;
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
        if ($unit[8]=~/ID=(.*?);/){
            $id=$1;
        }
        $record="$_\n";
        $hash{$id}=$record;
        $pos{$id}="$unit[0]_$unit[3]";
    }elsif($unit[0] eq $seq and $unit[8] =~ /Parent=$id/){
        $hash{$id}.="$_\n";
    }

}
close IN;

return (\%hash,\%pos);
}

