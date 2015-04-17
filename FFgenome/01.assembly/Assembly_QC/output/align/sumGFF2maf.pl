#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"--genegff:s","--tegff:s","--maf:s","verbose","help");


my $help=<<USAGE;
Given gene, te gff file and maf alignment. Summary the number and length of gene that not present or present
as gaps in compared sequence.
perl $0 --genegff --tegff --maf
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my @maf=glob("$opt{maf}/*.maf");
open OUT, ">summary.txt" or die "$!";
print OUT "BAC\tDNA TE\tLOSS\tPercent%\tRNA TE\tLOSS\tPercent%\tGENE NUM\tLOSS\tPercent%\tGENE LEN\tLOSS\tPercent%\n";
for (my $i=0;$i<@maf;$i++){
   if ($maf[$i]=~/maf\/(.*?)\.fas\.maf/){
      my $title=$1;
      print "$title\t$maf[$i]\n";
      my $gene="$opt{genegff}/$title".".gff";
      my $te  ="$opt{tegff}/$title".".gff";
      if (-f $gene and -f $te){
          print "$gene\n$te\n";
          my $geneinf=sumgene($gene,$maf[$i]);
          my $teinf=sumte($te,$maf[$i]);
          my $p1=sprintf("%.1f%",($teinf->{DNA}->[1]/$teinf->{DNA}->[0])*100);
          my $p2=sprintf("%.1f%",($teinf->{RNA}->[1]/$teinf->{RNA}->[0])*100);
          my $p3=sprintf("%.1f%",($geneinf->{number}->[1]/$geneinf->{number}->[0])*100);
          my $p4=sprintf("%.1f%",($geneinf->{len}->[1]/$geneinf->{len}->[0])*100);
          print OUT "$title\t$teinf->{DNA}->[0]\t$teinf->{DNA}->[1]\t$p1\t$teinf->{RNA}->[0]\t$teinf->{RNA}->[1]\t$p2\t";
          print OUT "$geneinf->{number}->[0]\t$geneinf->{number}->[1]\t$p3\t$geneinf->{len}->[0]\t$geneinf->{len}->[1]\t$p4\n";
      }
   }
}
close OUT;

sub sumgene
{
my ($gff,$maf)=@_;
my $refgff=GENEGFF($gff);

my %hash;
foreach my $id (keys %$refgff){
   print ">>>Start a $refgff->{$id}->[2] TE\n" if $opt{verbose};
   my ($align,$gap)=summaf($refgff->{$id},$maf);
   my $len=$refgff->{$id}->[1]-$refgff->{$id}->[0]+1;
   print "TE(id,start,end): $id\t$refgff->{$id}->[0]\t$refgff->{$id}->[1]\t" if $opt{verbose};
   print "RESULT(type,len,align,gap): $refgff->{$id}->[2]\t$len\t$align\t$gap\n" if $opt{verbose};
   $hash{$refgff->{$id}->[2]}->[0]+=$len;
   $hash{$refgff->{$id}->[2]}->[1]+=$align;
   $hash{$refgff->{$id}->[2]}->[2]+=$gap;
}
my %inf;
my @type=keys %hash; ## gene id
my ($gnum,$gnumloss,$glen,$galign,$glenloss);
$glenloss=0;
$gnumloss=0;
print "Type\tLength\tAlign\tGap\tLoss\n";
$gnum=@type;
for(my $i=0;$i<@type;$i++){
   $glen+=$hash{$type[$i]}->[0];
   $galign+=$hash{$type[$i]}->[1];
   my $temp=$hash{$type[$i]}->[0]-$hash{$type[$i]}->[1]+1;
   print "GENE ALIGN (gene,len,align,loss): $type[$i]\t$hash{$type[$i]}->[0]\t$hash{$type[$i]}->[1]\t$temp\n" if $opt{verbose};
   if ($hash{$type[$i]}->[1] < ($hash{$type[$i]}->[0]) * 0.5){
      print "GENE LOSS (gene,len,align): $type[$i]\t$hash{$type[$i]}->[0]\t$hash{$type[$i]}->[1]\n" if $opt{verbose};
      $gnumloss++;
   }
}
$glenloss=$glen-$galign+1;
$inf{number}=[$gnum,$gnumloss];
$inf{len}=[$glen,$glenloss];
return \%inf;   
}


sub sumte
{
my ($gff,$maf)=@_;
my $refgff=TEGFF($gff);
my @type=("En-Spm","MuDR","MITE","hAT","DNA","Gypsy","Copia","LTR","Helitron","OtherTE");

my %hash;
foreach my $id (keys %$refgff){
   print ">>>Start a $refgff->{$id}->[2] TE\n" if $opt{verbose};
   my ($align,$gap)=summaf($refgff->{$id},$maf);
   my $len=$refgff->{$id}->[1]-$refgff->{$id}->[0]+1;
   print "TE(id,start,end): $id\t$refgff->{$id}->[0]\t$refgff->{$id}->[1]\t" if $opt{verbose};
   print "RESULT(type,len,align,gap): $refgff->{$id}->[2]\t$len\t$align\t$gap\n" if $opt{verbose};
   $hash{$refgff->{$id}->[2]}->[0]+=$len;
   $hash{$refgff->{$id}->[2]}->[1]+=$align;
   $hash{$refgff->{$id}->[2]}->[2]+=$gap;
}
my %inf;
print "Type\tLength\tAlign\tGap\tLoss\n";
for(my $i=0;$i<@type;$i++){
   my $lost = $hash{$type[$i]}->[0]-$hash{$type[$i]}->[1]-$hash{$type[$i]}->[2];
   print "$type[$i]\t$hash{$type[$i]}->[0]\t$hash{$type[$i]}->[1]\t$hash{$type[$i]}->[2]\t$lost\n";
   if ($type[$i]=~/LTR/ or $type[$i]=~/Gypsy/ or $type[$i]=~/Copia/){
      $inf{"RNA"}->[0]+=$hash{$type[$i]}->[0];
      #$inf{"RNA"}->[1]+=$hash{$type[$i]}->[1];
      $inf{"RNA"}->[1]+=$hash{$type[$i]}->[0]-$hash{$type[$i]}->[1]+1;
   }else{
      $inf{"DNA"}->[0]+=$hash{$type[$i]}->[0];
      #$inf{"DNA"}->[1]+=$hash{$type[$i]}->[1];
      $inf{"DNA"}->[1]+=$hash{$type[$i]}->[0]-$hash{$type[$i]}->[1]+1;
   }
}
return \%inf;
}

#a score=72767
#s 61n13       923 911 + 125279 AAGCTTGAGAAAAAGATATGGGTGAGAGTTTCTGCTGCACTGGAAAAAAAATAGTATCTCCAGAGATCTCTGGAAATCAAGACAATCTCTCTCT
#s scaffoldFF  999 893 + 125926 AAGCTTGAGAAAAAGATATGGGTGAGAGTTTCTGCTGCACTGGAAAAAAAATAGTATCTCCAGAGATCTCTGGAAATCAAGACAATCTCTCTCT
sub summaf
{
my ($element,$maf)=@_; # $element->[0],[1],[2]: start, end, class
my $alignlen=0;
my $gaplen=0; # alignlen should be gap free alignmeng len
my @array;
open IN, "$maf" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    next if ($_=~/^#/);
    my $line=$_;
    if ($line=~/^a/){ # every aligned segment
       my $refseq=<IN>;
       my $ref=mafline($refseq);
       my $qryseq=<IN>;
       my $qry=mafline($qryseq);
       my $start=$element->[0];
       my $end  =$element->[1];
       my $class=$element->[2];
       my $length=$end-$start+1;
       print "$start\t$end\t$ref->[0]\t$ref->[1]\t$qry->[0]\t$qry->[1]\n" if $opt{verbose};
       my ($seqstart,$seqend,$type); ## start and end in alignment, not including gap, need to tranform to include gap
       unless ($end < $ref->[0] or $start > $ref->[1]){ # overlap here
           #print "OVERLAP: $start\t$end\t$ref->[0]\t$ref->[1]\t$qry->[0]\t$qry->[1]\n";
           if ($end > $ref->[0] and $end <= $ref->[1] and $start <= $ref->[0]){  ###  overlap in left of ref
              $seqstart=$ref->[0];
              $seqend  =$end;
              $type = "left";
           }elsif($start < $ref->[1] and $start > $ref->[0] and $end >= $ref->[1]){ ## overlap in right of ref
              $seqstart=$start;
              $seqend  =$ref->[1];
              $type ="right";
           }elsif($start >= $ref->[0] and $end <= $ref->[1]){  ## element in ref
              $seqstart=$start;
              $seqend  =$end;
              $type="in";
           }elsif($start < $ref->[0] and $end > $ref->[1]){  ## element cover ref
              $seqstart=$ref->[0];
              $seqend  =$ref->[1];
              $type="cover";
           }
           print "MAF: Type, $type\tSeqlength, $length\tStart in align, $seqstart\tEnd in align, $seqend\n" if $opt{verbose};
           my ($align,$gap)=checkalign($seqstart,$seqend,$ref,$qry,$type);##transform postion to include gap, get gap free align length and gap length
           push (@array,[$seqstart,$seqend,$align,$gap]);
           #$alignlen+=$align;
           #$gaplen+=$gap;
       }
    }
}
    #deal with overlap bewtten align
    for (my $i=1;$i<@array;$i++){
        if ($array[$i]->[0] < $array[$i-1]->[1] and $array[$i]->[0] > $array[$i-1]->[0] and $array[$i]->[1] > $array[$i-1]->[1]){
          
           my $overlap=$array[$i-1]->[1]-$array[$i]->[0]+1;
           $array[$i]->[2]-=$overlap;
           print "overlap: $array[$i]->[0]\t$array[$i-1]->[1]\t$overlap\t$array[$i]->[2]\n" if $opt{verbose};
           
        }elsif($array[$i]->[0] < $array[$i-1]->[1] and $array[$i]->[1] > $array[$i-1]->[1]){ # not overlap, possible redudency
           $array[$i]->[2]=0;
           $array[$i]->[3]=0;
        }
    }
    for (my $i=0;$i<@array;$i++){
        print "Adding: align, $array[$i]->[2]\tgap, $array[$i]->[3]\n if $opt{verbose}";
        $alignlen+=$array[$i]->[2];
        $gaplen+=$array[$i]->[3];
    }
close IN;
return ($alignlen,$gaplen);
}

sub checkalign
{
my ($start,$end,$ref,$qry,$type)=@_;
my ($align,$gap);
my @base=split("",$ref->[2]);
my ($count,$pos1,$pos2); ## $count: count of none - base, $pos1 and $pos2 postion in alignment
$count=$ref->[0];# start from the start coordinate of sequence
for(my $i=0;$i<@base;$i++){
   if ($base[$i]=~/\w+/){
      if ($count == $start){
         $pos1=$i;
      }
      if ($count == $end){
         $pos2=$i;
         last;
      }
      $count++;
   }
}
my $len=$pos2-$pos1+1;
my $copy =substr($ref->[2],$pos1,$len); ### reference sequence
my $check=substr($qry->[2],$pos1,$len); ### check sequence
my ($align,$gap)=align($copy,$check);
#my $align=$check=~tr/atcgnATCGN/atcgnATCGN/;
#my $gap=$len-$align;
print "ALIGN: $align\tGAP: $gap\n" if $opt{verbose};
print "Overlap type: $type\n" if $opt{verbose};
print "Transformed position: P1, $pos1\tP2, $pos2\n" if $opt{verbose};
print "$copy\n$check\n" if $opt{verbose};
return ($align,$gap);
}

sub align
{
my ($copy,$check)=@_;
my $align=0;
my $gap=0;
my @ref=split("",$copy);
my @qry=split("",$check);
for (my $i=0;$i<@ref;$i++){
   if ($ref[$i]=~/[atcg]/i){
      $ref[$i]=~tr/atcg/ATCG/;
      $qry[$i]=~tr/atcg/ATCG/;
      if ($qry[$i] eq $ref[$i]){
         $align++;
      }else{
         $gap++;
         print "$i\t$qry[$i]\t$ref[$i]\n" if $opt{verbose};
      }
   }
}
return ($align,$gap);
}

sub overlap_size {
        my $block1_p = shift;
        my $block2_p = shift;
        
        my $combine_start = ($block1_p->[0] < $block2_p->[0]) ?  $block1_p->[0] : $block2_p->[0];
        my $combine_end   = ($block1_p->[1] > $block2_p->[1]) ?  $block1_p->[1] : $block2_p->[1];
        
        my $overlap_size = ($block1_p->[1]-$block1_p->[0]+1) + ($block2_p->[1]-$block2_p->[0]+1) - ($combine_end-$combine_start+1);

        return $overlap_size;
}



#s scaffoldFF  999 893 + 125926 AAGCTTGAGAAAAAGATATGGGTGAGAGTTTCTGCTGCACTGGAAAAAAAATAGTATCTCCAGAGATCTCTGGAAATCAAGACAATCTCTCTCT
sub mafline
{
my ($line)=@_;
my @unit=split(" ",$line);
my $start=$unit[2];
my $end  =$unit[2]+$unit[3]-1;
my $seq  =$unit[6]; #have gap
#$refseq=~s/\-//g;
my @inf=($start,$end,$seq);
return \@inf;
}

sub GENEGFF
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $id;    ##ID for element
my $class;
my $count;
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/CDS/){
        if ($unit[8]=~/Parent=(.*?);/){
            $count++;
            $id="exon".$count;
            $class=$1;
            $hash{$id}=[$unit[3],$unit[4],$class];
        }
    }
}
close IN;
return \%hash;
}



sub TEGFF
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $seq;
my $id;    ##ID for element
my $class;
my $record;##all line for this record
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/Trans/){
        $seq=$unit[0];
        if ($unit[8]=~/ID=(.*?);.*Class=(.*?);/){
            $id=$1;
            $class=$2;
            my $len=$unit[4]-$unit[3]+1;
            if ($class=~/DNA/){
               if ($class=~/En-Spm/ or $class=~/CACTA/){
                  $class="En-Spm";
               }elsif($class=~/MuDR/ or $class=~/MULE/){
                  $class="MuDR";
               }elsif($class=~/Stowaway/ or $class=~/Tourist/){
                  $class="MITE";
               }elsif($class=~/(hAT)/){
                  $class=$1;
               }else{
                  $class="DNA";
               }
            }elsif($class=~/LTR/){
               if ($class=~/(Gypsy)/i){
                  $class="Gypsy";
               }elsif($class=~/(Copia)/i){
                  $class="Copia";
               }else{
                  $class="LTR";
               }
            }elsif($class=~/Helitron/){
               $class="Helitron";
            }else{
               $class="OtherTE";
            }
            $hash{$id}=[$unit[3],$unit[4],$class];
        }
    }
}
close IN;
return \%hash;
}
 
