#!/usr/bin/perl
use Getopt::Long;
use Data::Dumper;

GetOptions (\%opt,"refgff:s","refpep:s","refgff3:s","reffasta:s","qrygff:s","qrypep:s","qrygff3:s","qryfasta:s","align:s","help");


my $help=<<USAGE;

perl $0 -refgff Os.gff -refpep Os.pep -refgff3 Os.gff3 -reffasta Os.fasta -qrygff Ob.gff -qrypep Ob.pep -qrygff3 Ob.gff3 -qryfasta Ob.fasta -align ob_os.align 
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}


$blastall="/home/biosoftware/blast-2.2.14/bin/blastall";
$formatdb="/home/biosoftware/blast-2.2.14/bin/formatdb";
$bestalign="/home/jfchen/FFproject/tools/bin/bestAlign.pl";
$protein2genome="/home/jfchen/FFproject/FFgenome/03.evolution/protein-map-genome/bin/protein_map_genome.pl";


my $refgff=gff($opt{refgff});
our $refpep=getfastaseq($opt{refpep});
my $qrygff=gff($opt{qrygff});
our $qrypep=getfastaseq($opt{qrypep});
my $align =parsealign($opt{align});

my $additionalpair;

open INTER, ">interval.txt" or die "$!";
foreach(sort {$a <=> $b} keys %$align){
	 my $block=$_; 
	 print "Block: $block\n";
	 my $reftemp=$align->{$_};
         my $chr=shift @$reftemp;
         my $qrychrgff=$qrygff->{$chr->[0]};
         my $refchrgff=$refgff->{$chr->[1]};
         print "$chr->[0]\t$chr->[1]\n";
         for(my $i=1;$i<@$reftemp;$i++){
             $counter++;
             my $last=$reftemp->[$i-1];
             my $present=$reftemp->[$i];
             print "Pair: $counter\n";
             print "Last: $last->[0]\t$last->[1]\n";
             print "Present: $present->[0]\t$present->[1]\n";
             ####deal with query####
             my ($qrystart,$qryend);
             if ($qrychrgff->{$last->[0]}->[0] > $qrychrgff->{$present->[0]}->[0]){
                $qryend   =$qrychrgff->{$last->[0]}->[1]+10;
	        $qrystart =$qrychrgff->{$present->[0]}->[0]-10;
             }else{
                $qrystart =$qrychrgff->{$last->[0]}->[0]-10;
                $qryend   =$qrychrgff->{$present->[0]}->[1]+10;
             }
	     my $qrygene =countgene($qrychrgff,$qrystart,$qryend);
             my $tempn1=@$qrygene;
             print "$tempn1\n@$qrygene\n";
             print INTER "$tempn1\t";
             my $qrychr=$chr->[0];
             $qrychr=~s/Ob/chr/;
             my $qryhead="OB"."_".$qrychr."_".$qrystart."_".$qryend;
             ####deal with refer#### 
             my ($refstart,$refend);
             if ($refchrgff->{$last->[1]}->[0] > $refchrgff->{$present->[1]}->[0]){
                $refend   =$refchrgff->{$last->[1]}->[1]+10;
                $refstart =$refchrgff->{$present->[1]}->[0]-10;
             }else{
                $refstart =$refchrgff->{$last->[1]}->[0]-10;
                $refend   =$refchrgff->{$present->[1]}->[1]+10;
             }
             my $refgene =countgene($refchrgff,$refstart,$refend);
             my $tempn2=@$refgene;
             print "$tempn2\n@$refgene\n";
             print INTER "$tempn2\n";
             my $refchr=$chr->[1];
             $refchr=~s/Os/chr/;
             my $refhead="OS"."_".$refchr."_".$refstart."_".$refend;
             ####ortholog pair
             #my $orthpair=ortholog ($qrygene,$refgene);
             #foreach (keys %$orthpair){
             #      print "$_\t$orthpair->{$_}\n";
             #}
             my $temp=dealtwolist ($qrygene,$refgene,$qryhead,$refhead,$counter);
             $additionalpair+=$temp;
         }
}
close INTER;

print "Get additional pair of ortholog: $additionalpair\n";

###############sub functions####################

sub dealtwolist
{
my ($qrygene,$refgene,$qryhead,$refhead,$pairnum)=@_;
my $additionalpair=-2;
my $pairnum="Pair".$pairnum;
& getlistseq ($qrygene,$qrypep,"qry.pep");
& getlistseq ($refgene,$refpep,"ref.pep");
system ("$formatdb -i qry.pep -p T");
system ("$formatdb -i ref.pep -p T");
system ("$blastall -p blastp -i qry.pep -d qry.pep -e 1e-5 -o qry_qry.blastm8 -m 8 > qry.log 2> qry.log2");
system ("$blastall -p blastp -i ref.pep -d ref.pep -e 1e-5 -o ref_ref.blastm8 -m 8 > ref.log 2> ref.log2");
my ($qrymost)=mosthitm8("qry_qry.blastm8");
my ($refmost)=mosthitm8("ref_ref.blastm8");
if ($qrymost >= 0){
#if ($qrymost >= 2 or $refmost >= 2){
#   print "Have tandem duplicate\n";
   #print "$opt{qrygff3}\t$opt{refgff3}\n";
   getsubgff3($opt{qrygff3},$qryhead,"+");
   getsubgff3($opt{refgff3},$refhead,"+");
   getsubfasta($opt{qryfasta},$qryhead,"+");
   getsubfasta($opt{reffasta},$refhead,"+");
   `perl GFF2embl.pl -gff $qryhead.gff -embl $qryhead.embl -fasta $qryhead.fasta`;
   `perl GFF2embl.pl -gff $refhead.gff -embl $refhead.embl -fasta $refhead.fasta`;
   `perl runblast2seq.pl`;
   `perl run2act.pl`;
   `rm *.fasta.n* *.blast *.temp`;
   `mkdir ../output/$pairnum`;
   `mv $qryhead* ../output/$pairnum/`;
   `mv $refhead* ../output/$pairnum/`;
   #`rm qry* ref*`;
#}else{
   my $orthpair=ortholog();
   my %orthqry;###id of qry gene in ortholog pairs
   my %orthref;###id of ref gene in ortholog pairs 
   foreach (keys %$orthpair){
      $additionalpair++;
      $orthqry{$_}=1;
      $orthref{$orthpair->{$_}}=1;
      print "$_\t$orthpair->{$_}\n";
   }
   ####deal non-collinarity gene in qry
   if (@$qrygene > keys %orthqry){
      foreach(@$qrygene){
         unless (exists $orthqry{$_}){
            my $testgene=$_;
            my @tempgene;
            push (@tempgene,$testgene);
            & getlistseq (\@tempgene, $qrypep, "temp.pep");
            system ("$blastall -p blastp -i temp.pep -d ref.pep -e 1e-5 -o temp_ref.blastm8 -m 8 > temp.blast.log 2> temp.blast.log2");
            my ($tempmost)=mosthitm8("temp_ref.blastm8");
            & getsubgff3($opt{refgff3},$refhead,"+");
            & getsubfasta($opt{reffasta},$refhead,"+");
            if ($tempmost == 0){ ### if this gene have no homolog gene in ref region
                `perl $protein2genome -cpu 1 -run multi -step 12 temp.pep $refhead.fasta > temp.log 2> temp.log2`;
                my $temphit=parsesolornr("temp.pep.blast.solar.filter.nr");
                if ($temphit == 0){
                   & write2file("$testgene\t$pairnum","FF.specific.id");
                }else{
                   `perl $protein2genome -cpu 1 -run multi -step 34 temp.pep $refhead.fasta > temp.log 2> temp.log2`;
                   ## add this gene in Rice genome and add this pair of ortholog to orthlog list
                   my $reftempgff=parseGFF("temp.pep.solar.genewise.gff");
                   if (keys %$reftempgff == 0){
                       & write2file("$testgene\t$pairnum","FF.specific.id"); 
                   }else{
                       foreach (keys %$reftempgff){
                           #& write2file($testgene, "FF.lostinRice.id");
                           & write2file($reftempgff->{$_}, "FF.add2Rice.gff"); 
                       }
                       & write2file("$testgene\t$pairnum", "FF.SingleLostRice.id");
                   } 
                }
            }else{
                `perl $protein2genome -cpu 1 -run multi -step 1234 temp.pep $refhead.fasta > temp.log 2> temp.log2`;
                ## filter out overlap gene model with qrygff and add the additional gene into list.
                my $reftempgff=parseGFF("temp.pep.solar.genewise.gff");
                my $newgff=0;
                foreach(keys %$reftempgff){
                     my $overlap=overlapgff("$refhead.gff",$reftempgff->{$_});
                     unless ($overlap){
                         $newgff++;
                         #& write2file($testgene, "FF.lostinRice.id");
                         & write2file($reftempgff->{$_}, "FF.add2Rice.gff");
                     }
                }
                if ($newgff == 0){
                     & write2file("$testgene\t$pairnum", "FF.duplcated.id");
                }else{
                     & write2file("$testgene\t$pairnum", "FF.FamilyLostRice.id");
                }
            }
            `rm -R temp* $refhead*`;  
         }
      }
   }

   ####deal non-collinarity gene in ref



   `rm qry* ref*`;   
}
return $additionalpair;
}


sub ortholog
{
my %best_pair1;
my %best_pair2;
my %best_pair;
system ("$blastall -p blastp -i qry.pep -d ref.pep -e 1e-5 -o qry_ref.blastm8 -m 8 > qry.log 2> qry.log2");
system ("$blastall -p blastp -i ref.pep -d qry.pep -e 1e-5 -o ref_qry.blastm8 -m 8 > ref.log 2> ref.log2");
system ("perl $bestalign -f m8 qry_ref.blastm8 > qry_ref.blastm8.best");
system ("perl $bestalign -f m8 ref_qry.blastm8 > ref_qry.blastm8.best");
& read_m8_best_pair ("qry_ref.blastm8.best",\%best_pair1);
& read_m8_best_pair ("ref_qry.blastm8.best",\%best_pair2);
#`rm qry* ref*`;
foreach my $pep1 (sort keys %best_pair1) {
        my $pep2 = $best_pair1{$pep1};
        if (exists $best_pair2{$pep2} && $best_pair2{$pep2} eq $pep1) {
           $best_pair{$pep1} = $pep2;
        }
}
return \%best_pair;
}


sub write2file
{
my ($content,$file)=@_;
open OUTFILE, ">>$file" or die "$!";
    print OUTFILE "$content\n";
close OUTFILE;
}

sub mosthitm8
{
my ($m8)=@_;
my $most=0;
my %hash;
open IN, "$m8" or die "$!";
while(<IN>){
   chomp $_;
   my @unit=split("\t",$_);
   if (exists $hash{$unit[0]}){
       $hash{$unit[0]}++;
   }else{
       $hash{$unit[0]}=1; ### self hit is included
   }
}
close IN;
foreach(keys %hash){
   if ($hash{$_} > $most){
      $most=$hash{$_};
   }
}
return ($most);
}

sub parsesolornr
{
my ($solornr)=@_;
my $hit=0;
open SOLOR, "$solornr" or die "$!";
while(<SOLOR>){
    chomp $_;
    $hit++ unless ($_ eq "");
}
close SOLOR;
return $hit;
}


sub getlistseq
{
my ($refarray,$refseq,$outfile)=@_;
open TEMP, ">$outfile" or die "$!";
foreach(@$refarray){
    if (exists $refseq->{$_}){
        print TEMP ">$_\n$refseq->{$_}\n";
    }
}
close TEMP;
}

sub read_m8_best_pair {
        my $file = shift;
        my $hash_p = shift;
        open IN,$file || die "fail $file";
        while (<IN>) {
                my @t = split /\t/;
                $hash_p->{$t[0]} = $t[1];
        }
        close IN;
}



sub getfastaseq
{
$/=">";
my %hash;
my ($file)=@_;
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my @temp1=split(" ",$temp);
    my $head=$temp1[0];
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    #print "$head\n";
    $hash{$head}=$seq;
}
$/="\n";
return \%hash;
}


sub countgene
{
my ($chrgff,$start,$end)=@_;
#if ($start > $end){
#   my $temp=$start;
#   $start= $end;
#   $end  = $temp;
#}
print "$start\t$end\n";
my $number=0;
my @gene=sort {$chrgff->{$a}->[0] <=> $chrgff->{$b}->[0]} keys %$chrgff;
my @count;
foreach(@gene){
    if ($chrgff->{$_}->[0] >= $start and $chrgff->{$_}->[0] <= $end) {
        $number++;
        push (@count,$_);
    }elsif($chrgff->{$_}->[0] > $end){
	return \@count;
    }
}
return \@count;
}





sub parsealign
{
my ($align)=@_;
my %hash;
my $block;
my $qry;
my $ref;
open IN, "$align" or die "$!";
while(<IN>){
    chomp $_;
    if ($_=~/\#\# Alignment (\d+)\: score=.* e_value=.* N=\d+ (\w+)\&(\w+) (\w+)/) {
	$block=$1;
	$qry=$2;
	$ref=$3; 
	#print "$block\t$qry\t$ref\n";
    }elsif($_=~/\s*$block\-\s*\d+\:/){
	my @unit=split("\t",$_);
	#print "1:$unit[1]\t2:$unit[2]\n";
	if (exists $hash{$block}) {
           my $reftemp=$hash{$block};
	   push (@$reftemp,[$unit[1],$unit[2]]);
	   $hash{$block}=$reftemp;		
	}else{
           my @temp;
           my @unit=split("\t",$_);
           push (@temp,[$qry,$ref]);
	   push (@temp,[$unit[1],$unit[2]]);
           $hash{$block}=\@temp;
        }
     }	
}
close IN;
return \%hash;
}

sub getsubfasta
{
my ($fasta,$head,$strand)=@_;
my @temp=split("_",$head);
my $chr=$temp[1];
my $start=$temp[2];
my $end  =$temp[3];
my $len=$end-$start+1;
my $refseq=getfastaseq($fasta);
if ($strand eq "+"){
   if (exists $refseq->{$chr}){
       my $subseq=substr($refseq->{$chr},$start,$len);
       open FAS, ">$head.fasta" or die "$!"; 
             print FAS ">$head\n$subseq\n";
       close FAS;  
   }else{
       print "$chr can not found in $fasta\n";
   }
}else{
   if (exists $refseq->{$chr}){
       my $subseq=substr($refseq->{$chr},$start,$len);
       $subseqrec=revcom($subseq);
       open FAS, ">$head.fasta" or die "$!";
             print FAS ">$head rec\n$subseqrec\n";
       close FAS;
   }else{
       print "$chr can not found in $fasta\n";
   }   
}
}


sub revcom
{
my ($seq)=@_;
my $rev=reverse $seq;
$rev=~tr/ATGCatgc/TACGtacg/;
return $rev;
}



sub getsubgff3
{
my ($gff,$head,$strand)=@_;
my @temp=split("_",$head);
my $chr=$temp[1];
my $start=$temp[2];
my $end  =$temp[3];
my $len=$end-$start+1;
#print "$chr\t$start\t$end\n";
if ($strand eq "+"){
   open IN1, "$gff" or die "$!";
   open OUT1, ">>$head.gff" or die "$!";
   while (<IN1>){
     my @unit=split("\t",$_);
     if ($unit[0] eq $chr and $unit[3] >= $start and $unit[4] <= $end){ 
       $unit[3]=$unit[3]-$start;
       $unit[4]=$unit[4]-$start;
       my $line=join("\t",@unit);
       print OUT1 "$line";
     } 

   }
   close OUT1;
   close IN1;
}else{
   open IN1, "$gff" or die "$!";
   open OUT1, ">>$head.gff" or die "$!";
   while (<IN1>){
     my @unit=split("\t",$_);
     if ($unit[0] eq $chr and $unit[3] >= $start and $unit[4] <= $end){ 
        my $tempend   =$len-($unit[3]-$start);
        my $tempstart =$len-($unit[4]-$start); 
        $unit[3]=$tempstart;
        $unit[4]=$tempend;
        if ($unit[6] eq "+"){
            $unit[6] = "-";
        }else{
            $unit[6] = "+";
        }
        my $line=join("\t",@unit);
        print OUT1 "$line";
     }
   }
   close OUT1;
   close IN1;
}
}

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

sub overlapgff
{
my ($gff,$gffrecord)=@_;
my @array1=split("\n",$gffrecord);
my @array2=split("\t",$array1[0]);
my $start=$array2[3];
my $end  =$array2[4];
my $overlap=0;
open GFF, "$gff" or die "$!";
while(<GFF>){
   chomp $_;
   next if ($_ eq "");
   my @unit=split("\t",$_);
   next if ($unit[2] ne "mRNA");
   unless ($start > $unit[4] or $end < $unit[3]){
      $overlap = 1;
      return $overlap;
   }
}
close GFF;
return $overlap;
}


sub gff
{
my ($gff)=@_;
my %hash;
open IN, "$gff" or die "$!";
while(<IN>){
        chomp $_;
	next if ($_ eq "");
	my @unit=split("\t",$_);
	if (exists $hash{$unit[0]}) {
        my $reftemp=$hash{$unit[0]};
	    $reftemp->{$unit[1]}=[$unit[2],$unit[3]];
            $hash{$unit[0]}=$reftemp;
	}else{
	    my %temp;
	    $temp{$unit[1]}=[$unit[2],$unit[3]];
	    $hash{$unit[0]}=\%temp;
	}
}
close IN;
return \%hash;
}
