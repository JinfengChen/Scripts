#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"fasta:s","link:s","gff:s","cut:s","project:s","help");


my $help=<<USAGE;
Link big-super-scaffold from fpc based manual check table.
The scaffold sequence that mapped to the fpc physical map is super-scaffold, which construted by BES scaffolding. 
After link, we create two new sequence file, one for big-super-scaffold and the other for chromosome.
We also create two gff gene annotation file for each new fasta file.
-fasta:   scaffold sequence
-link:    link table
-gff:     gff format gene annotation, optional.
-cut:     contain split point table for scaffold which are wrong assembled, optional.
-project: project name used to name output sequence file

Example:
Chr	Super	Scaffold	Strand	Evidence
chr1	1	Scaffold000152	plus	rice/fpc
chr1	1	Scaffold000414	minus	fpc

Run: perl ../superscafV2.pl -fasta /share/raid12/chenjinfeng/tools/Lastz/input/mask/Scaffold20k.fa.RepeatMasker.masked -link ../../input/joinScaffold.txt -cut ../../input/split.txt -project aMasked_OBa > log
USAGE

if ($opt{help} or keys %opt < 1){
   print "$help\n";
   exit ();
}

my $scafseq=getfastaseq($opt{fasta});
our $totallen;
our %split;
our %bstart; ### b fregment start position in original scaffold
##########split scaffold according the cut point, which are identified to be misassembly###############################
if (exists $opt{cut}){
   open IN, "$opt{cut}" or die "$!";
      while (<IN>){
         next if (length $_ < 2);
         my @unit=split("\t",$_);
         if (exists $scafseq->{$unit[0]}){
             my $aname=$unit[0]."a";
             my $bname=$unit[0]."b";
             my $aseq =substr($scafseq->{$unit[0]},0,$unit[1]);
             my $bseq =substr($scafseq->{$unit[0]},$unit[1]);
             my $newaseq=trim3N($aseq);
             my $newbseq=trim5N($bseq); 
             $scafseq->{$aname}=$newaseq;
             $scafseq->{$bname}=$newbseq;
             $split{$unit[0]}=$unit[1];
             my $templen1=length $scafseq->{$unit[0]};
             my $templen2=length $newbseq;
             $bstart{$bname}=$templen1 - $templen2 + 1; 
             delete $scafseq->{$unit[0]};
         }
      }
   close IN;
}
################################################################################################################

###############read data from link table and join the scaffold################################################
my $n="NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
my %chr;
my %super;
my %scaf2chr;
my %scaf2super;
open IN, "$opt{link}" or die "$!";
<IN>;
while (<IN>){ 
    next if ($_ eq "");
    my @unit=split("\t",$_);
    my $supername=sprintf("Scaffold%06d",$unit[1]);
    my $chrname;
    if ($unit[0]=~/chr(\d+)/){
       $chrname=sprintf("chr%02d",$1);
    }
    if ($unit[3] eq "plus"){
       if (exists $scafseq->{$unit[2]}){ 
          if (exists $chr{$chrname}){
              $chr{$chrname}.=$n.$scafseq->{$unit[2]};
          }else{
              $chr{$chrname}=$scafseq->{$unit[2]};
          }
          my $chrlen=length $chr{$chrname};
          my $scaflen=length $scafseq->{$unit[2]};
          my $chrstart =$chrlen-$scaflen+1;
          my $strand="+";
          my @temp1=[$chrname,$chrstart,$scaflen,$strand];
          $scaf2chr{$unit[2]}=\@temp1;          
          
          #print "$unit[2]\t$chrname\t$chrstart\n";

          if (exists $super{$unit[1]}){
              $super{$unit[1]}.=$n.$scafseq->{$unit[2]};
          }else{
              $super{$unit[1]}=$scafseq->{$unit[2]};
          }
          my $suplen=length $super{$unit[1]};
          my $supstart=$suplen-$scaflen+1;
          my @temp2=[$supername,$supstart,$scaflen,$strand];
          $scaf2super{$unit[2]}=\@temp2;

          #print "$unit[2]\t$supername\t$supstart\t+\n";
          delete $scafseq->{$unit[2]};
       }else{
          print "Can not found sequence of $unit[2]\n";
       }
    }else{
       if (exists $scafseq->{$unit[2]}){
          my $recseq=revcom($scafseq->{$unit[2]});
          if (exists $chr{$chrname}){
              $chr{$chrname}.=$n.$recseq;
          }else{
              $chr{$chrname}=$recseq;
          }
          my $chrlen=length $chr{$chrname};
          my $scaflen=length $scafseq->{$unit[2]};
          my $chrstart =$chrlen-$scaflen+1;
          my $strand="-";
          my @temp1=[$chrname,$chrstart,$scaflen,$strand];
          $scaf2chr{$unit[2]}=\@temp1;

          if (exists $super{$unit[1]}){
              $super{$unit[1]}.=$n.$recseq;
          }else{
              $super{$unit[1]}=$recseq;
          }
          my $suplen=length $super{$unit[1]};
          my $supstart=$suplen-$scaflen+1;
          my @temp2=[$supername,$supstart,$scaflen,$strand];
          $scaf2super{$unit[2]}=\@temp2;

          #print "$unit[2]\t$supername\t$supstart\t-\n";
          delete $scafseq->{$unit[2]};
       }else{
          print "Can not found sequence of $unit[2]\n";
       }
    }

}
close IN;


################# write chromosome and super into files###########################################
open OUT, ">$opt{project}.chr.fa" or die "$!";
foreach (sort keys %chr){
    my $fseq=formatseq($chr{$_},50);
    print OUT ">$_\n$fseq\n";
}
close OUT;

open OUT, ">$opt{project}.super.fa" or die "$!";
foreach (sort {$a <=> $b} keys %super){
    my $fseq=formatseq($super{$_},50);
    printf OUT (">Scaffold%06d\n",$_);
    print OUT "$fseq\n";
}
close OUT;

my $leftlen;
my $leftscafn=36;
our %unmapscaf;
open OUT1, ">>$opt{project}.chr.fa" or die "$!";
open OUT2, ">>$opt{project}.super.fa" or die "$!";
foreach (sort keys %$scafseq){
   $leftscafn++;
   $leftlen+=length $scafseq->{$_};
   my $fseq=formatseq($scafseq->{$_},50);
   #my $leftscaf2=sprintf("Scaffold%06d",$leftscafn);
   #my $leftscaf1=sprintf("super%04d",$leftscafn);
   my $leftscaf2=$_;
   my $leftscaf1=$_;
   if (exists $unmapscaf{$_}){
      print "Multi Scaffold in Unmaped Scaffold\n";
   }else{
      $unmapscaf{$_}=[$leftscaf1,$leftscaf2];
      print OUT1 ">$leftscaf1\n$fseq\n";
      print OUT2 ">$leftscaf2\n$fseq\n";
   }
}
close OUT1;
close OUT2;


################# write gff into files###########################################################

if (-f $opt{gff}){

& mergegff($opt{gff},\%scaf2chr,"$opt{project}.chr");
& mergegff($opt{gff},\%scaf2super,"$opt{project}.super");

}

my $rate=$leftlen/$totallen;
print "Left length: $leftlen\n";
print "Percent: $rate\n";
##############################################################################################################


=pod
foreach (keys %$scafseq){
    my $seq=formatseq($scafseq->{$_},100);
    print "$seq\n";
}
for (my $i=1;$i<=10;$i++){
    printf ("Scaffold%08d\n",$i);
    
}
=cut


sub trim5N
{
my ($seq)=@_;
my $revseq=reverse $seq;
my $trimseq=trim3N($revseq);
my $newseq=reverse $trimseq;
return $newseq;
}

sub trim3N
{
my ($seq)=@_;
my $do=1;
while($do==1){
  my $last=substr($seq,-1);
  if ($last =~/N/i){
     chop $seq
  }else{
     $do=0;
  }
}
return $seq;
}



sub formatseq
{
### format a single line sequence into lines with user specific length
my ($seq,$step)=@_;
my $length=length $seq;
my $run=int ($length/$step);
my $newseq;
for(my $i=0;$i<=$run;$i++){
   my $start=$i*$step;
   my $line=substr($seq,$start,$step);
   $newseq.="$line\n";
}
return $newseq;
}


sub revcom
{
my ($seq)=@_;
my $rev=reverse $seq;
$rev=~tr/ATGCatgc/TACGtacg/;
return $rev;
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
    $totallen+=length $seq;
    $hash{$head}=$seq;
}
$/="\n";
print "Total length: $totallen\n";
return \%hash;
}


############gff merge#####################################################
### $scaf2chr is a ref of hash scaffold->[ref][start][scaflen][strand]####
sub mergegff {
my ($gff,$scaf2chr,$name)=@_;
open IN, "$gff" or die "$!";
open OUT, ">$name.gff" or die "$!";
#open UN, ">$name.unmaped.gff" or die "$!";
while (<IN>){
    chomp $_;
    if ($_ eq ""){
       print OUT "$_\n";
       next;
    }
    if ($_ =~/^##/){
       print OUT "$_\n";
       next;
    }
    my @unit=split("\t",$_);
    my $scaf=$unit[0];
    my $scafa=$unit[0]."a";
    my $scafb=$unit[0]."b";
    if (exists $scaf2chr->{$unit[0]}){
     my $array=$scaf2chr->{$unit[0]};
     foreach (@$array){
       my @temp=@unit;
       my $scafstrand=$_->[3];
       my $scafstart =$_->[1];
       my $scaflen   =$_->[2];
       my $scafref   =$_->[0];
       #print "$scafref\t$scaf\t$scaflen\t$scafstart\n";
       if ($scafstrand eq "+"){
          $temp[3]+=$scafstart-1;
          $temp[4]+=$scafstart-1;
       }else{
          my $end   =$scafstart+$scaflen-$temp[3];
          my $start =$scafstart+$scaflen-$temp[4];
          $temp[3]=$start;
          $temp[4]=$end;
          if ($temp[6] eq "-"){
                   $temp[6] = "+";
          }elsif($temp[6] eq "+"){
                   $temp[6] = "-";
          }
       }
       #print "$scafref\t$scaf\t$scaflen\t$scafstart\t$temp[3]\n";
       $temp[0]=$scafref;
       my $line=join("\t",@temp);
       print OUT "$line\n";
     }
    }elsif(exists $split{$unit[0]}){
     #print "@unit\n$split{$unit[0]}\n";

     my $tempscaf=$unit[0];
     my $array;
     if ($unit[3] < $split{$unit[0]}){
        $array=$scaf2chr->{$scafa};
        $unit[0]=$scafa;
     }else{
        $array=$scaf2chr->{$scafb};
        $unit[0]=$scafb;
        $unit[3]=$unit[3]-$bstart{$scafb} + 1;
        $unit[4]=$unit[4]-$bstart{$scafb} + 1;
        #print "bstart\t$bstart{$scafb}\n";
     }

     foreach (@$array){
       my @temp=@unit;
       my $scafstrand=$_->[3];
       my $scafstart =$_->[1];
       my $scaflen   =$_->[2];
       my $scafref   =$_->[0];
       #print "$tempscaf\t$scafref\t$scafstart\t$scaflen\t$scafstrand\n";
       #print "@temp\n";
       if ($scafstrand eq "+"){
          $temp[3]+=$scafstart-1;
          $temp[4]+=$scafstart-1;
       }else{
          my $end   =$scafstart+$scaflen-$temp[3];
          my $start =$scafstart+$scaflen-$temp[4];
          $temp[3]=$start; 
          $temp[4]=$end;
          if ($temp[6] eq "-"){
                   $temp[6] = "+";
          }elsif($temp[6] eq "+"){
                   $temp[6] = "-";
          }
       }
       #print "$scafref\t$scaf\t$scaflen\t$scafstart\t$temp[3]\n";
       $temp[0]=$scafref;
       my $line=join("\t",@temp);
       #print "$line\n";
       print OUT "$line\n";
     }
  
    }else{
       #print "$unit[0]\t$unmapscaf{$unit[0]}\n";
       if (exists $unmapscaf{$unit[0]}){
           $unit[0]=$unmapscaf{$unit[0]}->[0] if ($name=~/chr/);
           $unit[0]=$unmapscaf{$unit[0]}->[1] if ($name=~/super/);
       }else{
           print "Unknown Unmap scaffold\n";
       }
       #$unit[0]="chrUN";
       my $line=join("\t",@unit);
       print OUT "$line\n";
    }
}
close OUT;
#close UN;
close IN;
}

