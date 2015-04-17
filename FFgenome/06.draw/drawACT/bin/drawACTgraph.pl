#!/usr/bin/perl
=header
read two embl and 4ACT, draw comparision figure as well as annotation if present in embl
perl drawACT.pl -q query -t target -a 4ACT > log & 
Usage: perl drawACT.pl -q query.embl -t target.embl -a 4ACT -d rice2ff > log
=cut

use warnings;
use strict;
use SVG;
use Getopt::Long;

our %opt;
GetOptions(\%opt,"query:s","target:s","act:s","desc:s","qRNAseq:s","qmC:s","qH3K4:s","qH3K9:s","tRNAseq:s","tmC:s","tH3K4:s","tH3K9:s","help");
my $query=$opt{query};
my $target=$opt{target};
my $act=$opt{act};
my $desc=$opt{desc};

if (keys %opt  < 1 or $opt{help}){
    print "Usage: perl $0 -query ../data/Chr7_8901000_9404000.embl -target ../data/ghd7_FF169400_455000.embl -act ../data/Chr7_8901000_9404000VSghd7_FF169400_4550004ACT.4550004ACT -d graphtest -qRNAseq ../data/test.bed -qmC ../data/test.bed -qH3K4 ../data/test.bed -qH3K9 ../data/test.bed -tRNAseq ../data/test.bed -tmC ../data/test.bed -tH3K4 ../data/test.bed -tH3K9 ../data/test.bed\n";
    exit();
}


#################################### readin input files##################3
####parse embl file store position in reference of array
my ($qmin,$qmax,$qlength,$qgene,$qdnate,$qrnate)=parse($query);
my ($tmin,$tmax,$tlength,$tgene,$tdnate,$trnate)=parse($target);
####
my $identity=80;
my @links;
open IN, "$act" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_ eq "");
    my @array=split(" ",$_);
    if ($array[1] >= $identity){
        push (@links,"$array[2]\t$array[3]\t$array[5]\t$array[6]");
    }
}
close IN;

########################### drawing    ###########################
our $width=1200;
our $height=1000;
our $querypos=400;
our $targetpos=590;
my $rate1=$qlength/($width-200);
my $rate2=$tlength/($width-200);
my $svg=SVG->new(width=>$width,height=>$height);


#####title, description
my $title=$svg->text(
       x=>100,y=>50,
       style=>{
           fontsize=>'2'
       }
)->cdata($desc);

### draw query line;
my $qwidth=$qlength/$rate1;
my $qseq=$svg->rectangle(
            x=>100,y=>$querypos,
            width=>$qwidth,height=>10,  
            style=>{
              stroke=>'black'
            }
);
my $qanno=$svg->text(
           x=>90,y=>$querypos+10,
           style=>{'stroke-opacity'=>'0',
               fontsize=>'0.01','text-anchor'=>'end','font-weight'=>100,'stroke-width'=>0.01
           }
)->cdata("OB");


### draw target line;
my $twidth=$tlength/$rate2;
my $tseq=$svg->rectangle(
            x=>100,y=>$targetpos,
            width=>$twidth,height=>10,
            style=>{
              stroke=>'black'
            }
);
my $tanno=$svg->text(
           x=>90,y=>$targetpos+10,
           style=>{'stroke-opacity'=>'0',
               fontsize=>'0.01','text-anchor'=>'end','font-weight'=>100,'stroke-width'=>0.01
           }
)->cdata("OS");



my $blue="blue";
my $red="red";
my $black="black";
my $pos1="query";
my $pos2="target";
##draw query box
$svg=drawbox($qgene,$black,$svg,$rate1,$pos1,"gene");
$svg=drawbox($qdnate,$blue,$svg,$rate1,$pos1,"DNA");
$svg=drawbox($qrnate,$red,$svg,$rate1,$pos1,"RT");
##draw target box
$svg=drawbox($tgene,$black,$svg,$rate2,$pos2,"gene");
$svg=drawbox($tdnate,$blue,$svg,$rate2,$pos2,"DNA");
$svg=drawbox($trnate,$red,$svg,$rate2,$pos2,"RT");

##draw x axis for ref and qry
my $qaxisend=100+$qlength/$rate1;
my $qaxispos=$querypos-15;
print "$qmin\t$qmax\n$tmin\t$tmax\n";
$svg=drawxaxis($svg,"100",$qaxisend,$qaxispos,$qmin,$qmax,"200000",$rate1);
my $taxisend=100+$tlength/$rate2;
my $taxispos=$targetpos+20;
$svg=drawxaxis($svg,"100",$taxisend,$taxispos,$tmin,$tmax,"200000",$rate2);
#### draw links between query and target
foreach (@links){
    my @unit=split("\t",$_);
    my $qleft=$unit[0]/$rate1+100;
    my $qright=$unit[1]/$rate1+100;
    my $tleft=$unit[2]/$rate2+100;
    my $tright=$unit[3]/$rate2+100;
    my $color;
    if ($tright < $tleft){
         $color='red';
    }else{
         $color='#778899';
    }
   
    my $qheight=$querypos+11;
    my $theight=$targetpos-1;
    my $xv=[$qleft,$qright,$tright,$tleft];
    my $yv=[$qheight,$qheight,$theight,$theight];
    my $points =$svg->get_path(
                     x=>$xv,y=>$yv,
                     -type=>'polyline',
                     -closed=>'true'
              );
    my $tag=$svg->polyline(
                     %$points,
                     style=>{
                        fill=>$color
                        #fill=>'#FFDAB9' 
                     }
              );

}
######### draw graph of signal 
##query RNAseq
if (exists $opt{qRNAseq}){
  $svg=graph($svg,"100","300",$opt{qRNAseq},"50",$rate1,"RNAseq","Purple");
}
##query mC
if (exists $opt{qmC}){
  $svg=graph($svg,"100","250",$opt{qmC},"50",$rate1,"mC","brown");
}
##query H3K4
if (exists $opt{qH3K4}){
  $svg=graph($svg,"100","200",$opt{qH3K4},"50",$rate1,"H3K4","cyan");
}
##query H3K9
if (exists $opt{qH3K9}){
  $svg=graph($svg,"100","150",$opt{qH3K9},"50",$rate1,"H3K9","Orange");
}

##target RNAseq
if (exists $opt{tRNAseq}){
  $svg=graph($svg,"100","740",$opt{tRNAseq},"50",$rate2,"RNAseq","Purple");
}
##target mC
if (exists $opt{tmC}){
  $svg=graph($svg,"100","790",$opt{tmC},"50",$rate2,"mC","brown");
}
##target H3K4
if (exists $opt{tH3K4}){
  $svg=graph($svg,"100","840",$opt{tH3K4},"50",$rate2,"H3K4","cyan");
}
##target H3K9
if (exists $opt{tH3K9}){
  $svg=graph($svg,"100","890",$opt{tH3K9},"50",$rate2,"H3K9","Orange");
}
######### draw legend
$svg=legend($svg,100,"gene","black");
$svg=legend($svg,170,"LTR","red");
$svg=legend($svg,220,"DNA TE","blue");




my $outfile="$desc.svg";
writesvg($outfile,$svg);



################ sub functions#############################################
sub graph
{
my ($svg,$x,$y,$bed,$win,$rate,$title,$color)=@_;

my $axisx=$x-20;
my $axisy1=$y-40;
my $axisy2=$y;
my $rate4yaxis;
#($svg,$rate4yaxis)=drawyaxis($svg,$axisx,$axisy1,$axisy2,"0","40","10",$title);

my %hash;
open IN, "$bed" or die "$!";
while (<IN>){
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my $index=int $unit[1]/$win;
    my $pos  =$index*$win+1;
    if (exists $hash{$pos}){
       $hash{$pos}++;
    }else{
       $hash{$pos}=1;
    }
}
close IN;
my @signal=values %hash;
my $tops=topnum(@signal);
my $step=10;
if ($tops <= 40){
   $tops=40;
}elsif($tops > 40 and $tops < 60){
   $step=10;
}elsif($tops >= 60 and $tops <= 100){
   $step=20;
}else{
   $step=40;
}
($svg,$rate4yaxis)=drawyaxis($svg,$axisx,$axisy1,$axisy2,"0",$tops,$step,$title);

foreach (sort {$a <=> $b} keys %hash){
       #print "$_\t$hash{$_}\t$rate4yaxis\n";
    if ($hash{$_} > $tops){
       $hash{$_}=$tops;
       my $tleft=$_/$rate+100;
       my $tright=($_+$win)/$rate+100;
       my $bright=$tright;
       my $bleft=$tleft;
       my $theight=$y-$hash{$_}*$rate4yaxis+1;
       my $bheight=$y;
       #print "$tleft\t$tright\t$bright\t$bleft\t$theight\t$bheight\n";
       my $xv=[$tleft,$tright,$bright,$bleft];
       my $yv=[$theight,$theight,$bheight,$bheight];
       my $points =$svg->get_path(
                     x=>$xv,y=>$yv,
                     -type=>'polyline',
                     -closed=>'true'
              );
       my $tag=$svg->polyline(
                     %$points,
                     style=>{
                        #fill=>'red'
                        fill=>'#FFDAB9' 
                     }
              );
    }else{
       my $tleft=$_/$rate+100;
       my $tright=($_+$win)/$rate+100;
       my $bright=$tright;
       my $bleft=$tleft;
       my $theight=$y-$hash{$_}*$rate4yaxis+1;
       my $bheight=$y;
       #print "$tleft\t$tright\t$bright\t$bleft\t$theight\t$bheight\n";
       my $xv=[$tleft,$tright,$bright,$bleft];
       my $yv=[$theight,$theight,$bheight,$bheight];
       my $points =$svg->get_path(
                     x=>$xv,y=>$yv,
                     -type=>'polyline',
                     -closed=>'true'
              );
       my $tag=$svg->polyline(
                     %$points,
                     style=>{
                        fill=>$color
                        #fill=>'#FFDAB9' 
                     }
              );
    }
}
       
return $svg;
}

##########
sub topnum
{
####read in an array, the function will return a number, bellow which the numbers are consist of 90% of total number.
my (@array)=@_;
@array =sort {$a <=> $b} @array;
my $total=@array;
my $index=int ($total*0.95);
return $array[$index];
}
#########
sub drawxaxis
{
my ($svg,$x1,$x2,$y,$min,$max,$step,$rate)=@_;
#print "$x1\t$x2\t$y\n";
my $xaxis=$svg->line(
     x1=>$x1,y1=>$y,
     x2=>$x2,y2=>$y,
     style=>{stroke=>'black'}
);
my $topline=$svg->line(
     x1=>$x1,y1=>$y,
     x2=>$x1,y2=>$y-10,
     style=>{stroke=>'black'}
);
my $toptext=$svg->text(
     x=>$x1,y=>$y+13,
     style=>{'stroke-opacity'=>'0',
         fontsize=>'0.01','text-anchor'=>'start','font-weight'=>100,'stroke-width'=>0.01
     }
)->cdata("$min");
my $bottomline=$svg->line(
     x1=>$x2,y1=>$y,
     x2=>$x2,y2=>$y-10,
     style=>{stroke=>'black'}
);
my $bottomtext=$svg->text(
     x=>$x2,y=>$y+13,
     style=>{'stroke-opacity'=>'0',
         fontsize=>'0.01','text-anchor'=>'start','font-weight'=>100,'stroke-width'=>0.01
     }
)->cdata($max);

for(my $i=$min+$step;$i<$max;$i=$i+$step){
     my $tempx=$x1+($i-$min)/$rate;
     #print "$min\t$step\t$i\t$rate\n";
     #print "$tempx\t$y\n";
     my $line=$svg->line(
         x1=>$tempx,y1=>$y,
         x2=>$tempx,y2=>$y-5,
         style=>{stroke=>'black'}
     );
     my $text=$svg->text(
         x=>$tempx,y=>$y+13,
         style=>{'stroke-opacity'=>'0',
             fontsize=>'0.01','text-anchor'=>'end','font-weight'=>100,'stroke-width'=>0.01
         }
     )->cdata($i);
}
return $svg;
}


##########
sub drawyaxis
{
my ($svg,$x,$y1,$y2,$min,$max,$step,$title)=@_;
my $yaxis=$svg->line(
     x1=>$x,y1=>$y1,
     x2=>$x,y2=>$y2,
     style=>{stroke=>'black'}
);
my $topline=$svg->line(
     x1=>$x,y1=>$y1,
     x2=>$x+10,y2=>$y1,
     style=>{stroke=>'black'}
);
my $toptext=$svg->text(
     x=>$x,y=>$y1+5,
     style=>{'stroke-opacity'=>'0',
         fontsize=>'0.01','text-anchor'=>'end','font-weight'=>100,'stroke-width'=>0.01
     }
)->cdata($max);
my $bottomline=$svg->line(
     x1=>$x,y1=>$y2,
     x2=>$x+10,y2=>$y2,
     style=>{stroke=>'black'}
);
my $bottomtext=$svg->text(
     x=>$x,y=>$y2,
     style=>{'stroke-opacity'=>'0',
         fontsize=>'0.01','text-anchor'=>'end','font-weight'=>100,'stroke-width'=>0.01
     }
)->cdata("0");

my $rate=($y2-$y1+1)/($max-$min+1);
for(my $i=$min+$step;$i<$max;$i=$i+$step){
     my $tempy=$y2-$i*$rate;
     my $line=$svg->line(
         x1=>$x,y1=>$tempy,
         x2=>$x+5,y2=>$tempy,
         style=>{stroke=>'black'}
     );
}
my $textx=$x-50;
my $texty=$y1+($y2-$y1)/2;
my $anno=$svg->text(
     x=>$textx,y=>$texty,
     style=>{stroke=>'black',
         fontsize=>'0.4','text-anchor'=>'start','font-weight'=>100,'stroke-width'=>0.1
     }
)->cdata($title);
return ($svg,$rate);
}

#########
sub legend{
my ($svg,$x,$name,$color)=@_;
my $xtext=$x+20; 
 $svg->rectangle(
            x=>$x,y=>$height-50,
            width=>15,height=>10,
            style=>{
                fill=>$color
            }
 );
 $svg->text(
            x=>$xtext,y=>$height-40,
            style=>{stroke=>'black',
                   fontsize=>'1','text-anchor'=>'start','font-weight'=>100,'stroke-width'=>0.1
            }
 )->cdata($name);
 
return $svg;
}
########
################################### sub for write svg to file
sub writesvg {
my ($file,$svg)=@_;
#print "$file\n";
open OUT, ">$file" or die "can not open my file";
       print OUT $svg->xmlify;
close OUT;
       system "/home/jfchen/FFproject/tools/draw/svg2xxx_release/svg2xxx $file -t pdf";
}


################################ sub for draw feature box

sub drawbox{
my ($refarray,$color,$svg,$rate,$pos,$type)=@_;
my $y;
my $ytext;
if ($pos eq "query"){
   if ($type eq "gene"){
      $y=$querypos-35;
      $ytext=$y-10;
   }elsif($type eq "RT"){
      $y=$querypos-50;
      $ytext=$y-10;
   }elsif($type eq "DNA"){
      $y=$querypos-65;
      $ytext=$y-10;
   }
}else{
   if ($type eq "gene"){
      $y=$targetpos+35;
      $ytext=$y+10;
   }elsif($type eq "RT"){
      $y=$targetpos+50;
      $ytext=$y+10;
   }elsif($type eq "DNA"){
      $y=$targetpos+65;
      $ytext=$y+10;
   }
}
unless ($type eq "gene"){
   foreach (@$refarray){
      my @unit=split("\t",$_);
      my $length;
      my $start;
      my $name=$unit[0];
      if ($unit[1] > $unit[2]){
           $start=$unit[2]/$rate+100;
           $length=($unit[1]-$unit[2])/$rate;  ## we stored the end in $unit[1] if it was on the minus strand
      }else{
           $start=$unit[1]/$rate+100;
           $length=($unit[2]-$unit[1])/$rate;
      }
      $svg->rectangle(
            x=>$start,y=>$y,
            width=>$length,height=>8,
            style=>{
                fill=>$color
            }
      );
   }
}else{
   foreach (@$refarray){
      my @unit=split("\t",$_);
      my $name=shift @unit;
      my $start;
      my $length;
      #print "$unit[0]\n$unit[@unit-1]\n";
      my @firstexon=split(/\.\./,$unit[0]);
      my @lastexon=split(/\.\./,$unit[@unit-1]);
      #print "$firstexon[0]\t$firstexon[1]\n";
      my $genestart=$firstexon[0]/$rate+100;
      my $genelen=($lastexon[1]-$firstexon[0]+1)/$rate;
      $svg->rectangle(
           x=>$genestart,y=>$y+3,
           width=>$genelen,height=>2,
           style=>{
                fill=>$color
           } 
      );
      foreach(@unit){
         my ($exons,$exone)=split(/\.\./,$_);
         if ($exons > $exone){
             $start=$exone/$rate+100;
             $length=($exons-$exone)/$rate;
         }else{
             $start=$exons/$rate+100;
             $length=($exone-$exons)/$rate;
         }
         $svg->rectangle(
             x=>$start,y=>$y,
             width=>$length,height=>8,
             style=>{
                  fill=>$color
             }
         );
      }
   }

}
=pod
      if ($name =~/(15770)/ or $name =~/g(15640)/ or $name =~/(15680)/ or $name =~/(15880)/){
           $name=$1;
           $svg->text(
                x=>$start,y=>$ytext,
                style=>{stroke=>'black',
                   fontsize=>'2','text-anchor'=>'start','font-weight'=>100,'stroke-width'=>0.1
                } 
           )->cdata($name);  

      }
=cut
return $svg;
}

#######################sub for parse embl
sub parse {
my ($file)=@_;
print "parsing embl file: $file\n";
my @gene;
my @dnate;
my @rnate;
my $length;
my %hash;
my $counter;
my $smin;
my $smax;
if ($file=~/\_(\d+)\_(\d+)\.embl/){
   $smin=$1;
   $smax=$2;
}
open IN, "$file" or die "$!";
  while (<IN>){
      if ($_=~/FT\s+CDS\s+complement/){
             $_=~s/[A-Za-z]//g;
             $_=~s/\s+//g;
             $_=~s/\(|\)//g;
             $_=~/(\d+\..*\.\d+)/;
             #print "$_\n$1\t$2\n";
             #my $start=$2;
             #my $end=$1;
             my @position=split(",",$1);
             my $pos=join("\t",@position);
             my $loc=<IN>;
             $loc=~/FT\s+\/gene\=\"(\w+)\"/;
             my $name=$1;
             if (exists $hash{$name}){
                 next;
             }else{
                 $counter++;
                 $hash{$name}=$counter;
             }
             push (@gene,"$name\t$pos");
      }elsif($_=~/FT\s+CDS\s+join/){
             $_=~s/[A-Za-z]//g;
             $_=~s/\s+//g;
             $_=~s/\(|\)//g;
             $_=~/(\d+\..*\.\d+)/;
             #$_=~/(\d+)\..*\.(\d+)/;
             #print "$_\n$1\t$2\n";
             #my $start=$1;
             #my $end=$2;
             my @position=split(",",$1);
             my $pos=join("\t",@position);
             my $loc=<IN>;
             $loc=~/FT\s+\/gene\=\"(\w+)\"/;
             my $name=$1;
             if (exists $hash{$name}){
                 next;
             }else{
                 $counter++;
                 $hash{$name}=$counter;
             }
             #push (@gene,"$name\t$start\t$end");
             push (@gene,"$name\t$pos");
      }elsif($_=~/FT\s+repeat_region\s+\d+/){
             $_=~s/[A-Za-z]//g;
             $_=~s/\s+//g;
             $_=~s/\_//g;
             $_=~/(\d+)\..*\.(\d+)/;
             my $start=$1;
             my $end=$2;
             my $loc=<IN>;
             $loc=~/FT\s+\/note\=\"(\w+)\"/;
             my $name=$1;
             push (@dnate,"$name\t$start\t$end");
      }elsif($_=~/FT\s+repeat_region\s+complement/){
             $_=~s/[A-Za-z]//g;
             $_=~s/\s+//g;
             $_=~s/\_//g;
             $_=~s/\(|\)//g;
             $_=~/(\d+)\..*\.(\d+)/;
             #print "$_\n$1\t$2\n";
             my $start=$2;
             my $end=$1;
             my $loc=<IN>;
             $loc=~/FT\s+\/note\=\"(\w+)\"/;
             my $name=$1;
             push (@dnate,"$name\t$start\t$end");
      }elsif($_=~/FT\s+LTR\s+complement/){
             $_=~s/[A-Za-z]//g;
             $_=~s/\s+//g;
             $_=~s/\(|\)//g;
             $_=~/(\d+)\..*\.(\d+)/;
             #print "$_\n$1\t$2\n";
             my $start=$2;
             my $end=$1;
             my $loc=<IN>;
             $loc=~/FT\s+\/note\=\"(\w+)\"/;
             my $name=$1;
             push (@rnate,"$name\t$start\t$end");
      }elsif($_=~/FT\s+LTR\s+\d+/){
             $_=~s/[A-Za-z]//g;
             $_=~s/\s+//g;
             $_=~s/\(|\)//g;
             $_=~/(\d+)\..*\.(\d+)/;
             #print "$_\n$1\t$2\n";
             my $start=$1;
             my $end=$2; 
             my $loc=<IN>; 
             $loc=~/FT\s+\/note\=\"(\w+)\"/;
             my $name=$1; 
             push (@rnate,"$name\t$start\t$end");
      }elsif($_=~/SQ\s+Sequence/){
             $_=~s/[A-Za-z]//g;
             $_=~s/\s+//g;
             $_=~s/\;//g; 
             $length=$_;
             print "$length\n";  
      }

  }
close IN;
my $refgene=\@gene;
my $refdnate=\@dnate;
my $refrnate=\@rnate;
return ($smin,$smax,$length,$refgene,$refdnate,$refrnate);
}



