#!/usr/bin/perl

use strict;
use FindBin qw ($Bin);
use Getopt::Long;
use SVG;

my $output="$Bin/../output";
#################################Usage ##############################################
my %opt;
GetOptions(\%opt,"dir:s","heat","qrytegff:s","qrygenegff:s","qryroot:s","qryshoot:s","reftegff:s","refgenegff:s","refseq:s","queryseq:s","type:s","help:s");
my $help=<<USAGE;

perl drawFeature.pl --dir ../data/axt --qrytegff ../data/gff/all.gff.chr --qrygenegff ../data/gff/ass.scafSeq.gapfill3.glean.gff.chr --qryroot ../data/gff/root.soap.chr --qryshoot ../data/gff/shoot.soap.chr --reftegff ../data/gff/all.con.RepeatMasker.out.gff.chr --refgenegff ../data/gff/all.con.gff3.clear.chr --refseq ../data/seq/all.con --queryseq ../data/seq/ass.scafSeq.gapfill3 --type lastz > log 2> log2 &

--dir       : directory of filter.net.axt for all chromosomes from Lastz.
--qrytegff  : dir contain gff of repeat annotations for chr. 
--qrygenegff: dir contain gff of gene annotations for chr.
--qryroot   : dir contain root.bed of transcriptome for chr.
--qryshoot  : dir contain shoot.bed of transcriptome for chr.
--reftegff  : dir contain te annotion of reference for chr.
--refgenegff: dir contain gene annoation of reference for chr.
--refshoot  : dir contain shoot.bed of transcriptome for chr.
--refseq    : reference sequence of rice genome.
--queryseq  : queryseq of FF genome scaffold.
--type      : connet line type, lastz hit region or ortholog
--heat      : draw headmap
--help      : print this table of USAGE.

USAGE
if(defined $opt{help} or !-f $opt{refseq}){
        die  $help ;
}
####################################################################################
my $reflen =fastalen($opt{refseq});  ###ref of hash that store length of each seq
my $scaflen=fastalen($opt{queryseq});

my @axt=<$opt{dir}/*.axt.check>; ## read all axt file in dir to @axt
open OUT, ">chr.scaffold" or die "$!";
print OUT "Chr\tScaffold\tScafLen\tScaffold Start ON Chr\tScaffold End On Chr\tScaffold strand\n";
foreach (@axt){
print "$_\n";
############################# for each axt file ,read and store position informations####
my $scafchrlen=0;
my %scafstart;
my @position;
my %reverse;
my %hit;
my $ref;
open IN, "$_" or die "$!";
while (<IN>){
    next if ($_ eq "");
    next if ($_ =~/#/);
    if ($_ =~/^\d+/){
        my @unit=split(" ",$_);
        $ref =$unit[1];
        my $scaf=$unit[4];
        my $refstart =$unit[2];
        my $refend   =$unit[3];
        unless (exists $scafstart{$scaf}){
           $scafstart{$scaf}=$scafchrlen;
           $scafchrlen+=$scaflen->{$scaf};
           print OUT "$ref\t$scaf\t$scaflen->{$scaf}\t$scafstart{$scaf}\t$scafchrlen\t$unit[7]\n";
        }
        my $qrystart =$unit[5]+$scafstart{$scaf};
        my $qryend   =$unit[6]+$scafstart{$scaf};
        push (@position,"$refstart\t$refend\t$qrystart\t$qryend");
        $hit{$scaf}++;
        if ($unit[7] eq "-"){
           $reverse{$scaf}++; ### the direction is considered in the axt so we only need to read it here
        }
    }
}
close IN;
#########################################################################
#####################draw chromosome feature of two species##############
my $svg=SVG->new(width=>800,height=>600);
my $startw=100; 
my $endw  =750;
my $starth=100;
my $endh  =550;
my $scale=$reflen->{$ref}/($endw-$startw);

#### chromosome line
my $refstarth=$starth+200;
my $refwidth =$endw-$startw;
my $refline  =$svg->rectangle(
              x=>$startw,y=>$refstarth,
              width=>$refwidth,height=>10,
              style=>{
                  stroke=>'black'
              }
);
my $refnote  =$svg->text(
              x=>60, y=>$refstarth+8,
              style=>{
                   fontsize=>'2','text-anchor'=>'start','stroke-width'=>0.1
              }
)->cdata($ref);


my $qrystarth=$starth+280;
my $qrywidth =$scafchrlen/$scale;
my $qryline  =$svg->rectangle(
              x=>$startw,y=>$qrystarth, 
              width=>$qrywidth,height=>10
);
my $qrynote  =$svg->text( 
              x=>60, y=>$qrystarth+8,
              style=>{
                   fontsize=>'2','text-anchor'=>'start','stroke-width'=>0.1
              }
)->cdata("O.bra");
#################### draw centremere
my $refcentfile="$Bin/../data/cent/IRGSPCent.txt";
my $qrycentfile="$Bin/../data/cent/OBaCent.txt";
my $refcent=cent($refcentfile);
my $qrycent=cent($qrycentfile); #### draw qry centremere in draw scaffold line foreach loop
if (exists $refcent->{$ref}){
        my $centstart=shift @{$refcent->{$ref}};
        my $centend  =shift @{$refcent->{$ref}};
        my $centlen  =($centend-$centstart)/$scale;
        my $centx    =100+$centstart/$scale;
        my $cent=$svg->rectangle(
                  x=>$centx,y=>$refstarth-5,
                  width=>$centlen,height=>5,
                  style=>{
                      stroke=>'red'  
                  }
        );  
}

#############################################
############################draw scale bar for figure
$svg=barline($svg,$reflen->{$ref},$scale,"180");

sub barline {
my ($svg,$len,$scale,$y)=@_;
my $end=100+$len/$scale;
my $line=$svg->line(
         x1=>100,y1=>$y,	
         x2=>$end,y2=>$y,
         style=>{
               stroke=>'black'
         }
);
my $bin=1000000;
my $num=int ($len/$bin);
my $l1=$svg->line(
         x1=>100,y1=>$y,	
         x2=>100,y2=>$y-4,
         style=>{
               stroke=>'black'
         }
);
my $t1=$svg->text(
         x=>100,y=>$y+15,
         style=>{
             fontsize=>'0.7','text-anchor'=>'start','stroke-width'=>0.1
         }	
)->cdata("0");
my $l2=$svg->line(
         x1=>$end,y1=>$y,	
         x2=>$end,y2=>$y-4,
         style=>{
               stroke=>'black'
         }
);
my $noteend=sprintf ("%.1f",$len/1000000);
my $notee=$noteend." Mb";
my $t2=$svg->text(
         x=>$end,y=>$y+15,
         style=>{
               fontsize=>'0.7','text-anchor'=>'start','stroke-width'=>0.1
         }	
)->cdata($notee);
for (my $i=1;$i<=$num;$i++) {
         my $x=$i*$bin/$scale+100;
	 my $l=$svg->line(
		 x1=>$x,y1=>$y,
		 x2=>$x,y2=>$y-2,
                 style=>{
                     stroke=>'black'
                 } 
	 );
}
for (my $i=5;$i<=$num;$i+=5){
         my $x=$i*$bin/$scale+100;
	 my $t=$svg->text(
		 x=>$x,y=>$y+15,
                 style=>{
                        fontsize=>'0.7','text-anchor'=>'end','stroke-width'=>0.1
                 }
	 )->cdata($i);
}
return $svg;
}


############################draw scaffold line
foreach (sort {$scafstart{$a} <=> $scafstart{$b}} keys %scafstart) {
         my $scafsh=$qrystarth+20;
	 my $scafwh=$scaflen->{$_}/$scale;
	 my $scafsw=100+$scafstart{$_}/$scale;
         my $scaf=$_;
	 $_=~s/[a-z]//g;
	 my $scafline=$svg->rectangle(
		    x=>$scafsw, y=>$scafsh, 
	            width=>$scafwh,height=>10,
		    style=>{
	              'fill-opacity'=>0,
		       stroke=>'black'
	            }
	 );
         ##########draw centermere line for scaffold 
         if (exists $qrycent->{$scaf}){
                my $qrycent3=pop @{$qrycent->{$scaf}};
                my @qrycentstart=split(",",$qrycent3);
                foreach (@qrycentstart){
                     my $qrycw;
                     if (exists $reverse{$scaf}){
                         $qrycw=100+($scafstart{$scaf}+$scaflen->{$scaf}-$_)/$scale;
                     }else{
                         $qrycw=100+($scafstart{$scaf}+$_)/$scale;
                     }
                     my $qrycl=$svg->line(
                         x1=>$qrycw,y1=>$scafsh,
                         x2=>$qrycw,y2=>$scafsh+10,
                         style=>{
                              stroke=>'red'
                         }
                     );
                }     
         }
	 #my $scafnote=$svg->text(
	#	    x=>$scafsw+5,y=>$scafsh+8,
	#	    style=>{
        #                fontsize=>'2','text-anchor'=>'start','stroke-width'=>0.1
        #            }
	# )->cdata($_);
}



####Draw lastz or ortholog connetion line
my $up=$refstarth+10;
my $down=$qrystarth;
if ($opt{type} eq "lastz"){
   $svg=drawlastz($svg,$up,$down,$scale,\@position);
}
####
#if ($opt{heat}){
#####################draw feature density bin###############################
my $win=1000000; ##1Mb
my $step=200000; ##200kb
####ref gene density
my $binsetr1=280; ## position to set the bin draw for gene density
my $minr1=0;  ## 0% base pair in every 1Mb bin is CDS, can be adjusted by view data
my $maxr1=20; ## 20% base pair in every 1Mb bin is CDS, can be adjusted by view data
$svg=drawrefgene($svg,$opt{refgenegff},$ref,$reflen->{$ref},$scale,$win,$step,$minr1,$maxr1,$binsetr1);
####ref TE density
#DNA
my $binsetr2=260;
my $minr2=0;
my $maxr2=25;
$svg=drawrefte($svg,$opt{reftegff},$ref,$reflen->{$ref},$scale,$win,$step,$minr2,$maxr2,$binsetr2,"DNA");
#LTR
my $binsetr3=240;
my $minr3=0;
my $maxr3=40;
$svg=drawrefte($svg,$opt{reftegff},$ref,$reflen->{$ref},$scale,$win,$step,$minr3,$maxr3,$binsetr3,"LTR");
#Other
my $binsetr4=220;
my $minr4=0;
my $maxr4=10;
$svg=drawrefte($svg,$opt{reftegff},$ref,$reflen->{$ref},$scale,$win,$step,$minr4,$maxr4,$binsetr4,"Other");
### Rice transcriptome Shoot
my $binsetr5=200;
my $minr5=1000;
my $maxr5=6000;
$svg=drawrefsoap($svg,$opt{refshoot},$ref,$reflen->{$ref},$scale,$win,$step,$minr5,$maxr5,$binsetr5,"Shoot");
##############################################################################################################
####qry gene density
my $binsetr1a=430; ## position to set the bin draw for gene density 
my $minr1a=0;  ## 0% base pair in every 1Mb bin is CDS, can be adjusted by view data 
my $maxr1a=20; ## 20% base pair in every 1Mb bin is CDS, can be adjusted by view data 
$svg=drawrefgene($svg,$opt{qrygenegff},$ref,$scafchrlen,$scale,$win,$step,$minr1a,$maxr1a,$binsetr1a); 
####qry TE density 
#DNA 
my $binsetr2a=450; 
my $minr2a=0; 
my $maxr2a=25; 
$svg=drawrefte($svg,$opt{qrytegff},$ref,$scafchrlen,$scale,$win,$step,$minr2a,$maxr2a,$binsetr2a,"DNA"); 
#LTR 
my $binsetr3a=470; 
my $minr3a=0; 
my $maxr3a=40; 
$svg=drawrefte($svg,$opt{qrytegff},$ref,$scafchrlen,$scale,$win,$step,$minr3a,$maxr3a,$binsetr3a,"LTR"); 
#Other 
my $binsetr4a=490; 
my $minr4a=0; 
my $maxr4a=10; 
$svg=drawrefte($svg,$opt{qrytegff},$ref,$scafchrlen,$scale,$win,$step,$minr4a,$maxr4a,$binsetr4a,"Other"); 
####qry transcriptome
#Root
my $binsetr5a=510;
my $minr5a=1000;
my $maxr5a=6000;
$svg=drawrefsoap($svg,$opt{qryroot},$ref,$scafchrlen,$scale,$win,$step,$minr5a,$maxr5a,$binsetr5a,"Root");
#Shoot
my $binsetr6a=530;
my $minr6a=1000; 
my $maxr6a=6000; 
$svg=drawrefsoap($svg,$opt{qryshoot},$ref,$scafchrlen,$scale,$win,$step,$minr6a,$maxr6a,$binsetr6a,"Shoot");

########################################################################
#} ### end of if heat 

###Read Cent File
### return a hash with chr or scaffold as key and last three colume array as values
sub cent {
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while (<IN>){
   chomp $_;
   next if ($_ eq "");
   my @unit=split("\t",$_);
   my $head =shift @unit;
   $hash{$head}=\@unit;
}
close IN;
return \%hash;
}
########################################################################

####write file
my $svgfile="$ref.svg";
writesvg($svgfile,$svg);
system "mv *.png *.svg *.pdf $output";
#########################################################################
} ### end of foreach of axt file
close OUT;
#####################################################################

sub drawrefte {
my ($svg,$gff,$chr,$chrlen,$scale,$win,$step,$min,$max,$binset,$te)=@_;
chomp $gff;
my %start;
open IN, "$gff/$chr" or die "$!";
print "$gff\n";
while (<IN>){
     chomp $_;
     next if ($_ =~/^#/);
     next if ($_ eq "" );
     next if ($_ =~/^chrUN/);
     my @unit=split("\t", $_);
     $chr=~/(\d+)/;
     $chr=$1;
     $unit[0]=~/(\d+)/;
     $unit[0]=$1;
     #print "$chr\t$unit[0]\n";
     if ($unit[0] == $chr and $unit[2] eq "Transposon" or $unit[2] eq "TEprotein" or $unit[2] eq "TandemRepeat"){
        if ($te eq "DNA"){
             if ($unit[8] =~/Class=DNA/){
                 my $len=$unit[4]-$unit[3]+1;
                 $start{$unit[3]}=$len;
             }
        }elsif($te eq "LTR"){
             if ($unit[8] =~/Class=LTR/){
                 my $len=$unit[4]-$unit[3]+1;
                 $start{$unit[3]}=$len;
             }
        }elsif($te eq "Other"){
             unless ($unit[8] =~/Class=LTR/ or $unit[8] =~/Class=DNA/){
                 my $len=$unit[4]-$unit[3]+1;
                 $start{$unit[3]}=$len;
             } 

        }else{
             print "Need DNA,LTR,other Type specified\n";
             exit;
        }
     }

}
close IN;
my $basename =`basename $gff`;
chomp $basename;
open OUT1, ">>$basename.$te.density" or die "$!";
my @color=colortab();
my $binstart=100;
my $binend=100;
my $run=int(($chrlen-$win)/$step)+1;
for (my $i=0;$i<=$run;$i++){
    my $start=$i*$step;
    my $end  =$start+$win;
    my $totallen=1;
    #my @length;
    foreach (sort { $a <=> $b } keys %start){
        if ($_ >= $start and $_ <= $end ){
           $totallen+=$start{$_};
           if ($start == 0 and $end == 1000000){
           #print "Start $start\t End $end\n";
           print "$_\t$start{$_}\t$totallen\n";
           }
        }
    }
    my $density=100*$totallen/$win;
    print OUT1 "$density\n"; ## print out density data, $min and $max can be found using the data figure
    #my $min=10;  ## 0% base pair in every 1Mb bin is CDS, can be adjusted by view data
    #my $max=70; ## 20% base pair in every 1Mb bin is CDS, can be adjusted by view data
    my $interval=($max-$min)/26;
    my $index;
    my $number=0;
    for (my $j=$min;$j<=$max;$j+=$interval){
        my $limit=$j+$interval; 
        if ($density >= $j and $density <= $limit){
           $index=$number;
        }
        $number++;   
    }
    if ($density < $min){
       $index=0;
    }elsif($density > $max){
       $index=25; 
    }
    if ($index > 25){
       $index=25;
    }
    #print "$index\n";
    my $color=$color[$index]; ## red color stand for high gene density
    my ($red,$yellow,$blue)=split("\t",$color);
    $binstart=$binend;
    if ($i == $run){ 
      $binend=$chrlen/$scale+100;
    }else{
      $binend=($win/2+$start+$step/2)/$scale+100;
    }
    my $binwidth=$binend-$binstart;
    #print "$binstart\t$density\n";
    my $rec=$svg->rectangle(
                x=>$binstart,y=>$binset,           
                width=>$binwidth,height=>10,
                style=>{
                     fill=>"rgb($red,$yellow,$blue)"
                }
    );
}
    my $text=$svg->text(
                x=>60,y=>$binset+8,
                style=>{
                   fontsize=>'2','text-anchor'=>'start','stroke-width'=>0.1
                }
    )->cdata($te);
close OUT1;
return $svg;
}
#####################################################################
###########give sub a array return mean min max######################
sub sumarray {
my (@array)=@_;
my @sort=sort {$a <=> $b} @array;
my $sum;
my $mean;
foreach (@array){
   $sum+=$_;
}
$mean=$sum/scalar @_;
my $min=shift @sort;
my $max=pop @sort;   
return ($mean,$min,$max);
}
#####################################################################
#####################sub for draw refgene density####################

sub drawrefgene {
my ($svg,$gff,$chr,$chrlen,$scale,$win,$step,$min,$max,$binset)=@_;
my %start;
chomp $gff;
print "GFF $gff";
open IN, "$gff/$chr" or die "$!";
#print "$gff\n";
while (<IN>){
     chomp $_;
     next if ($_ =~/^#/);
     next if ($_ eq "" );
     next if ($_ =~/^chrUN/);
     my @unit=split("\t", $_);
     $chr=~/(\d+)/;
     $chr=$1;
     $unit[0]=~/(\d+)/;
     $unit[0]=$1;
     #print "$chr\t$unit[0]\n";
     if ($unit[0] == $chr and $unit[2] eq "CDS"){
        #print "$unit[3]\n";
        my $len=$unit[4]-$unit[3]+1;
        $start{$unit[3]}=$len;
     }
}
close IN;
my $basename =`basename $gff`;
chomp $basename;
open OUT1, ">>$basename.density" or die "$!";
my @color=colortab();
my $binstart=100;
my $binend=100;
my $run=int(($chrlen-$win)/$step)+1;
for (my $i=0;$i<=$run;$i++){
    my $start=$i*$step;
    my $end  =$start+$win;
    my $totallen=1;
    foreach (sort {$a <=> $b} keys %start){
        if ($_ >= $start and $_ <= $end ){
           $totallen+=$start{$_};
           #if ($start==0 and $end == 1000000){
           #print "Start $start\t End $end\n";
           #print "$_\t$start{$_}\t$totallen\n";
           #}
        }
    }
    my $density=100*$totallen/$win;
    print OUT1 "$density\n"; ## print out density data, $min and $max can be found using the data figure
    #my $min=0;  ## 0% base pair in every 1Mb bin is CDS, can be adjusted by view data
    #my $max=20; ## 20% base pair in every 1Mb bin is CDS, can be adjusted by view data
    my $interval=($max-$min)/26;
    my $index;
    my $number=0;
    for (my $j=$min;$j<=$max;$j+=$interval){
        my $limit=$j+$interval; 
        if ($density >= $j and $density <= $limit){
           $index=$number;
        }
        $number++;   
    }
    if ($density < $min){
       $index=0;
    }elsif($density > $max){
       $index=25; 
    }
    if ($index > 25){
       $index=25;
    }
    #print "$index\n";
    my $color=$color[$index]; ## red color stand for high gene density
    my ($red,$yellow,$blue)=split("\t",$color);
    $binstart=$binend;
    if ($i == $run){ 
      $binend=$chrlen/$scale+100;
    }else{
      $binend=($win/2+$start+$step/2)/$scale+100;
    }
    my $binwidth=$binend-$binstart;
    #print "$binstart\t$density\n";
    my $rec=$svg->rectangle(
                x=>$binstart,y=>$binset,           
                width=>$binwidth,height=>10,
                style=>{
                     fill=>"rgb($red,$yellow,$blue)"
                }
    );
}
    my $text=$svg->text(
                x=>60,y=>$binset+8,
                style=>{
                   fontsize=>'2','text-anchor'=>'start','stroke-width'=>0.1
                }
    )->cdata("Gene");
close OUT1;
return $svg;
}
####################draw transcriptome##########################################
sub drawrefsoap {
my ($svg,$gff,$chr,$chrlen,$scale,$win,$step,$min,$max,$binset,$tissue)=@_;
my %start;
chomp $gff;
my $readn;
my @all=<$gff/*>;
foreach (@all){
   my @read=split(" ",`wc -l $_`);
   $readn+=$read[0];
}
open IN, "$gff/$chr" or die "$!";
print "$gff\n";
while (<IN>){
     chomp $_;
     next if ($_ eq "" );
     next if ($_ =~/^chrUN/);
     my @unit=split("\t", $_);
     $chr=~/(\d+)/;
     $chr=$1;
     $unit[0]=~/(\d+)/;
     $unit[0]=$1;
     #print "$chr\t$unit[0]\n";
     if ($unit[0] == $chr){
        if (exists $start{$unit[1]}){  ## in BED format, $unit[1] is start position
           $start{$unit[1]}+=1;
        }else{
           $start{$unit[1]}=1;
        }
     }
}
close IN;
my $basename =`basename $gff`;
chomp $basename;
open OUT1, ">>$basename.density" or die "$!";
my @color=colortab();
my $binstart=100;
my $binend=100;
my $run=int(($chrlen-$win)/$step)+1;
for (my $i=0;$i<=$run;$i++){
    my $start=$i*$step;
    my $end  =$start+$win;
    my $totalread=1;
    foreach (sort {$a <=> $b} keys %start){
        if ($_ >= $start and $_ <= $end ){
           $totalread+=$start{$_};
        }
    }
    my $density=$totalread*1000000/$readn; ###similar as RPKM in transcriptome analysis,total read per Mb per Million read
    print OUT1 "$density\n"; ## print out density data, $min and $max can be found using the data figure
    #my $min=0;  ## 0% base pair in every 1Mb bin is CDS, can be adjusted by view data
    #my $max=20; ## 20% base pair in every 1Mb bin is CDS, can be adjusted by view data
    my $interval=($max-$min)/26;
    my $index;
    my $number=0;
    for (my $j=$min;$j<=$max;$j+=$interval){
        my $limit=$j+$interval; 
        if ($density >= $j and $density <= $limit){
           $index=$number;
        }
        $number++;   
    }
    if ($density < $min){
       $index=0;
    }elsif($density > $max){
       $index=25; 
    }
    if ($index > 25){
       $index=25;
    }
    #print "$index\n";
    my $color=$color[$index]; ## red color stand for high gene density
    my ($red,$yellow,$blue)=split("\t",$color);
    $binstart=$binend;
    if ($i == $run){ 
      $binend=$chrlen/$scale+100;
    }else{
      $binend=($win/2+$start+$step/2)/$scale+100;
    }
    my $binwidth=$binend-$binstart;
    #print "$binstart\t$density\n";
    my $rec=$svg->rectangle(
                x=>$binstart,y=>$binset,           
                width=>$binwidth,height=>10,
                style=>{
                     fill=>"rgb($red,$yellow,$blue)"
                }
    );
}
    my $text=$svg->text(
                x=>60,y=>$binset+8,
                style=>{
                   fontsize=>'2','text-anchor'=>'start','stroke-width'=>0.1
                }
    )->cdata($tissue);

close OUT1;
return $svg;
}
#######################################################################
sub colortab{
my @color;
for(my $red=0;$red<=255;$red+=10){
       my $blue;
       my $yellow;
       if ($red <= 125){
           $blue=255-$red;
           $yellow=$red;
       }elsif($red >125){
           $blue=0;
           $yellow=255-$red;
       }
       push (@color,"$red\t$yellow\t$blue");
}
return @color;
}
#######################################################################
#####################sub for draw lastz hit regions####################
sub drawlastz{
my ($svg,$up,$down,$scale,$position)=@_;
foreach (@$position){
    my @unit=split("\t",$_);
    my $upleft=$unit[0]/$scale+100;
    my $upright=$unit[1]/$scale+100;
    my $downleft=$unit[2]/$scale+100;
    my $downright=$unit[3]/$scale+100;
    my $xv=[$upleft,$upright,$downright,$downleft];
    my $yv=[$up,$up,$down,$down];
    my $points=$svg->get_path(
               x=>$xv,y=>$yv,
               -type=>'polyline',
               -closed=>'true'
    );
    my $tag=$svg->polyline(
               %$points,
               style=>{
                    fill=>'#FF6347'
               }
    );
}
return $svg;
}
######################################################################

##################### sub to get len hash of fasta file ###############
sub fastalen {
my ($file)=@_;
my %len;
$/=">";
print "$file\n";
open IN, "$file" or die "$!";
      while (<IN>){
          chomp $_;
          next if (length $_ < 2);
          my @unit=split("\n",$_);
          my $head=shift @unit;
          my @array=split(" ",$head);
          my $name=$array[0];
          my $seq=join("",@unit);
          $seq=~s/\s//g;
          $seq=~s/\r//g;
          $seq=~s/\>//g;
          $len{$name}=length $seq;
      }
close IN;
$/="\n";
return \%len;
}
######################################################################

########################write svg######################################
sub writesvg {
my ($file,$svg)=@_;
#print "$file\n";
open OUT1, ">$file" or die "can not open my file";
       print OUT1 $svg->xmlify;
close OUT1;
      system "/share/raid12/chenjinfeng/tools/draw/svg2xxx_release/svg2xxx $file -m 500 -t png";
      #system "/share/raid12/chenjinfeng/tools/draw/svg2xxx_release/svg2xxx $file -t pdf";
}
########################################################################
