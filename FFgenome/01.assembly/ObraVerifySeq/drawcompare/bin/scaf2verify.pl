### parse blastout and map scaffold to reference, draw a figure to view the homolog regions.
### perl scaf2ref.pl --query query --target target --blastout blastout --identity identity --lengthcut hitlencut
### perl scaf2ref.pl --query chr04.scarSeq --target chr04.txt --blastout ../data/scaf2chr04blastm8 --identity 99 --lenthcut 500
### version1, the target sequence is only one and queries may be hidden if they are overlaped.
### 20091210, by chen jinfeng
### version2, can do multi target and the queries are order on the chromosome
### 20091211, by chen jinfeng
### Modified from scaf2refV2.pl flip query if they are reverse hitted on ref. 
### 20091216, by chen jinfeng

use strict;
use warnings;
use Getopt::Long;
use SVG;
### get in the options 
my ($query,$target,$blastm8,$identity,$lencut,$help);
GetOptions(                ## --query or other could be writen in short, such as --que or --q or -q
    "query:s"=>\$query,    ## blast query                                ## s stand for string and i stand for number
    "target:s"=>\$target,  ## blast target database
    "blastout:s"=>\$blastm8,  ## blastout file -m 8 format
    "identity:i"=>\$identity, ## cut off identify of hit
    "lengthcut:i"=>\$lencut,  ## cut off  of hit length
    "help:s"=> \$help  ###
);
print "$query\t$target\t$blastm8\t$identity\t$lencut\n";
#die `pod2text $0` if (@ARGV == 0 || $Help);
die "Usage: perl scaf2ref.pl --query chr04.scarSeq --target chr04.txt --blastout ../data/scaf2chr04blastm8 --identity 99 --lengthcut 500\n" unless (-f $blastm8 and -f $target);
die "Usage: perl scaf2ref.pl --query chr04.scarSeq --target chr04.txt --blastout ../data/scaf2chr04blastm8 --identity 99 --lengthcut 500\n" if ($help);
#print "$blastm8\n";

### read scaffold sequence and store scaffold length in hash as to draw line
my %scaflen;
$/=">";
open IN, "$query" or die "can not open my scaffold file";;
      while (<IN>){
          my @unit=split("\n",$_);
          my $head=shift @unit;
          my $seq=join("",@unit);
          $scaflen{$head}=length $seq;
      }
close IN;
$/="\n";
#open LOG, ">runing.txt";


#### draw chromosome and return figure in reference not in file
my ($refsvg,$refrate)=drawref($target);
#print LOG "drawref done!\n";
###
### parse blastout and store position information in hash
$refsvg=map2ref($blastm8,$identity,$lencut,$refsvg,$refrate,\%scaflen);
#print LOG "map2ref done!\n";

### write to svg file
foreach (keys %$refsvg){
     print "$_\n"; 
     my $svgfile=$_.".svg";
     my $svgref=$refsvg->{$_};
     writesvg($svgfile,$svgref);
}

#print LOG "writefile done!\n";
#close LOG;
### map scaf to ref

sub map2ref {
my ($blastm8,$identity,$lencut,$refsvg,$refrate,$scaflen)=@_;
my %scaf2ref; ## stroe hit chromosome for scaffold, value is ref of hash chr->@position
my %scaf2num; ## store number that a chromosome that contain hit for each scaffold
my %scaf2max; ## store anchor chromosome for each scaffold
my $counter;
open BLAST, "$blastm8" or die "can not open my blastout";
     while (<BLAST>){
           chomp $_;
           my @unit=split("\t",$_);
           if ($unit[2] < $identity){next};
           my $hitlength=$unit[7]-$unit[6];
           #print "$hitlength\t$lencut\n";
           if ($hitlength < $lencut){next};

           my @position;
           my %chr2num;
           my %chr2pos;
           if (exists $scaf2ref{$unit[0]}){
                if (exists $scaf2ref{$unit[0]}->{$unit[1]}){
                        $scaf2num{$unit[0]}->{$unit[1]}+=1;
                        my $refarray=$scaf2ref{$unit[0]}->{$unit[1]};
                        #print "Exist: @$refarray";
                        push (@$refarray,"$unit[6]\t$unit[7]\t$unit[8]\t$unit[9]");  
                        $scaf2ref{$unit[0]}->{$unit[1]}=$refarray;
                }else{
                        $scaf2num{$unit[0]}->{$unit[1]}+=1;
                        push (@position,"$unit[6]\t$unit[7]\t$unit[8]\t$unit[9]");
                        my $refhash=$scaf2ref{$unit[0]};
                        $refhash->{$unit[1]}=\@position;
                        $scaf2ref{$unit[0]}=$refhash;
                } 
           }else{
                $counter++;
                print "$counter\n"; 
                push (@position,"$unit[6]\t$unit[7]\t$unit[8]\t$unit[9]");
                $chr2pos{$unit[1]}=\@position; ## store position of query and hit in array and give them to hash chr->@
                $scaf2ref{$unit[0]}=\%chr2pos; ## chr->@ as value and scaffold as key, need to update when new chr come 
                $chr2num{$unit[1]}+=1;
                $scaf2num{$unit[0]}=\%chr2num; ## chr->num as value and scaffold as key, need to update when new chr
           }
           #ckarray(\%scaf2ref);
          # print "Run\n"; 
     }
close BLAST;
print "$counter\n";
#ckarray(\%scaf2ref);
### check @position array
sub ckarray {
my ($scaf2ref)=@_;
foreach (keys %$scaf2ref){
        print "$_\n";
        my $ref=$scaf2ref->{$_};
        foreach (keys %$ref){
             print "$_\n";
             my $refarray=$ref->{$_};
             foreach (@$refarray){
                     print "$_\n";
             };
        }
}
}
###

my %scafonchr; ## store scaf->midpoint of scaf on psudo chr(query scaffold added chromosome)
my %scafchrlen; ## store chr->scaf added length;
#open CK, ">check.txt";
foreach (sort keys %scaf2ref){ ## in each scaffold, find the chromosome that have maximum number of hit for this scaf
          my $scaf=$_;
          my $refnum=$scaf2num{$scaf};
          my ($maxchr,$maxnum)=maxkey($refnum);
          $scaf2max{$scaf}=$maxchr;
          #print CK "$scaf\t$maxchr\t$maxnum\n";
}
#close CK;

foreach (values %scaf2max){  ### in each hit chr, we sort the scaf on the chr and save midpoint of scaf on chr
      my $chr=$_;
      my @scaf4chr; ### array store the scaffold anchored on this chr
      my %scafmid; ### hash store scaf mean target start point on chr;
      my $endonchr=0;
      my $midonchr; 
      foreach (keys %scaf2max){
             my $scaf=$_;
                
             if ($scaf2max{$scaf} eq $chr){
                 
                 push (@scaf4chr,$scaf);
                 my $refarray=$scaf2ref{$scaf}->{$chr};
                 my @hitstart;
                 my $addlen=0;
                 foreach (@$refarray){
                       my @pos=split("\t",$_);
                       $addlen+=$pos[1]-$pos[0]; ## record overall hit length of a scaffold to ref
                       push (@hitstart,$pos[2]);
                 }
                 if ($addlen < 1000){
                     delete $scaf2ref{$scaf};
                     next;
                 } 
                 my $mid=mean(@hitstart); 
                 $scafmid{$scaf}=$mid;         
             }
      }
      my @scaffold=sort {$scafmid{$a} <=> $scafmid{$b}} keys %scafmid; ## sorted scaffold for a chr
      foreach (@scaffold){
            my $length=$scaflen->{$_};
           # print "Check scaflen: $length\t$_\n";
            $midonchr=$length/2+$endonchr;
            $endonchr+=$length;
            $scafonchr{$_}=$midonchr;
      }
      $scafchrlen{$chr}=$endonchr;
}

foreach (keys %scaf2ref){
        my $scaf=$_;
        my $refscaf=$scaf2ref{$scaf}; ### position array ref
        my $scafchr=$scaf2max{$scaf}; ### anchored chromosome
        my $scaflen=$scaflen->{$scaf};
        $refsvg=drawscaf(\%scafchrlen,\%scafonchr,$scaflen,$scafchr,$refsvg,$scaf,$refrate,$refscaf);  
}

return $refsvg;
} ## sub map2ref end

#my ($refsvg,$scalar,$head,$file,$starth,$endh,$startw,$endw)=drawref($target);
### draw scaffold
sub drawscaf{
my ($scafchrlen,$scafonchr,$scaflen,$scafchr,$refsvg,$scaf,$refrate,$refscaf)=@_;
my $svg=$refsvg->{$scafchr};
my $refarray=$refscaf->{$scafchr};
my (@qstart,@qend,@tstart,@tend);
my $reverse=0;
#my $addlen=0;
my $scafrate=$scafchrlen->{$scafchr}/800;
#my $scafrate=40000000/800; ### change to small number if figure is two samll
my $chrrate=$refrate->{$scafchr};
foreach (@$refarray){
     # print "$scaf\nRefarray\t $_\n";
      my @pos=split("\t",$_);
      if ($pos[2] > $pos[3]){ 
          $reverse+=1;  ## record the hit direction if it was reverse hit
      }
      #$addlen+=$pos[1]-$pos[0]; ## record overall hit length of a scaffold to ref
      push (@qstart,$pos[0]);
      push (@qend,$pos[1]); 
      push (@tstart,$pos[2]);
      push (@tend,$pos[3]);
}
      #if ($addlen < 5000){ ## if overall hit length smaller than 5kb, exit sub function;
      #     return $refsvg;
      #}
      my $refqs=\@qstart;
      my $refqe=\@qend;
      my $refts=\@tstart;
      my $refte=\@tend;
      #print "array: @qstart\n@qend\n@tstart\n@tend\n";
      my $mean=$scafonchr->{$scaf};
      #print "$scaf\t$mean\n";
      my $start=($mean-$scaflen/2)/$scafrate+100;
      my $end  =($mean+$scaflen/2)/$scafrate+100;
      #print "scaf $scaf chr $scafchr start $start end $end mean $mean scaflen $scaflen\n";
      my $scafheight=$end-$start;
      #print "scafheight $scafheight\n";
      my $rec=$svg->rectangle(
                   x=>400,y=>$start,
                   width=>5,height=>$scafheight,
                   style=>{
                          stroke=>'black',
                          fill  =>'red'
                   }
      );
      my $textw=400+15;
      my $texth=$mean/$scafrate+100;
      my $text=$svg->text(
                   x=>$textw, y=>$texth,
                   style=>{stroke=>'black',
                   fontsize=>'7','text-anchor'=>'start','font-weight'=>100,'stroke-width'=>0.1
                   }
      )->cdata($scaf);
      for (my $i=0;$i<@$refqs;$i++){
              my ($scafstart,$scafend,$refstart,$refend);    
              if ($reverse > @$refqs/2){  ## if reverse hit beyond half we reverse the scaffold in the figure
                  $scafstart=($scaflen-$$refqs[$i])/$scafrate+$start;
                  $scafend  =($scaflen-$$refqe[$i])/$scafrate+$start;
                  $refstart =$$refts[$i]/$chrrate+100;
                  $refend   =$$refte[$i]/$chrrate+100;
              }else{
                  $scafstart=$$refqs[$i]/$scafrate+$start;
                  $scafend  =$$refqe[$i]/$scafrate+$start;
                  $refstart =$$refts[$i]/$chrrate+100;
                  $refend   =$$refte[$i]/$chrrate+100;
              }
            #  print "Refrate: $chrrate\n";
             # print "$scaf position: $scafstart\t$scafend\t$refstart\t$refend\n";
              my $left=50+10;
              my $color;
              my $xv=[$left,400,400,$left];
              my $yv=[$refstart,$scafstart,$scafend,$refend];
              if ($refend > $refstart){
                 $color="#FFDAB9";
              }else{
                 #$color="#778899";
                 $color="#FFDAB9";
              }
                                        
              ## filled polyline, the xv, yv should be in the order of clock.
              #my $xv=[$left,$endw,$leftw,$endw];
              #my $yv=[$refstart,$scafstart,$refend,$scafend];
              my $points =$svg->get_path(
                     x=>$xv,y=>$yv,
                     -type=>'polyline',
                     -closed=>'true'
              );
              my $tag=$svg->polyline(
                     %$points,
                     style=>{
                          fill=>$color
                     }
              );
=pod                
              my $line=$svg->line(
                      x1=>500,y1=>$scafstart,
                      x2=>$left,y2=>$refstart,
                      style=>{
                            stroke=>$color
                      }
              );
              my $line=$svg->line(
                      x1=>500,y1=>$scafend,
                      x2=>$left,y2=>$refend,
                      style=>{
                            stroke=>$color
                      }
              );
=cut
              
              
      }
$refsvg->{$scafchr}=$svg;
return $refsvg;
} ## sub draw scaf
#writesvg($file,$svg);


sub maxkey { ##give the sub a hash, it will return you a key for max value of hash ###
my ($hash)=@_;   
my $maxvalue=0;
my $maxkey;
foreach (keys %$hash){
      if ($hash->{$_} > $maxvalue){
          $maxvalue=$hash->{$_};
          $maxkey=$_;
      }
}
return ($maxkey,$maxvalue);
}

sub mean{
my (@array)=@_;
my $sum;
my $mean;
foreach (@array){
   $sum+=$_;
}
   $mean=$sum/scalar @_;
return $mean;
}

sub drawref { ## one sequence in a fasta file

my ($target)=@_;
my %svghash;
my %ratehash;
my %hash;
$/=">";
open IN, "$target" or die "can not open my reference file";   
     while (<IN>){
          next if (length $_ < 2);
          my @unit=split("\n",$_);
          my $head=shift @unit;
          #print "$head\n";
          my $seqlen=length join("",@unit);
          $hash{$head}=$seqlen;
     }
close IN;     
$/="\n";

foreach (keys %hash){ ## for each target
my $head=$_;
my $seqlen=$hash{$head};
my $svg=SVG->new(width=>600,height=>1000);
my $starth=100;
my $endh=900;
my $startw=50;
my $endw=400;
my $height=$endh-$starth;
my $scalar=$seqlen/($endh-$starth);
my $rec=$svg->rectangle(
        x=>$startw, y=>$starth,
        width=>10, height=>$height,
        style=>{
            #stroke=>'black',
            fill=> 'red'
        }
);
my $note=$svg->text(
        x=>150, y=>90,
        style=>{
             stroke=>'black'
        }
)->cdata($head);
my $file=$head.".svg";
#writesvg($file,$svg);
$svghash{$head}=$svg;
$ratehash{$head}=$scalar;
}  ## end of foreach
my $refrate=\%ratehash;
my $refsvg=\%svghash;
#print "%$refsvg\n";
#foreach (keys %$refsvg){
#      print "Check svg: $_\t$refrate->{$_}\n";
#}
return ($refsvg,$refrate);   
}

sub writesvg {
my ($file,$svg)=@_;
#print "$file\n";
open OUT, ">$file" or die "can not open my file";
       print OUT $svg->xmlify;
close OUT;
       system "/share/raid12/chenjinfeng/tools/draw/svg2xxx_release/svg2xxx $file -t png";
}
