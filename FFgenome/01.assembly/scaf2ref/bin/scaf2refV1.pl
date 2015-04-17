### parse blastout and map scaffold to reference, draw a figure to view the homolog regions.
### perl scaf2ref.pl --query query --target target --blastout blastout --identity identity --hitlength hitlencut
### perl scaf2ref.pl --query chr04.scarSeq --target chr04.txt --blastout ../data/scaf2chr04blastm8 --identity 99 --hitlenth 500
### version1, the target sequence is only one and queries may be hidden if they are overlaped.
### 20091210, by chen jinfeng

use Getopt::Long;
use SVG;
### get in the options 
my ($query,$target,$blastm8,$identity,$lencut,$help);
GetOptions(                ## --query or other could be writen in short, such as --que or --q or -q
    "query:s"=>\$query,    ## s stand for string and i stand for number
    "target:s"=>\$target,
    "blastout:s"=>\$blastm8,
    "identity:i"=>\$identity,
    "lengthcut:i"=>\$lencut,
    "help:s"=> \$help
);
print "$query\t$target\t$blastm8\t$identity\t$lencut\n";
#die `pod2text $0` if (@ARGV == 0 || $Help);
die "Usage: perl scaf2ref.pl --query chr04.scarSeq --target chr04.txt --blastout ../data/scaf2chr04blastm8 --identity 99 --hitlength 500\n" unless (-f $blastm8 and -f $target);
die "Usage: perl scaf2ref.pl --query chr04.scarSeq --target chr04.txt --blastout ../data/scaf2chr04blastm8 --identity 99 --hitlength 500\n" if ($help);
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


### parse blastout and store position information in hash
my %scaf2ref;
my (%qstart,%qend,%tstart,%tend); ## store position information of each scaffold in hash, the key is array for all the hits.
open BLAST, "$blastm8" or die "can not open my blastout";
     while (<BLAST>){
           my @unit=split("\t",$_);
           if ($unit[2] < $identity){next};
           my $hitlength=$unit[7]-$unit[6];
           print "$hitlength\t$lencut\n";
           if ($hitlength < $lencut){next};
           my (@qstart,@qend,@tstart,@tend);
           if (exists $scaf2ref{$unit[0]}){
                my $refqs=$qstart{$unit[0]};
                my $refqe=$qend{$unit[0]};
                my $refts=$tstart{$unit[0]};
                my $refte=$tend{$unit[0]};
                push (@$refqs,$unit[6]); 
                push (@$refqe,  $unit[7]);
                push (@$refts,$unit[8]);
                push (@$refte,  $unit[9]);
           }else{
                $counter++;
                $scaf2ref{$unit[0]}=$counter;
                push (@qstart,$unit[6]);
                push (@qend,  $unit[7]);
                push (@tstart,$unit[8]);
                push (@tend,  $unit[9]);
                $qstart{$unit[0]}=\@qstart;
                $qend{$unit[0]}  =\@qend;
                $tstart{$unit[0]}=\@tstart;
                $tend{$unit[0]}  =\@tend;
           }
     }
close BLAST;


my ($refsvg,$scalar,$head,$file,$starth,$endh,$startw,$endw)=drawref($target);
my $svg=$refsvg->{$head};
foreach (sort keys %scaf2ref){
      #print "$_\n";
      my $scaf=$_;
      #my $refqs=$qstart{$_};
      #print "@$refqs,\n";
      my $refqs=$qstart{$_};
      my $refqe=$qend{$_};
      my $refts=$tstart{$_};
      my $refte=$tend{$_};
      my $mean=mean(@$refts);
      my $start=($mean-$scaflen{$_}/2)/$scalar+100;
      my $end  =($mean+$scaflen{$_}/2)/$scalar+100;
      my $scafheight=$end-$start;
      my $rec=$svg->rectangle(
                   x=>$endw,y=>$start,
                   width=>5,height=>$scafheight,
                   style=>{
                          stroke=>'black',
                          fill  =>'red'
                   }
      );
      my $textw=$endw+15;
      my $texth=$mean/$scalar+100;
      my $text=$svg->text(
                   x=>$textw, y=>$texth,
                   style=>{stroke=>'black',
                   fontsize=>'7','text-anchor'=>'start','font-weight'=>100,'stroke-width'=>0.1
                   }
      )->cdata($scaf);
      for (my $i=0;$i<@$refqs;$i++){
              my $scafstart=$$refqs[$i]/$scalar+$start;
              my $scafend  =$$refqe[$i]/$scalar+$start;
              my $refstart =$$refts[$i]/$scalar+100;
              my $refend   =$$refte[$i]/$scalar+100;
              my $left=$startw+10;
              my $color;
              my $xv=[$left,$endw,$endw,$left];
              my $yv=[$refstart,$scafstart,$scafend,$refend];
              if ($refend > $refstart){
                 $xv=[$left,$endw,$endw,$left];
                 $yv=[$refstart,$scafstart,$scafend,$refend];
                 $color="#FFDAB9";
              }else{
                 $color="#778899";
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
  ## use line to connect start and end site
              my $line=$svg->line(
                      x1=>$endw,y1=>$scafstart,
                      x2=>$left,y2=>$refstart,
                      style=>{
                            stroke=>$color
                      }
              );
              my $line=$svg->line(
                      x1=>$endw,y1=>$scafend,
                      x2=>$left,y2=>$refend,
                      style=>{
                            stroke=>$color
                      }
              );
=cut
      }
}
writesvg($file,$svg);




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
my $head;
my $seqlen;
$/=">";
open IN, "$target" or die "can not open my reference file";   
     while (<IN>){
          next if ($_ eq "");
          my @unit=split("\n",$_);
          $head=shift @unit;
          #print "$head\n";
          $seqlen=length join("",@unit);
     }
close IN;     
$/="\n";
my $svg=SVG->new(width=>600,height=>1000);
my $starth=100;
my $endh=900;
my $startw=50;
my $endw=500;
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
my $refsvg=\%svghash;
return ($refsvg,$scalar,$head,$file,$starth,$endh,$startw,$endw);   
}

sub writesvg {
my ($file,$svg)=@_;
#print "$file\n";
open OUT, ">$file" or die "can not open my file";
       print OUT $svg->xmlify;
close OUT;
       #system "/home/chenjinfeng/tools/draw/svg2xxx_release/svg2xxx $file -t png";
}
