### read blast m8 file of verify sequence aganist scaffold.
### draw scaffold compared with verify sequence
### 
### Usage: perl verify2scaf.pl -q verify -t scaffold -b blastm8 -i 99 -l 500 > log &
### Version1.0, add gap on the figure. But solexa scaffolds were too large, we need to cut out on
### these hitted region to view the gap distribution and gene or TE annotation.


use SVG;
use Getopt::Long;
use warnings;
use strict;



my ($query,$target,$blastm8,$identity,$lencut,$help);
GetOptions(
         "query:s" => \$query, 
         "target:s"=> \$target, 
         "blastm8:s"=> \$blastm8,
         "identity:i"=> \$identity,
         "lencut:i" => \$lencut,
         "help:s"   => \$help
);
die "perl verify2scaf.pl -q verify -t scaffold -b blastm8 -i 99 -l 500 > log &" if ($help);

################### get sequence length and store them in hash
my %verifylen;
my %verifygap;
my %scaflen;
my %scafgap;
$/=">";
open IN, "$query" or die "$!";
while (<IN>){
     chomp $_;
     next if ($_ eq "");
     my @unit=split("\n",$_);
     my $name=shift @unit;
     my $seq=join("",@unit);
     $seq=~s/\s//g;
     $seq=~s/\r//g;
     $seq=~s/\n//g;
     my $seqlen=length $seq;
     $verifylen{$name}=$seqlen;
     my @gaparray;
     while ($seq=~/(N+)/g) {
          my $gaplen=length $1;
          my $gappos=pos($seq);
          my $gapstart=$gappos-$gaplen;
          my $gapend  =$gappos-1;
          push (@gaparray,"$gapstart\t$gapend");
     }
     $verifygap{$name}=\@gaparray;
}
close IN; 
open IN, "$target" or die "$!";
while (<IN>){
     chomp $_;
     next if ($_ eq "");
     my @unit=split("\n",$_);
     my $name=shift @unit;
     my $seq=join("",@unit);
     $seq=~s/\s//g;
     $seq=~s/\r//g;
     $seq=~s/\n//g;
     my $seqlen=length $seq;
     $scaflen{$name}=$seqlen;
     my @gaparray;
     while($seq=~/(N+)/g){
          my $gaplen=length $1;
          my $gappos=pos($seq);
          my $gapstart=$gappos-$gaplen;
          my $gapend  =$gappos-1;
          push (@gaparray,"$gapstart\t$gapend");
     }
     $scafgap{$name}=\@gaparray;

}
close IN;
$/="\n";
####################


############################## read blastm8 store hit information in structured data 
my %verify2scaf;
my %verify2num; ## record scaffold hit number of each verify seq
open IN, "$blastm8" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    if ($unit[2] < $identity){next};
    my $hitlength=$unit[7]-$unit[6];
    if ($hitlength < $lencut){next};    
    my @position;
    my %scaf2pos;
    my @hitscaf;
    if (exists $verify2scaf{$unit[0]}){
        if (exists $verify2scaf{$unit[0]}->{$unit[1]}){
             push (@{$verify2scaf{$unit[0]}->{$unit[1]}},"$unit[6]\t$unit[7]\t$unit[8]\t$unit[9]");
        }else{
             push (@position,"$unit[6]\t$unit[7]\t$unit[8]\t$unit[9]");
             my $refhash=$verify2scaf{$unit[0]};
             $refhash->{$unit[1]}=\@position;
             $verify2scaf{$unit[0]}=$refhash;  
        }   

    }else{
        push (@position,"$unit[6]\t$unit[7]\t$unit[8]\t$unit[9]");  
        $scaf2pos{$unit[1]}=\@position;
        $verify2scaf{$unit[0]}=\%scaf2pos;
    }
}
close IN;
################################

###################### ordered scaffold by midpoint on verify sequence

foreach(sort keys %verify2scaf) {
       my $verify=$_;
       print "$verify\n";
       my $verilen=$verifylen{$verify};
       my $refgap=$verifygap{$verify};
       my ($svghash,$ratehash)=drawref($verify,$verilen,$refgap);
       my $rate4veri=$ratehash->{$verify}; 
       my %midpoint;
       my %scafrev; ### record is scaffold reverse hitted on ref
       foreach (sort keys %{$verify2scaf{$verify}}){
           my $scaffold=$_;
           print "$scaffold\n";
           my $refarray=$verify2scaf{$verify}->{$scaffold};
           #print "@$refarray\n";
           my @startonref;
           my $addlen; ### add hit length for a scaffold, if < 1000, next
           foreach (@$refarray){
                  my @unit=split("\t",$_);
                  $addlen+=$unit[1]-$unit[0];
                  push (@startonref,$unit[0]); 
                  if ($unit[2] > $unit[3]){
                     $scafrev{$scaffold}+=1;
                  }
           }
           if ($addlen < 1000){
                  next;
                  delete $verify2scaf{$verify}->{$scaffold};
           } 
           my $mean=mean(@startonref);
           $midpoint{$scaffold}=$mean;
       }
       my @scaffold=sort {$midpoint{$a} <=> $midpoint{$b}} keys %midpoint;
       my $scafchrlen;
       my $midonchr=0;
       my $endonchr=0;
       my %scafonchr;
       foreach (@scaffold){
             #print "$_\t$midpoint{$_}\n";
             $scafchrlen+=$scaflen{$_};
             $midonchr=$scaflen{$_}/2+$endonchr;
             $endonchr+=$scaflen{$_};
             $scafonchr{$_}=$midonchr; 
       }
       my $refsvg=$svghash->{$verify};
       my $rate4scaf=$scafchrlen/600;
       foreach (@scaffold){
             my $scaf=$_;
             my $start=($scafonchr{$scaf}-$scaflen{$scaf}/2)/$rate4scaf+150;
             my $end  =($scafonchr{$scaf}+$scaflen{$scaf}/2)/$rate4scaf+150;
             my $width=$end-$start;
             my $rec  =$refsvg->rectangle(
                      x=>$start,y=>500,
                      width=>$width,height=>10,
                      style=>{
                           stroke=>'black',
                           fill=>'red'
                      }
             );
             my $textx=$start+$width/2;
             my $text =$refsvg->text(
                      x=>$textx,y=>540,
                      style=>{
                           #stroke=>'black',
                           fontsize=>'0.7',
                           'text-anchor'=>'end'
                      }
             )->cdata($scaf);
=pod
             my $refgap2=$scafgap{$scaf};
             foreach (@$refgap2){
                    my @unit=split("\t",$_);
                    my $gapstart=$unit[0]/$scale+;
                    my $gapend  =$unit[1]/$scale+150;
                    my $gaplen=$gapend-$gapstart+1;
                    my $gap=$svg->rectangle(
                                    x=>$gapstart,y=>95,
                                    width=>$gaplen,height=>5,
                                    style=>{
                                        fill=>'black'
                                    }
                    );
             }
=cut
             my $refarray=$verify2scaf{$verify}->{$scaf}; 
             my $reverse=0;
             foreach (@$refarray){
                   my @unit=split("\t",$_);
                   my ($veristart,$veriend,$scafstart,$scafend);
                   if (exists $scafrev{$scaf} and $scafrev{$scaf} > @$refarray/2){
                      $reverse=1;
                      $veristart=$unit[0]/$rate4veri+150;
                      $veriend  =$unit[1]/$rate4veri+150;
                      $scafstart=($scaflen{$scaf}-$unit[2])/$rate4scaf+$start;
                      $scafend  =($scaflen{$scaf}-$unit[3])/$rate4scaf+$start;
                   }else{ 
                      $veristart=$unit[0]/$rate4veri+150;
                      $veriend  =$unit[1]/$rate4veri+150;
                      $scafstart=$unit[2]/$rate4scaf+$start;
                      $scafend  =$unit[3]/$rate4scaf+$start;
                   }
                   my $xv=[$veristart,$veriend,$scafend,$scafstart];
                   my $yv=[110,110,500,500];
                   my $points=$refsvg->get_path(
                              x=>$xv,y=>$yv,
                              -type=>'polyline',
                              -closed=>'true'
                   ); 
                   my $tag=$refsvg->polyline(
                              %$points,
                              style=>{
                                  fill=>'#FFDAB9'
                              }
                   );
             }
             my $refgap2=$scafgap{$scaf};
             foreach (@$refgap2){
                    my @unit=split("\t",$_);
                    my ($gapstart,$gapend,$gaplen);
                    if ($reverse){
                          $gapend=($scaflen{$scaf}-$unit[0])/$rate4scaf+$start;
                          $gapstart  =($scaflen{$scaf}-$unit[1])/$rate4scaf+$start;
                          $gaplen=abs($gapend-$gapstart)+1;
                    }else{
                          $gapstart=$unit[0]/$rate4scaf+$start;
                          $gapend  =$unit[1]/$rate4scaf+$start;
                          $gaplen=abs($gapend-$gapstart)+1;
                    }
                    my $gap=$refsvg->rectangle(
                                    x=>$gapstart,y=>515,
                                    width=>$gaplen,height=>5,
                                    style=>{
                                        fill=>'black'
                                    }
                    );
             }
 
       }  

       my $svgfile="../output/$verify.svg";
       writesvg($svgfile,$refsvg);
}

#### create svg obj and draw verify sequence, return a obj of svg; 
sub drawref {

my ($verify,$seqlen,$refgap)=@_;
my $svg=SVG->new(width=>800,height=>600);
my $starth=100;
my $endh=500;
my $startw=150;
my $endw=750;
my $width=$endw-$startw;
my $scale=$seqlen/($endw-$startw);
my $rec=$svg->rectangle(
        x=>$startw, y=>$starth,
        width=>$width, height=>10,
        style=>{
            #stroke=>'black',
            fill=> 'red'
        }
);
my $note=$svg->text(
        x=>400, y=>90,
        style=>{
             fontsize=>'0.7', 'text-anchor'=>'end',
             strok=>'black'
        }
)->cdata($verify);
my $note1=$svg->text(
        x=>410, y=>90,
        style=>{
             fontsize=>'0.7', 'text-anchor'=>'start',
             strok=>'black'
        }
)->cdata($seqlen);

foreach (@$refgap){
     my @unit=split("\t",$_);
     my $gapstart=$unit[0]/$scale+150;
     my $gapend  =$unit[1]/$scale+150;
     my $gaplen=$gapend-$gapstart+1;
     my $gap=$svg->rectangle(
                  x=>$gapstart,y=>95,
                  width=>$gaplen,height=>5,
                  style=>{
                      fill=>'black'  
                  }
     );
}

my %svghash;
my %ratehash;
$svghash{$verify}=$svg;
$ratehash{$verify}=$scale;
return (\%svghash,\%ratehash);
}

#######sub mean, return a mean for a array
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

#######sub writesvg, write svg to file
sub writesvg {
my ($file,$svg)=@_;
#print "$file\n";
open OUT, ">$file" or die "can not open my file";
       print OUT $svg->xmlify;
close OUT;
       system "/share/raid12/chenjinfeng/tools/draw/svg2xxx_release/svg2xxx $file -t png";
}








