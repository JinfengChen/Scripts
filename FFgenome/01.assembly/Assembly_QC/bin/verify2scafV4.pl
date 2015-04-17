### read blast m8 file of verify sequence aganist scaffold.
### draw scaffold compared with verify sequence
### 
### Usage: perl verify2scaf.pl -q verify -t scaffold -b blastm8 -g genegff -r repeatgff -i 99 -l 500 > log &
### Version3.0, cut only hitted region for comparision and store sequence in file so we can do alignment.
### modify on 20091220
### Version4.0, modify the figure height.
### Fix a bug in drawing gap and gene annotation. $gaplen  =abs($unit[1]-$unit[0]+1)/$rate4scaf;
### Write the postion of hit region into a file HitPosition4Verify.txt

use SVG;
use Getopt::Long;
#use warnings;
#use strict;

our $y=300;


my ($query,$target,$blastm8,$genegff,$repeatgff,$identity,$lencut,$help);
GetOptions(
         "query:s" => \$query, 
         "target:s"=> \$target, 
         "blastm8:s"=> \$blastm8,
         "genegff:s"=> \$genegff,
         "repeatgff:s"=> \$repeatgff,
         "identity:i"=> \$identity,
         "lencut:i" => \$lencut,
         "help:s"   => \$help
);
die "perl verify2scaf.pl -q verify -t scaffold -b blastm8 -g genegff -r repeatgff -i 99 -l 500 > log &" if ($help);
############### parse gene annotation and repeat annotation gff file
my $scaf2gene;
if (-e $genegff) {
     $scaf2gene=parsegene($genegff);
}
my $scaf2te;
if (-e $repeatgff) {
	 $scaf2te  =parsete($repeatgff);
}


##############
################### get sequence length and store them in hash
my %verifylen;
my %verifygap;
my %scaflen;
my %scafgap;
my %verifyseq;
my %scafseq;
$/=">";
open IN, "$query" or die "$!";
while (<IN>){
     chomp $_;
     next if ($_ eq "");
     my @unit=split("\n",$_);
     my $name=shift @unit;
     my @tempname=split(" ",$name);
     $name=shift @tempname;
     my $seq=join("",@unit);
     $seq=~s/\s//g;
     $seq=~s/\r//g;
     $seq=~s/\n//g;
     my $seqlen=length $seq;
     $verifylen{$name}=$seqlen;
     $verifyseq{$name}=$seq;
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
     my @tempname=split(" ",$name);
     $name=shift @tempname;
     my $seq=join("",@unit);
     $seq=~s/\s//g;
     $seq=~s/\r//g;
     $seq=~s/\n//g;
     my $seqlen=length $seq;
     $scaflen{$name}=$seqlen;
     $scafseq{$name}=$seq;
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
system `mkdir ../output/align` unless (-d "../output/align");
foreach(sort keys %verify2scaf) {
open OUT, ">../output/align/$_.fas" or die "$!";
       my $verify=$_;
       print "$verify\n";
       print OUT ">$verify\n$verifyseq{$verify}\n";
       my $verilen=$verifylen{$verify};
       my $refgap=$verifygap{$verify};
       my ($svghash,$ratehash)=drawref($verify,$verilen,$refgap);
       my $rate4veri=$ratehash->{$verify}; 
       my %midpoint;
       my %hitscaflen; ### hit scaffold length for each scaffold 
       my %hitstart; ###hit start for each scaffold
       my %hitend;   ###hit end for each scaffold
       my %scafrev; ### record is scaffold reverse hitted on ref
       my $compareseq; ### all scaffold seq mapped to on verify seq
       my %hitseq; ## hit seq for each scaffold 
       foreach (sort keys %{$verify2scaf{$verify}}){
           my $scaffold=$_;
           #print "$scaffold\n";
           my $refarray=$verify2scaf{$verify}->{$scaffold};
           #print "@$refarray\n";
           my @startonref;
           my @startonscaf; ### record start and end of hit on scaffold and find out hitted region by min and max sub
           my @endonscaf;
           my $addlen; ### add hit length for a scaffold, if < 1000, next
           foreach (@$refarray){
                  my @unit=split("\t",$_);
                  $addlen+=$unit[1]-$unit[0];
                  push (@startonref,$unit[0]); 
                  if ($unit[2] > $unit[3]){
                     $scafrev{$scaffold}+=1;
                     push (@startonscaf,$unit[3]);
                     push (@endonscaf,$unit[2]);
                  }else{
                     push (@startonscaf,$unit[2]);
                     push (@endonscaf,$unit[3]);
                  }
           }
           my @array1=sort {$a <=> $b} @startonscaf;
           my $hstart;
           for(my $i=0;$i<@array1-2;$i++){
               my $hstart1=$array1[$i];
               my $hstart2=$array1[$i+1];
               my $hstart3=$array1[$i+2];
               my $inter1=abs ($hstart2-$hstart1);
               my $inter2=abs ($hstart3-$hstart2);
               if ($inter1 < 10000 and $inter2 < 10000){
                   $hstart=$array1[$i];
                   last;
               }
           }
           my @array2=sort {$a <=> $b} @endonscaf;
           my $hend;
           for(my $i=@array1-1;$i>=0;$i--){
               my $hend1=$array1[$i];
               my $hend2=$array1[$i-1];
               my $hend3=$array1[$i-2];
               my $inter1=abs ($hend1-$hend2);
               my $inter2=abs ($hend2-$hend3);
               if ($inter1 < 10000 and $inter2 < 10000){
                   $hend=$array1[$i];
                   last;
               }
           } 
           if ($hstart-1000 > 0){ #### add 10k to the start and end of hit region
               $hstart=$hstart-1000;  
           }else{
               $hstart=0;
           }
           if ($hend+1000 < $scaflen{$scaffold}){
               $hend=$hend+1000;
           }else{
               $hend=$scaflen{$scaffold};
           }
           $hitstart{$scaffold}=$hstart;
           $hitend{$scaffold}  =$hend;
           my $hitlength=$hend -$hstart;
           $hitseq{$scaffold}=substr($scafseq{$scaffold},$hstart,$hitlength);
           $hitscaflen{$scaffold}=$hend -$hstart;
           if ($addlen < 5000){
                  delete $verify2scaf{$verify}->{$scaffold};
                  next;
           }else{
                  ####Write verify name and start/end postions of compared region on the chr.
                  open OUT1, ">>HitPosition4Verify.txt" or die "$!";
                      print OUT1 "$verify\t$scaffold\t$hstart\t$hend\t$addlen\n";
                  close OUT1;
                  #####
           } 
           my $mean=mean(@startonref);
           $midpoint{$scaffold}=$mean;
       }
       my @scaffold=sort {$midpoint{$a} <=> $midpoint{$b}} keys %midpoint;
       my $scafchrlen;
       my $hitchrlen;
       my $midonchr=0;
       my $endonchr=0;
       my %scafonchr;
       foreach (@scaffold){
             print "$_\t$midpoint{$_}\n";
             #$scafchrlen+=$scaflen{$_};
             $hitchrlen+=$hitscaflen{$_};
             $midonchr=$hitscaflen{$_}/2+$endonchr;
             $endonchr+=$hitscaflen{$_};
             $scafonchr{$_}=$midonchr; 
       }
       my $refsvg=$svghash->{$verify};
       my $rate4scaf=$hitchrlen/600;
       foreach (@scaffold){
             my $scaf=$_;
             my $start=($scafonchr{$scaf}-$hitscaflen{$scaf}/2)/$rate4scaf+150;
             my $end  =($scafonchr{$scaf}+$hitscaflen{$scaf}/2)/$rate4scaf+150;
             my $width=$end-$start;
             my $rec  =$refsvg->rectangle(
                      x=>$start,y=>$y-100,
                      width=>$width,height=>10,
                      style=>{
                           stroke=>'black',
                           fill=>'red'
                      }
             );
             my $textx=$start+$width/2;
             my $text =$refsvg->text(
                      x=>$textx,y=>$y-45,
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
                      if ($unit[3] >= $hitstart{$scaf} and $unit[2] <= $hitend{$scaf}){
                         $scafstart=($hitscaflen{$scaf}-$unit[2]+$hitstart{$scaf})/$rate4scaf+$start;
                         $scafend  =($hitscaflen{$scaf}-$unit[3]+$hitstart{$scaf})/$rate4scaf+$start;
                      }else{
                         next;
                      }
                   }else{ 
                      $veristart=$unit[0]/$rate4veri+150;
                      $veriend  =$unit[1]/$rate4veri+150;
                      if ($unit[2] >= $hitstart{$scaf} and $unit[3] <= $hitend{$scaf}){
                         $scafstart=($unit[2]-$hitstart{$scaf})/$rate4scaf+$start;
                         $scafend  =($unit[3]-$hitstart{$scaf})/$rate4scaf+$start;
                      }{
                         next;
                      }
                   }
                   my $xv=[$veristart,$veriend,$scafend,$scafstart];
                   my $yv=[60,60,$y-100,$y-100];
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


             ########################## add scaffold sequence to compareseq
             if ($reverse){
                 my $rev=reverse $hitseq{$scaf};
                 $rev=~tr/ATCGatcg/TAGCtagc/; 
                 $compareseq.=$rev; 
             }else{
                 $compareseq.=$hitseq{$scaf};
             }
             ##########################
             ################draw gap for scaffold
             my $refgap2=$scafgap{$scaf};
             foreach (@$refgap2){
                    my @unit=split("\t",$_);
                    my ($gapstart,$gapend,$gaplen);
                    if ($reverse){
                        if ($unit[0] >= $hitstart{$scaf} and $unit[1] <= $hitend{$scaf}){
                          $gapend    =($hitscaflen{$scaf}-$unit[0]+$hitstart{$scaf})/$rate4scaf+$start;
                          $gapstart  =($hitscaflen{$scaf}-$unit[1]+$hitstart{$scaf})/$rate4scaf+$start;
                          #$gaplen    =abs($gapend-$gapstart)+1;
                          $gaplen  =abs($unit[1]-$unit[0]+1)/$rate4scaf;
                        }else{
                          next;
                        }
                    }else{
                        if ($unit[0] >= $hitstart{$scaf} and $unit[1] <= $hitend{$scaf}){  
                          $gapstart=($unit[0]-$hitstart{$scaf})/$rate4scaf+$start;
                          $gapend  =($unit[1]-$hitstart{$scaf})/$rate4scaf+$start;
                          #$gaplen  =abs($gapend-$gapstart)+1;
                          $gaplen  =abs($unit[1]-$unit[0]+1)/$rate4scaf;
                        }else{
                          next;
                        } 
                    }
                    my $gap=$refsvg->rectangle(
                                    x=>$gapstart,y=>$y-85,
                                    width=>$gaplen,height=>5,
                                    style=>{
                                        fill=>'black'
                                    }
                    );
             }
             #################draw gene annotation for scaffold
         if (-e $genegff){    
             print "Draw gene annotation\n"; 
             my $refgene=$scaf2gene->{$scaf};
             foreach (@$refgene){
                    my @unit=split("\t",$_);  
                    my ($genestart,$geneend,$genelen);
                    if ($reverse){
                         if ($unit[0] >= $hitstart{$scaf} and $unit[1] <= $hitend{$scaf}){    
                            $geneend   =($hitscaflen{$scaf}-$unit[0]+$hitstart{$scaf})/$rate4scaf+$start;
                            $genestart =($hitscaflen{$scaf}-$unit[1]+$hitstart{$scaf})/$rate4scaf+$start;
                            #$genelen   =abs($geneend-$genestart)+1;
                            $genelen  =abs($unit[1]-$unit[0]+1)/$rate4scaf;
                         }else{
                            next;
                         }
                    }else{
                         if ($unit[0] >= $hitstart{$scaf} and $unit[1] <= $hitend{$scaf}){
                            $genestart=($unit[0]-$hitstart{$scaf})/$rate4scaf+$start;
                            $geneend  =($unit[1]-$hitstart{$scaf})/$rate4scaf+$start;
                            #$genelen  =abs($geneend-$genestart)+1;
                            $genelen  =abs($unit[1]-$unit[0]+1)/$rate4scaf;
                         }else{
                            next;
                         }
                    }
                    my $gene=$refsvg->rectangle(
                                     x=>$genestart,y=>$y-75,
                                     width=>$genelen,height=>5,
                                     style=>{
                                         fill=>'red' 
                                     }
                    ); 
             }
           }
             #################draw te annotation for scaffold
           if (-e $repeatgff){
             print "Draw TE annotation\n";
             my $refte=$scaf2te->{$scaf};
             foreach (@$refte){
                    my @unit=split("\t",$_);
                    my ($testart,$teend,$telen);
                    if ($reverse){
                         if ($unit[0] >= $hitstart{$scaf} and $unit[1] <= $hitend{$scaf}){  
                            $teend=($hitscaflen{$scaf}-$unit[0]+$hitstart{$scaf})/$rate4scaf+$start;
                            $testart =($hitscaflen{$scaf}-$unit[1]+$hitstart{$scaf})/$rate4scaf+$start;
                            #$telen   =abs($teend-$testart)+1;
                            $telen  =abs($unit[1]-$unit[0]+1)/$rate4scaf;
                         }else{
                            next;
                         }
                    }else{ 
                         if ($unit[0] >= $hitstart{$scaf} and $unit[1] <= $hitend{$scaf}){
                            $testart=($unit[0]-$hitstart{$scaf})/$rate4scaf+$start; 
                            $teend  =($unit[1]-$hitstart{$scaf})/$rate4scaf+$start;
                            #$telen  =abs($teend-$testart)+1;
                            $telen  =abs($unit[1]-$unit[0]+1)/$rate4scaf;
                         }else{
                            next;
                         }
                    }
                    my $te=$refsvg->rectangle(
                                     x=>$testart,y=>$y-65,
                                     width=>$telen,height=>5,
                                     style=>{
                                         fill=>'blue'
                                     }
                    ); 
             }
           } 
             ############################ draw legend
             $refsvg=legend($refsvg,150,"gap","black");
             $refsvg=legend($refsvg,220,"gene","red");
             $refsvg=legend($refsvg,290,"TE","blue"); 
 
             ###########################
       } 
       print OUT ">scaffoldFF\n$compareseq\n"; 
       close OUT;
       my $svgfile="../output/$verify.svg";
       writesvg($svgfile,$refsvg);
}

####################draw legend for gap, gene and te

sub legend{
my ($svg,$x,$name,$color)=@_;
my $xtext=$x+20; 
 $svg->rectangle(
            x=>$x,y=>$y-30,
            width=>15,height=>5,
            style=>{
                fill=>$color
            }
 );
 $svg->text(
            x=>$xtext,y=>$y-20,
            style=>{stroke=>'black',
                   fontsize=>'1','text-anchor'=>'start','font-weight'=>100,'stroke-width'=>0.1
            }
 )->cdata($name);
 
return $svg;
}


#### create svg obj and draw verify sequence, return a obj of svg; 
sub drawref {

my ($verify,$seqlen,$refgap)=@_;
my $svg=SVG->new(width=>800,height=>$y);
my $starth=50;
my $endh=$y-100;
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
        x=>400, y=>$starth-10,
        style=>{
             fontsize=>'0.7', 'text-anchor'=>'end',
             strok=>'black'
        }
)->cdata($verify);
my $note1=$svg->text(
        x=>410, y=>$starth-10,
        style=>{
             fontsize=>'0.7', 'text-anchor'=>'start',
             strok=>'black'
        }
)->cdata($seqlen);

foreach (@$refgap){
     my @unit=split("\t",$_);
     my $gapstart=$unit[0]/$scale+150;
     my $gapend  =$unit[1]/$scale+150;
     my $gaplen=($unit[1]-$unit[0]+1)/$scale;
     my $gap=$svg->rectangle(
                  x=>$gapstart,y=>$starth-5,
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

############# sub parsegene, return a hash, scaf->\@position
sub parsegene{
my ($genegff)=@_;
my %scaf2gene;
open IN, "$genegff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/##gff-version 3/);
    my @unit=split("\t",$_);
    my @position;
    if ($unit[2] eq "CDS"){
        if (exists $scaf2gene{$unit[0]}){
            my $refarray=$scaf2gene{$unit[0]};
            push (@$refarray,"$unit[3]\t$unit[4]");
            $scaf2gene{$unit[0]}=$refarray;
        }else{
            push (@position,"$unit[3]\t$unit[4]");
            $scaf2gene{$unit[0]}=\@position;
        }
    }
}
close IN;
return \%scaf2gene;
}
##############################################################

############ sub parsete, return a hash, scaf->\@position
sub parsete{
my ($repeatgff)=@_;
my %scaf2te;
open IN, "$repeatgff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/##gff-version 3/);
    my @unit=split("\t",$_);
    my @position;
    if (exists $scaf2te{$unit[0]}){
            my $refarray=$scaf2te{$unit[0]};
            push (@$refarray,"$unit[3]\t$unit[4]");
            $scaf2te{$unit[0]}=$refarray;
    }else{
            push (@position,"$unit[3]\t$unit[4]");
            $scaf2te{$unit[0]}=\@position;
    }
}
close IN;
return \%scaf2te;
}
##############################################################

#######sub writesvg, write svg to file
sub writesvg {
my ($file,$svg)=@_;
#print "$file\n";
open OUT, ">$file" or die "can not open my file";
       print OUT $svg->xmlify;
close OUT;
       system "/home/jfchen/FFproject/tools/draw/svg2xxx_release/svg2xxx $file -t pdf -m 200";
}








