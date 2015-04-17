### map scaffold sequence to fpc map,draw figure to view a good look 
use SVG;
## version 1, scaffolds are anchored to the mid position of each contig, but if multi scaffold hit on one ctg it is unable 
## to sort them out. Need to sort before draw 
## It works well when the sequence if very few and large.


#die "The scripts need contig file and sort.scafmap file\nUsage: perl map2fpc.pl > log &\n program should be run in dir where the files are";
##The scripts need contig file and sort.scafmap file.
##Usage: perl map2fpc.pl > log &

###
 my %scaflenhash;
 open SCAF, "scaffold" or die "can not open my scaffold file";
        while (<SCAF>){
             chomp $_;
             my @unit=split("\t",$_);
             $scaflenhash{$unit[0]}=$unit[1];
        }
close SCAF;

####draw fpc map for each chromosome and store in gif##
my %ctglength;
my %chrlength;
my %ctgchr;
##read contig information from contig and store in hash
open CONTIG, "contig" or die "can not open my contig\nThe scripts need contig file and sort.scafmap file\nUsage: perl map2fpc.pl > log &\n program should be run in dir where the files are";
while (<CONTIG>){
      chomp $_;
      my @unit=split("\t",$_);
      if (exists $chrlength{$unit[2]}){
         $chrlength{$unit[2]}=$chrlength{$unit[2]}+$unit[1];
         #print "$_\n$chrlength{$unit[2]}\n";
         $ctglength{$unit[0]}=$unit[1];
         $ctgchr{$unit[0]}=$unit[2];  
      }else{ 
         $chrlength{$unit[2]}=$unit[1];
         $ctglength{$unit[0]}=$unit[1];
         $ctgchr{$unit[0]}=$unit[2];
      }
}
close CONTIG;    
my $refsvg=drawchr(\%ctglength,\%chrlength,\%ctgchr); ## if the chromosome figure has been drawn, this line can be killed in the later.

$refsvg=mapscaf(\%ctglength,\%chrlength,\%ctgchr,$refsvg);

foreach (keys %$refsvg){
     print "$_\n"; 
     my $svgfile=$_.".svg";
     my $svgref=$refsvg->{$_};
     writesvg($svgfile,$svgref);
}

sub mapscaf{
my ($ctglength,$chrlength,$ctgchr,$refsvg)=@_;
my %scaf2ctg;##store all the ctg for a scaffold, the value is a ref of hash that contain $ctg->@position;
my %scaf2num;##store the ctg the contain maximum hit for a scaffold, the value the ctg number that has maximum number;
open IN, "sort.scafmap" or die "can not open my sort.scafmap";
       while (<IN>){
              chomp $_;
              my @unit=split("\t",$_);
              if (length $unit[2] == 0) {next};
              my %ctg2pos;
              my %ctg2num;
              my @position;
              if (exists $scaf2ctg{$unit[0]}){
                   if (exists $scaf2ctg{$unit[0]}->{$unit[2]}){
                          $scaf2num{$unit[0]}->{$unit[2]}+=1;
                          my $refarray=$scaf2ctg{$unit[0]}->{$unit[2]};
                          push (@$refarray,"$unit[3]\t$unit[4]"); 
                          $scaf2ctg{$unit[0]}->{$unit[2]}=$refarray;
                   }else{ ## if while into a new ctg for this scaffold 
                          $scaf2num{$unit[0]}->{$unit[2]}+=1; 
                          push (@position,"$unit[3]\t$unit[4]");
                          my $refhash= $scaf2ctg{$unit[0]};
                          $refhash->{$unit[2]}=\@position;
                          $scaf2ctg{$unit[0]}=$refhash;
                   }
                   
              }else{
                   push (@position,"$unit[3]\t$unit[4]");
                   $ctg2pos{$unit[2]}=\@position;  ## store postion in @ and stored as hash ctg->@pos
                   $scaf2ctg{$unit[0]}=\%ctg2pos;  ## ctg->@pos as hash and stored as value in scaf->hash
                   $ctg2num{$unit[2]}+=1;
                   $scaf2num{$unit[0]}=\%ctg2num;  ## ctg->num as hash and stored as value in scaf->hash
              }
       }  
close IN;
foreach (sort keys %scaf2ctg){  ## in each scaffold
       print "Scaf\t$_\n";
       my $scaf=$_;
       my $refscaf=$scaf2ctg{$scaf};
       my $refnum=$scaf2num{$scaf};
       my ($maxctg,$maxnum)=maxkey($refnum);
       ## para for sub drawscaf
       my $scaflen=$scaflenhash{$scaf};
       
       $scafchr=$ctgchr->{$maxctg};
       my $scafctg=$maxctg;
       ##
       print "Maxhitctg\t$maxctg\tMaxnum\t$maxnum\n";
       $refsvg=drawscaf($scaflen,$scafchr,$scafctg,$refsvg,$scaf,$refscaf,$chrlength,$ctglength,$ctgchr);
       foreach (sort {$a <=> $b} keys %$refscaf){ ## for each contig
              print "Ctg\t$_\n";
              my $refhash=$scaf2ctg{$scaf};
              my $refarray=$refhash->{$_};
              print "Position\t@$refarray\n";
              print "HitNum\t$scaf2num{$scaf}->{$_}\n";
       }
}
return $refsvg;
}
###get scaf length from it name chr04.con10360001:11694001
sub len {
my ($head)=@_;
$head=~/chr04.con(\d+)\:(\d+)/;
my $length=abs($2-$1);
print "$head\t$length\n";
return $length;
}
###
####draw every scaffold into the png of each chromosome
sub drawscaf{
my ($scaflen,$scafchr,$scafctg,$refsvg,$scaf,$refscaf,$chrlength,$ctglength,$ctgchr)=@_;
my $svg=$refsvg->{$scafchr};
my $notefile=$scafchr.".note";
my %ctgpos;
open IN, "$notefile" or die "can not open my $scafchr note file";
     while (<IN>){
        chomp $_;  
        my @unit=split("\t",$_);
        $ctgpos{$unit[0]}=$unit[1];## mid point of ctg in svg figure
     }
close IN;   
my $scafpos=$ctgpos{$scafctg};
open ER, ">check.txt";

my $rate=20000000/1000; # rate for scaffold drawing
my $ratio=$chrlength->{$scafchr}/800; # rate for ctg drawing

my $conpoint1;
my $conpoint2;
my @twopoint;
###get two point on scaf and ctg to be connected
   #open ER, ">check.txt";
   foreach (sort {$a <=> $b} keys %$refscaf){ ## for each contig, the position here is clone position on ctg as CB unit. CB unit/ratio=length in figure
              my $refarray=$refscaf->{$_};## each element contain (the position on ctg and position on scaffold) for one hit, the whole array is for one ctg
              my $halfctglen=$ctglength->{$_}/2;
              my $halfscaflen=$scaflen/2;
              my $ctgin=$_;
              if ($ctgchr->{$ctgin} ne $scafchr){
                    next;
              }
              foreach (@$refarray){ ## each time for one hit
                    print ER "$_\n";         
                    my @unit=split("\t",$_);
                    my $ctgtest=$unit[0];
                    if ($ctgtest < $halfctglen ){
                       $conpoint1=$ctgpos{$ctgin}-($halfctglen-$ctgtest)/$ratio; ## scaffold mid points is the same with cresponding ctg mid point 
                    }else{
                       $conpoint1=$ctgpos{$ctgin}+($ctgtest-$halfctglen)/$ratio;
                    }
                    my $scaftest=$unit[1];
                    if ($scaftest < $halfscaflen){
                       $conpoint2=$scafpos-($halfscaflen-$scaftest)/$rate;
                    }else{
                       $conpoint2=$scafpos+($scaftest-$halfscaflen)/$rate;
                    }
              print ER "$conpoint1\t$conpoint2\n";
              push (@twopoint, "$conpoint1\t$conpoint2");    
            } 
              
   }
  # close ER;

###

my $x=300;
my $y1=$scafpos-$scaflen/($rate*2);
my $y2=$scafpos+$scaflen/($rate*2);
print ER "$scaf\t$y1\t$y2\n";
my $line=$svg->line( ## scaf line
         x1=>$x, y1=> $y1,      
         x2=>$x, y2=> $y2,
         style=>{stroke=>'red'}
   );

foreach (@twopoint){ ## connection line
     my @unit=split("\t",$_);
     my $line=$svg->line(
            x1=>130,y1=>$unit[0],
            x2=>300,y2=>$unit[1],
            style=>{stroke=>'red'}
     );
}

my $notex=320;
my $note=$svg->text( ## note for scaf
         x=>$notex,y=>$scafpos,
         style=>{stroke=>'black',fontsize=>'7','stroke-width'=>0.1,'font-weight'=>100}
   )->cdata($scaf);
$refsvg->{$scafchr}=$svg;
close ER;
return $refsvg;
}
####
###give the sub a hash, it will return you a key for max value of hash ###
sub maxkey {
my ($hash)=@_;   
my $maxvalue;
my $maxkey;
foreach (keys %$hash){
      if ($hash->{$_} > $maxvalue){
          $maxvalue=$hash->{$_};
          $maxkey=$_;
      }
}
return ($maxkey,$maxvalue);
}

#### draw chromosome ###
sub drawchr{
my ($ctglength,$chrlength,$ctgchr)=@_;
my %svghash;
foreach (sort keys %$chrlength){
     print "$_\n";
     my $chr=$_;
     my $svg=SVG->new(width=>600,height=>1000);
     my $file=$_.".svg";
     my $starth=100;
     my $endh=900;
     my $startw=110;
     my $endw=500;
     my $ratio=$chrlength->{$_}/($endh-$starth);
     my @ctg;
     print "$chr\t$chrlength->{$_}\t$ratio\n";
     foreach (keys %$ctgchr){ ## sort contig in a chromosome, so as to draw them in the order of physical map
          #print "$_\t";
          if ($ctgchr->{$_} eq $chr){push (@ctg, $_)}
     }
     @ctg= sort {$a <=> $b} @ctg;
     #print "@ctg", "\n";
     my $midfile=$chr.".note";## file has the infomation for midpoint of each ctg, used to draw scaffold into the figure   
     open NN, ">$midfile" or die "can not open my out nn";
     foreach (@ctg){
         print "ctg: $_\n";
         my $ctgwidth=20;
         my $ctgheight=$ctglength->{$_}/$ratio;
         print "$startw\t$starth\t$ctgwidth\t$ctgheight\n";
         my $midheight=$ctgheight/2+$starth;
         print NN "$_\t$midheight\n";
         my $rec=$svg->rectangle(
                       x=>$startw, y=>$starth,
                       width=>$ctgwidth,height=>$ctgheight,
                       style=>{
                             stroke=>'black',
                             fill=>'#FFDAB9' 
                       }
         );
         my $txtw=$startw-20;
         my $note=$svg->text(
                 x=>$txtw, y=>$midheight,
                 style=>{stroke=>'black',
                 fontsize=>'7','text-anchor'=>'end','font-weight'=>100,'stroke-width'=>0.1
                 }
         )->cdata($_);
         $starth+=$ctgheight;
         #print "$ctgheight\n";
     }
     close NN;
=pod     
     my $line=$svg->line(
               x1=>100, y1=>100,
               x2=>500, y2=>100,
               style=>{stroke=>'red'}
     );
     my $rec=$svg-> rectangle (
                  x=>100, y=> 100,
                  width=>20,height=>800,
                  style=> {
                       stroke=>'black',
                       fill=>'red'
                  }
     );
     open OUT, ">$file" or die "can not open my file";
            print OUT $svg->xmlify;
     close OUT;
`     system "/home/chenjinfeng/tools/draw/svg2xxx_release/svg2xxx $file -t png";
=cut     
     $svghash{$chr}=$svg;
}
my $refsvg=\%svghash;
return $refsvg;
}
####

sub writesvg {
my ($file,$svg)=@_;
open OUT, ">$file" or die "can not open my file";
       print OUT $svg->xmlify;
       close OUT;
       #system "/home/chenjinfeng/tools/draw/svg2xxx_release/svg2xxx $file -t png";
}


