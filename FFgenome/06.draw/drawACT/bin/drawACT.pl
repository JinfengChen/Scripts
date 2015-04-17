#!/usr/bin/perl

#### read two embl and 4ACT, draw comparision figure as well as annotation if present in embl
#### perl drawACT.pl -q query -t target -a 4ACT > log & 
use warnings;
use strict;
use SVG;
use Getopt::Long;

my ($query,$target,$act,$desc,$help);

GetOptions(
        "query:s" => \$query, 
        "target:s"=> \$target,
        "act:s"   => \$act,
        "desc:s" => \$desc, ### title of figure
        "help:s"  => \$help
);
#print "$query\n$target\n$act\n";
die "Usage: perl drawACT.pl -q query.embl -t target.embl -a 4ACT -d rice2ff > log &" if ($help);
die "Usage: perl drawACT.pl -q query.embl -t target.embl -a 4ACT -d rice2ff> log &" unless (-f $query and -f $target and -f $act);

#################################### readin input files##################3
####parse emble file store position in reference of array
my ($qlength,$qgene,$qdnate,$qrnate)=parse($query);
my ($tlength,$tgene,$tdnate,$trnate)=parse($target);
####
my $identity=50;
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
my $rate; ### scale for draw sequence, 200kb/800.
if ($qlength > $tlength){
    $rate=$qlength/600;
}else{
    $rate=$tlength/600;
}
my $svg=SVG->new(width=>800,height=>600);


#####title, description
my $title=$svg->text(
       x=>300,y=>50,
       style=>{
           fontsize=>'2'
       }
)->cdata($desc);

### draw query line;
my $qwidth=$qlength/$rate;
my $qseq=$svg->rectangle(
            x=>100,y=>100,
            width=>$qwidth,height=>10,  
            style=>{
            stroke=>'black'
            }
);


### draw target line;
my $twidth=$tlength/$rate;
my $tseq=$svg->rectangle(
            x=>100,y=>500,
            width=>$twidth,height=>10,
            style=>{
            stroke=>'black'
            }
);

my $blue="blue";
my $red="red";
my $black="black";
my $pos1="query";
my $pos2="target";
##draw query box
$svg=drawbox($qgene,$black,$svg,$rate,$pos1);
$svg=drawbox($qdnate,$blue,$svg,$rate,$pos1);
$svg=drawbox($qrnate,$red,$svg,$rate,$pos1);
##draw target box
$svg=drawbox($tgene,$black,$svg,$rate,$pos2);
$svg=drawbox($tdnate,$blue,$svg,$rate,$pos2);
$svg=drawbox($trnate,$red,$svg,$rate,$pos2);

#### draw links between query and target
foreach (@links){
    my @unit=split("\t",$_);
    my $qleft=$unit[0]/$rate+100;
    my $qright=$unit[1]/$rate+100;
    my $tleft=$unit[2]/$rate+100;
    my $tright=$unit[3]/$rate+100;
    my $color;
    if ($tright < $tleft){
         $color='red';
    }else{
         $color='#778899';
    }
   
    my $qheight=110;
    my $theight=500;
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

######### draw legend
 
$svg=legend($svg,100,"gene","black");
$svg=legend($svg,170,"LTR","red");
$svg=legend($svg,220,"DNA TE","blue");
sub legend{
my ($svg,$x,$name,$color)=@_;
my $xtext=$x+20; 
 $svg->rectangle(
            x=>$x,y=>550,
            width=>15,height=>10,
            style=>{
                fill=>$color
            }
 );
 $svg->text(
            x=>$xtext,y=>560,
            style=>{stroke=>'black',
                   fontsize=>'1','text-anchor'=>'start','font-weight'=>100,'stroke-width'=>0.1
            }
 )->cdata($name);
 
return $svg;
}

########



my $outfile="$desc.svg";
writesvg($outfile,$svg);
##############################################################


################################### sub for write svg to file
sub writesvg {
my ($file,$svg)=@_;
#print "$file\n";
open OUT, ">$file" or die "can not open my file";
       print OUT $svg->xmlify;
close OUT;
       system "/share/raid12/chenjinfeng/tools/draw/svg2xxx_release/svg2xxx $file -t png";
}


################################ sub for draw feature box

sub drawbox{
my ($refarray,$color,$svg,$rate,$pos)=@_;

my $y;
my $ytext;
if ($pos eq "query"){
   $y=85;
   $ytext=$y-10;
}else{
   $y=515;
   $ytext=$y+20;
}

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
            width=>$length,height=>10,
            style=>{
                fill=>$color
            }
      );

      if ($name =~/(15770)/ or $name =~/(15460)/ or $name =~/(15640)/ or $name =~/(15680)/ or $name =~/(15880)/){
           $name=$1;
           $svg->text(
                x=>$start,y=>$ytext,
                style=>{stroke=>'black',
                   fontsize=>'2','text-anchor'=>'start','font-weight'=>100,'stroke-width'=>0.1
                } 
           )->cdata($name);  

      }
}
return $svg;
}

#######################sub for parse embl
sub parse {
my ($file)=@_;
my @gene;
my @dnate;
my @rnate;
my $length;
my %hash;
my $counter;
open IN, "$file" or die "$!";
  while (<IN>){
      if ($_=~/FT\s+CDS\s+complement/){
             $_=~s/[A-Za-z]//g;
             $_=~s/\s+//g;
             $_=~s/\(|\)//g;
             $_=~/(\d+)\..*\.(\d+)/;
             #print "$_\n$1\t$2\n";
             my $start=$2;
             my $end=$1;
             my $loc=<IN>;
             $loc=~/FT\s+\/gene\=\"(\w+)\"/;
             my $name=$1;
             if (exists $hash{$name}){
                 next;
             }else{
                 $counter++;     
                 $hash{$name}=$counter;
             }
             push (@gene,"$name\t$start\t$end");
      }elsif($_=~/FT\s+CDS\s+join/){
             $_=~s/[A-Za-z]//g;
             $_=~s/\s+//g;
             $_=~s/\(|\)//g;
             $_=~/(\d+)\..*\.(\d+)/;
             #print "$_\n$1\t$2\n";
             my $start=$1;
             my $end=$2;
             my $loc=<IN>;
             $loc=~/FT\s+\/gene\=\"(\w+)\"/;
             my $name=$1;
             if (exists $hash{$name}){
                 next;
             }else{
                 $counter++;
                 $hash{$name}=$counter;
             }
             push (@gene,"$name\t$start\t$end");
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
return ($length,$refgene,$refdnate,$refrnate);
}



