#!/usr/bin/perl
=header
Draw feature distribution along chromosome.
Usage: perl drawChrFeature.pl --chrlen ../input/OBa.chrlen --TEfeature ../input/OBa_feature_wave --methylation ../input/OBa_methylation_wave --project OBa > log 2> log2 &
--chrlen chromosome length
--TEfeature: TE feature bar.txt file
--methylation: methylation CG/CHG/CHH bar.txt file
--project name
=cut

use warnings;
use strict;
use SVG;
use Getopt::Long;

our %opt;
GetOptions(\%opt,"chrlen:s","TEfeature:s","GENEfeature:s","methylation:s","project:s","help");

if (keys %opt  < 1 or $opt{help}){
    print "perl drawChrFeature.pl --chrlen ../input/OBa.chrlen --TEfeature ../input/OBa_feature_wave --methylation ../input/OBa_methylation_wave --project OBa > log 2> log2 &\n";
    exit();
}

#########chr length###############
my %chrlen;
open IN, "$opt{chrlen}" or die "$!";
while(<IN>){
    chomp $_;
    my @unit=split("\t",$_);
    $chrlen{$unit[0]}=$unit[1];
}
close IN;
#################################

########################### drawing    ###########################
our $width=1000;
our $height=800;

foreach my $chr (keys %chrlen){
     ###svg
     my $svg=SVG->new(width=>$width,height=>$height);
     ###titile
     my $desc=$opt{project}.".".$chr;
     $svg=title($svg,$desc);
     ###x axis
     my $rate=$chrlen{$chr}/800;
     my $xpos=$height-80;
     $svg=drawxaxis($svg,100,900,$xpos,0,$chrlen{$chr},1000000,$rate);
     
     ###draw methlation
     $svg=drawmethylation($svg,$opt{methylation},$chr,$rate);

     ###draw TE feature
     $svg=drawTEfeature($svg,$opt{TEfeature},$chr,$rate);
     
     ###write svg 
     my $outfile="$desc.svg";
     writesvg($outfile,$svg);
} ## foreach chromosome


##draw x axis for ref and qry
#my $qaxisend=100+$qlength/$rate1;
#my $qaxispos=$querypos-15;
#print "$qmin\t$qmax\n$tmin\t$tmax\n";
#$svg=drawxaxis($svg,"100",$qaxisend,$qaxispos,$qmin,$qmax,"200000",$rate1);
#my $taxisend=100+$tlength/$rate2;
#my $taxispos=$targetpos+20;
#$svg=drawxaxis($svg,"100",$taxisend,$taxispos,$tmin,$tmax,"200000",$rate2);
######### draw graph of signal 
##query H3K9
#if (exists $opt{qH3K9}){
#  $svg=graph($svg,"100","150",$opt{qH3K9},"50",$rate1,"H3K9","Orange");
#}


##########draw methylation ###################
##chr01.CG.wave.fill.wig.0.smooth.bar.txt
sub drawmethylation
{
my ($svg,$dir,$chr,$rate)=@_;
my @file=glob("$dir/$chr.*.bar.txt");
foreach my $file (@file){
    if ($file=~/$chr.*.0.smooth.bar.txt/){
        my $tempy=$height-200;
        print "$file\n$tempy\n";
        $svg=graph($svg,100,$tempy,$file,"50000",$rate,"CG","Red");
    }elsif($file=~/$chr.*.1.smooth.bar.txt/){
        my $tempy=$height-150;
        $svg=graph($svg,100,$tempy,$file,"50000",$rate,"CHG","blue");
    }elsif($file=~/$chr.*.2.smooth.bar.txt/){
        my $tempy=$height-100;
        $svg=graph($svg,100,$tempy,$file,"50000",$rate,"CHH","green");
    }

}
return $svg;
}

#########draw TE feature############################
##chr01.MITE.wave.fill.smooth.bar.txt
sub drawTEfeature
{
my ($svg,$dir,$chr,$rate)=@_;
my @file=glob("$dir/$chr.*.bar.txt");
foreach my $file (@file){
    if ($file=~/$chr.(GENE).*.smooth.bar.txt/){
        my $note=$1;
        my $tempy=$height-600;
        $svg=graph($svg,100,$tempy,$file,"50000",$rate,$note,"red");
    }elsif ($file=~/$chr.(RT).*.smooth.bar.txt/){
        my $note=$1;
        my $tempy=$height-550;
        $svg=graph($svg,100,$tempy,$file,"50000",$rate,$note,"brown");
    }elsif($file=~/$chr.(COPIA).*.smooth.bar.txt/){
        my $note=$1;
        my $tempy=$height-500;
        $svg=graph($svg,100,$tempy,$file,"50000",$rate,$note,"brown");
    }elsif($file=~/$chr.(GYPSY).*.smooth.bar.txt/){
        my $note=$1;
        my $tempy=$height-450;
        $svg=graph($svg,100,$tempy,$file,"50000",$rate,$note,"brown");
    }elsif($file=~/$chr.(DNA).*.smooth.bar.txt/){
        my $note=$1;
        my $tempy=$height-400;
        $svg=graph($svg,100,$tempy,$file,"50000",$rate,$note,"purple");
    }elsif($file=~/$chr.(CACTA).*.smooth.bar.txt/){
        my $note=$1;
        my $tempy=$height-350;
        $svg=graph($svg,100,$tempy,$file,"50000",$rate,$note,"purple");
    }elsif($file=~/$chr.(MUDR).*.smooth.bar.txt/){
        my $note=$1;
        my $tempy=$height-300;
        $svg=graph($svg,100,$tempy,$file,"50000",$rate,$note,"purple");
    }elsif($file=~/$chr.(MITE).*.smooth.bar.txt/){
        my $note=$1;
        my $tempy=$height-250;
        $svg=graph($svg,100,$tempy,$file,"50000",$rate,$note,"purple");
    }
}
return $svg;
}


##########title ################
sub title
{
my ($svg,$desc)=@_;
my $title=$svg->text(
       x=>50,y=>50,
       style=>{
           fontsize=>'2'
       }
)->cdata($desc);
return $svg;
}


################ sub functions#############################################
sub graph
{
my ($svg,$x,$y,$bar,$win,$rate,$title,$color)=@_;

my $axisx=$x-20;
my $axisy1=$y-40;
my $axisy2=$y;
#my $rate4yaxis=1;
#($svg,$rate4yaxis)=drawyaxis($svg,$axisx,$axisy1,$axisy2,"0","40","10",$title);

my %hash;
open IN, "$bar" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my $pos  =$unit[1];
    $hash{$pos}=$unit[2];
}
close IN;
my @signal=values %hash;
my $tops=topnum(@signal);
chomp $tops;
my $step=10;

if ($tops <= 20){
   $step=10;
}elsif($tops > 20 and $tops < 60){
   $step=10;
}elsif($tops >= 60 and $tops <= 100){
   $step=20;
}else{
   $step=40;
}



print "$tops\t$step\n";
my ($svg1,$rate4yaxis)=drawyaxis($svg,$axisx,$axisy1,$axisy2,"0",$tops,$step,$title);
$svg=$svg1;
#print "$rate4yaxis\n";

foreach (sort {$a <=> $b} keys %hash){
       #print "$_\t$hash{$_}\t$rate4yaxis\n";
    #if ($hash{$_} > 100){
    if ($hash{$_} > $tops){
       $hash{$_}=$tops;
       my $tleft=$_/$rate+100;
       my $tright=($_+$win)/$rate+100;
       my $bright=$tright;
       my $bleft=$tleft;
       my $theight=$y-$hash{$_}*$rate4yaxis;
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
       my $theight=$y-$hash{$_}*$rate4yaxis;
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
#my $index=int ($total*0.95);
my $top=pop @array;
$top=10*(int ($top/10))+10;
print "$top\n";
return $top;
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
my $tail =int ($max/100000)/10;

my $bottomtext=$svg->text(
     x=>$x2,y=>$y+13,
     style=>{'stroke-opacity'=>'0',
         fontsize=>'0.01','text-anchor'=>'start','font-weight'=>100,'stroke-width'=>0.01
     }
)->cdata("$tail Mb");

for(my $i=$min+$step;$i<$max;$i=$i+$step){
     my $tempx=$x1+($i-$min)/$rate;
     #print "$min\t$step\t$i\t$rate\n";
     #print "$tempx\t$y\n";
     my $line=$svg->line(
         x1=>$tempx,y1=>$y,
         x2=>$tempx,y2=>$y-5,
         style=>{stroke=>'black'}
     );
     my $tempi=int ($i/1000000);
     my $text=$svg->text(
         x=>$tempx+3,y=>$y+13,
         style=>{'stroke-opacity'=>'0',
             fontsize=>'0.01','text-anchor'=>'end','font-weight'=>100,'stroke-width'=>0.01
         }
     )->cdata($tempi);
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
     x2=>$x+5,y2=>$y1,
     style=>{stroke=>'black'}
);
my $toptext=$svg->text(
     x=>$x-5,y=>$y1+5,
     style=>{'stroke-opacity'=>'0',
         fontsize=>'0.01','text-anchor'=>'end','font-weight'=>100,'stroke-width'=>0.01
     }
)->cdata($max);
my $bottomline=$svg->line(
     x1=>$x,y1=>$y2,
     x2=>$x+5,y2=>$y2,
     style=>{stroke=>'black'}
);
my $bottomtext=$svg->text(
     x=>$x-5,y=>$y2,
     style=>{'stroke-opacity'=>'0',
         fontsize=>'0.01','text-anchor'=>'end','font-weight'=>100,'stroke-width'=>0.01
     }
)->cdata("0");

my $rate=($y2-$y1+1)/($max-$min+1);
for(my $i=$min+$step;$i<$max;$i=$i+$step){
     my $tempy=$y2-$i*$rate;
     my $line=$svg->line(
         x1=>$x,y1=>$tempy,
         x2=>$x+3,y2=>$tempy,
         style=>{stroke=>'black'}
     );
}
my $textx=$x-60;
my $texty=$y1+($y2-$y1)/2;
my $anno=$svg->text(
     x=>$textx,y=>$texty,
     style=>{'stroke-opacity'=>'0',
         fontsize=>'0.01','text-anchor'=>'start','font-weight'=>100,'stroke-width'=>0.01
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
            style=>{'stroke-opacity'=>'0',
                   fontsize=>'0.01','text-anchor'=>'start','font-weight'=>10,'stroke-width'=>0.1
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


