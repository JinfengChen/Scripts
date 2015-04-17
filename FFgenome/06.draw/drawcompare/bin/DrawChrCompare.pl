#!/usr/bin/perl
use SVG;
use Getopt::Long;

GetOptions(\%opt,"chrlen:s","table:s","suffix:s","title:s","log2","project:s","help");

my $help=<<USAGE;
perl drawChrFeature.pl --chrlen ../input/OBa.chrlen --table --suffix --title  --project OBa > log 2> log2
--table: directory contains table for each chromosome
--suffix: suffix for table, blocks for chr01.blocks.table
--table: yaxis annotation
--log2: if draw log ratio or absolute value
USAGE


if (keys %opt  < 1 or $opt{help}){
    print "$help\n";
    exit();
}

$opt{title} ||= "log2(Os/Ob)";

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
our $height=400;

foreach my $chr (keys %chrlen){
     print "$chr\n";
     my $table="$opt{table}/$chr.$opt{suffix}.table";
     print "$table\n";
     next unless (-f $table);
     print "Drawing\n"; 
     ###svg
     my $svg=SVG->new(width=>$width,height=>$height);
     ###titile
     my $desc=$opt{project}.".$opt{suffix}".".".$chr;
     #$svg=title($svg,$desc);
     ###x axis
     #my $rate=$chrlen{$chr}/800;
     my $rate=45000000/800;
     my $end =$chrlen{$chr}/$rate+100;
     print "$end\n";
     my $xpos=$height-80;
     $svg=drawxaxis($svg,100,$end,$xpos,0,$chrlen{$chr},2000000,$rate);
     ###y axis
     my $up=-10;
     my $down=10;
     my $step=2;
     my $ypos=($xpos-100)/2+100;
     my $title=$opt{title} ? $opt{title} : "Number";
     my ($svg,$rate4yaxis)=drawyaxis($svg,100,100,$xpos,$up,$down,$step,$title);
     ### draw Os, upward
     my $svg=drawbar($svg,$table,$ypos,$rate,$rate4yaxis);  
    
     ###write svg 
     my $outfile="$desc.svg";
     writesvg($outfile,$svg);
} ## foreach chromosome

#####
sub log2 {
    my ($n) = shift;
    return log($n)/log(2);
}


########## draw bar ######################
sub drawbar
{
my ($svg,$table,$ypos,$xrate,$yrate)=@_;
open IN, "$table" or die "$!";
while(<IN>){
 chomp $_;
 my @unit=split("\t",$_);
 if ($opt{log2}){
   my $x1=$unit[2]/$xrate+100;
   my $log=log2($unit[3]/$unit[6]);
   my $y1=$log > 0 ? $ypos-$log*$yrate : $ypos;
   my $h1=abs $log*$yrate;
   my $ratiobar=$svg->rectangle(
       x=>$x1, y=>$y1,
       width=>'0.2',height=>$h1,
       style=>{
           fill=>'brown'
       }
   );
 }else{
   ###up
   my $x1=$unit[2]/$xrate+100;
   my $y1=$ypos-$unit[3]*$yrate;
   my $h1=$unit[3]*$yrate;
   my $osbar=$svg->rectangle(
       x=>$x1, y=>$y1,
       width=>'0.2',height=>$h1,
       style=>{
           fill=>'brown'
       }
   ); 
   ###down
   my $x2=$unit[2]/$xrate+100;
   my $y2=$ypos;
   my $h2=$unit[6]*$yrate;
   my $osbar=$svg->rectangle(
       x=>$x2, y=>$y2,
       width=>'0.2',height=>$h2,
       style=>{
           fill=>'brown'
       }
   );
 }
}
close IN;
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

sub drawchr
{
#drawchr($svg,$opt{segment},100,$cpos,$chrlen{$chr},$chr,$rate);
my ($svg,$dir,$x,$y,$chrlen,$chr,$rate)=@_;
=pod
my $chrline=$svg->rectangle(
         x=>$x, y=>$y,
         width=>800,height=>10,
         style=>{
             'stroke-opacity'=>'0'
         }
);
=cut
my @file=glob("$dir/$chr.*.wig.segment.txt");
foreach my $file (@file){
    if ($file=~/$chr.*.wig.segment.txt/){
        open IN, "$file" or die "$!";
        while(<IN>){
           chomp $_;
           next if ($_=~/^$/);
           my @unit=split("\t",$_);
           if ($unit[3] == 1){
              my $hx=$unit[1]/$rate+100;
              my $hw=($unit[2]-$unit[1]+1)/$rate;
              my $chr=$svg->rectangle(
                    x=>$hx, y=>$y,
                    width=>$hw,height=>10,
                    style=>{
                        fill=>'brown'  
                    }
              );
           }
        }
        close IN;
    }
}
return $svg;
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

my $bottomline=$svg->line(
     x1=>$x2,y1=>$y,
     x2=>$x2,y2=>$y+5,
     style=>{stroke=>'black'}
);
my $tail =int ($max/100000)/10;

my $bottomtext=$svg->text(
     x=>$x2,y=>$y+20,
     style=>{'stroke-opacity'=>'0',
         fontsize=>'0.01','text-anchor'=>'start','font-weight'=>100,'stroke-width'=>0.01
     }
)->cdata("$tail Mb");

for(my $i=$min;$i<$max;$i=$i+$step){
     my $tempx=$x1+($i-$min+1)/$rate;
     print "$tempx\t$min\t$step\t$i\t$rate\n";
     #print "$tempx\t$y\n";
     my $line=$svg->line(
         x1=>$tempx,y1=>$y,
         x2=>$tempx,y2=>$y+5,
         style=>{stroke=>'black'}
     );
     my $tempi=int ($i/1000000);
     my $text=$svg->text(
         x=>$tempx+3,y=>$y+20,
         style=>{
             'stroke-opacity'=>'0','fontsize'=>'0.001','text-anchor'=>'end'
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
my $rate=($y2-$y1)/($max-$min);
my $len=$max-$min;
for(my $i=0;$i<=$len;$i=$i+$step){
     my $tempy=$y2-$i*$rate;
     my $line=$svg->line(
         x1=>$x,y1=>$tempy,
         x2=>$x-5,y2=>$tempy,
         style=>{stroke=>'black'}
     );
     my $tempi=$opt{log2} ? $min+$i : abs ($min+$i);
     my $text=$svg->text(
         x=>$x-15,y=>$tempy,
         style=>{
             'stroke-opacity'=>'0','fontsize'=>'0.001','text-anchor'=>'end'
         }
     )->cdata($tempi);
 
}
my $textx=$x-80;
my $texty=$y1+($y2-$y1)/2;
my $rx=$textx+30;
my $ry=$texty;
my $anno=$svg->text(
     x=>$textx,y=>$texty,
     style=>{'stroke-opacity'=>'0',
         fontsize=>'0.01','text-anchor'=>'start','font-weight'=>100,'stroke-width'=>0.01
     },
     transform=> "rotate(270,$rx,$ry)"
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


