use Getopt::Long;
use SVG;

GetOptions (\%opt,"chrlen:s","bar:s","pointtype:s","type:s","project:s","help");


my $help=<<USAGE;
Plot features along chromosomes.
Run: perl plotFeature2chr.pl --chrlen ../input/IRGSP.chrlen --bar ../input/rice_feature_bar --type all --project rice > log 2> log2 &
--chrlen file contains chr length for each chromosome
--pointtype: line/point/rectangle
--bar directory contains bar file for all chromosome
--type draw this type of Feature only, all means draw all in one figure
--project title
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

$opt{pointtype} ||="line";
my $refchr=chrlen($opt{chrlen});
my $type;
if ($opt{type} eq "all"){
   $type="GENE";
}else{
   $type=$opt{type};
}

while(glob("$opt{bar}/*.$type.bar.txt")){
print "$_\n";
my $infile=$_;
my $name;
if ($infile=~/(chr\d+)\.(\w+)\.bar\.txt/){
   $name=$1;
}
my $head=$opt{project}."_".$name."_".$type;
my $svg=SVG-> new (width=>500,heith=>400);
our $xstart=100;
our $xend=400;
our $ystart=100;
our $yend=300;
our $yheight=$yend-$ystart;
our $xlength=$xend-$xstart+4;

my %pos2value;
my @chr;
open IN, "$infile" or die "can not open my gc file";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    $pos2value{$unit[1]}=$unit[2];
    push (@chr,$unit[0]);   
}
close IN;
my $chrlen=$refchr->{$chr[0]};
undef @chr;

my $filename="$head.svg";
my $xcut=10;

my $ystartvalue=0; ##y axis start,0
my $ylength=150; ##y axis length,5
my $ycut=$ylength/10; #/1
my $xnote="Position (Mb)";
my $ynote="$type counts";
my $title=$head;
my $color="red";
my $legandn=1;
my $legandnote=$type;
print "Drawing axis\n";
$svg=axis($svg,$xcut,$ycut,$ystartvalue,$ylength,$chrlen);
print "Drawing note\n";
$svg=note($svg,$xnote,$ynote,$title);
#drawpoint(\@gc);
print "Drawing line\n";
$svg=drawline($svg,\%pos2value,$color,$legandn,$legandnote,$ystartvalue,$ylength,$chrlen) if ($opt{pointtype}=~/line/);
$svg=drawpoint($svg,\%pos2value,$color,$legandn,$legandnote,$ystartvalue,$ylength,$chrlen) if ($opt{pointtype}=~/rec/);
if ($opt{type} eq "all"){
   my $chg="$opt{bar}/$name.TE.bar.txt";
   if (-f $chg){
      my $chgvalue=value($chg);
      $svg=drawline($svg,$chgvalue,"blue","2","TE",$ystartvalue,$ylength,$chrlen);
   }
}

print "Writing file\n";
writefile($svg,$filename);
}#while end

###sub get value
sub value
{
my ($file)=@_;
my %pos2value;
open IN, "$file" or die "can not open my gc file";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    $pos2value{$unit[1]}=$unit[2];
}
close IN;
return \%pos2value;
}


###sub note for figure##
sub note{
    my ($svg,$xnote,$ynote)=@_;
    my $x1=$xstart+($xend-$xstart)/2-30;## think of length of word used for title 
    my $y1=$yend+30;
    my $x2=$xstart-30;
    my $y2=$ystart+($yend-$ystart)/2+30;
    my $xaxis=$svg->text(
       x=>$x1,y=>$y1,
       style=>{
          #stroke=>'black',
          'font-size'=>7
       }
     
    )->cdata($xnote);
    my $yaxis=$svg->text(
       x=>$x2,y=>$y2,
       style=>{
           #stroke=>'black',
           'font-size'=>7
       },
       transform=>"rotate(-90,$x2,$y2)"
    ) ->cdata($ynote);
    return $svg;
}
### sub draw axis#####
sub axis{
my ($svg,$xcut,$ycut,$ystartvalue,$ylength,$chrlen) = @_;
my $xinterval=($xend-$xstart)/$chrlen;
my $yinterval=($yend-$ystart)/$ylength;


###draw x and y lines
my $rec=$svg->rectangle(
         x=>$xstart,y=>$ystart,
         width=>$xlength,height=>$yheight,
         style=>{
             stroke=>'black',
             fill=>'white'
         } 
);

### x annotation##
my $interx=int($chrlen/5000000);
for(my $i=0;$i<=$interx;$i=$i+1){
   my $x=$i*5000000*$xinterval+$xstart;
   my $y=$yend-4;
   my $line=$svg->line(    
       x1 => $x,y1=>$y,
       x2 => $x,y2=>$yend,
       style=>{stroke=>'black'}
   );
   my $xano=$x-5;
   my $yano=$yend+15;
   my $mb=$i*5;
   my $ano=$svg->text(
       x=>$xano,y=>$yano,
       style => { 
          color=>'black',
          'font-size'=>7
          },
       transform=> "rotate(30 $xano $yano)"
   )->cdata($mb); 
}
####y annotation###
my $intery=$ylength/$ycut; 
for (my $i=0;$i<$ylength;$i=$i+$intery){
    my $y=$yend-$i*$yinterval;
    my $x=$xstart+4;
    my $line=$svg->line(
       x1 => $xstart,y1=>$y,
       x2 => $x,y2=>$y,
       style=>{stroke=>'black'}
    );
    my $xano=$xstart-20;
    my $yano=$y+2;
    my $yword=$ystartvalue+$i;
    my $ano=$svg->text(
       x=>$xano,y=>$yano,
       style => { 
           color=>'black',
           'font-size'=>7
           },
       #transform => 'rotate(-5)'
    )->cdata($yword);
}
return $svg;
}


### get chr length hash
sub chrlen
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
   chomp $_;
   my @unit=split("\t",$_);
   $hash{$unit[0]}=$unit[1];
}
close IN;
return \%hash;
}

=pod
###sub drawpoint#### 
sub drawpoint{

my ($gc,$chrlen)=@_;
my $counter;
my $xinterval=($xend-$xstart)/$chrlen;
my $yinterval=($yend-$ystart)/100;
#print "$interval\n";
foreach(@$gc){
     $counter++;
     my $x=$counter*$xinterval+100;
     my $y=$yend-$_*$yinterval*100;
     #print "$x\t$y\n";
     my $point=$svg->circle(cx=>$x,cy=>$y,r=>1); 
}
}
=cut

###sub drawline ####
sub drawpoint{
my ($svg,$value,$color,$legandn,$legandnote,$ystartvalue,$ylength,$chrlen)=@_;
my $xinterval=($xend-$xstart)/$chrlen;
my $yinterval=($yend-$ystart)/$ylength;
my $xlegand1=$xend-90;
my $ylegand1=$ystart+$legandn*20;
my $legand=$svg->rectangle(
     x=>$xlegand1,y=>$ylegand1,
     width=>10,height=>2,
     style=>{'fill'=>$color}

);
my $xlegandt=$xend-70;
my $ylegandt=$ylegand1+2;
my $legandnote=$svg->text(
       x=>$xlegandt,y=>$ylegandt,
       style=>{
             color=>'black',
             'font-size'=>7
          }
)->cdata($legandnote);

my @pos=sort {$a <=> $b} keys %$value;
for(my $i=0;$i<@pos;$i++){
     my $x1=$pos[$i]*$xinterval+$xstart;
     my $y1=$yend-($value->{$pos[$i]}-$ystartvalue)*$yinterval;
     my $yheight=($value->{$pos[$i]}-$ystartvalue)*$yinterval; 
     my $xwidth =50000*$xinterval;
     #my $point=$svg->circle(cx=>$x1,cy=>$y1,r=>1);
     my $rec  =$svg->rectangle(
         x=>$x1,y=>$y1,
         width=>$xwidth,height=>$yheight,
         style=>{
             fill=>$color
         }
     );
}
return $svg;
}

####sub drawline ####
sub drawline{
my ($svg,$value,$color,$legandn,$legandnote,$ystartvalue,$ylength,$chrlen)=@_;
my $xinterval=($xend-$xstart)/$chrlen;
my $yinterval=($yend-$ystart)/$ylength;
my $xlegand1=$xend-90;
my $xlegand2=$xend-80;
my $ylegand1=$ystart+$legandn*20;
my $legand=$svg->line(
     x1=>$xlegand1,y1=>$ylegand1,
     x2=>$xlegand2,y2=>$ylegand1,
     style=>{'stroke'=>$color}

);
my $xlegandt=$xend-70;
my $ylegandt=$ylegand1+2;
my $legandnote=$svg->text(
       x=>$xlegandt,y=>$ylegandt,
       style=>{
             color=>'black',
             'font-size'=>7
          }
)->cdata($legandnote);

my @pos=sort {$a <=> $b} keys %$value;
for(my $i=0;$i<@pos-1;$i++){
     
     my $x1=$pos[$i]*$xinterval+$xstart;
     my $x2=$pos[$i+1]*$xinterval+$xstart;
     my $y1=$yend-($value->{$pos[$i]}-$ystartvalue)*$yinterval;
     my $y2=$yend-($value->{$pos[$i+1]}-$ystartvalue)*$yinterval;
     #print "x1:$x1\ty1:$y1\tx2:$x2\ty2:$y2\n";
     my $line=$svg->line(
          x1=>$x1,y1=>$y1,
          x2=>$x2,y2=>$y2, 
          style=>{'stroke'=>$color}
     );
}
return $svg;
}

####sub writefile#####
sub writefile{
my ($svg,$file)=@_;
open OUT, ">$file";
print OUT $svg->xmlify;
close OUT;
system "/home/jfchen/FFproject/tools/draw/svg2xxx_release/svg2xxx $file -m 1000 -t pdf";
#system "/home/jfchen/bgitraining/draw/svg2xxx_release/svg2xxx $file -t png";
}
