use SVG;
##draw gc content for chromosome 
##write by chenjinfeng on 2009.9.25
##y axis from 0-100, version 1.

while(glob("*.gc")=~/(.*)\.gc/){
my $infile="$1\.gc";
my $name=$1;

our $svg=SVG-> new (width=>1000,heith=>800);
our $xstart=100;
our $xend=800;
our $ystart=100;
our $yend=300;
our $yheight=$yend-$ystart;
our $xlength=$xend-$xstart+4;

my @gaps;
my $gapfile="$1\.gaps";
open IN, "$gapfile" or die "$!";
our ($head,$seqlen)=split("\t",<IN>);
while (<IN>){
    chomp $_;
    next if ($_=~/^$/ );
    my @unit=split("\t",$_);
    push (@gaps,[$unit[0],$unit[1]]);
}
close IN;



my @gc;
my @position;
open IN, "$infile" or die "can not open my gc file";
while(<IN>){
    my @unit=split("\t",$_);
	#print "$_\n";
    push (@gc,$unit[1]);
    push (@position,$unit[0]);     
}
my $filename="$name.svg";
my $xcut=10;
my $ycut=5;
my $xnote="Position";
my $ynote="GC Content";
my $title="ChrN";
my $color="red";
my $legandn=1;
my $legandnote="Rice";
axis(\@position,\@gc,$xcut,$ycut);
note($xnote,$ynote,$title);
#drawpoint(\@gc);
drawline(\@gc,$color,$legandn,$legandnote);
drawgaps(@gaps);
writefile($filename);
close IN;
}##while end


sub drawgaps{
my (@gaps)=@_;
my $xinterval=($xend-$xstart)/$seqlen;
for(my $i=0;$i<@gaps;$i++){
    my $x1=$gaps[$i][0]*$xinterval+100;
    my $x2=$gaps[$i][1]*$xinterval+100;
    my $y =$yend-20;
    print "$gaps[$i][0]\t$gaps[$i][1]\t$y\n";
    my $line=$svg->line(
       x1 => $x1,y1=>$ystart,
       x2 => $x1,y2=>$yend,
       style=>{stroke=>'blue'}
    );
    my $line=$svg->line(
       x1 => $x2,y1=>$ystart,
       x2 => $x2,y2=>$yend,
       style=>{stroke=>'blue'}
    );

}
}

###sub note for figure##
sub note{
    my ($xnote,$ynote)=@_;
    my $x1=$xstart+($xend-$xstart)/2-50;## think of length of word used for title 
    my $y1=$yend+50;
    my $x2=$xstart-50;
    my $y2=$ystart+($yend-$ystart)/2+50;
    my $xaxis=$svg->text(
       x=>$x1,y=>$y1,
       style=>{stroke=>'black'}
     
    )->cdata($xnote);
    my $yaxis=$svg->text(
       x=>$x2,y=>$y2,
       style=>{stroke=>'black'},
       transform=>"rotate(-90,$x2,$y2)"
    ) ->cdata($ynote);
}
### sub draw axis#####
sub axis{
my ($position,$gc,$xcut,$ycut) = @_;
my $xinterval=($xend-$xstart)/@$gc;
my $yinterval=($yend-$ystart)/100;


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
my $interx=int(@$gc/$xcut);
for(my $i=1;$i<@$gc;$i=$i+$interx){
   my $x=$i*$xinterval+$xstart;
   my $y=$yend-4;
   my $line=$svg->line(    
       x1 => $x,y1=>$y,
       x2 => $x,y2=>$yend,
       style=>{stroke=>'black'}
   );
   my $xano=$x-10;
   my $yano=$yend+15;

   my $ano=$svg->text(
       x=>$xano,y=>$yano,
       style => { stroke=>'black'},
       transform=> "rotate(30 $xano $yano)"
   )->cdata($$position[$i-1]); 
}
####y annotation###
my $intery=100/$ycut; 
for (my $i=1;$i<100;$i=$i+$intery){
    my $y=$yend-$i*$yinterval;
    my $x=$xstart+4;
    my $line=$svg->line(
       x1 => $xstart,y1=>$y,
       x2 => $x,y2=>$y,
       style=>{stroke=>'black'}
    );
    my $xano=$xstart-20;
    my $yano=$y+2;
    my $yword=$i-1;
    my $ano=$svg->text(
       x=>$xano,y=>$yano,
       style => { stroke=>'black'},
       #transform => 'rotate(-5)'
    )->cdata($yword);

}

}


###sub drawpoint#### 
sub drawpoint{

my ($gc)=@_;
my $counter;
my $xinterval=($xend-$xstart)/@$gc;
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

####sub drawline ####
sub drawline{
#print "in drawline\n";
my ($gc,$color,$legandn,$legandnote)=@_;
my $xinterval=($xend-$xstart)/@$gc;
my $yinterval=($yend-$ystart)/100;
my $xlegand1=$xend-30;
my $xlegand2=$xend-10;
my $ylegand1=$ystart+$legandn*20;

my $legand=$svg->line(
     x1=>$xlegand1,y1=>$ylegand1,
     x2=>$xlegand2,y2=>$ylegand1,
     style=>{'stroke'=>$color}

);
my $xlegandt=$xend-70;
my $legandnote=$svg->text(
       x=>$xlegandt,y=>$ylegand1,
       style=>{'stroke'=>'black'}
)->cdata($legandnote);

for(my $i=1;$i<@$gc;$i++){
     my $x1=$i*$xinterval+100;
     my $x2=($i+1)*$xinterval+100;
     my $y1=$yend-$$gc[$i-1]*$yinterval*100;
     my $y2=$yend-$$gc[$i]*$yinterval*100;
	 #print "$x1\t$y1\t$x2\t$y2\n";
     my $line=$svg->line(
          x1=>$x1,y1=>$y1,
          x2=>$x2,y2=>$y2, 
          style=>{'stroke'=>$color}
     );
}
}

####sub writefile#####
sub writefile{
my ($file)=@_;
open OUT, ">$file";
print OUT $svg->xmlify;
close OUT;
system "/share/raid12/chenjinfeng/tools/draw/svg2xxx_release/svg2xxx $file -t pdf";
system "/share/raid12/chenjinfeng/tools/draw/svg2xxx_release/svg2xxx $file -t png";
}
