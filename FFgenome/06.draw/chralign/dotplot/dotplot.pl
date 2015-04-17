##Draw dot plot figure using the result fo DAGalign or Mcscan 
##Write by Chen jinfeng on 20090929
##modified by jfchen on 20091009,modify the selfhit for chromosome length 
##Uesage: perl dotplot.pl os_sb, need 4dotplot .chr 

use strict;
use SVG;

our $svg=SVG-> new (width=>500,heith=>1200);
our $xstart=50;
our $xend=450;
our $ystart=50;
our $yend=350;
our $yheight=$yend-$ystart;
our $xlength=$xend-$xstart;
my ($spec1,$spec2)=split("_",$ARGV[0]);
my $infile1="$ARGV[0]"."4dotplot";
my $infile2="$ARGV[0].chr";
######draw x,y lines for figure and write the note

my $rec=$svg->rectangle(
         x=>$xstart,y=>$ystart,
         width=>$xlength,height=>$yheight,
         style=>{
             stroke=>'black',
             fill=>'white',
             'stroke-width'=>0.7
         } 
 
);######draw the lines

my %os2sb;###hash for inter species chromosome relation
my %sb2os;
my %oschr;##store chr length inhash
my %sbchr; 
my $oslength;
my $sblength;
open CHR, "$infile2" or die "can not open my chr file";
my @head=split("\t",<CHR>);
my $yname=$head[0];
my $xname=$head[2];
print "$yname\t$xname\n";
while (<CHR>){
         chomp $_;
         my @unit=split("\t",$_);
         if ($unit[0]=~/(\d+)/){$unit[0]=$1};
         if ($unit[2]=~/(\d+)/){$unit[2]=$1};
         unless(exists $os2sb{$unit[0]}){$os2sb{$unit[0]}=$unit[2] };
         unless(exists $sb2os{$unit[2]}){$sb2os{$unit[2]}=$unit[0] }; 
         unless(exists $oschr{$unit[0]}){
             $oschr{$unit[0]}=$unit[1];
             #print "$unit[0]\n";
             $oslength+=$unit[1];   
         }
         unless(exists $sbchr{$unit[2]}){
             $sbchr{$unit[2]}=$unit[3];
             #print "$unit[2]\n";
             $sblength+=$unit[3];
         }
         
}
close CHR;
our $osratio=$yheight/$oslength;
our $sbratio=$xlength/$sblength;

#print "OS\t$osratio\nSB\t$sbratio\n";
###read in chrinfor and store infro in hash and scalar

##draw lines seperate for chromosome
my @osstart; ## this array will store the start x of every chromomore for y axis
push (@osstart,$ystart); 
my $yv=$ystart;
my @array1=sort keys %oschr;
my $narray1=@array1;
my $lastkeyy=pop @array1;
my $lastyvy=$yend-$oschr{$lastkeyy}*$osratio/2+2;
my $lastyvx=$xstart-10;
my $lastnotey="C$narray1";
my $ylastnote =$svg->text(
           x=>$lastyvx, y=>$lastyvy,
           style=>{'text-anchor'=>'end','font-size'=>'xx-small','font-weight'=>100,'stroke-width'=>0.1}     
)->cdata($lastnotey); ###write note for the last chromosome

my $countery;
foreach (@array1) {
    $countery++;
    my $notey="C$countery";  
    #print "$_"; 
    $yv+=$oschr{$_}*$osratio; 
    push (@osstart,$yv); ## store start x for every chromosome, except last one 
    #print "$yv\n";
    my $yvy=$yv-$oschr{$_}*$osratio/2+2;
    my $yvx=$xstart-10;
    my $ynote =$svg->text(
           x=>$yvx, y=>$yvy,
           style=>{'text-anchor'=>'end','font-size'=>'xx-small','font-weight'=>100,'stroke-width'=>0.1}     
    )->cdata($notey); ###write note for each chromosome except for the last
    my $yy=$yheight/2+$ystart;
    my $yx=$xstart-35;
    my $yano=$svg->text(
       x=>$yx,y=>$yy,
       style=>{'font-size'=>'small','font-style'=>'italic','font-weight'=>100,'stroke-width'=>0.1},
       transform=>"rotate(-90,$yx,$yy)"
    ) ->cdata($yname);## write y anotation for y axis
    my $yline=$svg->line(
           x1=>$xstart, y1=>$yv,
           x2=>$xend,   y2=>$yv,  
           style=>{stroke=>'red','stroke-width'=>0.4}  
    );##draw internal line to seperate chromosome
}

###draw line for x and write note, similar with y
my @sbstart;
push (@sbstart,$xstart);
my $xv=$xstart;
my @array2=sort keys %sbchr;
my $narray2=@array2;
my $lastkeyx=pop @array2; 
my $lastxvy=$yend+20;
my $lastxvx=$xend-$sbchr{$lastkeyx}*$sbratio/2-2;
my $lastnotex="C$narray2";
my $xlastnote=$svg->text(
           x=>$lastxvx, y=>$lastxvy,
           style=>{'font-size'=>'xx-small','font-weight'=>100,'stroke-width'=>0.1}
)->cdata($lastnotex);
my $counterx;
foreach (@array2){
    $counterx++;
    my $notex="C$counterx";
    #print "$_\n"; 
    $xv+=$sbchr{$_}*$sbratio;  
    push (@sbstart,$xv); 
    #print "Sb\t$sbchr{$_}\n";
    #print "P\t$xv\n";
    my $xvx=$xv-$sbchr{$_}*$sbratio/2-2;
    my $xvy=$yend+20;
    my $xnote =$svg->text(
           x=>$xvx, y=> $xvy,
           style=>{'font-size'=>'xx-small','font-weight'=>100,'stroke-width'=>0.1}
    )->cdata($notex);
    my $xx=$xlength/2+$xstart-5;
    my $xy=$yend+40;
    my $xano=$svg->text(
         x=>$xx, y=>$xy,
         style=>{fill=>'black','font-size'=>'small','font-style'=>'italic','font-weight'=>100,'stroke-width'=>0.1},
    )->cdata($xname);
    my $xline=$svg->line(
           x1=>$xv,  y1=>$ystart,
           x2=>$xv,  y2=>$yend,
           style=>{stroke=>'red','stroke-width'=>0.4}
    );
}
#####

#######read in os_sb4dotplot and draw dot in each cell###
$/="##";
my $dotfile="$ARGV[0]"."4dotplot";
open DOT, "$dotfile" or die "can not open my dotfile";
while(<DOT>){
     my @unit=split("\n",$_);
     shift @unit;
     $unit[@unit]=~s/#//g;
     my @word=split("\t",$unit[0]);
     my $sub1=substr($word[1],0,2);
     my $sub2=substr($word[5],0,2);
     my $chromosome1=$word[0];
     my $chromosome2=$word[4];
     #print "$sub1\t$sub2\n"; 
     unless($sub1 eq $sub2){ ##if the blast result is from a self comparition,us if replace unless and delete the elsif steps 
           foreach(@unit){
			     unless($_=~/\w+/){next};
                 my @para=split("\t",$_);
                 my $dotx;
                 if ($para[6]>$para[7]){
                     $dotx=$sbratio*(($para[6]-$para[7])/2+$para[7]);
                 }else{
                     $dotx=$sbratio*(($para[7]-$para[6])/2+$para[6]);
                 } 
                 my $doty;
                 if ($para[2]>$para[3]){
                     $doty=$osratio*(($para[2]-$para[3])/2+$para[3]);
                 }else{
                     $doty=$osratio*(($para[3]-$para[2])/2+$para[2]);
                 }
                 $dotx+=$sbstart[$para[4]-1];
                 $doty+=$osstart[$para[0]-1];
                 my $dot=$svg->circle(cx=>$dotx,cy=>$doty,r=>0.6,'fill'=>'#DAA520');
           }

     }elsif($sub1=~/$spec1/i){
          
          foreach(@unit){
			     unless($_=~/\w+/){next};
                 my @para=split("\t",$_);
                 my $dotx;
                 if ($para[6]>$para[7]){
                     $dotx=$sbratio*(($para[6]-$para[7])/2+$para[7]);
                 }else{
                     $dotx=$sbratio*(($para[7]-$para[6])/2+$para[6]);
                 }
                 my $doty;
                 if ($para[2]>$para[3]){
                     $doty=$osratio*(($para[2]-$para[3])/2+$para[3]);
                 }else{
                     $doty=$osratio*(($para[3]-$para[2])/2+$para[2]);
                 }
				 #print "$_\n";
                 #print "PARA4\t$para[4]\t$oschr{$para[4]}\n";
                 
                 my $rx=$sbchr{$os2sb{$para[4]}}/$oschr{$para[4]}; 
                 $dotx=$dotx*$rx;
                 $dotx+=$sbstart[$os2sb{$para[4]}-1];
                 
                 $doty+=$osstart[$para[0]-1];
                 my $dot=$svg->circle(cx=>$dotx,cy=>$doty,r=>0.6,'fill'=>'black');
          } 
     }elsif($sub1=~/$spec2/i){
          foreach(@unit){
			     unless($_=~/\w+/){next};
                 my @para=split("\t",$_);
                 my $dotx;
                 if ($para[6]>$para[7]){
                     $dotx=$sbratio*(($para[6]-$para[7])/2+$para[7]);
                 }else{
                     $dotx=$sbratio*(($para[7]-$para[6])/2+$para[6]);
                 }
                 my $doty;
                 if ($para[2]>$para[3]){
                     $doty=$osratio*(($para[2]-$para[3])/2+$para[3]);
                 }else{
                     $doty=$osratio*(($para[3]-$para[2])/2+$para[2]);
                 }
                 
                 $dotx+=$sbstart[$para[4]-1];
                 
				 
                 $doty=$doty*$oschr{$sb2os{$para[0]}}/$sbchr{$para[0]};
                 $doty+=$osstart[$sb2os{$para[0]}-1];
                 
                 my $dot=$svg->circle(cx=>$dotx,cy=>$doty,r=>0.6,'fill'=>'red');
          }
     } 




}

close DOT;
$/="\n";
######
my $outfile="$ARGV[0].svg";
writefile($outfile);



####sub writefile#####
sub writefile{
my ($file)=@_;
open OUT, ">$file";
print OUT $svg->xmlify;
close OUT;
system "/home/jfchen/bgitraining/draw/svg2xxx_release/svg2xxx $file -t pdf";
#system "/home/jfchen/bgitraining/draw/svg2xxx_release/svg2xxx $file -t png";
}
