#!/usr/bin/perl
use Getopt::Long;
use SVG;

GetOptions (\%opt,"refseq:s","refgenegff:s","reftegff:s","qryseq:s","qrygenegff:s","qrytegff:s","compare:s","type:s","project:s","help");


my $help=<<USAGE;
Dotplot 2 region with gene and TE annotation.
perl $0 --refseq --qryseq --refgenegff --reftegff --qrygenegff --qrytegff --compare --type --project 
--compare: compare file, nummer coord/blastm8/4ACT
--type: coord/blastm8/4ACT
if --compare or --type is not specified, we use blast to compare the two sequence.
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

$opt{project} ||="nucmer";
#`nucmer --coords -p $opt{project} --nosimplify $opt{refseq} $opt{qryseq} > log 2> log2 &`;

unless (defined $opt{type} or defined $opt{compare}){
   `formatdb -i $opt{qryseq} -p F`;
   `blastall -p blastn -i $opt{refseq} -d $opt{qryseq} -e 1e-5 -m 8 -o compare.blast`;
   $opt{type}="blastm8";
   $opt{compare}="compare.blast";
   `rm *.nhr *.nin *.nsq`;
}

our $height=800;
our $width=800;
my $svg=SVG->new(width=>$height,height=>$width);

### common value
my $reflen=getfastalen($opt{refseq});
my $qrylen=getfastalen($opt{qryseq});
my $maxlen=$reflen > $qrylen ? $reflen : $qrylen;
my $rate  =$maxlen/(800-100-50);
my $step=10000; ## 10kb

print "Reflen:$reflen\nQrylen:$qrylen\n";

### X axis
my $xx1=100;
my $xx2=$width-50;
my $xy =$height-50;
$svg=drawxaxis($svg,$xx1,$xx2,$xy,0,$reflen,$step,$rate);

### Y axis
my $yx=$xx2;
my $yy1=$xy;
my $yy2=$xy-$qrylen/$rate;
$svg=drawyaxis($svg,$yy1,$yy2,$yx,0,$qrylen,$step,$rate);

### box
$svg=drawbox($svg,$xx1,$xx2,$yy1,$yy2);

### draw gene on X axis
my $refgenegff=parseGFF($opt{refgenegff});
$svg=drawXgene($svg,$refgenegff,$rate,$xx1,$xx2,$yy2);
### draw te on X axis
if (defined $opt{$opt{reftegff}}){
my $reftegff=parseTEGFF($opt{reftegff});
$svg=drawXTE($svg,$reftegff,$rate,$xx1,$yy2);
}
### draw gene on Y axis
my $qrygenegff=parseGFF($opt{qrygenegff});
$svg=drawYgene($svg,$qrygenegff,$rate,$yy1,$yy2,$xx1);
### draw te on Y axis
if (defined $opt{$opt{qrytegff}}){
my $qrytegff=parseTEGFF($opt{qrytegff});
$svg=drawYTE($svg,$qrytegff,$rate,$xx1,$yy1);
}

###
if ($opt{type} =~/ACT/){
   my $compare=parseACT($opt{compare});
   $svg=drawdotplot($svg,$compare,$rate,$xx1,$yy1);
}elsif($opt{type} =~/coord/){
   my $compare=parsecoord($opt{compare});
   $svg=drawdotplot($svg,$compare,$rate,$xx1,$yy1);
}elsif($opt{type} =~/blastm8/){
   my $compare=parseblastm8($opt{compare});
   $svg=drawdotplot($svg,$compare,$rate,$xx1,$yy1);
}else{
   print "The type of file $opt{compare} is unknown, use --type\n";
}

my $outfile="$opt{project}.svg";
writesvg($outfile,$svg);

#############
sub parseACT
{
my ($file)=@_;
my @com;
open IN, "$file" or die "$!";
while(<IN>){
   chomp $_;
   my @unit=split(" ",$_);
   push (@com,[$unit[2],$unit[3],$unit[5],$unit[6]]);
}
close IN;
return \@com;
}

#############
sub parsecoord
{
my ($file)=@_;
my @com;
open IN, "$file" or die "$!";
while(<IN>){
   chomp $_;
   next unless ($_=~/^\s*\d+/);
   my @unit=split(" ",$_);
   push (@com,[$unit[1],$unit[2],$unit[4],$unit[5]]);
}
close IN;
return \@com;
}


#############
sub parseblastm8
{
my ($file)=@_;
my @com;
open IN, "$file" or die "$!";
while(<IN>){
   chomp $_;
   my @unit=split("\t",$_);
   push (@com,[$unit[6],$unit[7],$unit[8],$unit[9]]);
}
close IN;
return \@com;
}


sub drawdotplot
{
my ($svg,$compare,$rate,$x,$y)=@_;
foreach my $match (@$compare){
   my $x1=$match->[0]/$rate+$x;
   my $y1=$y-$match->[2]/$rate;
   my $x2=$match->[1]/$rate+$x;
   my $y2=$y-$match->[3]/$rate;
   #print "$x1\t$x2\t$y1\t$y2\n";
   my $line=$svg->line(
     x1=>$x1,y1=>$y1,
     x2=>$x2,y2=>$y2,
     style=>{stroke=>'black'}
   );   
}
return $svg;
}


sub drawXgene
{
my ($svg,$refgenegff,$rate,$x1,$x2,$y)=@_;
my $strandline=$svg->line(
     x1=>$x1,y1=>$y-30,
     x2=>$x2,y2=>$y-30,
     style=>{stroke=>'black'}
);
foreach my $g (keys %$refgenegff){
    my @line=split("\n",$refgenegff->{$g});
    my @pos;
    my $strand;
    foreach my $e (@line){
        my @unit=split("\t",$e);
        if ($unit[2] eq "mRNA"){
           $strand=$unit[6];
        }else{
           push (@pos,[$unit[3],$unit[4]]);
        }
    }
    @pos=sort {$a->[0] <=> $b->[1]} @pos;
    my $gstart=$pos[0][0]/$rate+$x1;
    my $gend  =$pos[$#pos][1]/$rate+$x1; 
    if ($strand eq "+"){
       my $geneid=$svg->text(
          x=>$gstart,y=>$y-50,
          style=>{
             'font-size'=>'50%','text-anchor'=>'start','font-weight'=>'100'
          }
       )->cdata("$g");

       my $geneline=$svg->line(
          x1=>$gstart,y1=>$y-36,
          x2=>$gend,y2=>$y-36,
          style=>{stroke=>'brown'}
       );
       foreach my $e (sort {$a->[0] <=> $b->[1]} @pos){
           my $start=$e->[0]/$rate+$x1;
           my $elen =($e->[1]-$e->[0]+1)/$rate;
           my $exony=$y-40;
           my $exon=$svg->rectangle(
              x=>$start, y=>$exony,
              width=>$elen,height=>8,
              style=>{
                fill=>'brown'
              }
           );
       }   
    }else{
       my $geneid=$svg->text(
          x=>$gstart,y=>$y-10,
          style=>{
             'font-size'=>'50%','text-anchor'=>'start','font-weight'=>'100'
          }
       )->cdata("$g");

       my $geneline=$svg->line(
          x1=>$gstart,y1=>$y-24,
          x2=>$gend,y2=>$y-24,
          style=>{stroke=>'brown'}
       );
       foreach my $e (sort {$a->[0] <=> $b->[1]} @pos){
           my $start=$e->[0]/$rate+$x1;
           my $elen =($e->[1]-$e->[0]+1)/$rate;
           my $exony=$y-28;
           my $exon=$svg->rectangle(
              x=>$start, y=>$exony,
              width=>$elen,height=>8,
              style=>{
                fill=>'brown'
              }
           );
       }
    }
}
return $svg;
}
#####
sub drawXTE
{
my ($svg,$reftegff,$rate,$x1,$y)=@_;
foreach my $te (keys %$reftegff){
    my @line=split("\t",$reftegff->{$te});
    my $gstart=$line[3]/$rate+$x1;
    my $gend  =$line[4]/$rate+$x1;
    my $strand=$line[6];
    my $type=$1 if ($line[8]=~/Class=(.*?);/);
    $type=~s/DNA\///;
    $type=~s/LTR\///;
    #print "$te\t$y\t$gstart\t$gend\t$strand\t$type\n";
    if ($strand eq "+"){
       my $geneid=$svg->text(
          x=>$gstart,y=>$y-50,
          style=>{
             'font-size'=>'50%','text-anchor'=>'start','font-weight'=>'100'
          }
       )->cdata("$type");
              my $qleft =$gstart;
              my $qright=$gend;
              my $tright=$gend;
              my $tleft =$gstart;
              my $qheight=$y-40;
              my $theight=$y-32;
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
                        fill=>'gray'
                     }
              );

    }else{
       my $geneid=$svg->text(
          x=>$gstart,y=>$y-10,
          style=>{
             'font-size'=>'50%','text-anchor'=>'start','font-weight'=>'100'
          }
       )->cdata("$type");
              my $qleft =$gstart;
              my $qright=$gend;
              my $tright=$gend;
              my $tleft =$gstart;
              my $qheight=$y-28;
              my $theight=$y-20;
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
                        fill=>'gray'
                     }
              );      
    }
}
return $svg;
}


####
sub drawYgene
{
my ($svg,$refgenegff,$rate,$y1,$y2,$x)=@_;
my $strandline=$svg->line(
     x1=>$x-30,y1=>$y1,
     x2=>$x-30,y2=>$y2,
     style=>{stroke=>'black'}
);
foreach my $g (keys %$refgenegff){
    my @line=split("\n",$refgenegff->{$g});
    my @pos;
    my $strand;
    foreach my $e (@line){
        my @unit=split("\t",$e);
        if ($unit[2] eq "mRNA"){
           $strand=$unit[6];
        }else{
           push (@pos,[$unit[3],$unit[4]]);
        }
    }
    @pos=sort {$a->[0] <=> $b->[1]} @pos;
    my $gstart=$y1-$pos[0][0]/$rate;
    my $gend  =$y1-$pos[$#pos][1]/$rate;
    #print "$g\t$pos[0][0]\t$pos[$#pos][1]\t$gstart\t$gend\n";
    if ($strand eq "+"){
       my $rx=$x-50;
       my $ry=$gstart;
       my $geneid=$svg->text(
          y=>$gstart,x=>$x-50,
          style=>{
             'font-size'=>'50%','text-anchor'=>'start','font-weight'=>'100'
          },
          transform=> "rotate(270,$rx,$ry)"
       )->cdata("$g");
       my $geneline=$svg->line(
          y1=>$gstart,x1=>$x-36,
          y2=>$gend,x2=>$x-36,
          style=>{stroke=>'brown'}
       );
       foreach my $e (sort {$a->[0] <=> $b->[1]} @pos){
              my $qleft =$x-40;
              my $qright=$x-32;
              my $tright=$x-32;
              my $tleft =$x-40;
              my $qheight=$y1-$e->[1]/$rate;
              my $theight=$y1-$e->[0]/$rate;
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
                        fill=>'brown' 
                     }
              );

       }
    }else{
       my $rx=$x-10;
       my $ry=$gstart;
       my $geneid=$svg->text(
          y=>$gstart,x=>$x-10,
          style=>{
             'font-size'=>'50%','text-anchor'=>'start','font-weight'=>'100'
          },
          transform=> "rotate(270,$rx,$ry)"
       )->cdata("$g");
       my $geneline=$svg->line(
          y1=>$gstart,x1=>$x-24,
          y2=>$gend,x2=>$x-24,
          style=>{stroke=>'brown'}
       );
       foreach my $e (sort {$a->[0] <=> $b->[1]} @pos){
              my $qleft =$x-28;
              my $qright=$x-20;
              my $tright=$x-20;
              my $tleft =$x-28;
              my $qheight=$y1-$e->[1]/$rate;
              my $theight=$y1-$e->[0]/$rate;
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
                        fill=>'brown'
                     }
              );

       }
    }
}
return $svg;
}
#####
sub drawYTE
{
my ($svg,$reftegff,$rate,$x,$y)=@_;
foreach my $te (keys %$reftegff){
    my @line=split("\t",$reftegff->{$te});
    my $gstart=$y-$line[3]/$rate;
    my $gend  =$y-$line[4]/$rate;
    my $strand=$line[6];
    my $type=$1 if ($line[8]=~/Class=(.*?);/);
    $type=~s/DNA\///;
    $type=~s/LTR\///;
    #print "$te\t$y\t$gstart\t$gend\t$strand\t$type\n";
    if ($strand eq "+"){
       my $rx=$x-50;
       my $ry=$gstart;
       my $geneid=$svg->text(
          y=>$gstart,x=>$x-50,
          style=>{
             'font-size'=>'50%','text-anchor'=>'start','font-weight'=>'100'
          },
          transform=> "rotate(270,$rx,$ry)"
       )->cdata("$type");
              my $qleft =$x-40;
              my $qright=$x-32;
              my $tright=$x-32;
              my $tleft =$x-40;
              my $qheight=$gend;
              my $theight=$gstart;
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
                        fill=>'gray'
                     }
              );

    }else{
       my $rx=$x-10;
       my $ry=$gstart;
       my $geneid=$svg->text(
          y=>$gstart,x=>$x-10,
          style=>{
             'font-size'=>'50%','text-anchor'=>'start','font-weight'=>'100'
          },
          transform=> "rotate(270,$rx,$ry)"
       )->cdata("$type");
              my $qleft =$x-28;
              my $qright=$x-20;
              my $tright=$x-20;
              my $tleft =$x-28;
              my $qheight=$gend;
              my $theight=$gstart;
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
                        fill=>'gray'
                     }
              );
    }
}
return $svg;
}


###
sub getfastalen
{
$/=">";
my %hash;
my $len;
my ($file)=@_;
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my @temp1=split(" ",$temp);
    my $head=$temp1[0];
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    $len=length $seq;
}
$/="\n";
return $len;
}


####
sub drawbox
{
my ($svg,$x1,$x2,$y1,$y2)=@_;
my $hline=$svg->line(
     x1=>$x1,y1=>$y2,
     x2=>$x2,y2=>$y2,
     style=>{stroke=>'black'}
);
my $vline=$svg->line(
     x1=>$x1,y1=>$y1,
     x2=>$x1,y2=>$y2,
     style=>{stroke=>'black'}
);
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
my $tail =int ($max/1000);
=pod
my $bottomtext=$svg->text(
     x=>$x2,y=>$y+20,
     style=>{
         'font-size'=>'70%','text-anchor'=>'start','font-weight'=>'100'
     }
)->cdata("$tail kb");
=cut

for(my $i=$min;$i<$max;$i=$i+$step){
     my $tempx=$x1+($i-$min+1)/$rate;
     #print "$tempx\t$min\t$step\t$i\t$rate\n";
     #print "$tempx\t$y\n";
     my $line=$svg->line(
         x1=>$tempx,y1=>$y,
         x2=>$tempx,y2=>$y+5,
         style=>{stroke=>'black'}
     );
     my $tempi=int ($i/1000);
     my $text=$svg->text(
         x=>$tempx+3,y=>$y+20,
         style=>{
             'font-size'=>'70%','text-anchor'=>'end','font-weight'=>'100'
         }
     )->cdata($tempi);
}
return $svg;
}

#########
sub drawyaxis
{
my ($svg,$y1,$y2,$x,$min,$max,$step,$rate)=@_;
#print "$x1\t$x2\t$y\n";
my $yaxis=$svg->line(
     x1=>$x,y1=>$y1,
     x2=>$x,y2=>$y2,
     style=>{stroke=>'black'}
);

my $bottomline=$svg->line(
     x1=>$x,y1=>$y2,
     x2=>$x+5,y2=>$y2,
     style=>{stroke=>'black'}
);
my $tail =int ($max/1000);
=pod
my $bottomtext=$svg->text(
     x=>$x,y=>$y2+20,
     style=>{
         'font-size'=>'70%','text-anchor'=>'start','font-weight'=>'100'
     }
)->cdata("$tail kb");
=cut
for(my $i=$min;$i<=$max;$i=$i+$step){
     my $tempy=$y1-($i-$min+1)/$rate;
     #print "$tempy\t$min\t$step\t$i\t$rate\n";
     #print "$tempx\t$y\n";
     my $line=$svg->line(
         x1=>$x,y1=>$tempy,
         x2=>$x+5,y2=>$tempy,
         style=>{stroke=>'black'}
     );
     my $tempi=int ($i/1000);
     my $text=$svg->text(
         x=>$x+20,y=>$tempy+3,
         style=>{
             'font-size'=>'70%','text-anchor'=>'end','font-weight'=>'100'
         }
     )->cdata($tempi);
}
return $svg;
}
#####
sub parseGFF
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $seq;   ##Scaffold
my $id;    ##ID for element
my $record;##all line for this record
my $index; ##key, Seq_id
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/mRNA/){
        $seq=$unit[0];
        if ($unit[8]=~/ID=(.*?);/ or $unit[8] =~/ID=(.*)/){
            $id=$1;
        }
        $record="$_\n";
        $hash{$id}=$record;
    }elsif($unit[0] eq $seq and $unit[8] =~ /Parent=$id/){
        $hash{$id}.="$_\n";
    }

}
close IN;
return \%hash;
}

#####
sub parseTEGFF
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $seq;   ##Scaffold
my $id;    ##ID for element
my $record;##all line for this record
my $index; ##key, Seq_id
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[8]=~/ID=(.*?);/ or $unit[8] =~/ID=(.*)/){
        $id=$1;
        $hash{$id}="$_";
    }

}
close IN;
return \%hash;
}


################################### sub for write svg to file
sub writesvg {
my ($file,$svg)=@_;
#print "$file\n";
open OUT, ">$file" or die "can not open my file";
       print OUT $svg->xmlify;
close OUT;
       system "/home/jfchen/FFproject/tools/draw/svg2xxx_release/svg2xxx $file -t pdf";
}


