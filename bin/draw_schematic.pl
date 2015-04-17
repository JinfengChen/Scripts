#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);
use SVG;

GetOptions (\%opt,"project:s","help");


my $help=<<USAGE;
perl draw_schematic.pl --project test

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}



my $svg=SVG->new(width=>800,height=>600);


####Inversion within curve
my $width=30;
###exterior arc for inversion
my $x1=400;
my $y1=200;
my $x2=400;
my $y2=400;
my $rx1=($y2-$y1)/2;
my $ry1=$rx1;
###interior arc for inversion
my $x3=$x1;
my $x4=$x2;
my $y3=$y1+$width;
my $y4=$y2-$width;
my $rx2=($y4-$y3)/2;
my $ry2=$rx2;
###
my $string = "M$x1,$y1 A$rx1,$ry1 0 0,1 $x2,$y2 L$x4,$y4 A$rx2,$ry2 0 0,0 $x3,$y3 L$x1,$y1";
#my $string = "M$x1,$y1 A$rx1,$ry1 0 0,1 $x2,$y2 L$x4,$y4 ";
print $string,"\n";

my $tag = $svg->path(
        d => $string,
        style => {
            'fill' => 'gray',
            'stroke'=> 'black'
            #'fill'   => 'green',
        }
    );


my $outfile="$opt{project}.svg";
writesvg($outfile,$svg);



##################
sub readtable
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash{$unit[0]}=1;
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
       system "/rhome/cjinfeng/software/tools/draw/svg2xxx_release/svg2xxx $file -t pdf -m 700";
}


 
