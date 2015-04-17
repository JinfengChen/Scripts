#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin $Script);
use lib $Bin;
use File::Basename qw(basename dirname);
use Getopt::Long;
use SVG;
use SVG::Font;

my($X_step,$X_end,$Y_step,$Y_end,$help);
GetOptions(
	"X_step:i"=>\$X_step,
	"X_end:i"=>\$X_end,
	"Y_end:f"=>\$Y_end,
	"Y_step:f"=>\$Y_step,
	"help"=>\$help
);
die "Usage:perl $0 <.gff> <Genomo_size>" if(@ARGV!=2 || $help);
$X_step ||= 1;
$X_end ||= 40;
$Y_step ||= 0.5;
$Y_end ||= 2;

my (@len,$genome_size);
$genome_size=$ARGV[1];
&read_gff($ARGV[0],\@len);

&draw_svg(\@len);
sub read_gff{
	my $file=shift;
	my $p=shift;
	open IN,$file or die "$!";
	while(<IN>)
	{
#			$_=~s/^\s+//;
#			next if(/^\D/);
			next if(/^\#\#/);
			chomp;
			my @c=split/\t/,$_;
			my $div=0;
			if($c[8]=~/PercDiv=(\S+);PercDel/)
			{
					$div=int $1;
			}
			if($c[8]=~/LINE/)
			{
					$p->[$div]{'LINE'}+=abs($c[3]-$c[4]);
			}
			elsif($c[8]=~/SINE/)
			{
					$p->[$div]{'SINE'}+=abs($c[3]-$c[4]);
			}
			elsif($c[8]=~/LTR/)
			{
					$p->[$div]{'LTR'}+=abs($c[3]-$c[4]);
			}
			elsif($c[8]=~/DNA/)
			{
					$p->[$div]{'DNA'}+=abs($c[3]-$c[4]);
			}
			else
			{
					next;
			}
	}
	close IN;
}
sub draw_svg{
	my $p=shift;
	my $svg=SVG->new('width',1200,'height',600);
	my $font_family="ArialNarrow";
	my $font_size= 32;
	my $x=200;
	my $y=100;
	my $width = 16 * $X_end;
	my $y_resolution=200;
	my $height = $y_resolution * $Y_end;
	$svg->rect('x',$x,'y',$y,'width',$width,'height',$height,'fill','none','stroke-width',3,'stroke','black');
	my $X_tag="Sequence divergence rate(%)";
	$svg->text('x',$x+$width/2-textWidth($font_family,$font_size,$X_tag)/2,'y',100+$height+10+textHeight($font_size)*2,'-cdata',$X_tag,'font-family',$font_family,'font-size',$font_size);
	
	my $Y_tag="Percentage of genome";
	my $Y_x=200-15-textWidth($font_family,$font_size,'0.25');
	my $Y_y=100+$height/2+textWidth($font_family,$font_size,$Y_tag)/2;
	#$svg->path('id','left','d',"M $Y_x $Y_y L $Y_x 100",'fill','none');
	#my $path=$svg->text('font-size',$font_size, 'font-family',$font_family,'fill',"black");
	#my $pp=$path->textPath('xlink:href',"#left");
 	#$pp->tspan('-cdata',$Y_tag);
        my $ynote=$svg->text('x',$Y_x,'y',$Y_y,'font-size',$font_size, 'font-family',$font_family,'fill',"black",'-cdata',$Y_tag,'transform',"rotate(-90,$Y_x,$Y_y)");   
 
	for (my $i=0;$i<=$X_end;$i+=$X_step)
	{
			my $x= $i*12 +200;
			my $y= 100+$height;
			if($i%10==0){
			$svg->line('x1',$x,'y1',$y,'x2',$x,'y2',$y-5,'stroke-width',2,'stroke','black');
			$svg->line('x1',$x,'y1',100,'x2',$x,'y2',105,'stroke-width',2,'stroke','black');
			$svg->text('x',$x-textWidth($font_family,$font_size,$i)/2,'y',$y+5+textHeight($font_size),'-cdata',$i,'font-family',$font_family,'font-size',$font_size);
			}
	}
	for (my $j=0;$j<=$Y_end;$j+=$Y_step)
	{
			my $x=200;
			my $y=100+$height-$j*$y_resolution;
			$svg->line('x1',$x,'y1',$y,'x2',203,'y2',$y,'stroke-width',2,'stroke','black');
			$svg->line('x1',$x+$width,'y1',$y,'x2',$x+$width-3,'y2',$y,'stroke-width',2,'stroke','black');
			$svg->text('x',200-5-textWidth($font_family,$font_size,$j),'y',$y+textHeight($font_size)/2,'-cdata',$j,'font-family',$font_family,'font-size',$font_size);
	}
	my $h_width=textHeight($font_size);
	my $r_width=$h_width+textWidth($font_family,$font_size,"LINE")+15;
	my $r_hight=$h_width *4 + 20;

################ Tag ###########################################
	$svg->rect('x',200+$width-$r_width-5,'y',103,'height',$r_hight,'width',$r_width,'fill','none','stroke-width',2,'stroke','black');
	$svg->rect('x',200+$width-$r_width,'y',108,'width',$h_width,'height',$h_width,'fill','yellow','stroke-width',2,'stroke','black');
	$svg->text('x',200+$width-$r_width+5+$h_width,'y',105+$h_width,'-cdata',"DNA",'font-family',$font_family,'font-size',$font_size);
	$svg->rect('x',200+$width-$r_width,'y',108+$h_width+5,'width',$h_width,'height',$h_width,'fill','red','stroke-width',2,'stroke','black');
	$svg->text('x',200+$width-$r_width+5+$h_width,'y',105+$h_width*2+5,'-cdata',"LINE",'font-family',$font_family,'font-size',$font_size);
	$svg->rect('x',200+$width-$r_width,'y',108+$h_width*2+10,'width',$h_width,'height',$h_width,'fill','green','stroke-width',2,'stroke','black');
	$svg->text('x',200+$width-$r_width+5+$h_width,'y',105+$h_width*3+10,'-cdata',"LTR",'font-family',$font_family,'font-size',$font_size);
	$svg->rect('x',200+$width-$r_width,'y',108+$h_width*3+15,'width',$h_width,'height',$h_width,'fill','blue','stroke-width',2,'stroke','black');
	$svg->text('x',200+$width-$r_width+5+$h_width,'y',105+$h_width*4+15,'-cdata',"SINE",'font-family',$font_family,'font-size',$font_size);

	foreach my $i(0..$X_end-1)
	{
		my $per_DNA=$$p[$i]{'DNA'}/$genome_size*100;
		my $per_SINE=$$p[$i]{'SINE'}/$genome_size*100;
		my $per_LINE=$$p[$i]{'LINE'}/$genome_size*100;
		my $per_LTR=$$p[$i]{'LTR'}/$genome_size*100;
		my $x=$i*16+200;
		my $w=8;
		my $y=100+$height-$per_SINE*$y_resolution;
		$svg->rect('x',$x,'y',$y,'width',$w,'height',$per_SINE*$y_resolution,'fill','blue','stroke-width',2);
		$y-=$per_LTR*$y_resolution;
		$svg->rect('x',$x,'y',$y,'width',$w,'height',$per_LTR*$y_resolution,'fill','green','stroke-width',2);
		$y-=$per_LINE*$y_resolution;
		$svg->rect('x',$x,'y',$y,'width',$w,'height',$per_LINE*$y_resolution,'fill','red','stroke-width',2);
		$y-=$per_DNA*$y_resolution;
		$svg->rect('x',$x,'y',$y,'width',$w,'height',$per_DNA*$y_resolution,'fill','yellow','stroke-width',2);
			
	}
	open OUT,">$ARGV[0].TEdivergence.svg" or die "$!";
	print OUT $svg->xmlify();
	close OUT;
}
