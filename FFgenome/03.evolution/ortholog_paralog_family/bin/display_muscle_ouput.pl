#!/usr/bin/perl

=head1 Name

display_muscle_output.pl --draw figure to display multiple sequence alignment.

=head1 Description

To display alignment and use the svg2xxx.pl to convert figure format.

=head1 Version

Author: Li Jianwen, lijianwen@genomics.org.cn
Version: 1.0,  Date: 2009-3-14
 
				       
=head1 Usage
	--join <T/F>		set whether to join the tree figure, default=F
	--figure_width <num>    set the width of the figure, if join=T, default=1000, else default=500
	--left_margin	<num>	set the left margin, if join=T, default=510, if join=F, default=15
	--right_margin <num>	set the right margin, default=15
	--top_margin <num>	set the top margin, default=15
	--bottom_margin <num>   set the bottom_margin, default=15
	--sequence_width <num>	set the width of each sequence, default=14
	--vertical_skip <num>   set the vertical distance between neighbor sequence, default=35
	--match_color <str>     set the color for the matched region, default="green"
	--font <str>		set the font, default="ArialNarrow"			
	--fsize <num>		set the font size, default=12
	--annotation_file <str> set the annotation file
	--exon_mark <str>	set whether to mark the exon border: none, triangle or line, default=none
	--structure_file <str>	set the structure file, when exon_mark=triangle/line
	--id_file <str>		set the id file to set the sequence order
	--seq_file <str>	set the sequence file in fasta format, every sequence must have the same length. Necessary.
	--type <str>            set the result figure type: png, tiff, jpeg, pdf, default=svg
	--outdir <str>		set the output directory, default=./
	--coordinate		use the vertical coordiante in the id_file
	--help			output help information to screen
	
=head1 Example

perl display_muscle_output.pl --seq_file ../input/12.cds.muscle

=cut

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin";
use SVG;
use SVG::Font;
use Data::Dumper;
use File::Basename qw(dirname basename);

my ($join,$figure_width, $left_margin, $right_margin, $top_margin, $bottom_margin, $sequence_width, $vertical_skip, $match_color, $font, $fsize, $annotation_file, $exon_mark, $structure_file, $id_file, $seq_file, $type, $outdir, $help, $coordinate);
GetOptions(
	"join:s"=>\$join,
	"figure_width:i"=>\$figure_width,
	"left_margin:i"=>\$left_margin,
	"right_margin:i"=>\$right_margin,
	"top_margin:i"=>\$top_margin,
	"bottom_margin:i"=>\$bottom_margin,
	"sequence_width:i"=>\$sequence_width,
	"vertical_skip:i"=>\$vertical_skip,
	"match_color:s"=>\$match_color,
	"font:s"=>\$font,
	"fsize:i"=>\$fsize,
	"annotation_file:s"=>\$annotation_file,
	"exon_mark:s"=>\$exon_mark,
	"structure_file:s"=>\$structure_file,
	"id_file:s"=>\$id_file,
	"seq_file:s"=>\$seq_file,
	"type:s"=>\$type,
	"outdir:s"=>\$outdir,
	"coordinate"=>\$coordinate,
	"help"=>\$help,
);
die `pod2text $0` if (!$seq_file || $help);


$outdir||=".";
$outdir=~s/\/$//;
$join ||="F";
$figure_width ||=($join eq "T" ) ? 1000: 500;
$left_margin ||= ($join eq "T" ) ? 520 : 30;
$right_margin ||= 15;
$top_margin ||=40;
$bottom_margin ||=15;
$sequence_width ||=14;
$vertical_skip ||=34;
$match_color ||="green";
$font ||="ArialNarrow";
$fsize ||=12;
$exon_mark ||="line";
$structure_file ||="";
$id_file ||="";
$seq_file ||="";
$type ||="svg";

my (%genes, @id, $gene_number, $figure_height, $svg, %structure, %annotation, %coordinate);

ReadFasta($seq_file, \%genes);
ReadId($id_file, \@id, \%coordinate) if ($id_file);
ReadStructure($structure_file, \%structure) if ($structure_file);
@id=keys %genes if (!@id);
my $long_text="";
$long_text=ReadAnnotation($annotation_file, \%annotation) if ($annotation_file);

$gene_number=@id;
$figure_height=$vertical_skip*($gene_number-1)+$top_margin+$bottom_margin;
my $width=$figure_width+textWidth($font,$fsize,$long_text);
$svg=SVG->new('width',$width,'height',$figure_height);

my ($begin, $end, $vertical);
my $counter=0;
my $line_length=$figure_width-$left_margin-$right_margin;

foreach my $gene (@id){
	my $length=$genes{$gene}{len};
	my @seq=split//,$genes{$gene}{seq};
	my $unit=$length/$line_length;
	
##draw two lines, the frame of the gene################
	if ($coordinate){
		$vertical=$coordinate{$gene}-$sequence_width/2;
	}else{
		$vertical=$top_margin+($counter)*$vertical_skip-$sequence_width/2;
		print "!!!!!!!!\n";

	}
	$svg->line('x1',$left_margin,'y1',$vertical,'x2',$left_margin+$line_length-1,'y2',$vertical,'stroke',$match_color);
	$svg->line('x1',$left_margin,'y1',$vertical+$sequence_width-1,'x2',$left_margin+$line_length-1,'y2',$vertical+$sequence_width-1,'stroke',$match_color);
	$counter++;
	
##get the first bp of the first matched block##########
	my $i=0;
	for ($i=0; $i<$length; $i++){
		if ($seq[$i]=~/[^-]/){
			$begin=$i;
			last;
		}
	}
	$end=$begin;
	
##Find which block is matched##########################
	for ($i=$begin+1; $i<$length; $i++){
		if ($seq[$i]=~/[^-]/ ){
			if($end==($i-1)){
				$end=$i;
			}else{
				$begin=$i;
				$end=$i;
			}
		}elsif ($end==($i-1)){
			Draw($begin, $end, $unit)
		}
	}
	Draw($begin, $end, $unit) if ($end==$length-1);

## Mark borders of the exons###########################
	if (%structure){
		my $last=0;
		my $ture_length=0;
		my $n=@{$structure{$gene}};
		foreach my $point (@{$structure{$gene}}){
			for (my $j=$last; $j<@seq; $j++){
				$ture_length++ if ($seq[$j]!~/-/);
				if ($ture_length==$point){
					Mark($j, $exon_mark, $unit);
					$last=$j+1;
					$ture_length=0;
					last;
				}
			}
		}
	}

## Add annotation######################################
	if (%annotation){
		$svg->text('x',$figure_width,'y',$vertical+textHeight($fsize),'-cdata',$annotation{$gene},'font-family',$font,'fong-size',$fsize,'fill',"black");
	}
	
}
## Print and convert the figure format#################
my $name = (split/\//, $seq_file)[-1];
open (OUT, ">$outdir/$name.svg") || die $!;
print OUT $svg->xmlify();
close OUT;

if ($type && $type ne "svg"){
	`perl /share/raid1/self-software/biosoft/svg_kit/svg2xxx.pl $outdir/$name.svg --type $type`;
}


############################################################################################
######################### Sub Routines #####################################################
############################################################################################
sub ReadFasta{
	my ($file, $p)=@_;
	
	open (IN, $file) || die "Fail to open seq_file $file\n$!";
	$/=">"; <IN>; $/="\n";
	while (<IN>){
		chomp;
		my $id=$1 if (/^(\S+)/);
		my $title=$_;
		$/=">";
		my $seq=<IN>;
		chomp $seq;
		$/="\n";

		$seq=~s/\s//g;
		$p->{$id}{title}=$title;
		$p->{$id}{seq}=$seq;
		$p->{$id}{len}=length($seq);
	}
	close IN;
}

###########################################################################
sub ReadId{
	my ($file, $p, $c_p)=@_;
	open (IN, $file) || die $!;
	while (<IN>){
		chomp;
		my @c=split;
		push @{$p}, $c[0];
		$c_p->{$c[0]}=$c[1] if ($c[1]);
	}
	close IN;
}

###########################################################################
sub Draw{
	my ($begin, $end, $unit)=@_;
	my $block_width=int(($end-$begin+1)/$unit);
	$begin=int($begin/$unit);
	$svg->rect('x',$left_margin+$begin,'y',$vertical,'width',$block_width,'height',$sequence_width-1,'fill',$match_color,'stroke',$match_color);
}

###########################################################################
sub ReadStructure{
	my ($file, $p)=@_;
	open (IN, $file) || die $!;
	while (<IN>){
		chomp;
		my @c=split;
		my @d=split/;/, $c[4];
		foreach my $tmp (@d){
			my ($begin, $end)=split/,/, $tmp;
			push @{$p->{$c[0]}},$end-$begin+1;
		}
		@{$p->{$c[0]}}=reverse(@{$p->{$c[0]}}) if ($c[2] eq "-");
		pop @{$p->{$c[0]}};
	}
	close IN;
}

##########################################################################
sub Mark{
	my ($position, $mark, $unit)=@_;
	my $color="red";
	$position=int($position/$unit);
	if ($mark=~/triangle/){
		$svg->polygon('points',[$left_margin+$position,$vertical-1,$left_margin+$position-2,$vertical-4,$left_margin+$position+2,$vertical-4],'fill',$color,'stroke',$color);
	}elsif ($mark=~/line/ ){
		$svg->line('x1',$left_margin+$position,'y1',$vertical,'x2',$left_margin+$position,'y2',$vertical+$sequence_width-1,'stroke',$color);
	}
}

##########################################################################
sub ReadAnnotation{
	my ($file, $p)=@_;
	my $longest="";
	open (IN, $file) || die $!;
	while (<IN>){
		chomp;
		my $id=$1 if (/^(\S+)/);
		s/^$id\s+//;
		$p->{$id}=$_;
		$longest=$_ if (length($_)>length($longest));
	}
	close IN;
	return $longest;
}

