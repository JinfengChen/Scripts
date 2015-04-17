#!/usr/bin/perl

=head1 Program Description

Map all the annotions in a wanted region of a sequence, the input file format should be .psl or .gff.

Note that the number of prediction files can be varied (>=0). Genes are displayed with exon-intron 
structure. While the other features (such as TEs,ncRNAs) are displayed just in block. The number of feature file can also 
be varied (>=0).


=head1 Contact & Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 2.0,  Date: 2007-12-21

=head1 Command-line Option

  perl map_chromosome.pl <options>
  --sequence <str>      specify the sequence wanted
  --startpos <num>      specify the starting positon on sequence
  --endpos <num>        speicfy the ending positon on seqeunce
  
  --fasta <str>         specify the genome fasta file

  --gene <file>         specify predicted gene input file
  --genemark <str>      specify mark for prediction
  --feature <file>   specify transposon file
  --featmark <str>      specify mark for transposon
  
  --figurerate <float>  specify the figure rate, default=0.01, 1 dot for 100bp
  --figuredir <str>     specify path of the result figures
  --verbose             output verbose information to screen  
  --help                output help information to screen  

=head1 Usage Exmples

  perl ../bin/map_chromosome.pl --sequence chr1  --startpos 1000001  --endpos  2000000 --gene ../input/test_3chrs.fa.bgf.gff --genemark bgf  --gene ../input/test_3chrs.fa.genscan.gff --genemark genscan  --feature ../input/test_3chrs.fa.RepeatMasker.out.gff --featmark TE --feature ../input/test_3chrs.fa.trf.dat.gff --featmark TRF --figuredir ./chr1

=cut

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);
use File::Path;  ## function " mkpath" and "rmtree" deal with directory
use lib "$Bin/../lib";
use SVG;
use SVG::Font;


my $Left_margin = 30;
my $Right_margin = 30;
my $Vert_skip = 20;   ## vertical skip distance 
my $Top_margin = 30;
my $Bottom_margin = 30;
my $Font_family = 'ArialNarrow';
my $Font_size = 18;

my @Color = ('black','blue','green','red','orange','brown','yellow','black','blue','green','red','orange','brown','yellow');
my $Color_id = 0;	


my (@Gene_file,@Gene_mark,@Feature_file,@Feature_mark,$Max_mark,$Fasta_file);
my ($Sequence,$Startpos,$Endpos);
my ($Fig_rate,$Fig_dir,$Verbose,$Help);
GetOptions(
	"sequence:s"=>\$Sequence,
	"startpos:n"=>\$Startpos,
	"endpos:n"=>\$Endpos,
	
	"fasta:s"=>\$Fasta_file,

	"gene:s"=>\@Gene_file,
	"genemark:s"=>\@Gene_mark,
	"feature:s"=>\@Feature_file,
	"featmark:s"=>\@Feature_mark,
	
	"figurerate:f"=>\$Fig_rate,
	"figuredir:s"=>\$Fig_dir,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Fig_dir ||= "./map_chr_figs";
$Fig_dir =~ s/\/$//;
$Startpos ||= 1;
$Endpos ||= 1000000;
$Fig_rate ||= 0.02; ## 0.01 means 1 dot stand for 100bp
die `pod2text $0` if (! defined $Sequence || $Help);

for (my $i=0; $i<@Gene_file; $i++) {
	$Gene_mark[$i] ||= "prediction$i";
	my $pre_mark = textWidth($Font_family,$Font_size,$Gene_mark[$i]);
	$Max_mark = $pre_mark if($pre_mark > $Max_mark);
}
for (my $i=0; $i<@Feature_file; $i++) {
	$Feature_mark[$i] ||= "transposon$i";
	my $te_mark = textWidth($Font_family,$Font_size,$Feature_mark[$i]);
	$Max_mark = $te_mark if($te_mark > $Max_mark);
}
die `pod2text $0` if (!$Sequence || $Help);

my $Main_width = ($Endpos - $Startpos) * $Fig_rate ; 
my $Mark_width = $Max_mark + 10;
my $Fig_width = $Main_width + $Mark_width + $Left_margin + $Right_margin;
my $Fig_height = $Top_margin + $Bottom_margin + $Vert_skip * ( scalar(@Gene_file) + scalar(@Feature_file) + 1);


my @Gene;
my @Feature;

foreach my $pre_file (@Gene_file) {
	my %Gene;
	read_gene_psl($pre_file,\%Gene) if($pre_file=~/\.psl/);
	read_gene_gff($pre_file,\%Gene) if($pre_file=~/\.gff/);
	push @Gene,\%Gene;
}

foreach my $feature_file (@Feature_file) {
	my %Feature;
	read_feature_gff($feature_file,\%Feature) if($feature_file=~/\.gff/);
	push @Feature, \%Feature;
}

mkpath($Fig_dir) unless(-d $Fig_dir);

my $svg = SVG->new('width',$Fig_width,'height',$Fig_height);
my $draw_Y = $Fig_height - $Bottom_margin;

## draw ruler
plot_ruler($svg,$draw_Y,$Left_margin+$Mark_width,$Fig_width-$Right_margin,$Startpos,$Endpos,'Mb','left');
$draw_Y -= $Vert_skip;

## draw GC content curve
if (defined $Fasta_file) {
	my %seq_fasta;
	Read_fasta($Fasta_file,\%seq_fasta);
	my $Window = 200;
	my $Step = 200;
	my $seq = $seq_fasta{$Sequence}{seq};
	my %GC_content;
	for (my $i=$Window/2; $i+$Window/2<$seq_fasta{$Sequence}{len}; $i+=$Step) {
		my $str = substr($seq,$i-$Window/2,$Window);
		my $gc_rate = gc_content($str);
		$GC_content{$i} = $gc_rate;
	}
	
	plot_GC($svg,\%GC_content,$draw_Y,$Vert_skip);
	
	$svg->text('x',$Left_margin,'y',$draw_Y+textHeight($Font_size),'fill','#000000','-cdata',"GC%",'font-size',$Font_size, 'font-family',$Font_family);
	$draw_Y -= $Vert_skip;
}

## draw feature sets
for (my $i=0; $i<@Feature_file; $i++) {
	my $pre_p = $Feature[$i];
	if (exists $pre_p->{$Sequence}) {
		my $pre_chr_p = $pre_p->{$Sequence};
		foreach my $pre_feature (sort keys %$pre_chr_p) {
			my $pre_feature_p = $pre_chr_p->{$pre_feature};
			my $feature_start = $Left_margin + $Mark_width + ($pre_feature_p->{start} - $Startpos) * $Fig_rate;
			my $feature_end = $Left_margin + $Mark_width + ($pre_feature_p->{end} - $Startpos) * $Fig_rate;
		
			draw_feature($svg,"$Feature_mark[$i]: $pre_feature",$feature_start,$feature_end,$draw_Y,$pre_feature_p->{strand},$Color[$Color_id]);
		}
	}	
	$svg->text('x',$Left_margin,'y',$draw_Y+textHeight($Font_size),'fill','#000000','-cdata',$Feature_mark[$i],'font-size',$Font_size, 'font-family',$Font_family);
	$draw_Y -= $Vert_skip;
	$Color_id++;
}


## draw gene sets
for (my $i=0; $i<@Gene_file; $i++) {
	my $pre_p = $Gene[$i];
	if (exists $pre_p->{$Sequence}) {
		my $pre_chr_p = $pre_p->{$Sequence};
		foreach my $pre_gene (sort keys %$pre_chr_p) {
			my $pre_gene_p = $pre_chr_p->{$pre_gene};
			my @xary;
			foreach my $p (@{$pre_gene_p->{exon}}) {
				my $exon_start = $Left_margin + $Mark_width + ($p->[0] - $Startpos) * $Fig_rate;
				my $exon_end = $Left_margin + $Mark_width + ($p->[1] - $Startpos) * $Fig_rate;
				push @xary,$exon_start,$exon_end;
			}
			draw_gene($svg,"$Gene_mark[$i]: $pre_gene",\@xary,$draw_Y,$pre_gene_p->{strand},$Color[$Color_id]);
			
		}
	}	
	$svg->text('x',$Left_margin,'y',$draw_Y+textHeight($Font_size),'fill','#000000','-cdata',$Gene_mark[$i],'font-size',$Font_size, 'font-family',$Font_family);
	$draw_Y -= $Vert_skip;
	$Color_id++;
}



open OUT,">$Fig_dir/$Sequence\_$Startpos\_$Endpos.svg" || die "fail to creat figure";
print OUT $svg->xmlify();
close OUT;

`/share/raid1/self-software/biosoft/svg_kit/svg2xxx.pl $Fig_dir/$Sequence\_$Startpos\_$Endpos.svg`;

####################################################
################### Sub Routines ###################
####################################################


##return GC% of a given string
#############################################
sub gc_content{
	my $str = shift;
	my $GC = $str =~ tr/GgCc//;
	my $all = $str =~ tr/GgCcAaTt//;
	my $rate = $GC / $all;
	return $rate;
}

##plot GC% part
#############################################
sub plot_GC {
	my $svg = shift;
	my $hash_p = shift;
	my $y_coor = shift;
	my $gc_figure_heigth = shift;

	$y_coor -= $gc_figure_heigth/2;
	my $line_color = "red";
	my $gc_max_rate = 1; ## GC range: 0 to 1
	my $gc_resolution = $gc_figure_heigth / $gc_max_rate;
	
	##draw a frame for the gc part
	##$svg->rect('x',$Left_margin, 'y',$y_coor,'width',scalar(@$ary_p),'height',$gc_figure_heigth,'stroke','black','fill','none');

	my @pos = sort {$a <=> $b} keys %$hash_p;
	##$Left_margin + $Mark_width + ($p->[0] - $Startpos) * $Fig_rate;

	for (my $i=0; $i<@pos-1; $i++) {
		my $x1_coor = $Left_margin + $Mark_width + ($pos[$i] - $Startpos) * $Fig_rate;
		my $x2_coor = $Left_margin + $Mark_width + ($pos[$i+1] - $Startpos) * $Fig_rate;
		my $y1_coor = $y_coor + $gc_figure_heigth - $hash_p->{$pos[$i]} * $gc_resolution;
		my $y2_coor = $y_coor + $gc_figure_heigth - $hash_p->{$pos[$i+1]} * $gc_resolution;

		$svg->line('x1',$x1_coor,'y1',$y1_coor,'x2',$x2_coor,'y2',$y2_coor,'stroke',$line_color,'stroke-width',1);
		
	}

}


#read fasta file
#usage: Read_fasta($file,\%hash);
#############################################
sub Read_fasta{
	my $file=shift;
	my $hash_p=shift;
	
	my $total_num;
	open(IN, $file) || die ("can not open $file\n");
	$/=">"; <IN>; $/="\n";
	while (<IN>) {
		chomp;
		my $head = $_;
		my $name = $1 if($head =~ /^(\S+)/);
		
		$/=">";
		my $seq = <IN>;
		chomp $seq;
		$seq=~s/\s//g;
		$/="\n";
		
		if (exists $hash_p->{$name}) {
			warn "name $name is not uniq";
		}

		$hash_p->{$name}{head} =  $head;
		$hash_p->{$name}{len} = length($seq);
		$hash_p->{$name}{seq} = $seq;

		$total_num++;
	}
	close(IN);
	
	return $total_num;
}


## draw exon-intron structure
sub draw_gene{
	my ($svg, $id, $pp,$y,$strand,$color) = @_ ;
	my $exon_size = 12;
	my $g = $svg->group('id'=>$id);
	
	## draw intron lines
	for (my $i=0; $i<@$pp-2; $i+=2) {
		my ($x_start,$x_end) = ($pp->[$i+1],$pp->[$i+2]);
		$g->line('x1',$x_start,'y1',$y+$exon_size/2,'x2',$x_end,'y2',$y+$exon_size/2,'stroke','black','stroke-width',1);		
	}
	## draw exon rectanges
	for (my $i=0; $i<@$pp; $i+=2) {
		my ($x_start,$x_end) = ($pp->[$i],$pp->[$i+1]);
		$g->rect('x',$x_start, 'y',$y,'width',($x_end-$x_start || 2),'height',$exon_size,'fill',$color,'onclick',"alert('$id')");
	}

	
	## draw exon polygons to present strand
	if ($strand eq "+") {
		my $x_end = $pp->[-1];
		my $x_start = ($pp->[-1] - $pp->[-2] > $exon_size) ? ($x_end - $exon_size) : ($pp->[-2]) ;
		$g->rect('x',$x_start, 'y',$y,'width',$x_end-$x_start,'height',$exon_size,'fill','white');
		$g->polygon('points',[$x_start, $y,$x_start, $y+$exon_size,$x_end,$y+$exon_size/2],'fill',$color,"stroke",$color);


	}elsif ($strand eq "-") {
		my $x_start = $pp->[0];

		my $x_end = ($pp->[1] - $pp->[0] > $exon_size ) ? ($x_start + $exon_size) :  ($pp->[1]);
		$g->rect('x',$x_start, 'y',$y,'width',$x_end-$x_start,'height',$exon_size,'fill','white'); 
		$g->polygon('points',[$x_start,$y+$exon_size/2,$x_end,$y+$exon_size,$x_end,$y],'fill',$color,"stroke",$color);
	}

}

## draw exon-intron structure
sub draw_feature{
	my ($svg, $id, $x_start,$x_end,$y,$strand,$color) = @_ ;
	my $block_size = 12;

	$svg->rect('x',$x_start, 'y',$y,'width',($x_end-$x_start || 2),'height',$block_size,'fill',$color,'onclick',"alert('$id')");
}

sub plot_ruler{
	my ($svg,$Y,$X_start,$X_end,$bp_start,$bp_end,$type,$pos) = @_ ;
	
	my $scale_size = 6;
	my $divid = 50;

	my $bp_len = $bp_end - $bp_start;
	my $X_len = $X_end - $X_start;

	my $g = $svg->group('id'=>'scale ruler');
	
	## draw the main axis
	$g->line('x1',$X_start,'y1',$Y,'x2',$X_end,'y2',$Y,'stroke','#000000','stroke-width',1);		
	return if($bp_end - $bp_start  == 0);
	
	##draw ruler mark text at the specified postion(left or right of the ruler)
	if ($pos eq '' || $pos eq "left") {
		$g->text('x',$X_start-textWidth($Font_family,$Font_size,$type)-6,'y',$Y,'-cdata',$type,"font-family",$Font_family,"font-size",$Font_size,"fill",'#000000');
	}
	if ($pos eq "right") {
		$g->text('x',$X_end + 6,'y',$Y,'-cdata',$type,"font-family",$Font_family,"font-size",$Font_size,"fill",'#000000');
	}


	## draw small scale lines
	$bp_start--;
	for (my $i=$bp_start; $i<=$bp_end; $i+=1000) {
		my $X = $X_start + ($i - $bp_start) / $bp_len * $X_len;
		$g->line('x1',$X,'y1',$Y - $scale_size/2,'x2',$X,'y2',$Y,'stroke','#000000','stroke-width',1);
	}
	
	for (my $i=$bp_start; $i<=$bp_end; $i+=10000) {
		my $X = $X_start + ($i - $bp_start) / $bp_len * $X_len;
		$g->line('x1',$X,'y1',$Y - $scale_size,'x2',$X,'y2',$Y,'stroke','#000000','stroke-width',1);
		$svg->text('x',$X - textWidth($Font_family,$Font_size,$i/1000000) / 2,'y',$Y+textHeight($Font_size)+4,'fill','#000000','-cdata',$i/1000000,'font-size',$Font_size, 'font-family',$Font_family);
	}
}





sub read_gene_psl{
	my $file=shift;
	my $pos_hp=shift;
	
	open (IN,$file) || die ("fail open $file\n");
	while (<IN>) {
		my @temp=split(/\s+/,$_);
		my $chr=$temp[13];
		next if($chr ne $Sequence);
		my $gene=$temp[9];
		my $strand=$temp[8];
		my $start=$temp[15]+1;
		my $end=$temp[16];
		next if($start < $Startpos || $end > $Endpos);
		my @starts=split(/,/,$temp[20]);
		my @sizes=split(/,/,$temp[18]);
		my @exon;
		for (my $i=0; $i<@starts; $i++) {
			push @exon, [$starts[$i]+1, $starts[$i]+$sizes[$i]];
		}

		$$pos_hp{$chr}{$gene}{strand}=$strand;
		$$pos_hp{$chr}{$gene}{start}=$start;
		$$pos_hp{$chr}{$gene}{end}=$end;
		$$pos_hp{$chr}{$gene}{exon}=\@exon;
	}
	close(IN);
}


sub read_gene_gff{
	my $file=shift;
	my $ref=shift;

	open (IN,$file) || die ("fail open $file\n");
	while (<IN>) {
		chomp;
		s/^\s+//;
		next if(/^\#/);
		my @t = split(/\t/);
		my $tname = $t[0];
		next if($tname ne $Sequence);
		my $qname;
		if ($t[2] eq 'mRNA' || $t[2] eq 'CDS') {
			$qname = $1 if($t[8] =~ /^GenePrediction\s+(\S+)/ || $t[8] =~ /^ID=(\S+?);/ || $t[8] =~ /^Parent=(\S+?);/);
		}
		if ($t[2] eq 'match' || $t[2] eq 'HSP') {
			$qname = $1 if($t[8] =~ /Target\s+\"(\S+)\"/);
		}
		

		if ($t[2] eq 'mRNA' || $t[2] eq 'match') {
			my $start = $t[3];
			my $end = $t[4];
			my $strand = $t[6];
			next if($start < $Startpos || $end > $Endpos);
			$ref->{$tname}{$qname}{strand} = $strand;
			$ref->{$tname}{$qname}{start} = $start;
			$ref->{$tname}{$qname}{end} = $end;
		}
		if ($t[2] eq 'CDS' || $t[2] eq 'HSP') {
			push @{$ref->{$tname}{$qname}{exon}}, [$t[3],$t[4]] if(exists $ref->{$tname}{$qname});
		}
	}
	close(IN);

	foreach my $scaff (sort keys %$ref) {
		my $scaff_p = $ref->{$scaff};
		foreach my $gene (sort keys %$scaff_p) {
			my $gene_p = $scaff_p->{$gene};
			my @exon = @{$gene_p->{exon}};
			@exon = reverse @exon if($exon[0][0] > $exon[-1][0]);
			$ref->{$scaff}{$gene}{exon} = \@exon;
		}
	}

}

sub read_feature_gff {
	my $file=shift;
	my $ref=shift;

	open (IN,$file) || die ("fail open $file\n");
	while (<IN>) {
		chomp;
		s/^\s+//;
		my @t = split(/\t/);
		my $tname = $t[0];
		next if($tname ne $Sequence);
		my $qname = $1 if($t[8] =~ /^ID=(\S+);/ || $t[8] =~ /^Parent=(\S+);/);
		my $start = $t[3];
		my $end = $t[4];
		next if($start < $Startpos || $end > $Endpos);
		my $strand = $t[6];
		$ref->{$tname}{$qname}{strand} = $strand;
		$ref->{$tname}{$qname}{start} = $start;
		$ref->{$tname}{$qname}{end} = $end;

	}
	close(IN);

}