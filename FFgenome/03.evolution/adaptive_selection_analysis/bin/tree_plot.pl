#!/usr/bin/perl

=head1 Name

tree_plot.pl  --  a generic program for plotting tree figures

=head1 Description

invoke the Tree::nhx_svg package, and make most of the Tree::nhx_svg's method accessible from commond line.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2008-8-6
  Note:

=head1 Usage

  perl tree_plot.pl [option] <tree.nhx>
  --figure_width <num>   set the width of tree figure, default 640
  --vertical_skip <num>  set the vertical distance between neighbor nodes, default 20
  --line_color <str>     set the color for branch lines, default "black"
  --line_width <num>     set the width for branch lines, default 1
  --show_box             show the node box
  --box_color <str>      set the color for node boxes, default "blue"
  --box_width <num>      set the width for node boxes, default 3
  --show_inner           show the inner names
  --inner_color <str>    set the color for inner names, default "red"
  --font_family <str>    set the font family, default "Times-Bold"
  --font_size <num>      set the font size, default 12
  --equal_branch         use equal branch length
  --show_kaks       show ka/ks on the branch, default no
  --show_B          show Bootstrap values on the branch, default no
  --show_groupKaKs  show the group Ka/Ks rank sum test result
  --groupKaKs_significance <num>  set the significance level of groupKaKs
  --show_ruler      show the ruler, default no
  --sort_tree       sort the leaf nodes according to child number
  --png             convert svg to png figure, default no
  --outdir          set the output directory, default "./"
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

  perl tree_plot.pl  tree.nhx

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory
use lib $Bin;
use lib "$Bin/../lib";
use Tree::nhx_svg;

##get options from command line into variables and set default values
my ($figure_width,$vertical_skip,$line_color,$line_width,$show_box,$box_color,$box_width);
my ($show_inner,$inner_color,$font_family,$font_size,$is_real,$equal_branch,$show_kaks,$show_B,$show_ruler,$png);
my ($show_groupKaKs, $groupKaKs_significance,$Sort_tree);
my ($Outdir,$Verbose,$Help);
GetOptions(
	"figure_width:i"=>\$figure_width,
	"vertical_skip:i"=>\$vertical_skip,,
	"line_color:s"=>\$line_color,
	"line_width:i"=>\$line_width,
	"show_box"=>\$show_box,
	"box_color:s"=>\$box_color,
	"box_width:i"=>\$box_width,
	"show_inner"=>\$show_inner,
	"inner_color:s"=>\$inner_color,
	"font_size:i"=>\$font_size,
	"font_family:s"=>\$font_family,
	"equal_branch"=>\$equal_branch,
	"show_kaks"=>\$show_kaks,
	"show_B"=>\$show_B,
	"show_ruler"=>\$show_ruler,
	"show_groupKaKs"=>\$show_groupKaKs,
	"groupKaKs_significance:f"=>\$groupKaKs_significance,
	"sort_tree"=>\$Sort_tree,
	"png"=>\$png,
	"outdir:s"=>\$Outdir,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
die `pod2text $0` if (@ARGV == 0 || $Help);

my $nhx_file = shift;
my $nhx_core = basename($nhx_file);

my $svg2xxx = "/share/raid12/chenjinfeng/tools/draw/svg2xxx_release/svg2xxx";

$Outdir ||= ".";
$Outdir =~ s/\/$//;
mkdir($Outdir) unless(-d $Outdir);

$figure_width ||= 640;
$vertical_skip ||= 20;
$line_color ||= "black";
$line_width ||= 1;
$show_box  = (defined $show_box) ? 1 : 0;
$box_color ||= "blue";
$box_width ||= 3;
$show_inner  = (defined $show_inner) ? 1 : 0;
$inner_color ||= "red";
$font_family ||= 'Times-Bold';
$font_size ||= 12;
$is_real = (defined $equal_branch) ? 0 : 1;
$show_kaks  = (defined $show_kaks) ? 1 : 0;
$show_B  = (defined $show_B) ? 1 : 0;
$show_ruler  = (defined $show_ruler) ? 1 : 0;
$show_groupKaKs = (defined $show_groupKaKs)  ? 1 : 0;
$groupKaKs_significance ||= 0.05;

my $nhx = Tree::nhx_svg->new("width",$figure_width,"skip",$vertical_skip,"c_line",$line_color,"line_width",$line_width,"show_box",$show_box,"c_node",$box_color,"half_box",$box_width,"show_inter",$show_inner,"c_inter",$inner_color, "font",$font_family,"fsize",$font_size,"is_real",$is_real,"show_W",$show_kaks,"show_B",$show_B,"show_ruler",$show_ruler,"show_groupKaKs",$show_groupKaKs,"groupKaKs_significance",$groupKaKs_significance);
$nhx->parse($nhx_file,"file");
$nhx->sort_tree() if(defined $Sort_tree);

open OUT,">$Outdir/$nhx_core.svg" || die "fail svg";
print OUT $nhx->plot();
close OUT;

`$svg2xxx -t png $Outdir/$nhx_core.svg` if(defined $png);

####################################################
################### Sub Routines ###################
####################################################
