#!/usr/bin/perl

=head1 Name

draw the svg figure for ortholog_paralog_family pipeline

=head1 Description

Make this part to be a single program, in order to save memory.
在一个程序内建很多对象时，内存累积不释放.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage
  
  perl draw_tree_svg.pl <tree.nhx> [key_tag]
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple



=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory
use lib "$Bin";
use Tree::nhx_align_func_svg;


##get options from command line into variables and set default values
my ($Verbose,$Help);
GetOptions(
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
die `pod2text $0` if (@ARGV == 0 || $Help);

my $tree_nhx_file = shift;
my $Key_tag = shift;

my (%status);

my $nhx_anno_svg = Tree::nhx_align_func_svg->new("show_B","1");
$nhx_anno_svg->parse($tree_nhx_file,"file");
foreach my $p ($nhx_anno_svg->node()) {
	if (!$p->{C} && $p->{N}) {
		my $gene = $p->{N};
		$status{$gene} = 1 if($Key_tag && $gene =~ /_$Key_tag$/);
	}
}
$nhx_anno_svg->add_character("status",\%status);

my $tree_svg_file = $tree_nhx_file;
$tree_svg_file =~ s/\.nhx$/\.svg/;
open OUT,">$tree_svg_file" || die "fail";
print OUT $nhx_anno_svg->plot();
close OUT;

my $tree_list_file = $tree_nhx_file;
$tree_list_file =~ s/\.nhx$/\.list/;
open OUT,">$tree_list_file" || die "fail";
my %Leaf;
foreach my $p ($nhx_anno_svg->node()) {
	if (!$p->{C} && $p->{N}) {
		my $gene = $p->{N};
		my $y_pos = $p->{y};
		$Leaf{$gene} = $y_pos;
	}
}
foreach my $gene (sort {$Leaf{$a}<=>$Leaf{$b}} keys %Leaf) {
	my $y_pos = $Leaf{$gene};
	print OUT "$gene\t$y_pos\n";
}
close OUT;

####################################################
################### Sub Routines ###################
####################################################



