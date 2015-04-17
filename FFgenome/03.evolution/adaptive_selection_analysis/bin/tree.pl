#!/usr/bin/perl

=head1 Name

tree.pl  --  a generic program for tree applications

=head1 Description

invoke the Tree::nhx package, and make most of the Tree::nhx's method accessible from commond line.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2008-8-6
  Note:

=head1 Usage
  
  perl tree.pl [option] <tree.nhx>
  --leaf      output leaf names
  --leafnum   output leaf number
  --info      output details information for each nodes
  --mark      set mark for inter node which does not have name
  --sort      sort the tree according to children number of each node
  --reformat <type>  reformat the tree text: nh_low, nhx_low, nh_high, nhx_high, default nhx_high
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

  perl tree.pl --leaf tree.nhx

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory
use lib $Bin;
use lib "$Bin/../lib";
use Tree::nhx;

##get options from command line into variables and set default values
my ($Verbose,$Help);
my ($Leaf,$Leafnum,$Info,$Mark,$Sort,$Reformat);
GetOptions(
	"leaf"=>\$Leaf,
	"leafnum"=>\$Leafnum,
	"info"=>\$Info,
	"sort"=>\$Sort,
	"mark"=>\$Mark,
	"reformat:s"=>\$Reformat,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
die `pod2text $0` if (@ARGV == 0 || $Help);

my $nhx_file = shift;

my $nhx = Tree::nhx->new();
$nhx->parse($nhx_file,"file");

if (defined $Leaf || defined $Leafnum) {
	my @leaf = $nhx->leaf();
	my $leafnum = @leaf;
	print $leafnum."\n" if(defined $Leafnum);
	if (defined $Leaf) {
		foreach  (@leaf) {
			print "$_\n";
		}
	}
}

if (defined $Info) {
	print $nhx->info();
}

if (defined $Sort) {
	$nhx->sort_tree();
	print $nhx->string_nhx_format();
}

if (defined $Mark) {
	$nhx->mark_tree();
	print $nhx->string_nhx_format();
}

if (defined $Reformat) {
	print $nhx->string($nhx->root(),"nh") if($Reformat eq "nh_low");
	print $nhx->string($nhx->root(),"nhx") if($Reformat eq "nhx_low");
	print $nhx->string_nhx_format($nhx->root(),"nh") if($Reformat eq "nh_high");
	print $nhx->string_nhx_format($nhx->root(),"nhx") if($Reformat eq "nhx_high");
}


####################################################
################### Sub Routines ###################
####################################################
