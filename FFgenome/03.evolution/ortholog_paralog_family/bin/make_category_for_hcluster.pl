#!/usr/bin/perl

=head1 Name

make_category_for_hcluster.pl  --  make a category file for hcluster_sg

=head1 Description

this program generate a categroy file of genes, the format is like below:
At3g04520.1_ARATH       1
At3g04530.1_ARATH       1
Os03g0698800_ORYSA      1
Os03g0698900_ORYSA      1
W05B2.1_CAEEL   2
W05B2.2_CAEEL   2

Each gene has a category numeric id, where 1 stands for the main group,
and 2 stands for the outgroup.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage

  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

  perl make_category_for_hcluster.pl category.species.txt all.cds.mark > category.genes.txt 

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory

##get options from command line into variables and set default values
my ($Verbose,$Help);
GetOptions(
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
die `pod2text $0` if (@ARGV == 0 || $Help);

my $category_species_file = shift;
my $all_cds_file = shift;

my %Category;

open IN,$category_species_file || die $!;
while (<IN>) {
	chomp;
	my ($spec,$cate) = (split /\s+/)[0,1];
	$Category{$spec} = $cate;
}
close IN;

open IN,$all_cds_file || die $!;
while (<IN>) {
	chomp;
	if (/>(\S+)/) {
		my $gene = $1;
		my $spec = $1 if($gene =~ /_([^_]+)$/);
		print "$gene\t$Category{$spec}\n";
	}
}
close IN;


####################################################
################### Sub Routines ###################
####################################################
