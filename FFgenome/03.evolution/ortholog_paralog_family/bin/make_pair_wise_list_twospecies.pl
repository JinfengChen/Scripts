#!/usr/bin/perl

=head1 Name

make_pair_wise_list_twospecies.pl  --  the the pair wise gene list for synteny analysis

=head1 Description



=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2008-9-5
  Note:

=head1 Usage
  
  perl make_pair_wise_list_twospecies.pl <hcluster_file> <structure_file> <key_species> <key_species2>
  --outdir    set the output directory
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

##get options from command line into variables and set default values
my ($Verbose,$Help,$Outdir);
GetOptions(
	"verbose"=>\$Verbose,
	"help"=>\$Help,
	"outdir:s"=>\$Outdir
);
die `pod2text $0` if (@ARGV == 0 || $Help);

my $solar_file = shift;
my $structure_file = shift;
my $key_species = shift;
my $key_species2 = shift;

my %Structure;

$Outdir ||= ".";
mkdir($Outdir);

##read the structure file
##At1g01010.1_ARATH       CHR1v01212004  -  6  3760,3913;3996,4276;4486,4605;4706,5095;5174,5326;5439,5630;  3760  1
open IN,$structure_file || die $!;
while (<IN>) {
	chomp;
	my ($gene,$chr,$strand,$position,$index) = (split /\s+/)[0,1,2,5,6];
	$Structure{$gene} = "$chr\t$strand\t$position\t$index"; ##chr1_id  start_position  strand  index_number
}
close IN;

##read hcluster file and generate result
open ALL, ">$Outdir/$key_species.$key_species2.pairwise.list" || die "fail";
open IN,$solar_file || die $!;
while (<IN>) {
	chomp;
	
	my @t = split /\t/;
	my $gene1 = $t[0];
	my $gene2 = $t[5];
	next if($gene1 eq $gene2);
	next if(!exists $Structure{$gene1} || !exists $Structure{$gene2});
	
	##solar结果已经为单份，此处可保证XY图为双份
	if ( $gene1 =~ /_$key_species$/ && $gene2 =~ /_$key_species2$/ ) {
		print ALL "$gene1\t$Structure{$gene1}\t$gene2\t$Structure{$gene2}\n";
	}
	if ( $gene1 =~ /_$key_species2$/ && $gene2 =~ /_$key_species$/ ) {
		print ALL "$gene2\t$Structure{$gene2}\t$gene1\t$Structure{$gene1}\n";
	}
}
close IN;
close ALL;





####################################################
################### Sub Routines ###################
####################################################
