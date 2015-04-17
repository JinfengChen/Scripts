#!/usr/bin/perl

=head1 Name

hcluster_stat.pl  --  stat the number of genes for each species in the familes

=head1 Description

this program read the hcluster_sg result, and do some statistics.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2008-6-3
  Note:

=head1 Usage

  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

  perl hcluster_stat.pl <category_species_file> <hcluster_result_file>

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
my $hcluster_result_file = shift;

my @Species;


open IN,$category_species_file || die $!;
while (<IN>) {
	chomp;
	my ($spec,$cate) = (split /\s+/)[0,1];
	push @Species,$spec if($spec);
}
close IN;

#print Dumper \@Species;
#exit;

open STAT, ">$hcluster_result_file.stat" || die "fail";
open SINGLE, ">$hcluster_result_file.stat.single-copy" || die "fail";

print STAT "fam_id\ttotal\t".join("\t",@Species)."\tspec_num\n";
print SINGLE "fam_id\ttotal\t".join("\t",@Species)."\tspec_num\n";

open IN,$hcluster_result_file || die $!;
while (<IN>) {
	chomp;
	my @t = split /\t/;
	$t[-1] =~ s/,$//;
	my @genes = split /,/, $t[-1];
	my $fam_id = $t[0];
	my $gene_num = scalar(@genes);
	my %sum;
	foreach  (@genes) {
		if (/_([^_]+)$/) {
			$sum{$1}++;
		}
	}
	my $output = "$fam_id\t$gene_num";
	my $is_single_copy = 1;
	foreach my $spec (@Species) {
		my $num = (exists $sum{$spec}) ? $sum{$spec} : 0;
		$output .= "\t$num";
		$is_single_copy = 0 if($num != 1);
	}
	my $spec_num = keys %sum;
	$output .= "\t$spec_num\n";
	
	print STAT $output;
	print SINGLE $output if($is_single_copy);
}
close IN;

close STAT;
close SINGLE;


####################################################
################### Sub Routines ###################
####################################################
