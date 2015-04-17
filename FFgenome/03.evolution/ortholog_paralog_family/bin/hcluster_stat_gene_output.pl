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

  perl hcluster_stat.pl  <hcluster_result_file>

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

my $hcluster_result_file = shift;


open IN,$hcluster_result_file || die $!;
while (<IN>) {
	next if(/^fam_id/);
	chomp;
	my @t = split /\t/;
	$t[-1] =~ s/,$//;
	my @genes = split /,/, $t[-1];
	my $fam_id = $t[0];
	my $gene_num = scalar(@genes);

	foreach my $this (@genes) {
		print $this."\n";
	}

}
close IN;




####################################################
################### Sub Routines ###################
####################################################
