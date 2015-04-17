#!/usr/bin/perl

=head1 Name

draw_4dtv_kaks.pl  -- draw distribution figure of 4dtv and Ks, Ka/Ks

=head1 Description

The input is a file with postfix .4dtv or .kaks;

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage

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
my ($Verbose,$Help);
GetOptions(
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
die `pod2text $0` if (@ARGV == 0 || $Help);

my $file = shift;

my $common_bin = "$Bin/../../../common_bin";

if ($file =~ /\.4dtv$/i) {
	##draw 4dtv_corrected distribution figure
	`awk '\$2!="NA" && \$1!="tag" {print \$2}' $file > $file.4dtv_corrected `;
	`perl $common_bin/distribute_pre.pl --frequence  --binsize 0.02 --minborder 0 --cut_right 3 --header rect --color blue --Note "Distribution of 4dtv values (corrected)" --Y "Number of orthologs" --X "4dtv (window=0.02)" $file.4dtv_corrected > $file.4dtv_corrected.lst`;
	`perl $common_bin/distribute_svg.pl $file.4dtv_corrected.lst $file.4dtv_corrected.svg`;

	##draw 4dtv_raw distribution figure
	`awk '\$2!="NA" && \$1!="tag" {print \$3}' $file > $file.4dtv_raw `;
	`perl $common_bin/distribute_pre.pl --frequence  --binsize 0.01 --minborder 0 --cut_right 3 --header rect --color blue --Note "Distribution of 4dtv values (raw)"  --Y "Number of orthologs"  --X "4dtv (window=0.01)" $file.4dtv_raw > $file.4dtv_raw.lst`;
	`perl $common_bin/distribute_svg.pl $file.4dtv_raw.lst $file.4dtv_raw.svg`;

}


if ($file =~ /\.kaks$/i) {
	##draw ks distribution figure
	`awk '\$1!="Sequence" && \$4<99 && \$5!="NA" && \$5!="nan" {print \$4}' $file > $file.ks`;
	`perl $common_bin/distribute_pre.pl --frequence  --binsize 0.02 --minborder 0 --cut_right 3 --header rect --color blue --Note "Distribution of Ks values"  --Y "Number of orthologs" --X "Ks (window=0.02)" $file.ks > $file.ks.lst`;
	`perl $common_bin/distribute_svg.pl $file.ks.lst $file.ks.svg`;

	##draw W(Ka/Ks) distribution figure
	`awk '\$1!="Sequence" && \$4<99 && \$5!="NA" && \$5!="nan" {print \$5}' $file > $file.W`;
	`perl $common_bin/distribute_pre.pl --frequence  --binsize 0.002 --minborder 0 --cut_right 3 --header rect --color blue --Note "Distribution of Ka/Ks values"  --Y "Number of orthologs" --X "Ka/Ks (window=0.002)" $file.W > $file.W.lst`;
	`perl $common_bin/distribute_svg.pl $file.W.lst $file.W.svg`;

}


####################################################
################### Sub Routines ###################
####################################################
