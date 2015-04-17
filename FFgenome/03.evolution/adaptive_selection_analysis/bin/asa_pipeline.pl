#!/usr/bin/perl

=head1 Name

asa_pipeline.pl  --  adaptive selection analysis pipeline

=head1 Description

The programs invoked: Kaks_Calculator , R, and SVG & Tree package.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2008-8-2

=head1 Usage
  
   perl asa_pipeline.pl [option] <cds.mfa> <tree.nhx>
   --outdir          set the output directory
   --verbose         output verbose information to screen  
   --help            output help information to screen  

=head1 Exmple

  perl ../bin/asa_pipeline.pl ../input/leptin.cds.fa.muscle ../input/leptin.cds.fa.muscle.nj.nhx -verbose &
  
=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory
use lib "$Bin/../lib";
use SVG;
use Tree::nhx_svg;

##global variables
my ($cds_mfa_file,$cds_core_file,$tree_nhx_file,$tree_core_file);
my ($Outdir,$Verbose,$Help);
my (%config,$kaks_calculator,$R);

##get options from command line
GetOptions(
	"outdir:s"=>\$Outdir,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
die `pod2text $0` if (@ARGV < 1 || $Help);

$cds_mfa_file = shift;
$tree_nhx_file = shift;

$cds_core_file = basename($cds_mfa_file);
$tree_core_file = basename($tree_nhx_file);

$Outdir ||= ".";
$Outdir =~ s/\/$//;
mkdir($Outdir) unless(-d $Outdir);

parse_config("$Bin/config.txt",\%config);
$kaks_calculator = $config{kaks_caculator};
$R = $config{R};

##calculate Kaks
`perl $Bin/calculate_kaks.pl --kaks_calculator $kaks_calculator -outdir $Outdir $cds_mfa_file`;
print STDERR "kaks calculation complete\n" if(defined $Verbose);

##wilcox rank-sum test
`perl $Bin/rank_sum_test.pl --R $R -outdir $Outdir $tree_nhx_file $Outdir/$cds_core_file.axt.kaks`;
print STDERR "R rank sum test complete\n" if(defined $Verbose);

##draw svg figure to display asa result
`perl $Bin/draw_asa_figure.pl -outdir $Outdir $tree_nhx_file  $Outdir/$cds_core_file.axt.kaks  $Outdir/$tree_core_file.split  $Outdir/$tree_core_file.R.test`;
print STDERR "draw asa figure complete\n" if(defined $Verbose);


####################################################
################### Sub Routines ###################
####################################################


##parse the config.txt file, and check the path of each software
####################################################
sub parse_config{
	my $conifg_file = shift;
	my $config_p = shift;
	
	die "$conifg_file not exist" unless(-f $conifg_file);

	my $error_status = 0;
	open IN,$conifg_file || die "fail open: $conifg_file";
	while (<IN>) {
		if (/(\S+)\s*=\s*(\S+)/) {
			my ($software_name,$software_address) = ($1,$2);
			$config_p->{$software_name} = $software_address;
			if (! -e $software_address){
				warn "Non-exist:  $software_name  $software_address\n"; 
				$error_status = 1;
			}
		}
	}
	close IN;
	die "\nExit due to error of software configuration\n" if($error_status);
}


