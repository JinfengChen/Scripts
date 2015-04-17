#!/usr/bin/perl

=head1 Name



=head1 Description



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

my $stat_file = shift;

##1       2       3       4       5       6       7       8       9       10      11      12      13
##fam_id  total   PANDA   CANFA   FELCA   HUMAN   PANTR   MOUSE   RATNO   BOVIN   HORSE   MONDO   ORNAN
open IN, $stat_file || die "fail $stat_file";
while (<IN>) {
	next if(/^fam_id/);
	chomp;
	my @t = split /\t/;
	my $ingroup = 0;
	my $outgroup = 0;
	my $species = 0;
	my $panda = 0;
	
	pop @t;
	my ($num_single,$num_lost,$num_dupl) = (0,0,0);
	for (my $i=2; $i<@t-2; $i++) {
		$num_single ++ if($t[$i] == 1);
		$num_lost ++ if($t[$i] == 0);
		$num_dupl ++ if($t[$i] >= 2);
		$ingroup++ if($t[$i] > 0);
		$species++ if($t[$i] > 0);
	}
	$outgroup++ if($t[-1] >0);
	$outgroup++ if($t[-2] >0);
	$species++ if($t[-1] > 0);
	$species++ if($t[-2] > 0);
	
	##具有内群和外群基因的家族
	#if ($ingroup && $outgroup ) {
	#	print $_."\n";
	#}

#	##最少有两个物种的家族
#	if ($species >= 2 ) {
#		print $_."\n";
#	}

#	##允许1个lost和1个duplicate
#	if ($num_lost <= 1 && $num_dupl <= 1) {
#		print $_."\n";
#	}
	
#	##允许1个lost和1个duplicate
	if ($num_lost <= 2 && $num_dupl <= 2) {
		print $_."\n";
	}

#	##允许0个lost和0个duplicate
#	if ($num_lost <= 0 && $num_dupl <= 0) {
#		print $_."\n";
#	}

}
