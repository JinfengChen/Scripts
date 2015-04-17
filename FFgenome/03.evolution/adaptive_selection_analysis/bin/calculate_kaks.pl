#!/usr/bin/perl

=head1 Name

calculate_kaks.pl  --  calculate pairwise kaks for genes in mfa file

=head1 Description

The programs invoked: Kaks_Calculator.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2008-8-2

=head1 Usage
   
   perl calculate_kaks.pl [option] <cds.mfa>
   --kaks_caculator <str>  set the path of kaks_caculator software
   --model <str>     set the model for kaks_caculator, defaul YN
   --outdir          set the output directory
   --verbose         output verbose information to screen  
   --help            output help information to screen  

=head1 Exmple

  perl ../bin/calculate_kaks.pl ./leptin.cds.fa.muscle

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory


##global variables
my ($Outdir,$Verbose,$Help);
my ($Kaks_Calculator,$Model);

##get options from command line
GetOptions(
	"kaks_calculator:s"=>\$Kaks_Calculator,
	"model:s"=>\$Model,
	"outdir:s"=>\$Outdir,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
die `pod2text $0` if (@ARGV < 1 || $Help);

my $cds_mfa_file = shift;
my $cds_core_file = basename($cds_mfa_file);

$Kaks_Calculator ||= "/home/jfchen/FFproject/FFgenome/03.evolution/KaKs_Calculator1.2/bin/KaKs_Calculator";
$Model ||= "YN";

$Outdir ||= ".";
$Outdir =~ s/\/$//;
mkdir($Outdir) unless(-d $Outdir);

mfa2axt($cds_mfa_file,"$Outdir/$cds_core_file.axt");
`$Kaks_Calculator -m $Model -i $Outdir/$cds_core_file.axt -o $Outdir/$cds_core_file.axt.kaks`;


####################################################
################### Sub Routines ###################
####################################################


##get matrix of Ka an Ks form KaKs_caculator's result
sub mfa2axt{
	my $mfa_file = shift;
	my $axt_file = shift;
	my (%name_seq,%pair,$output);

	open IN, $mfa_file || die "fail open $mfa_file\n";
	$/=">"; <IN>; $/="\n";
	while (<IN>) {
		my $name = $1 if(/^(\S+)/);
		$/=">";
		my $seq = <IN>;
		chomp $seq;
		$seq =~ s/\s//g;
		$/="\n";
		$name_seq{$name} = $seq;
	}
	close IN;

	foreach my $first (sort keys %name_seq) {
		foreach my $second (sort keys %name_seq) {
			next if($first eq $second || exists $pair{"$second&$first"});
			$pair{"$first&$second"} = 1;
		}
	}

	foreach (sort keys %pair) {
		if (/([^&]+)&([^&]+)/) {
			$output .= $_."\n".$name_seq{$1}."\n".$name_seq{$2}."\n\n";
		}
	}

	open OUT, ">$axt_file" || die "fail $axt_file";
	print OUT $output;
	close OUT;
}

