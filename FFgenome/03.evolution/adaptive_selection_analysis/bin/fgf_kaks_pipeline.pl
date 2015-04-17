#!/usr/bin/perl

=head1 Name

fgf_kaks_pipeline.pl  -- put pairwise kaks onto tree branches by FGF method

=head1 Description

Inputs: cds.mfa file and tree.nhx file.
Method: the same method as the FGF website: http://fgf.genomics.org.cn

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2008-8-2

=head1 Usage

  perl add_kaks_to_tree.pl [options] <cds.mfa> <tree.nhx>
  --kaks_caculator <str>  set the path of kaks_caculator software
  --treebest <str>  set the path of treebest software
  --outdir          set the output directory, default "./"
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

  perl ../bin/add_kaks_to_tree.pl ./leptin.cds.fa.muscle ./leptin.cds.fa.muscle.nj.nhx

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

##get options from command line into variables and set default values
my ($treebest,$kaks_calculator);
my ($Outdir,$Verbose,$Help);
GetOptions(
	"kaks_calculator:s"=>\$kaks_calculator,
	"treebest:s"=>\$treebest,
	"outdir:s"=>\$Outdir,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$treebest ||= "/share/raid6/fanw/GACP/GACP-5.0/software/treebest-1.9.2/treebest";
$kaks_calculator ||= "/share/raid6/fanw/GACP/GACP-5.0/software/KaKs_Calculator1.2/src/KaKs_Calculator";
die `pod2text $0` if (@ARGV == 0 || $Help);

my $cds_mfa_file = shift;
my $tree_nhx_file = shift;
my $cds_mfa_core = basename($cds_mfa_file);
my $tree_nhx_core = basename($tree_nhx_file);
my $axt_kaks_file;

$Outdir ||= ".";
$Outdir =~ s/\/$//;
mkdir($Outdir) unless(-d $Outdir);

`perl $Bin/calculate_kaks.pl $cds_mfa_file --kaks_calculator $kaks_calculator -outdir $Outdir`;
$axt_kaks_file = "$Outdir/$cds_mfa_core.axt.kaks";

get_KaKs_matrix($axt_kaks_file);

`$treebest estlen $tree_nhx_file $axt_kaks_file.Ks.matrix Ds > $Outdir/$tree_nhx_core.Ks`;
`$treebest estlen $Outdir/$tree_nhx_core.Ks $axt_kaks_file.Ka.matrix Dn > $Outdir/$tree_nhx_core.Ks.Ka`;

my $nhx = Tree::nhx_svg->new();
$nhx->parse("$Outdir/$tree_nhx_core.Ks.Ka","file");
foreach my $p ($nhx->node) {
	if(defined $p->{Dn} && defined $p->{Ds}){
		$p->{W} = sprintf("%.6f",$p->{Dn} / $p->{Ds}) if ($p->{Ds} > 0);
		if ($p->{Ds} == 0){
			$p->{W} = ($p->{Dn} > 0) ? "infinite" : 0;
		}
	}
}

##creat the final nhx file
open OUT,">$Outdir/$tree_nhx_core.Ks.Ka.fgf.nhx" || die "fail creat";
print OUT $nhx->string($nhx->root,"nhx");
close OUT;

##draw svg figure
my $nhx_svg = Tree::nhx_svg->new('show_W',2,'show_B',0,"show_ruler",1,"dist_type","dm", "width",640,"skip",20,"is_real",1);
$nhx_svg->parse("$Outdir/$tree_nhx_core.Ks.Ka.fgf.nhx","file");
open OUT,">$Outdir/$tree_nhx_core.Ks.Ka.fgf.nhx.svg" || die "fail creat svg";
print OUT $nhx_svg->plot;
close OUT;

####################################################
################### Sub Routines ###################
####################################################

##get matrix of Ka an Ks form KaKs_caculator's result
sub get_KaKs_matrix{
	my $KaKs_file = shift;
	my (%Ka,%Ks,$output_Ka,$output_Ks);

	open IN, $KaKs_file || die "fail open $KaKs_file\n";
	while (<IN>) {
		my @t = split /\s+/;
		if($t[0] =~ /([^&]+)&([^&]+)/){
			$Ka{$1}{$2} = $t[2];
			$Ka{$2}{$1} = $t[2];
			$Ks{$1}{$2} = $t[3];
			$Ks{$2}{$1} = $t[3];
			
			$Ka{$1}{$1} = 0.000000;
			$Ka{$2}{$2} = 0.000000;
			$Ks{$1}{$1} = 0.000000;
			$Ks{$2}{$2} = 0.000000;

		}
	}
	close IN;

	my $num = keys %Ka;
	$output_Ka = "   $num\n";
	foreach my $first (sort keys %Ka) {
		$output_Ka .= sprintf("%-10s  ",$first);
		my $pp = $Ka{$first};
		foreach my $second (sort keys %$pp) {
			$output_Ka .= sprintf("%.6f\t",$pp->{$second});
		}
		$output_Ka =~ s/\t$//;
		$output_Ka .= "\n";
	}

	my $num = keys %Ks;
	$output_Ks = "   $num\n";
	foreach my $first (sort keys %Ks) {
		$output_Ks .= sprintf("%-10s  ",$first);
		my $pp = $Ks{$first};
		foreach my $second (sort keys %$pp) {
			$output_Ks .= sprintf("%.6f\t",$pp->{$second});
		}
		$output_Ks =~ s/\t$//;
		$output_Ks .= "\n";
	}
	
	open OUT,">$KaKs_file.Ka.matrix" || die "fail open $KaKs_file.Ka.matrix";
	print OUT $output_Ka;
	close OUT;
	
	open OUT,">$KaKs_file.Ks.matrix" || die "fail open $KaKs_file.Ks.matrix";
	print OUT $output_Ks;
	close OUT;

}


##get matrix of Ka an Ks form KaKs_caculator's result
sub mfa2axt{
	my $mfa_file = shift;
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

	open OUT, ">$mfa_file.axt" || die "fail $mfa_file.axt";
	print OUT $output;
	close OUT;
}
