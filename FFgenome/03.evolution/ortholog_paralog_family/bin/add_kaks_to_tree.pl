#!/usr/bin/perl

=head1 Name

add_kaks_to_tree.pl  -- caculate Ka/Ks and add to tree

=head1 Description

read in cds.mfa and tree.nhx file, use treebest and KaKs_Caculator software.
use the same method as the FGF website: http://fgf.genomics.org.cn

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage

  perl add_kaks_to_tree.pl [options] <cds_mfa_file> <tree_nhx_file>
  --treebest <str>  set the path of treebest software
  --kaks_caculator <str>  set the path of kaks_caculator software
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

  perl ~/fanw/GACP/GACP-5.0/06.evolution-analysis/ortholog_paralog_family/bin/add_kaks_to_tree.pl 1.cds.muscle 1.cds.muscle.best.nhx

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory
use lib $Bin;
use Tree::nhx_align_func_svg;

##get options from command line into variables and set default values
my ($treebest,$Kaks_Calculator);
my ($Verbose,$Help);
GetOptions(
	"treebest:s"=>\$treebest,
	"kaks_caculator:s"=>\$Kaks_Calculator,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$treebest ||= "/share/raid6/fanw/GACP/GACP-5.0/software/treebest-1.9.2/treebest";
$Kaks_Calculator ||= "/share/raid6/fanw/GACP/GACP-5.0/software/KaKs_Calculator1.2/src/KaKs_Calculator";
die `pod2text $0` if (@ARGV == 0 || $Help);

my $cds_mfa_file = shift;
my $tree_nhx_file = shift;

mfa2axt("$cds_mfa_file");
`$Kaks_Calculator -m YN -i $cds_mfa_file.axt -o $cds_mfa_file.axt.KaKs`;
get_KaKs_matrix("$cds_mfa_file.axt.KaKs");

`$treebest estlen $tree_nhx_file $cds_mfa_file.axt.KaKs.Ks.matrix Ds > $tree_nhx_file.Ks`;
`$treebest estlen $tree_nhx_file.Ks $cds_mfa_file.axt.KaKs.Ka.matrix Dn > $tree_nhx_file.Ks.Ka`;

my $nhx = Tree::nhx_align_func_svg->new();
$nhx->parse("$tree_nhx_file.Ks.Ka","file");
foreach my $p ($nhx->node) {
	if(defined $p->{Dn} && defined $p->{Ds}){
		$p->{W} = ($p->{Ds} > 0) ? sprintf("%.6f",$p->{Dn} / $p->{Ds}) :  "-";
	}
}

open OUT,">$tree_nhx_file.Ks.Ka.KaKs.nhx" || die "fail creat";
print OUT $nhx->string($nhx->root,"nhx");
close OUT;

`rm $tree_nhx_file.Ks $tree_nhx_file.Ks.Ka $cds_mfa_file.axt.KaKs.Ks.matrix $cds_mfa_file.axt.KaKs.Ka.matrix`;

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


##get matrix of Ka and Ks from KaKs_caculator's result
sub mfa2axt{
	my $mfa_file = shift;
	my (%name_seq,%pair,$output);

	open OUT, ">$mfa_file.axt" || die "fail $mfa_file.axt";

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
			print OUT $_."\n".$name_seq{$1}."\n".$name_seq{$2}."\n\n";
		}
	}

	close OUT;
}
