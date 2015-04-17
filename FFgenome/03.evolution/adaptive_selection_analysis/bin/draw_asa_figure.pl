#!/usr/bin/perl

=head1 Name

draw_asa_figure.pl  --  draw figure to display the ASA result

=head1 Description

Calculate average Ka/Ks between two subtrees.
Show average Ka/Ks and P-value on the tree figure.

=head1 Version
  
  Author: Yujie Hu, huyj@genomics.org.cn
  Modify: Fan Wei, fanw@genomics.org.cn
  Version: 1.1,  Date: 2008-8-2

=head1 Usage
  
   perl draw_asa_figure.pl [option] <tree_nhx_file> <kaks_result_file> <tree_split_file> <R_result_file>
   --width <int>     set the width of figure, default = 1000
   --outdir          set the output directory, default "./"
   --verbose         output verbose information to screen  
   --help            output help information to screen  

=head1 Exmple
  
  perl ../bin/draw_asa_figure.pl ../input/leptin.cds.fa.muscle.nj.nhx ./leptin.cds.fa.muscle.axt.kaks ./leptin.cds.fa.muscle.nj.nhx.split ./leptin.cds.fa.muscle.nj.nhx.R.test
  
=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use lib "$Bin/../lib";
use SVG;
use Tree::nhx_svg;

my $Width;
my ($Outdir,$Verbose,$Help);
GetOptions(
	"width:i"=>\$Width,
	"outdir:s"=>\$Outdir,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Width ||= 1000;
die `pod2text $0` if ($Help);

my $tree_nhx_file = shift;
my $kaks_result_file = shift;
my $tree_split_file = shift;
my $R_result_file = shift;
my $tree_core_file = basename($tree_nhx_file);

my %Attribution;
my %kaks;
my (%root,%left,%right);

$Outdir ||= ".";
$Outdir =~ s/\/$//;
mkdir($Outdir) unless(-d $Outdir);

##extract P-value from the R result file
open (IN, "$R_result_file")||die "fail open $R_result_file";
while (<IN>) {
	chomp;
	my ($name,$pvalue);
	if (/data\:\s+(\S+)\_RootLeft_group1\s+and\s+(\S+)\_group2/) {
		$name = $1;
		my $next = <IN>;
		if ($next=~/W\s+\S+\s+\S+\,\s+p-value\s+\S+\s+(\S+)/) {
			$pvalue = $1;
			$Attribution{$name}{RootLeftP}=$pvalue;
		}
	}
	if (/data\:\s+(\S+)\_RootRight_group1\s+and\s+(\S+)\_group2/) {
		$name = $1;
		my $next = <IN>;
		if ($next=~/W\s+\S+\s+\S+\,\s+p-value\s+\S+\s+(\S+)/) {
			$pvalue = $1;
			$Attribution{$name}{RootRightP}=$pvalue;
		}
	}
	if (/data\:\s+(\S+)\_LeftRight_group1\s+and\s+(\S+)\_group2/) {
		$name = $1;
		my $next = <IN>;
		if ($next=~/W\s+\S+\s+\S+\,\s+p-value\s+\S+\s+(\S+)/) {
			$pvalue = $1;
			$Attribution{$name}{LeftRightP}=$pvalue;
		}
	}
}
close IN;


##extract kaks information for kaks result file
open (IN, $kaks_result_file)||die "fail open $kaks_result_file\n";
<IN>;
while (<IN>) {
	chomp;
	my @array = split;
	if ($array[0]=~/(\S+)&(\S+)/){ 
		my $key = ($1 lt $2)?"$1&$2":"$2&$1";
		if ($array[4]=~/[-\d\.eE]+/) {
			$kaks{$key} = $array[4];
		}
	}
}
close IN;

##extract node information from tree_split file
open (IN, $tree_split_file)||die "fail open $tree_split_file\n";
$/=">";<IN>;$/="\n";
while (<IN>) {
	my @temp = split (/:/,$_);
	my $name = $temp[1];
	$name =~s/\s+$//;$name =~s/^\s+//;;
	$/=">";
	my $seq = <IN>;
	chomp $seq;
	$/="\n";
	my @array = split (/\n/,$seq);
	my @a = split (/:/,$array[0]);my $root = $a[1];
	my @b = split (/:/,$array[1]);my $left = $b[1];
	my @c = split (/:/,$array[2]);my $right = $c[1];
	$root =~s/\s+$//;$root =~s/^\s+//; 
	$left =~s/\s+$//;$left =~s/^\s+//; 
	$right =~s/\s+$//;$right =~s/^\s+//; 
	
	$root{$name}=$root;
	$left{$name}=$left;
	$right{$name}=$right;
}
close IN;

##calculate average Ka/Ks between two subtrees
foreach my $n  (keys %left) {
	my @l = split (/\s+/,$left{$n});
	my @ri = split (/\s+/,$right{$n});
	my @ro = split (/\s+/,$root{$n});
	if (@ro>0 && @l>0) {
		my @rol_kaks;
		for (my $i=0;$i<@ro ;$i++) {
			for (my $j=0;$j<@l ;$j++) {
				my $title = ($ro[$i] lt $l[$j])?"$ro[$i]&$l[$j]":"$l[$j]&$ro[$i]";
				push @rol_kaks, $kaks{$title} if ($kaks{$title} ne '');
			}
		}
		my $sum;
		for (my $i=0;$i<@rol_kaks ;$i++) {
			 $sum += $rol_kaks[$i];
		}
		if (@rol_kaks > 0) {
			my $average = $sum/@rol_kaks;
			my $Wvalue = sprintf("%.4f",$average);
			$Attribution{$n}{RootLeftW}=$Wvalue;
		}
		else{
			$Attribution{$n}{RootLeftW}="nan";
		}	
	}
	if (@ro>0 && @ri>0) {
		my @rori_kaks;
		for (my $i=0;$i<@ro ;$i++) {
			for (my $j=0;$j<@ri ;$j++) {
				my $title = ($ro[$i] lt $ri[$j])?"$ro[$i]&$ri[$j]":"$ri[$j]&$ro[$i]";
				push @rori_kaks, $kaks{$title} if ($kaks{$title} ne '');
			}
		}	
		my $sum;
		for (my $i=0;$i<@rori_kaks ;$i++) {
			 $sum += $rori_kaks[$i];
		}
		if(@rori_kaks > 0){
			my $average = $sum/@rori_kaks;
			my $Wvalue = sprintf("%.4f",$average);
			$Attribution{$n}{RootRightW}=$Wvalue;
		}
		else{
			$Attribution{$n}{RootRightW}="nan";
		}
	}
	if (@l>0 && @ri>0) {
		my @lri_kaks;
		for (my $i=0;$i<@ri ;$i++) {
			for (my $j=0;$j<@l ;$j++) {
				my $title = ($ri[$i] lt $l[$j])?"$ri[$i]&$l[$j]":"$l[$j]&$ri[$i]";
				push @lri_kaks, $kaks{$title} if ($kaks{$title} ne '');
			}
		}	
		my $sum;
		for (my $i=0;$i<@lri_kaks ;$i++) {
			 $sum += $lri_kaks[$i];
		}
		if(@lri_kaks > 0){
			my $average = $sum/@lri_kaks;
			my $Wvalue = sprintf("%.4f",$average);
			$Attribution{$n}{LeftRightW}=$Wvalue;
		}
		else {
			$Attribution{$n}{LeftRightW}="nan";
		}
	}
}

##generate the final asa files: .asa.nhx and .asa.svg
my $tree = new Tree::nhx_svg("exter_include",1,"show_groupKaKs",1,"groupKaKs_significance",0.05,"width",$Width,"skip",40);
$tree->parse($tree_nhx_file,"file");
$tree->mark_tree();
$tree->add_attribution(\%Attribution);

open (OUT,">$Outdir/$tree_core_file.asa.nhx")||"fail open $Outdir/$tree_core_file.asa.nhx\n";
print OUT $tree->string_nhx_format();
close OUT;

open (OUT,">$Outdir/$tree_core_file.asa.nhx.svg")||die "fail open $Outdir/$tree_core_file.asa.nhx.svg\n";
print OUT $tree->plot();
close OUT;

