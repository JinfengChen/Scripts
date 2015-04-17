#!/usr/bin/perl

=head1 Name

bitScore_to_hclusterScore.pl  -- convert the blast bit score to hcluster Score

=head1 Description

The relation between blast bit score (S') and blast raw score (S) is that: S' = 0.385 * S + 4.6 ;
As blast raw score (S) is linear, so S' is nealy linear, because 4.6 is small enough that can be 
ignored. The blast bit score is affected mainly by two factors: query length and aa constitution.
Here we convert bit score to raw score, and convert raw score to hcluster score(HCS).
When two proteins are aligned, we assummed their potential best bit score(S_Best) should be the larger
value of the two self-aligning bit-scores. And this potential best bit was used as the denominator
to caculate hcluster score.

The two formulas are:
S = (S' - 4.6) /  0.385;
HCS = S / S_Best * 100;

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 2.0,  Date: 2006-12-6
  Note:

=head1 Usage

  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

  perl bitScore_to_hclusterScore.pl all_vs_all.blast.m8.solar > all_vs_all.blast.m8.solar.forHC 2> all_vs_all.blast.m8.solar.forHC.details

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

my $solar_file = shift;

my %Self; ##score the potential best bit score by self alignment independantly

##read self alignment bit score
open IN,$solar_file || die "fail open $solar_file";
while (<IN>) {
	my ($query,$target,$bitscore) = (split /\t/)[0,5,10];
	$Self{$query} = $bitscore if($query eq $target);
}
close IN;

##generate new hcluster score 
open IN,$solar_file || die "fail open $solar_file";
while (<IN>) {
	my ($query,$query_len,$target,$target_len,$bitscore,$query_pos) = (split /\t/)[0,1,5,6,10,11];
	print STDERR  "$query self alignment not exists!" if(!exists $Self{$query});
	print STDERR  "$target self alignment not exists!" if(!exists $Self{$target});
	
	##两条序列比对时，潜在的最佳分值应该是两个自身最佳分值当中较小的那一个
	my $best_bitscore = ($Self{$query} > $Self{$target}) ? $Self{$query} : $Self{$target};
	
	my $hcluster_score = bit2raw($bitscore) / bit2raw($best_bitscore) * 100;
	print STDERR  "Error:$query,$target,$hcluster_score\n" if($hcluster_score > 100);
	$hcluster_score = 100 if($hcluster_score > 100);
	printf("%s\t%s\t%d\n",$query,$target,$hcluster_score);
	
	##output the detailed information for checking bugs
	my $align_len;
	while ($query_pos =~ /(\d+),(\d+);/g) {
		$align_len = $2 - $1 + 1;
	}
	$hcluster_score = sprintf("%.2f",$hcluster_score);

}
close IN;



####################################################
################### Sub Routines ###################
####################################################


sub bit2raw{
	my $bit_score = shift;
	my $raw_score = ($bit_score - 4.6) /  0.385;
	return $raw_score;
}