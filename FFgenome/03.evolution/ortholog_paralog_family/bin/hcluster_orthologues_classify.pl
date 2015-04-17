#!/usr/bin/perl

=head1 Name

hcluster_orthologues_classify.pl  --  classify the genes on orthologues

=head1 Description

read in the *.hcluster.stat file, and generate a statistic table.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2008-6-19
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

my $hcluster_stat_file = shift;

my @Num;

my $class_result;


#classes
#0#1:1:1 orthologues
#1#X:X:X orthologues
#2#Patchy orthologues
#3#insect-specific
#4#vertebrate-specific
#5#species-specific mutlti-copy
#6#species-specific single-copy

#species
#family	BOMMO	DROME	ANOGA	APIME	HUMAN	GALGA	BRARE
#0		1		2		3		4		5		6		7


open IN,$hcluster_stat_file || die "fail $hcluster_stat_file";
while (<IN>) {
	next if(/^fam_id/);
	chomp;
	my ($fam_id,$total,$BOMMO,$DROME,$ANOGA,$APIME,$HUMAN,$GALGA,$FUGRU,$CAEEL) = split /\t/;
	if (all_species($BOMMO,$DROME,$ANOGA,$APIME,$HUMAN,$GALGA,$FUGRU)) {
		if (single_copy($BOMMO,$DROME,$ANOGA,$APIME,$HUMAN,$GALGA,$FUGRU)) {
			$Num[0][0]++; $Num[0][1]+=$BOMMO; $Num[0][2]+=$DROME; $Num[0][3]+=$ANOGA; $Num[0][4]+=$APIME; $Num[0][5]+=$HUMAN; $Num[0][6]+=$GALGA; $Num[0][7]+=$FUGRU; 
			$class_result .= "$fam_id\t1:1:1 orthologues\n";
		}else{
			$Num[1][0]++; $Num[1][1]+=$BOMMO; $Num[1][2]+=$DROME; $Num[1][3]+=$ANOGA; $Num[1][4]+=$APIME; $Num[1][5]+=$HUMAN; $Num[1][6]+=$GALGA; $Num[1][7]+=$FUGRU; 
			$class_result .= "$fam_id\tX:X:X orthologues\n";
		}
	}elsif(patchy($BOMMO,$DROME,$ANOGA,$APIME,$HUMAN,$GALGA,$FUGRU)){
		$Num[2][0]++; $Num[2][1]+=$BOMMO; $Num[2][2]+=$DROME; $Num[2][3]+=$ANOGA; $Num[2][4]+=$APIME; $Num[2][5]+=$HUMAN; $Num[2][6]+=$GALGA; $Num[2][7]+=$FUGRU; 
		$class_result .= "$fam_id\tPatchy orthologues\n";
	}elsif(insect_specific($BOMMO,$DROME,$ANOGA,$APIME,$HUMAN,$GALGA,$FUGRU)){
		$Num[3][0]++;$Num[3][1]+=$BOMMO; $Num[3][2]+=$DROME; $Num[3][3]+=$ANOGA; $Num[3][4]+=$APIME; $Num[3][5]+=$HUMAN; $Num[3][6]+=$GALGA; $Num[3][7]+=$FUGRU;
		$class_result .= "$fam_id\tinsect-specific\n";
	}elsif(vertebrate_specific($BOMMO,$DROME,$ANOGA,$APIME,$HUMAN,$GALGA,$FUGRU)){
		$Num[4][0]++; $Num[4][1]+=$BOMMO; $Num[4][2]+=$DROME; $Num[4][3]+=$ANOGA; $Num[4][4]+=$APIME; $Num[4][5]+=$HUMAN; $Num[4][6]+=$GALGA; $Num[4][7]+=$FUGRU; 
		$class_result .= "$fam_id\tvertebrate-specific\n";
	}else{
		if ($BOMMO >= 2 || $DROME >= 2 || $ANOGA >= 2 || $APIME >= 2 || $HUMAN >= 2 || $GALGA >= 2 || $FUGRU >= 2) {
			$Num[5][0]++; $Num[5][1]+=$BOMMO; $Num[5][2]+=$DROME; $Num[5][3]+=$ANOGA; $Num[5][4]+=$APIME; $Num[5][5]+=$HUMAN; $Num[5][6]+=$GALGA; $Num[5][7]+=$FUGRU; 
			$class_result .= "$fam_id\tspecies-specific mutlti-copy\n";
		}
	}
}
close IN;


open OUT,">$hcluster_stat_file.class" || die "fail";
print OUT $class_result;
close OUT;

open OUT,">$hcluster_stat_file.class.stat" || die "fail";
print OUT "class\tBOMMO\tDROME\tANOGA\tAPIME\tHUMAN\tGALGA\tBRARE\tfamily\n";
my @Class = ("1:1:1 orthologues","X:X:X orthologues","Patchy orthologues","insect-specific","vertebrate-specific","species-specific mutlti-copy");
for (my $i=0; $i<@Num; $i++) {
	my $output = $Class[$i];
	for (my $j=1; $j<@{$Num[$i]}; $j++) {
		$output .= "\t$Num[$i][$j]";
	}
	$output .= "\t$Num[$i][0]\n"; ##family
	print OUT $output;
}
close OUT;

####################################################
################### Sub Routines ###################
####################################################

##至少有两个昆虫物种
sub insect_specific{
	my ($BOMMO,$DROME,$ANOGA,$APIME,$HUMAN,$GALGA,$FUGRU) = @_;
	my $insect_num = 0;
	$insect_num++ if($BOMMO);
	$insect_num++ if($DROME);
	$insect_num++ if($ANOGA);
	$insect_num++ if($APIME);
	my $is_specific = ($insect_num >= 2 && !$HUMAN && !$GALGA && !$FUGRU) ? 1 : 0;
	return $is_specific;
}

##至少有两个脊椎动物物种
sub vertebrate_specific{
	my ($BOMMO,$DROME,$ANOGA,$APIME,$HUMAN,$GALGA,$FUGRU) = @_;
	my $vertebrate_num = 0;
	$vertebrate_num++ if($HUMAN);
	$vertebrate_num++ if($GALGA);
	$vertebrate_num++ if($FUGRU);
	my $is_specific = ($vertebrate_num >= 2 && !$BOMMO && !$DROME && !$ANOGA && !$APIME) ? 1 : 0;
	return $is_specific;
}


##至少有一个insect和一个vertebrate物种
sub patchy{
	my ($BOMMO,$DROME,$ANOGA,$APIME,$HUMAN,$GALGA,$FUGRU) = @_;
	my $with_insect = ($BOMMO || $DROME || $ANOGA || $APIME) ? 1 : 0;
	my $with_vertebrate = ($HUMAN || $GALGA || $FUGRU) ? 1 : 0;
	my $is_patchy = ($with_insect && $with_vertebrate) ? 1 : 0;
	return $is_patchy;
}

##允许一个物种没有
sub all_species{
	my ($BOMMO,$DROME,$ANOGA,$APIME,$HUMAN,$GALGA,$FUGRU) = @_;
	my $species_num = 0;
	$species_num++ if($BOMMO);
	$species_num++ if($DROME);
	$species_num++ if($ANOGA);
	$species_num++ if($APIME);
	$species_num++ if($HUMAN);
	$species_num++ if($GALGA);
	$species_num++ if($FUGRU);
	my $is_all = ($species_num >= 6) ? 1 : 0; 
	return $is_all;
}

##允许一个物种没有
sub single_copy{
	my ($BOMMO,$DROME,$ANOGA,$APIME,$HUMAN,$GALGA,$FUGRU) = @_;
	my $duplicate_num = 0; ##存储2的个数
	$duplicate_num++ if($BOMMO == 2);
	$duplicate_num++ if($DROME == 2);
	$duplicate_num++ if($ANOGA == 2);
	$duplicate_num++ if($APIME == 2);
	$duplicate_num++ if($HUMAN == 2);
	$duplicate_num++ if($GALGA == 2);
	$duplicate_num++ if($FUGRU == 2);

	my $is_single = ($duplicate_num <= 1 && $BOMMO <= 2 && $DROME <= 2 && $ANOGA <= 2 && $APIME <= 2 && $HUMAN <= 2 && $GALGA <= 2 && $FUGRU <= 2) ? 1 : 0; 
	return $is_single;
}

