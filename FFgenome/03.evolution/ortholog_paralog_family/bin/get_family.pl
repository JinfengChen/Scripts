use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory


my $family_list = shift;
my $fam_dir = shift;
my $target_dir = shift;
my $fam_in_subdir_num = 100;
my @fam;

read_family($family_list,\@fam);

#print Dumper \@fam;
#exit;

foreach my $family_id (@fam) {
	my $sub_dir = cacul_sub_dir($family_id,$fam_in_subdir_num);
	`cp $fam_dir/$sub_dir/$family_id/$family_id.cds.muscle.best.nhx.Ks.Ka.KaKs.svg.info.svg $target_dir`;
	`cp $fam_dir/$sub_dir/$family_id/$family_id.annotation $target_dir`;
}



##caculate subdir id in a two ranks direcotry
##the format is: 1-100, 101-200, 201,300, .......
####################################################
sub cacul_sub_dir {
	my $family_id = shift;
	my $fam_in_subdir_num = shift;
	
	$family_id -= 1;
	my $sub_dir_num = ($family_id - $family_id % $fam_in_subdir_num) / $fam_in_subdir_num;
	my $start_fam_id = $sub_dir_num * $fam_in_subdir_num + 1;
	my $end_fam_id = ($sub_dir_num + 1) * $fam_in_subdir_num;

	my $sub_dir = "$start_fam_id-$end_fam_id";
	return $sub_dir;
}


sub read_family{
	my $file = shift;
	my $ary_p = shift;

	open IN,$file || die "fail open $file";
	while (<IN>) {
		chomp;
		next if(/^fam_id/);
		my @t = split /\t/;
		push @$ary_p,$t[0];
	}
	close IN;
	
}