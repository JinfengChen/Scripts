use strict;

die "$0 <hscore_file> <hcluster_file>\n" if(@ARGV != 2);

my $hscore_file = shift;
my $hcluster_file = shift;

my %Fam;
my %Score;

open IN,$hscore_file || die "fail $hscore_file";
while (<IN>) {
	chomp;
	my @t = split /\s+/;
	my $score = pop @t;
	@t = sort @t;
	$Score{"$t[0]$t[1]"} = $score;
}
close IN;


open IN,$hcluster_file || die $!;
while (<IN>) {
	chomp;
	next if(/^fam_id/);
	my @t = split /\t/;
	$t[-1] =~ s/,$//;
	my @genes = split /,/, $t[-1];
	my $fam_id = $t[0];
	##my $gene_num = scalar(@genes);
	$Fam{$fam_id} = \@genes;
}
close IN;


##find links between clusters
print "fam_A_ID\tfam_B_ID\tA_gene_count\tB_gene_count\tlink_num\tlink_rate\tavg_score\n";
my @ID = sort {$a<=>$b} keys %Fam;
for (my $i=0; $i<@ID; $i++) {
	for (my $j=$i+1; $j<@ID; $j++) {
		my $fam_A = $ID[$i];
		my $fam_B = $ID[$j];
		my $fam_Ap = $Fam{$fam_A};
		my $fam_Bp = $Fam{$fam_B};
		my $fam_Anum = @$fam_Ap;
		my $fam_Bnum = @$fam_Bp;
		my ($link_num,$avg_score) = compare_family($fam_Ap,$fam_Bp);
		my $link_rate = $link_num / ($fam_Anum*$fam_Bnum);
		if ($link_num) {
			print "$fam_A\t$fam_B\t$fam_Anum\t$fam_Bnum\t$link_num\t$link_rate\t$avg_score\n";
		}
	}
}


##compare two families
sub compare_family{
	my ($ap,$bp) = @_;
	my $link_num = 0;
	my $avg_score = 0;
	for (my $i=0; $i<@$ap; $i++) {
		for (my $j=0; $j<@$bp; $j++) {
			my @t = ($ap->[$i],$bp->[$j]);
			@t = sort @t;
			if (exists $Score{"$t[0]$t[1]"}) {
				$link_num ++;
				$avg_score += $Score{"$t[0]$t[1]"};
			}
		}
	}
	$avg_score = ($link_num) ? $avg_score / $link_num : 0;
	return ($link_num,$avg_score);
}