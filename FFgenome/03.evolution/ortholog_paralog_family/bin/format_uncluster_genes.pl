use strict;

##24122   29378   91      1.000   1       3       ENSRNOP00000044356_RATNO,ENSRNOP00000051547_RATNO,ENSRNOP00000041004_RATNO,
die "$0 <gene_list_file> <fam_id_start>\n" if(@ARGV < 2);

my $gene_list_file = shift;
my $fam_id_start = shift;

open IN,$gene_list_file || die "fail $gene_list_file";
while (<IN>) {
	my $gene = $1 if(/^(\S+)/);
	print "$fam_id_start\tNo\tNo\tNo\tNo\t1\t$gene\n";
	$fam_id_start++;
}
close IN;