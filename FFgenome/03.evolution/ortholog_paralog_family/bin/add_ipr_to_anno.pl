use strict;

my $ipr_file = shift;
my $anno_mark_file = shift;

my %IPR;

##ENSFCAP00000005027      IPR007087       Znf_C2H2        Zinc finger, C2H2-type
open IN,$ipr_file || die "fail $ipr_file";
while (<IN>) {
	next if(/^Ensembl Protein ID/);
	chomp;
	my @t = split /\t/;
	my $prot_id = $t[0];
	my $ipr_id = $t[1];
	my $ipr_desc = $t[3];
	next if($ipr_id !~ /IPR\d+/);
	$IPR{$prot_id}{$ipr_id} = $ipr_desc;
}
close IN;

open IN,$anno_mark_file || die "fail $anno_mark_file";
while (<IN>) {
	chomp;
	my $prot_id = $1 if(/^(\S+)/);
	$prot_id =~ s/_[^_]+$//;
	if (exists $IPR{$prot_id}) {
		my $ipr_str = "InterPro: ";
		my $p = $IPR{$prot_id};
		foreach my $ipr_id (sort keys %$p) {
			my $ipr_desc = $p->{$ipr_id};
			$ipr_str .= "$ipr_id, $ipr_desc; ";
		}
		print "$_  $ipr_str\n";
	}else{
		print "$_\n";
	}
}
close IN;
