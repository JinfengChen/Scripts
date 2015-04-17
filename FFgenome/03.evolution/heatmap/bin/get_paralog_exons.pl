#!/usr/bin/perl
# usage: perl $0 para.cabbage.gff&

my $file= shift;
my %chr;

my ( $mRNA_s, $mRNA_e, $chr, $strand, $cds_num );
my @arr;

open IN, $file or die "Can't open $file:$!\n";
#open OUT_I, ">./gene_introns.gff" or die "Can't open ./gene_introns.gff: $!\n";
open OUT_E, ">./paralog_exons.gff" or die "Can't open ./paralog_exons.gff: $!\n";
#open OUT_E, ">./gene_exons.gff" or die "Can't open ./gene_exons.gff: $!\n";

while(<IN>){
	chomp;
	next unless ( /^A\d+/ );
	my @t= split;
	if ( $t[2] eq "mRNA" ){
		( $chr, $mRNA_s, $mRNA_e, $strand ) = ( $t[0], $t[3], $t[4], $t[6] );
		$cds_num= 0;
		@arr= ();
	}
	if ( $t[2] eq "CDS" ){
		$cds_num++;
		print OUT_E "$chr\tGLEAN\tparalog\t$t[3]\t$t[4]\n";
		#print OUT_E "$chr\tGLEAN\texons\t$t[3]\t$t[4]\n";
		#if ( ($strand eq "+") && ($cds_num > 1) ){
		#	print OUT_I "$chr\tGLEAN\tintrons\t".($arr[-1]+1)."\t".($t[3]-1)."\n";
		#}elsif( ($strand eq "-") && ($cds_num > 1) ){
		#	print OUT_I "$chr\tGLEAN\tintrons\t".($t[4]+1)."\t".($arr[-2]-1)."\n";
		#}
		push @arr, ( $t[3], $t[4] );
	}
}
