#!/usr/bin/perl 
use strict;

##author: fanwei, fanw@genomics.org.cn
##date: 2009-3-29

die "perl $0 AXTfile > outfile\n" unless( @ARGV == 1);

my %CODE = (
			
				'GCA' => 'A', 'GCC' => 'A', 'GCG' => 'A', 'GCT' => 'A',                               # Alanine
				'TGC' => 'C', 'TGT' => 'C',                                                           # Cysteine
				'GAC' => 'D', 'GAT' => 'D',                                                           # Aspartic Acid
				'GAA' => 'E', 'GAG' => 'E',                                                           # Glutamic Acid
				'TTC' => 'F', 'TTT' => 'F',                                                           # Phenylalanine
				'GGA' => 'G', 'GGC' => 'G', 'GGG' => 'G', 'GGT' => 'G',                               # Glycine
				'CAC' => 'H', 'CAT' => 'H',                                                           # Histidine
				'ATA' => 'I', 'ATC' => 'I', 'ATT' => 'I',                                             # Isoleucine
				'AAA' => 'K', 'AAG' => 'K',                                                           # Lysine
				'CTA' => 'L', 'CTC' => 'L', 'CTG' => 'L', 'CTT' => 'L', 'TTA' => 'L', 'TTG' => 'L',   # Leucine
				'ATG' => 'M',                                                                         # Methionine
				'AAC' => 'N', 'AAT' => 'N',                                                           # Asparagine
				'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P',                               # Proline
				'CAA' => 'Q', 'CAG' => 'Q',                                                           # Glutamine
				'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R', 'CGT' => 'R', 'AGA' => 'R', 'AGG' => 'R',   # Arginine
				'TCA' => 'S', 'TCC' => 'S', 'TCG' => 'S', 'TCT' => 'S', 'AGC' => 'S', 'AGT' => 'S',   # Serine
				'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T', 'ACT' => 'T',                               # Threonine
				'GTA' => 'V', 'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V',                               # Valine
				'TGG' => 'W',                                                                         # Tryptophan
				'TAC' => 'Y', 'TAT' => 'Y',                                                           # Tyrosine
				'TAA' => 'U', 'TAG' => 'U', 'TGA' => 'U',                                              # Stop
				'---' => '-'
			
	);

my $axtFile = shift;

open(AXT,"$axtFile")||die"Cannot open $axtFile\n";
$/ = "\n\n";
my @seqs = <AXT>;
$/ ="\n";
close AXT;

print "tag\tcds_identity\taa_identity\n";
foreach my $line ( @seqs ){
        chomp $line;
        if( $line =~ /^(\S+)\n(\S+)\n(\S+)$/ ){
                my $tag = $1;
                my $seq1 =$2;
                my $seq2 =$3;
                my $cds_identity = calculate_identity($seq1, $seq2);
				my $pep1 = cds2aa($seq1);
				my $pep2 = cds2aa($seq2);
				my $aa_identity = calculate_identity($pep1, $pep2);
               ##print "$tag\n$seq1\n$seq2\n\n";
			   ##print "$tag\n$pep1\n$pep2\n\n";
				print "$tag\t$cds_identity\t$aa_identity\n";
        }
}


##this function can be applied to both cds and pep
sub calculate_identity {
	my ($seq1,$seq2) = @_;
	
	my $base_len = 0;
	my $iden_len = 0;
	my $len = length($seq1);
	for (my $i=0; $i<$len; $i++) {
		my $base1 = substr($seq1,$i,1);
		my $base2 = substr($seq2,$i,1);
		next if($base1 eq "-" || $base2 eq "-");
		$base_len++;
		$iden_len++ if($base1 eq $base2);
	}
	my $iden_rate = ($base_len > 0) ? $iden_len/$base_len : "None";
	return $iden_rate;

}


## translate CDS to pep
####################################################
sub cds2aa {
	my $seq = shift;
	
	my $len = length($seq);
	
	my $prot;
	for (my $i=0; $i<$len; $i+=3) {
		my $codon = substr($seq,$i,3);
		last if(length($codon) < 3);
		$prot .= (exists $CODE{$codon}) ? $CODE{$codon} : 'X';
	}
	return $prot;

}