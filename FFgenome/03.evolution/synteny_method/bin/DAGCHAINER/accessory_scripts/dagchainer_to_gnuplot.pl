#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case bundling);


my $usage = <<_EOUSAGE_;


############################################################
#
#  Required:  
#
#  --genomeA    fasta file for genome A
#  --genomeB    fasta file for genome B
#  --aligncoords   DAGchainer aligncoords output file
# 
# Optional:
#  --prefix     prefix for plotting files (default: 'gp')
#  --min_diag_genes     minumum number of genes of a diagonal to plot
#  --accession      plot only those chains involving molecule accession
#
########################################################################

_EOUSAGE_

	;


my ($genomeA, $genomeB, $aligncoords, $help, $accession);

my $min_diag_genes = 2; #default minimum number of genes that can constitute a chain.
my $prefix = "gp";

&GetOptions ('h' => \$help,
			 'genomeA=s' => \$genomeA,
			 'genomeB=s' => \$genomeB,
			 'aligncoords=s' => \$aligncoords, 
			 'prefix=s' => \$prefix,
			 'min_diag_genes=i' => \$min_diag_genes,
			 'accession=s' => \$accession,
			 );

if ($help) { die $usage; }

unless ($genomeA && $genomeB && $aligncoords) { 
	die $usage;
}

main: {

	my %seqLengthsA = &parse_seqlengths_from_fasta($genomeA);
	my %seqLengthsB = &parse_seqlengths_from_fasta($genomeB);

	my @DAGchain_structs = &parse_DAGchainer_output($aligncoords);  
	
	if ($genomeA eq $genomeB) {
		@DAGchain_structs = &reciprocate_chains(\@DAGchain_structs);
	}

	
	my (%contigsA, %contigsB); # store the accessions for those entries encountered.
	@DAGchain_structs = &filter_DAGchains(\@DAGchain_structs, \%contigsA, \%contigsB);
		
	## create fake scaffolding for a big XY plot.
	my %scaffoldsA = &create_scaffolding(\%seqLengthsA, \%contigsA);
	my %scaffoldsB = &create_scaffolding(\%seqLengthsB, \%contigsB);
	
		
	&plot_DAGchains($prefix, \@DAGchain_structs, \%scaffoldsA, \%scaffoldsB);
	
	&write_gnuplot_script($prefix, \%scaffoldsA, \%scaffoldsB);
		
	exit(0);
}


####
sub parse_seqlengths_from_fasta {
	my ($genome_file) = @_;

	my $acc = "";
	
	my %lengths;

	my $seq_lengths_file = "$genome_file.seq_lengths";
	if (-e $seq_lengths_file) {
		open (my $fh, $seq_lengths_file) or die "Error, cannot open file $seq_lengths_file";
		while (<$fh>) {
			chomp;
			my ($acc, $length) = split (/\t/);
			$lengths{$acc} = $length;
		}
		close $fh;
	}
	else {
		## No sequence lengths file exists
		#  parse the genome sequence and write the seq lengths file for next time's use.
		open (my $fh, $genome_file) or die "Error, cannot open file $genome_file";
		while (<$fh>) {
			if (/^>(\S+)/) {
				$acc = $1;
				
				print STDERR "\r$genome_file: computing length for $acc         ";
				
			}
			else {
				s/\s+//g;
				my $len = length($_);
				if ($len && ! $acc) {
					die "Error, no accession read yet, but have length $len for $_";
				}
				$lengths{$acc} += $len;
			}
		}
		
		close $fh;
		
		print STDERR "\n";
		
		# write the seqlengths file
		open ($fh, ">$seq_lengths_file") or die "Error, cannot write file $seq_lengths_file";
		foreach my $acc (keys %lengths) {
			print $fh "$acc\t$lengths{$acc}\n";
		}
		close $fh;
	}
	
	return (%lengths);
	
}


####
sub parse_DAGchainer_output {
	my ($dagchainer_outfile) = @_;

	my @DAGchain_structs;
	
	###########################################################################################
	#
	# atts:
	#
	#   accA => ""
	#   accB => ""
	#   match_pairs => [   
	#                        [ [matchA1_lend, matchA1_rend] , [matchB1_lend, matchB1_rend] ],
	#                        [  [matchA2_lend, matchA2_rend] , [matchB2_lend, matchB2_rend] ],
	#                       
	#                  ]
	#   rangeA => [ far_A_lend, far_A_rend]
	#   rangeB => [ far_B_lend, far_B_rend]
	#
	#   orient => +|-     # for rangeB, since rangeA is always +
	#
	##########################################################################################
	
	
	my $dagchain;

	open (my $fh, $dagchainer_outfile) or die "Error, cannot open $dagchainer_outfile";
	while (<$fh>) {
		chomp;
		
		if (/^\#/) {
			## start a new dagchain
			$dagchain = { 
				match_pairs => [] 
				};
			
			push (@DAGchain_structs, $dagchain);
			
			my $line = $_;

			# parse molecule accessions:
			$line =~ /(\S+) vs\. (\S+)/;
			
			my ($contigA, $contigB) = ($1, $2);
			
			$dagchain->{accA} = $contigA;
			$dagchain->{accB} = $contigB;
			
			if ($line =~ /\(reverse\)/) {
				# reverse strand for molecule B
				$dagchain->{orient} = '-';
			}
			else {
				$dagchain->{orient} = '+';
			}
		}
		elsif (/\w/) {
			# add entry to chain:
			my ($contigA, $geneA, $end5A, $end3A, 
				$contigB, $geneB, $end5B, $end3B,
				$score) = split (/\t/);
			
			
			my $match_pairs_aref = $dagchain->{match_pairs};
			
			push (@$match_pairs_aref, [ [$end5A, $end3A],
										[$end5B, $end3B] ] );
	
		}

		
	}



	## get chain spans:
	foreach my $dagchain (@DAGchain_structs) {

		my $contigA = $dagchain->{accA};
		my $contigB = $dagchain->{accB};

		my $match_pairs_aref = $dagchain->{match_pairs};

		my @mol_A_coords;
		my @mol_B_coords;

		foreach my $match_pair (@$match_pairs_aref) {
			my ($match_A, $match_B) = @$match_pair;

			push (@mol_A_coords, @$match_A);
			push (@mol_B_coords, @$match_B);

		}

		@mol_A_coords = sort {$a<=>$b} @mol_A_coords;
		@mol_B_coords = sort {$a<=>$b} @mol_B_coords;

		my $lend_A = shift @mol_A_coords;
		my $rend_A = pop @mol_A_coords;

		my $lend_B = shift @mol_B_coords;
		my $rend_B = pop @mol_B_coords;
		
		
		$dagchain->{rangeA} = [ $lend_A, $rend_A ];
		$dagchain->{rangeB} = [ $lend_B, $rend_B ];
		
	}


	return (@DAGchain_structs);
}


sub filter_DAGchains {
	my ($DAGchain_structs_aref, $contigsA_href, $contigsB_href) = @_;
	
	my @filtered_structs;
	foreach my $dagchain (@$DAGchain_structs_aref) {
		
		if ($accession && $dagchain->{accA} ne $accession) { next; }
		
		$contigsA_href->{ $dagchain->{accA} } = 1;
		$contigsB_href->{ $dagchain->{accB} } = 1;
		
		my $match_pairs_aref = $dagchain->{match_pairs};
		my $num_match_pairs = scalar (@$match_pairs_aref);
		if ($num_match_pairs >= $min_diag_genes) {
			push (@filtered_structs, $dagchain);
		}
	}
	

	return (@filtered_structs);
}


####
sub create_scaffolding {
	my ($seqLengths_href, $contigs_href) = @_;
	
	my %scaffolding;

	my $seq_length = 1;
	
	foreach my $contig (reverse sort {$seqLengths_href->{$a}<=>$seqLengths_href->{$b}} keys %$contigs_href) {
		
		$scaffolding{$contig} = $seq_length;

		my $contig_len = $seqLengths_href->{$contig} or die "Error, no length for contig $contig";
		
		$seq_length += $contig_len;
	}
	
	return (%scaffolding);
}




####
sub plot_DAGchains {
	my ($prefix, $DAGchain_structs_aref,
		$scaffoldsA_href, $scaffoldsB_href) = @_;
	
	open (my $dots_forward_fh, ">$prefix.dots_F") or die $!;
	open (my $dots_reverse_fh, ">$prefix.dots_R") or die $!;
	open (my $chains_forward_fh, ">$prefix.chains_F") or die $!;
	open (my $chains_reverse_fh, ">$prefix.chains_R") or die $!;

	foreach my $chain (@$DAGchain_structs_aref) {
		my $accA = $chain->{accA};
		my $accB = $chain->{accB};

		my $scaff_pos_A = $scaffoldsA_href->{$accA} or die "Error, no scaff position for A: $accA";
		my $scaff_pos_B = $scaffoldsB_href->{$accB} or die "Error, no scaff position for B: $accB";
		
		my $match_pairs_aref = $chain->{match_pairs};
		my $rangeA = $chain->{rangeA};
		my $rangeB = $chain->{rangeB};
		
		my ($range_A_lend, $range_A_rend) = @$rangeA;
		my ($range_B_lend, $range_B_rend) = @$rangeB;
		
		$range_A_lend += $scaff_pos_A;
		$range_A_rend += $scaff_pos_A;

		$range_B_lend += $scaff_pos_B;
		$range_B_rend += $scaff_pos_B;
		
		my $orient = $chain->{orient};
		
		## make the chain bounds line:
		my $chain_fh = $chains_forward_fh;;
		my $dots_fh = $dots_forward_fh;
		if ($orient eq '-') {
			$chain_fh = $chains_reverse_fh;
			$dots_fh = $dots_reverse_fh;
			($range_B_lend, $range_B_rend) = ($range_B_rend, $range_B_lend);
		}
		
		print $chain_fh "$range_A_lend\t$range_B_lend\n"
			. "$range_A_rend\t$range_B_rend\n\n";
		
		## draw the dots
		foreach my $match (@$match_pairs_aref) {
			my ($geneA_coords, $geneB_coords) = @$match;
			my ($gene_A_end5, $gene_A_end3) = @$geneA_coords;
			my ($gene_B_end5, $gene_B_end3) = @$geneB_coords;

			$gene_A_end5 += $scaff_pos_A;
			$gene_A_end3 += $scaff_pos_A;

			$gene_B_end5 += $scaff_pos_B;
			$gene_B_end3 += $scaff_pos_B;
			
			print $dots_fh "$gene_A_end5\t$gene_B_end5\n"
				. "$gene_A_end3\t$gene_B_end3\n\n";
		}
	}
			
	return;
}


####
sub write_gnuplot_script {
	my ($prefix, $scaffoldsA_href, $scaffoldsB_href) = @_;

	my $gnuplot_script_name = "$prefix.gnuplot_script";

	open (my $fh, ">$gnuplot_script_name") or die $!;
	
	my @xtics = &build_tics($scaffoldsA_href);
	my @ytics = &build_tics($scaffoldsB_href);

	print $fh "set size 1,1\n"
		. "set grid\n"
		. "unset key\n"
		. "set border 10\n"
		. "set ticscale 0 0\n"
		. "set xrange [0:]\n"
		. "set yrange [0:]\n"
		. "set xtics rotate\n";
	

	print $fh "set xtics (\\\n" . join (", \\\n", @xtics) . ")\n";
	print $fh "set ytics (\\\n" . join (", \\\n", @ytics) . ")\n";

	print $fh "set style line 1  lt 1 lw 2 pt 6 ps 1\n"
		. "set style line 2  lt 2 lw 2 pt 6 ps 1\n"
		. "set style line 3  lt 3 lw 2 pt 6 ps 1\n"
		. "set style line 4  lt 4 lw 2 pt 6 ps 1\n";

	my @files_to_plot;
	push (@files_to_plot, " \"$prefix.chains_F\" w lp ls 1 ") if (-s "$prefix.chains_F");
	push (@files_to_plot, " \"$prefix.chains_R\" w lp ls 2 ") if (-s "$prefix.chains_R");
	push (@files_to_plot, " \"$prefix.dots_F\" w lp ls 3 ") if (-s "$prefix.dots_F");
	push (@files_to_plot, " \"$prefix.dots_R\" w lp ls 4 ") if (-s "$prefix.dots_R");
	
	if (@files_to_plot) {
		print $fh "plot " . join (", \\\n", @files_to_plot) . "\n\n";
	}
	else {
		die "Error, no data points to plot.\n\n";
	}

		
	close $fh;


	#print "\n\nrun:\n\t% gnuplot $gnuplot_script_name -\n\n\n";
	
	system ("gnuplot $gnuplot_script_name - ");
	
	return;
}

####
sub build_tics {
	my ($scaffolds_href) = @_;

	my @tics;

	foreach my $acc (sort {$scaffolds_href->{$a}<=>$scaffolds_href->{$b}} keys %$scaffolds_href) {
		push (@tics, "\"$acc\" $scaffolds_href->{$acc}");
	}

	return (@tics);
}
		
####
sub reciprocate_chains {
	my ($DAGchain_structs_aref) = @_;

	my @ret_structs;

	foreach my $DAGchain (@$DAGchain_structs_aref) {
		
		my $newchain = {
			accA => $DAGchain->{accB},
			accB => $DAGchain->{accA},

			rangeA => $DAGchain->{rangeB},
			rangeB => $DAGchain->{rangeA},

			orient => $DAGchain->{orient},

			match_pairs => [],

		};
		
		my $match_pairs_aref = $DAGchain->{match_pairs};
		foreach my $match_pair (@$match_pairs_aref) {
			my ($match_pairA, $match_pairB) = @$match_pair;
			
			push (@{$newchain->{match_pairs}}, [ $match_pairB, $match_pairA ] );


		}
		
		push (@ret_structs, $DAGchain, $newchain);

	}

	return (@ret_structs);
}

