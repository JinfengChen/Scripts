#!/usr/local/bin/perl

use Getopt::Std;
use strict;



my $usage = <<__EOUSAGE;

usage: $0 [-m max_intervening_genes] < match_file

match file is in dagchainer input format:

scaffold_id_A <tab> gene_acc_A <tab> end5_A <tab> end3_A <tab> scaffold_id_B <tab> gene_acc_B <tab> end5_B <tab> end3_B <tab> E-value



__EOUSAGE

    ;

our ($opt_m, $opt_v);

&getopts ('m:v');

my $MAX_NUM_INTERVENING = $opt_m || 1;
my $DEBUG = $opt_v;

## Read in Coordfile:
my %data; # structure: data{chromo}->[list of genes]

my %seen;

my %blast;

while (<STDIN>) {
    #print;
    chomp;
    my ($chromo_A, $gene_id_A, $c1_A, $c2_A, $chromo_B, $gene_id_B, $c1_B, $c2_B, $e_value) = split (/\t/);
    
    foreach my $entry ( [$chromo_A, $gene_id_A, $c1_A, $c2_A],
                        [$chromo_B, $gene_id_B, $c1_B, $c2_B] ) {
        
        my ($chromo, $gene_id, $c1, $c2) = @$entry;

        if ($seen{$gene_id}) { next; } # already seen it
        
        $seen{$gene_id} = 1;
        
        my $list_ref = $data{$chromo};
        unless (ref $list_ref) {
            $list_ref = $data{$chromo} = [];
        }
        my $mid = ($c1+$c2)/2;
        my $gene_struct = {
            gene_id=>$gene_id,
            c1=>$c1,
            c2=>$c2,
            mid=>$mid,
            chr=>$chromo};
        push (@$list_ref, $gene_struct);
        
    }

    my $existing_E_value = $blast{$gene_id_A}->{$gene_id_B};
    if (  (! defined($existing_E_value)) || ($existing_E_value > $e_value) ) {
        $blast{$gene_id_A}->{$gene_id_B} = $e_value;
        $blast{$gene_id_B}->{$gene_id_A} = $e_value;
    }
}



## Look for tandemly duplicated genes:
# allow one non-homolog to intervene a cluster of duplicates.
foreach my $chromo (sort keys %data) {
    print "Analyzing chromo: $chromo\n" if $DEBUG;
    my @all_tandem_arrays;
    ## sort genes by midpoint:
    my $gene_list_aref = $data{$chromo};
    @$gene_list_aref = sort {$a->{mid}<=>$b->{mid}} @$gene_list_aref;
    ## assign gene their indices:
    for (my $i = 0; $i <= $#{$gene_list_aref}; $i++) {
	$gene_list_aref->[$i]->{index} = $i;
    }
    
    my @curr_tandem_array;
    my $curr_gene = $gene_list_aref->[0];
    my $i = 1;
    my $num_intervening = 0;
    while ($i <= $#{$gene_list_aref}) {
	my $next_gene = $gene_list_aref->[$i];
	
	my $curr_gene_id = $curr_gene->{gene_id};
	my $next_gene_id = $next_gene->{gene_id};
	print "comparing $curr_gene_id to $next_gene_id\n" if $DEBUG;
	## check if next gene is homolog:
	my $evalue = $blast{$curr_gene_id}->{$next_gene_id};
	if ((defined ($evalue)) && ($evalue <= 1e-20) && ($num_intervening <= $MAX_NUM_INTERVENING)) {
	    ## got tandem duplicate:
	    print "* got tandem duplicate.\n" if $DEBUG;
	    push (@curr_tandem_array, [$curr_gene, $next_gene]);
	    $curr_gene = $next_gene;
	    $i++;
	} elsif
	    # See if next gene is a homolog:
	     ($num_intervening < $MAX_NUM_INTERVENING) {
		$num_intervening++;
		$i++;
		print "Adding a gap.  num_intervening = $num_intervening\n" if $DEBUG;
		## increment position, keep curr_gene the same.
	    } else {
		print "-resetting the automoton.\n" if $DEBUG;
		## reset the automoton:
		if (@curr_tandem_array) {
		    push (@all_tandem_arrays, [@curr_tandem_array]);
		    @curr_tandem_array = ();
		}
		
		## reset position 
		$num_intervening = 0;
		# move curr_gene to curr_gene + 1, reset i to next gene.
		my $curr_gene_pos = $curr_gene->{index};
		$curr_gene = $gene_list_aref->[$curr_gene_pos+1];
		$i = $curr_gene_pos + 2;
	    }
    }
    
    if (@curr_tandem_array) {
	push (@all_tandem_arrays, [@curr_tandem_array]);
    }
    
    &process_tandem_dups(@all_tandem_arrays);
}



####
sub process_tandem_dups {
    my @dup_listings = @_;
    foreach my $dup_listing (@dup_listings) {
        my @gene_pair_list = @$dup_listing;
        my $num_genes_in_tandem_dup = $#gene_pair_list + 2;
        print "## tandem duplication (n=$num_genes_in_tandem_dup)\n";
        foreach my $gene_pair (@gene_pair_list) {
            my ($g1, $g2) = @$gene_pair;
            my ($gene1_id, $gene2_id) = ($g1->{gene_id}, $g2->{gene_id});
            my ($g1_chromo, $g1_c1, $g1_c2, $g2_chromo, $g2_c1, $g2_c2, $e_value) = 
                ($g1->{chr},
                 $g1->{c1},
                 $g1->{c2},
                 $g2->{chr},
                 $g2->{c1},
                 $g2->{c2},
             $blast{$gene1_id}->{$gene2_id});
            
            printf ("%s\t%s\t%i\t%i\t%s\t%s\t%i\t%i\t%e\n", $g1_chromo, $gene1_id, $g1_c1, $g1_c2,
                    $g2_chromo, $gene2_id, $g2_c1, $g2_c2, $e_value);
        }
    }
}

	



