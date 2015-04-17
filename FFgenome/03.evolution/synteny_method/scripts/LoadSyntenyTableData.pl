#!/usr/bin/perl -w

=head1 NAME

LoadSyntenyTableData.pl - Populates a EnsEMBL compara database with synteny table data from Josh Stein's synteny analysis

=head1 SYNOPSIS

perl LoadSyntenyTableData.pl [options]

Options:
 
 -h|--help
 -r|--registry <file>
 -q|--query <registry_file_species>
 -t|--target <registry_file_species>
 -c|--compara <registry_file_compara_db_species> [optional]
 -o|--ort <ort_file> [optional]
 -d|--dist <dist_thresh> [optional]
 -l|--load

=head1 OPTIONS

B<-h|--help>
  Print a brief help message and exits.

B<-r|--registry_file> I<file>
  Use this Ensembl registry file for database connection info.

B<-q>|B<--query_species> I<registry_file_species>
  the species name in the registry_file for one of the synteny species, should correspond to the 1st species in the ort file

B<-t>|B<--target_species> I<registry_file_species>
  the species name in the registry_file for the other synteny species, should correspond to the 2nd species in the ort file

B<-c>|B<--compara> I<registry_file_species>
  the species name in the registry_file for the compara db

B<-o>|B<--ort> I<ort_file>
  the ort file 
  
B<-d>|B<--dist> I<dist_thresh>
  the distance threshold
  
=head1 DESCRIPTION

B<This program> 

  Populates the compara database with syntenic genes calculated from
  Josh Stein's synteny analysis. 

B<The ort file>

  The format of the file is (whitespace separated);
gene_stable_id sequence_name gene1_start gene1_end specie_id_1 homolog2_member_id sequence_name2 gene2_start gene2_end homology_relationship specie_id_2 

B<The Ensembl Registry>

  The database connection details for the Ensembl compara database
  must be configured in an ensembl registry file.

  An example registry file would be;
  ---
  use Bio::EnsEMBL::DBSQL::DBAdaptor; 
  use Bio::EnsEMBL::Registry;
  # The Ensembl core database
  Bio::EnsEMBL::DBSQL::DBAdaptor->new
  ( '-species' => 'rice', 
    '-group'   => 'core', 
    '-dbname'  => 'oryza_sativa_core_28_17', );
  ---

  The script would then be called using
  shell> perl load_protein_features.pl -s=rice -l=whatever


B<Restoring the database>

  If the script bombs out half-way through, your database will be in a
  partial loaded state (i.e. in a bit of a mess). Here is some SQL to
  help put it right;

  TODO: Add the SQL!

Created by Zhenyuan Lu <luj@cshl.edu>

=cut
use strict;
#use FindBin;
#use lib qq($FindBin::Bin/../modules);
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Analysis::RunnableDB::SyntenyTable;
use Bio::EnsEMBL::Analysis;
use Getopt::Long;
use Pod::Usage;

my ($HELP, $ORT, $QUERY, $TARGET, $COMPARA, $REGISTRY, $DIST, $LOAD);
my @ort;
GetOptions(
	'h|help'		=> \$HELP,
	'o|ort=s'     	=> \$ORT,
	'q|query=s'		=> \$QUERY,
	't|target=s'    => \$TARGET,
	'c|compara=s'	=> \$COMPARA,
	'r|registry=s'	=> \$REGISTRY,
	'd|dist=i'		=> \$DIST,
	'l|load'		=> \$LOAD
) or pod2usage();;

pod2usage(1) if $HELP;
$QUERY or print "Need query species\n" && pod2usage();
$TARGET or print "Need target species\n" && pod2usage();
$REGISTRY or print "Need registry file\n" && pod2usage();

$COMPARA ||= 'compara';

if ( $ORT ) {
	open ORT, $ORT or print "[*DIE] can't oprn $ORT\n" && pod2usage();
	while (<ORT>) {
		chomp;
		push @ort, $_;
	}
}

Bio::EnsEMBL::Registry->load_all($REGISTRY);
my $db = Bio::EnsEMBL::Registry->get_DBAdaptor( $QUERY, 'core' );
my $compara_db = Bio::EnsEMBL::Registry->get_DBAdaptor($COMPARA, 'compara');
my $target_db = Bio::EnsEMBL::Registry->get_DBAdaptor($TARGET,'core');

my $analysis=Bio::EnsEMBL::Analysis->new(
	-program         => "/data/ware/gramene-pipeline/bin/synteny_table.pl",
);
my $synteny_table = Bio::EnsEMBL::Analysis::RunnableDB::SyntenyTable->new(
	-analysis => $analysis,
	-input_id => 'synteny_table_analysis',
	-db	=> $db,
	-compara_db => $compara_db,
	-target_db => $target_db,
	-dist_thresh  => $DIST,
	-ort	=> scalar(@ort) ? \@ort : undef,
);

$synteny_table->fetch_input;
$synteny_table->run;
$LOAD and $synteny_table->write_output;

1;
