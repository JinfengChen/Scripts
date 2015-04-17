#!/usr/bin/perl

=head1 Name

initialize_taxonomy.pl  --  initialize the taxdmp data for later application

=head1 Description

This program is used to to initialize the taxdmp data, to creat
some useful files and stored in subdirectory "dat/".
When want to update the taxonomy datase, you should re-run this program.

Current database version is 2007-12-24. The input files are names.dmp, nodes.dmp,
division.dmp, gencode.dmp, gencode.dmp;  The output files are id_names,
id_parents, id_information.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage
  
  perl initialize_taxonomy.pl


=head1 Exmple



=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;


##get options from command line into variables and set default values
my ($Verbose,$Help);
GetOptions(
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
die `pod2text $0` if ($Help);



####################################################
################# Main  Function ###################
####################################################

print STDERR "\n\n\tcreat id_names.....\n\n";
&getId_Names();

print STDERR "\tcreat id_parents.....\n\n";
&getId_parents();

print STDERR "\tcreat id_information.....\n\n";
&getId_Information();

print STDERR "\tInitialization finished.\n\n";

####################################################
################# Main  Function ###################
####################################################



#Get tax_ids and all their names 
##########################################
sub getId_Names{

	my $infile="$Bin/../dat/names.dmp";
	my $outfile=">$Bin/../dat/id_names";
	my %id_names;
	open (IN,$infile) || die "fail open $infile\n";
		
	#find all nodes 
	while (<IN>) {
		if (/^(\d+)/) {
			$id_names{$1}={};
		}
	}

	#find all the names of each nodes
	seek(IN,0,0);
	while (<IN>) {
		if (/^(\d+)\s+\|\s+(.+?)\s+\|[^\|]+\|\s+(.+?)\s+\|/) {
			$id_names{$1}{$2}=$3; #$1:id, $2:name $3:type
		}
	}

	close(IN);


	#print out the result
	open (OUT,$outfile) || die ("can not open $outfile\n");
	select(OUT);

	foreach  (sort {$a<=>$b} keys %id_names) {
		print $_." | ";
		my $pp=$id_names{$_};
		foreach  (sort keys %$pp) {
			print $_." * ".$$pp{$_}." | ";
		}
		print "\n";

	}

	close(OUT);

}
##########################################



#Get tax_ids and all their parent tax_ids
##########################################
sub getId_parents{
	my $infile="$Bin/../dat/nodes.dmp";
	my $outfile=">$Bin/../dat/id_parents";
	my %forward; #correlation of each two nodes
	my %nodes; #all nodes and their parent nodes

	open (IN,$infile) || die ("can not open $infile\n");
	while(<IN>) {
		if (/^(\d+)\D+(\d+)/) {
			$forward{$1}=$2;
			$nodes{$1}=[$1];
		}	
	}
	close(IN);


	foreach my $id (sort keys %nodes) {
		my $i=0;
		while($nodes{$id}[$i]!=1) {
			$i++;	
			if ($forward{$nodes{$id}[$i-1]}) { #exsist
				$nodes{$id}[$i]=$forward{$nodes{$id}[$i-1]};
			}else{
				last; #if can not find, then last
			}
		}
	}


	#print out the result
	open (OUT,$outfile) || die ("can not open $outfile\n");
	select(OUT);

	foreach  my $key (sort {$a<=>$b} keys %nodes ) {
		my $pp=$nodes{$key};
		my $test=pop @$pp;
		push @$pp,$test;
		if ($test!=1) {
			print STDERR "root not 1\n"; #test whether there is nodes doesn't have root
		}
		foreach  (@$pp) {
			print $_."\t";
		}
		print "\n";
	}

	close(OUT);




}
##########################################




##########################################
sub getId_Information{
	my $names_file="$Bin/../dat/names.dmp";
	my $nodes_file="$Bin/../dat/nodes.dmp";
	my $division_file="$Bin/../dat/division.dmp";
	my $gcode_file="$Bin/../dat/gencode.dmp";
	my $out_file=">$Bin/../dat/id_information";
	my %id_info;

	#read names.dmp
	open (IN,$names_file) || die "fail open $names_file\n";

	#find all nodes 
	while (<IN>) {
		if (/^(\d+)/) {
			$id_info{$1}=[];
		}
	}
	
	#find scientific name of all nodes
	seek(IN,0,0);
	while (<IN>) {
		my $line=$_;
		if (/scientific name\s+\|$/) {
			if ($line=~/^(\d+)\s+\|\s+(.+?)\s+\|/) {
				$id_info{$1}[0]=$2;
			}

		}
	}

	close(IN);
	#print STDERR "\n\tnames.dmp read\n";
	
	#read nodes.dmp
	open (IN,$nodes_file) || die "fail open $nodes_file\n";
	
	while (<IN>) {
		if (/^(\d+)\s+\|[^\|]+\|\s+(.+?)\s+\|[^\|]+\|\s+(.+?)\s+\|[^\|]+\|\s+(.+?)\s+\|/) {
			if (exists $id_info{$1}) {
				$id_info{$1}[1]=$2;
				$id_info{$1}[2]=$3;
				$id_info{$1}[3]=$4;
			}
		}
	}
	
	close(IN);

	#print STDERR "\n\tnodes.dmp read\n";
	
	#read division.dmp 
	my %divid;
	open (IN,$division_file) || die "fail open $division_file\n";
	while (<IN>) {
		if (/^(\d+)\s+\|[^\|]+\|\s+(.+?)\s+\|/) {
			$divid{$1}=$2;
		}
	}
	close(IN);
	#print STDERR "\n\tdivision.dmp  read\n";

	#read gencode.dmp
	my %gc;
	open (IN,$gcode_file) || die "fail open $gcode_file\n";
	while (<IN>) {
		if (/^(\d+)\s+\|[^\|]+\|\s+(.+?)\s+\|/) {
			$gc{$1}=$2;
		}
	}
	
	close(IN);
	#print STDERR "\n\tgencode.dmp read\n";

	#change division_id and gcode_id to names
	foreach  (keys %id_info) {
		$id_info{$_}[2]=$divid{$id_info{$_}[2]};
		$id_info{$_}[3]=$gc{$id_info{$_}[3]};
	}

	#print out the result
	open (OUT,$out_file) || die "fail open $out_file\n";
	select(OUT);
	foreach  (sort {$a<=>$b} keys %id_info) {
		print $_." | ".$id_info{$_}[0]." | ".$id_info{$_}[1]." | ".$id_info{$_}[2]." | ".$id_info{$_}[3]." |\n";
	}
	
	close(OUT);

	#print STDERR "\n\tFinished\n";

}
##########################################


__END__
