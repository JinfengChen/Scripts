#!/usr/bin/perl

=head1 Name

taxonomy_dictionary.pl  --  get detailed information for a taxonomy node

=head1 Description

This is a dictionary program, a name word or id number is required as input of the program.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage
  
  perl taxonomy_dictionary.pl <input>
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

  perl bin/taxonomy_dictionary.pl human
  perl bin/taxonomy_dictionary.pl 9606
  perl bin/taxonomy_dictionary.pl Homo sapiens

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
die `pod2text $0` if (@ARGV == 0 || $Help);


####################################################
################# Main  Function ###################
####################################################
my $input=join(" ", @ARGV);
my $tax_id; #global variable

if ($input=~/^\d+/) {
	$tax_id=$input;
	query_id($tax_id);
} else{
	$tax_id=query_name($input);	
}

&query_information($tax_id);
&query_lineage($tax_id);

####################################################
################# Main  Function ###################
####################################################





#find the given name in id_names and return the tax_id
#Two ranks search, if first search fail, then do the 
#second  search
##########################################
sub query_name{
	my $query=shift; #name word
	print "\nYour species is: ".$query."\n";
	$query=lc($query);
	my $id=0; 
	my %pp_lines; #store the pp_lines which can partially pattern $query
	my $infile="$Bin/../dat/id_names";
	open (IN, $infile) || die ("can not find /dat/id_names\n");
	
	while (<IN>) {
		my $line=$_;
		$line=lc($line);
		if ($line=~/$query/) {
			$pp_lines{$line}=1;
			
			#first search
			if ($line=~/\|\s$query\s\*/) {    
				if ($line=~/^(\d+)/) {
					$id=$1;
					my @ary=split(/\|/,$line);
					print "\t tax_id: $ary[0]\n";
					for (my $i=1; $i<@ary; $i++) {
						print "\t".$ary[$i]."\n";
					}

				}
			}

		}
	}

	
	#second search
	if ($id==0) {
		print  "\nNot found, try the listed name again\n\n";
	}
	
	if ($id==0) {
		foreach  (keys %pp_lines) {
			my $line=$_;
			if (/\|\s([^\|]*$query[^\|]*)\s\|/) {  #second search
				print "\t".$1."\n";
				delete($pp_lines{$_});

			}
		}
		print "\n\n";
		close(IN);
		exit();

	}
	

	close(IN);

	return $id;

}
##########################################



#find the given tax_id in id_names and return names
##########################################
sub query_id{
	my $query=shift;
	print "\nYour tax_id is: ".$query."\n\n";
	my $infile="$Bin/../dat/id_names";
	open (IN, $infile) || die ("can not find /dat/id_names\n");
	my $id=0;
	while (<IN>) {
		my $line=$_;
		if (/^(\d+)/) {
			if ($1 eq $query) {
				$id=$1;
				my @ary=split(/\|/,$line);
				shift @ary; #rm tax_id
				foreach  (@ary) {
					print "\t".$_."\n";
				}
				last;
			}
		}
	}
	
	if (!$id) {
		print "\tError: invalid tax_id\n\n";
		exit();
	}
	
	
	close(IN);

}
##########################################



##########################################
sub query_information{
	my $id=shift;
	my $infile="$Bin/../dat/id_information";
	open (IN, $infile) || die ("can not find /dat/id_information\n");
	while (<IN>) {
		my $line=$_;
		if (/^$id\s/) {
			my @ary=split(/\|/,$line);
			print "\tScientific name:\t".$ary[1]."\n";
			print "\tClassification rank:\t".$ary[2]."\n";
			print "\tMain  division:\t\t".$ary[3]."\n";
			print "\tGenetic code:\t\t".$ary[4]."\n";
			last;
		}
	}
	close(IN);

}
##########################################



##########################################
sub query_lineage{
	my $id=shift;
	my @ary; #store parent id chain
	my $infile="$Bin/../dat/id_parents";
	open (IN, $infile) || die ("can not find /dat/id_parents\n");
	while (<IN>) {
		my $line=$_;
		if (/^$id\s/) {
			@ary=split(/\s+/,$line);
			@ary=reverse @ary;
			print "\nTax_id lineages:\n";
			foreach  (@ary) {
				print "\t\*".$_."\*";
			}
			print "\n\n";

			last;
		}
	}
	close(IN);
	
	
	my %id_name; #store id and scientific name
	foreach  (@ary) {
		$id_name{$_}=0;
	}
	my $infile="$Bin/../dat/id_names";
	open (IN, $infile) || die ("can not find /dat/id_names\n");
	while (<IN>) {
		my $line=$_;
		if (/^(\d+)/) {
			my $id=$1;
			if (exists $id_name{$1}) {
				if ($line=~/\|\s([^\|\*]+)\s\*\sscientific name/) {
					$id_name{$id}=$1;
				}
			}
		}
	}

	close(IN);

	print "\nScientific name lineages:\n";
	foreach  (@ary) {
		print "\t\*".$id_name{$_}."\*";
	}
	print "\n\n";


}
##########################################


__END__









