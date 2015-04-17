#!/usr/bin/perl

=head1 Name

add_treefam_species_mark.pl  -- convert gene_id to the treefam format 

=head1 Description

The treefam gene_id format is like ENSANGT00000000024.3_ANOGA, the last five 
capital characters is a shorcut of a species name. This program is used to
both cds and pep files.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2008-6-3

=head1 Usage

  $0 <fasta_file> <species_mark>
  --format    set input file format, fasta or table, default fasta
  --verbose   output verbose information to screen  
  --help      output help information to screen  

=head1 Exmple

 perl add_treefam_species_mark.pl Populus_trichocarpa.cds POPTR > Populus_trichocarpa.cds.mark
 perl add_treefam_species_mark.pl Populus_trichocarpa.pep POPTR > Populus_trichocarpa.pep.mark

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;

my ($Format);
my ($Verbose,$Help);
GetOptions(
	"format:s"=>\$Format,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Format ||= "fasta";
die `pod2text $0` if (@ARGV != 2 || $Help);

my $file = shift;
my $species_mark = shift;

if ($Format eq "fasta") {

	my %data;
	Read_fasta($file,\%data);

	foreach my $name (sort keys %data) {
		my $head = $data{$name}{head};
		my $seq = $data{$name}{seq};
		$head =~ s/^$name/$name\_$species_mark/;
		print ">$head\n$seq";
	}

}

if ($Format eq "table") {

	open IN, $file || die "fail $file";
	while (<IN>) {
		my @t = split /\t/;
		$t[0] = "$t[0]_$species_mark";
		my $line = join("\t",@t);
		print $line;
	}
	close IN;

}


####################################################
################### Sub Routines ###################
####################################################




#read fasta file
#usage: Read_fasta($file,\%hash);
#############################################
sub Read_fasta{
	my $file=shift;
	my $hash_p=shift;
	
	my $total_num;
	open(IN, $file) || die ("can not open $file\n");
	$/=">"; <IN>; $/="\n";
	while (<IN>) {
		chomp;
		my $head = $_;
		my $name = $1 if($head =~ /^(\S+)/);
		
		$/=">";
		my $seq = <IN>;
		chomp $seq;
		#$seq=~s/\s//g;
		$/="\n";
		
		if (exists $hash_p->{$name}) {
			warn "name $name is not uniq";
		}

		$hash_p->{$name}{head} =  $head;
		$hash_p->{$name}{len} = length($seq);
		$hash_p->{$name}{seq} = $seq;

		$total_num++;
	}
	close(IN);
	
	return $total_num;
}
