#!/usr/bin/perl

=head1 Name

reform the glean gff to structure format 

=head1 Description

add the index function

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2008-9-7
  Note:

=head1 Usage

  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple



=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory

##get options from command line into variables and set default values
my ($Verbose,$Help);
GetOptions(
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
die `pod2text $0` if (@ARGV == 0 || $Help);

my %Gene;
my %Index;

##glean gff format
while (<>) {
	chomp;
	my @t = split /\t/;
	my $type = $t[2];
	next if($type ne "CDS");
	next if(/^\#/);
	my $gene_id = $1 if($t[8] =~ /GenePrediction (\S+)/ || $t[8] =~ /Parent=([^;]+);?/ || $t[8] =~ /Gene (\S+)/ || $t[8] =~ /CDS (\w+)/);
	my $chr = $t[0];
	my $strand = $t[6];
	my $cds_start = $t[3];
	my $cds_end = $t[4];
	next if($cds_start <= 0 || $cds_end <= 0);
	$Gene{$gene_id}{"chr"} = $chr;
	$Gene{$gene_id}{"strand"} = $strand;
	push @{$Gene{$gene_id}{"cds"}},[$cds_start,$cds_end] if($cds_start && $cds_end);

}

#print Dumper \%Gene;
#exit;

##sort the exon order and prepare the index
foreach my $gene_id (sort keys %Gene) {
	my $gene_p = $Gene{$gene_id};
	
	##sort to regular exon order whith coordinate on the + strand
	foreach my $p (@{$gene_p->{cds}}) {
		if($p->[0] > $p->[1]){
			($p->[0],$p->[1]) = ($p->[1],$p->[0]);
		}
		
	}	
	@{$gene_p->{cds}} = reverse @{$gene_p->{cds}} if($gene_p->{cds}[0][0] > $gene_p->{cds}[-1][0]);
	
	my $exon_num = @{$gene_p->{cds}};
	
	my $output = "$gene_p->{chr}  $gene_p->{strand}  $exon_num  ";
	
	foreach my $p (@{$gene_p->{cds}}) {
		$output .= "$p->[0],$p->[1];";
	}

	$gene_p->{output} = $output;
	
	push @{$Index{$gene_p->{"chr"}}}, [$gene_p->{cds}[0][0],$gene_id];
}

##sort and make the index
foreach my $chr (sort keys %Index) {
	my $chr_p = $Index{$chr};
	my $index_num = 1;
	foreach my $p (sort {$a->[0] <=> $b->[0]} @$chr_p) {
		my $start_position = $p->[0];
		my $gene_id = $p->[1];
		$Gene{$gene_id}{start_position} = $start_position;
		$Gene{$gene_id}{index_num} = $index_num;
		$index_num++;
	}
}

##output the final result
foreach my $gene_id (sort keys %Gene) {
	my $gene_p = $Gene{$gene_id};
	
	print "$gene_id\t$gene_p->{output}  $gene_p->{start_position}  $gene_p->{index_num}\n";
}


####################################################
################### Sub Routines ###################
####################################################
