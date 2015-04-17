#!/usr/bin/perl

=head1 Name

findOverlap.pl  --  find overlap relations between two block sets on chromosome

=head1 Description

This program is designed to find overlap relations of blocks, the input file
can be in psl gff or table format. If not psl or gff format, then the program will
take it as table format, ie, tab delimitated text format.The program is originally 
desinged for finding overlap between gene sets, so it can recognize psl and gff 
formats. but it has a more common usage to find overalp between two sets of blocks, 
no matter what they are, so it is expanded to accept table format.

Inside this program, you will see two marks "ref" and "pre".  The "ref" means reference
genes file, the "pre" means predicted gene file, because the program is first used to 
find overlap between reference genes and predicted genes, to see how much difference 
between them. You can just imagine "ref" and "pre" as two block sets. If only one input 
file is given, the program will find overlap inside "ref" itself.
 
The algorithm is that: (1) sort the two block sets based on start positions, seperately;
(2) walk along the chromosome, and find out overlap between the two block sets.
(3) report who and who overlapped, as well as their own size and the overlapped size.

The format of output file: [example]
Variation_6833  25552   chr1    2       Variation_8315,46898,25552      Variation_6833,25552,25552
  Column 1: the query ID; 
  Column 2: the query size; 
  Column 3: chromosome ID; 
  Column 4: number of blocks overlapped with Variation_6833; 
  Column 5: the first subject block Variation_8315 overlapped with Variation_6833, 46898 is its own size, 25552 is the overlapped size; 
  Column 6: the second subject block Variation_6833 overlapped with Variation_6833, numbers has the same meaning as last column; 


=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2007-2-4

=head1 Usage
  
  perl findOverlap.pl <ref> [pre]
  --verbose   output verbose information to screen  
  --help      output help information to screen  

=head1 Exmple

  perl findOverlap.pl  gene.gff  mRNA.psl >  gene.gff.mRNA.overlap
  perl findOverlap.pl variation.txt > variation.txt.overlap

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;

my ($Verbose,$Help);
GetOptions(
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
die `pod2text $0` if (@ARGV < 1 || $Help);

my $ref_file=shift;
my $pre_file=shift || $ref_file;

my ( %Ref, %Pre );

if($ref_file=~/\.psl$/){
	read_psl($ref_file,\%Ref);
}elsif($ref_file=~/\.gff$/){
	read_gff($ref_file,\%Ref);
}else{
	read_table($ref_file,\%Ref);
}


print STDERR "read ref_file done\n" if($Verbose);


if($pre_file=~/\.psl$/){
	read_psl($pre_file,\%Pre);
}elsif($pre_file=~/\.gff$/){
	read_gff($pre_file,\%Pre);
}else{
	read_table($pre_file,\%Pre);
}

print STDERR "read pre_file done\n" if($Verbose);


find_overlap(\%Ref,\%Pre);

print STDERR "find overlap done\n" if($Verbose);


####################################################
################### Sub Routines ###################
####################################################



sub find_overlap{
	my $Ref_hp=shift;
	my $Pre_hp=shift;
	
	foreach  my $chr (sort keys %$Ref_hp) {
		
		my $output;
		my @ref_chr = (exists $Ref_hp->{$chr})  ? (sort {$a->[0] <=> $b->[0]} @{$Ref_hp->{$chr}}) : ();
		my @pre_chr = (exists $Pre_hp->{$chr})  ? (sort {$a->[0] <=> $b->[0]} @{$Pre_hp->{$chr}}) : ();
		
		print STDERR "find overlap on $chr\n" if($Verbose);
		
		my $pre_pos = 0;
		for (my $i=0; $i<@ref_chr; $i++) {
			my $ref_gene = $ref_chr[$i][2];
			my $ref_size = $ref_chr[$i][1] - $ref_chr[$i][0] + 1;
			my @overlap;
			
			for (my $j=$pre_pos; $j<@pre_chr; $j++) {
				if ($pre_chr[$j][1] < $ref_chr[$i][0]) {
					$pre_pos++;
					next;
				}
				if ($pre_chr[$j][0] > $ref_chr[$i][1]) {
					last;
				}
				
				my $pre_size = $pre_chr[$j][1] - $pre_chr[$j][0] + 1;
				my $overlap_size = overlap_size($pre_chr[$j],$ref_chr[$i]);
				
				push @overlap,"$pre_chr[$j][2],$pre_size,$overlap_size";
			}
			
			$output .= $ref_gene."\t".$ref_size."\t".$chr."\t".scalar(@overlap)."\t".join("\t",@overlap)."\n";
		}

		print $output;
	}

}


sub overlap_size {
	my $block1_p = shift;
	my $block2_p = shift;
	
	my $combine_start = ($block1_p->[0] < $block2_p->[0]) ?  $block1_p->[0] : $block2_p->[0];
	my $combine_end   = ($block1_p->[1] > $block2_p->[1]) ?  $block1_p->[1] : $block2_p->[1];
	
	my $overlap_size = ($block1_p->[1]-$block1_p->[0]+1) + ($block2_p->[1]-$block2_p->[0]+1) - ($combine_end-$combine_start+1);

	return $overlap_size;
}


sub read_gff{
	my $file=shift;
	my $ref=shift;

	open (IN,$file) || die ("fail open $file\n");
	while (<IN>) {
		chomp;
		s/^\s+//;
		my @t = split(/\t/);
		my $tname = $t[0];

		my $qname;
		if ($t[2] eq 'mRNA' || $t[2] eq 'CDS') {
			$qname = $1 if($t[8] =~ /GenePrediction\s+(\S+)/);
		}
		if ($t[2] eq 'match' || $t[2] eq 'HSP') {
			$qname = $1 if($t[8] =~ /Target\s+\"(\S+)\"/);
		}

		if ($t[2] eq 'mRNA' || $t[2] eq 'match') {
			my $start = $t[3];
			my $end = $t[4];
			push @{$ref->{$tname}},[$start,$end,$qname];
		}
	}
	close(IN);
	
}

sub read_psl{
	my $file=shift;
	my $ref=shift;
	open(REF,$file)||die("fail to open $file\n");
	while (<REF>) {
		chomp;
		my @temp=split(/\s+/,$_);
		my $chr=$temp[13];
		my $gene=$temp[9];
		my $start=$temp[15]+1;
		my $end=$temp[16];

		push @{$ref->{$chr}},[$start,$end,$gene];
	}
	close(REF);
}

sub read_table{
	my $file=shift;
	my $ref=shift;
	open(REF,$file)||die("fail to open $file\n");
	while (<REF>) {
		chomp;
		my @temp=split(/\t/,$_);
		my $chr=$temp[0];
		my $gene=$temp[1];
		my $start=$temp[2];
		my $end=$temp[3];
		
		next if($start !~ /\d/);
		push @{$ref->{$chr}},[$start,$end,$gene];
	}
	close(REF);
}

