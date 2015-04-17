#!/usr/bin/perl

=head1 Name

sim4_to_gff3.pl  --  convert blat result from sim4 format to gff3 format

=head1 Description

The score of mRNA is cacluted by identity * alignrate, you can choose the one 
with the best score as the best locus of this mRNA. However, in fact, some genes
have many duplications, and they may be quite similar to each other. So it is 
better to keep the duplicated genes, you can use the --tophit option to cut down
the number of duplicated genes.

This program is used for both EST and cDNA alignment results. It will automatically 
remove small introns less than 10bp, which is mostly caused by sequencing indels.
For full length cDNAs, if the  file [mrna.source] is give, then the CDS will be inferred.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6

=head1 Usage
  
  perl sim4_to_gff3.pl [options] <file.sim4> [mrna.source]

  --identity <num>  set identity cutoff, default 0.98
  --alignrate <num> set alignrate cutoff, default 0.98
  --tophit <num>    set the number of top hits, default no limitation
  --prefix <str>    set a prefix name for the gene ID in gff3
  --verbose         output running progress information to screen  
  --help            output help information to screen  

=head1 Exmple

 perl ../bin/sim4_to_gff3.pl rice_BGI_indicaEST_10000.fa.blat.sim4 > rice_BGI_indicaEST_10000.fa.blat.sim4.gff 

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;


my ($Identity_cutoff,$Alignrate_cutoff,$Prefix_tag,$Best_hit);
my ($Verbose,$Help);
GetOptions(
	"tophit:i"=>\$Best_hit,
	"identity:s"=>\$Identity_cutoff,
	"alignrate:s"=>\$Alignrate_cutoff,
	"prefix:s"=>\$Prefix_tag,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Prefix_tag = $Prefix_tag."_" if(defined $Prefix_tag);
$Identity_cutoff ||= 0.98;
$Alignrate_cutoff ||= 0.98;
die `pod2text $0` if (@ARGV == 0 || $Help);

my $infile = shift;
my $mrna_source_file = shift;

my %Data;  ## store the main data
my %Count; ## store matching number
my %mRNA;  ## store coding and annotation information for each mRNA

my $mini_intron_size = 10;

##read mRNA source file
Read_mRNA($mrna_source_file,\%mRNA) if(-f $mrna_source_file);

parse_sim4($infile,\%Data);

##print Dumper \%Data;
#
##output result in gff3 format
print "##gff-version 3\n";
foreach my $cDNA (sort keys %Data) {
	my $cDNA_p = $Data{$cDNA};
	my $hit_num;
	foreach my $match_p (sort {$b->[0] <=> $a->[0]} @$cDNA_p) {
		$hit_num++;
		my ($mRNA_score,$match_id,$seq2_id,$strand,$gene_start,$gene_end,$exon_p,$cds_p) = @$match_p;
		my $gene_id = $Prefix_tag.$cDNA;
		$gene_id .= "_D".$hit_num if($hit_num >= 2);
		
		my $output;
		$output .= "$seq2_id\tblat\tmRNA\t$gene_start\t$gene_end\t$mRNA_score\t$strand\t.\tID=$gene_id;Target=$cDNA $exon_p->[0][0] $exon_p->[-1][1];";
		$output .= "function=\"".$mRNA{$cDNA}{function}."\";" if(defined $mRNA{$cDNA}{function});
		$output .= "\n";
		
		foreach my $p (@$exon_p) {
			my ($exon_start,$exon_end,$exon_score) = ($p->[2],$p->[3],$p->[4]); 
			$output .= "$seq2_id\tblat\tCDS\t$exon_start\t$exon_end\t$exon_score\t$strand\t.\tParent=$gene_id;Target=$cDNA $p->[0] $p->[1];\n";
		}

		foreach my $p (@$cds_p) {
			my ($cds_start,$cds_end,$cds_score) = ($p->[2],$p->[3],$p->[4]); 
			$output .= "$seq2_id\tblat\tCDS\t$cds_start\t$cds_end\t$cds_score\t$strand\t.\tParent=$gene_id;Target=$cDNA $p->[0] $p->[1]; \n";
		}
		
		print $output;
		last if(defined $Best_hit && $hit_num >= $Best_hit);
	}
}

####################################################
################### Sub Routines ###################
####################################################



##alert that this subroutine contains many globle variables
####################################################
sub parse_sim4{
	my $infile = shift;
	my $data_p = shift;

	$/ = "seq1 = ";
	open IN,$infile || die "fail $infile"; 
	<IN>;
	while (<IN>) { 
		chomp;
		my $unit = "seq1 = ".$_;
		##print $unit."\n#################\n";

		my ($seq1_id,$seq2_id,$strand,$seq1_len,$seq2_len);
		my @exon;
		my ($math_len,$identity,$align_len,$align_rate);
		my ($gene_start,$gene_end);

		($seq1_id,$seq1_len) = ($1,$2) if($unit =~ /seq1 = (\S+), (\d+) bp/); 
		($seq2_id,$seq2_len) = ($1,$2) if($unit =~ /seq2 = (\S+), (\d+) bp/); 
		$strand = ($unit =~ /\(complement\)/) ? "-" : "+";
		
		while ($unit =~ /(\d+)-(\d+)  \((\d+)-(\d+)\)   (\d+)\%/g) {
			push @exon,[$1,$2,$3,$4,$5]; ## query_start, query_end, subject_start, subject_end, identity
			$math_len += ($2 - $1 + 1) * $5 / 100;
			$align_len += $2 - $1 + 1;
		}
		
		##remove intron <= 10 bp
		remove_small_intron(\@exon,$mini_intron_size);

		$identity = $math_len / $align_len;
		$align_rate = $align_len / $seq1_len;
		($gene_start,$gene_end) = ($exon[0][2],$exon[-1][3]);
				
		##find hits better than cutoff and infer cds
		if ($identity >= $Identity_cutoff && $align_rate >= $Alignrate_cutoff) {
			$Count{$seq1_id}++; 
			
			my $gene_id = $seq1_id;
			$gene_id .= "_match".$Count{$seq1_id} if($Count{$seq1_id} >= 2); ## add _match if have more than one locus

			my $mRNA_score = sprintf("%.4f",$identity * $align_rate);

			## if the mRNA source file is provided, infer cds structure
			my @cds;
			infer_cds(\@exon,\@cds,$strand,$seq1_len,$mRNA{$seq1_id}{orf_start},$mRNA{$seq1_id}{orf_end}) if (defined $mRNA{$seq1_id}{orf_start} && defined $mRNA{$seq1_id}{orf_end});	
		
			push @{$data_p->{$seq1_id}}, [$mRNA_score,$gene_id,$seq2_id,$strand,$gene_start,$gene_end,\@exon,\@cds];
		        #print STDERR "$seq1_id\t$mRNA_score\t$gene_id\t$seq2_id\t$gene_start\t$gene_end\n";
                }
		

	}
	close IN;

}



## remove small itrons and rebuild the exons
## my $exon = [ [$1,$2,$3,$4,$5],[$1,$2,$3,$4,$5],[$1,$2,$3,$4,$5],[$1,$2,$3,$4,$5] ]; ## query_start, query_end, subject_start, subject_end, identity
## remove_small_intron($exon,10);
####################################################
sub remove_small_intron {
	my $exon_p = shift;
	my $mimi_intron = shift || 10;

	my $do_status = 0;
	for (my $i=1; $i<@$exon_p; $i++) {
		if ($exon_p->[$i][2] - $exon_p->[$i-1][3] < $mimi_intron + 1) {
			$exon_p->[$i-1][3] = $exon_p->[$i][3];
			$exon_p->[$i-1][1] = $exon_p->[$i][1];
			$exon_p->[$i-1][4] = ($exon_p->[$i-1][4] + $exon_p->[$i][4]) / 2; ## a simple but not accurate caculation
			splice(@$exon_p,$i,1);
			$do_status = 1;
			last; ##删除一个就退出循环，每次只删除一个
		}
			
	}

	remove_small_intron($exon_p,$mimi_intron) if($do_status);
}




## query_start, query_end, subject_start, subject_end, identity
sub infer_cds {
	my ($exon_p,$cds_p,$strand,$mRNA_len,$orf_start1,$orf_end1) = @_;

	my ($orf_start,$orf_end);
	($orf_start,$orf_end) = ($mRNA_len - $orf_end1 + 1,$mRNA_len - $orf_start1 + 1) if($strand eq '-');
	($orf_start,$orf_end) = ($orf_start1,$orf_end1) if($strand eq '+');

	my ($mRNA_cds_start,$mRNA_cds_end,$genome_cds_start,$genome_cds_end);
	foreach my $p (@$exon_p) {
		next if($orf_start > $p->[1]);
		last if($orf_end < $p->[0]);
		$mRNA_cds_start = ( $p->[0] <= $orf_start && $orf_start <= $p->[1] ) ? $orf_start : $p->[0];
		$mRNA_cds_end = ( $p->[0] <= $orf_end && $orf_end <= $p->[1] ) ? $orf_end : $p->[1];
		
		$genome_cds_start = $p->[2] - $p->[0] + $mRNA_cds_start;
		$genome_cds_end = $p->[3] - $p->[1] + $mRNA_cds_end;
		
		push @$cds_p, [$mRNA_cds_start,$mRNA_cds_end,$genome_cds_start,$genome_cds_end,$p->[4]];


	}

}

#
#read fasta file
#usage: Read_fasta($file,\%hash);
#############################################

sub Read_mRNA{
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
		#chomp $seq;
		#$seq=~s/\s//g;
		$/="\n";

		$hash_p->{$name}{function} = $1 if($head =~ /\[.+?\]\s+(.+)/);
		if ($head =~ /\s+(\d+)\.\.(\d+)\s*/) {
			$hash_p->{$name}{orf_start} = $1;
			$hash_p->{$name}{orf_end} =   $2;
		}elsif($head =~ /,(\d+),(\d+),/){
			$hash_p->{$name}{orf_start} = $1;
			$hash_p->{$name}{orf_end} =   $2;
		}

		$total_num++;
	}
	close(IN);
	
}





__END__

1. gff3 format:
silk_chr1       genscan gene    24678   29209   .       +       .       ID=BMS000003; status=novel;
silk_chr1       genscan mRNA    24678   29209   .       +       .       ID=BMS000003-TA; Parent=BMS000003; status=novel;
silk_chr1       genscan CDS     24678   24695   -1.75   +       0       Parent=BMS000003-TA;
silk_chr1       genscan CDS     25189   25321   11.13   +       0       Parent=BMS000003-TA;
silk_chr1       genscan CDS     25393   25460   -2.09   +       1       Parent=BMS000003-TA;
silk_chr1       genscan CDS     25874   26064   14.78   +       0       Parent=BMS000003-TA;
silk_chr1       genscan CDS     26331   26389   -7.72   +       2       Parent=BMS000003-TA;
silk_chr1       genscan CDS     27430   27527   1.84    +       1       Parent=BMS000003-TA;
silk_chr1       genscan CDS     28158   28364   29.53   +       0       Parent=BMS000003-TA;
silk_chr1       genscan CDS     28444   28856   87.64   +       0       Parent=BMS000003-TA;
silk_chr1       genscan CDS     28972   29209   23.07   +       2       Parent=BMS000003-TA;



2. blat result in sim4 format:

seq1 = BMB000004-TA, 2412 bp
seq2 = silk_chr1, 19754638 bp

1-141  (30905-31045)   100% ->
142-820  (31723-32401)   100% ->
821-1002  (32482-32663)   100% ->
1003-1274  (33034-33305)   100% ->
1275-1447  (33409-33581)   100% ->
1448-1742  (33916-34210)   100% ->
1743-2166  (34497-34920)   100% ->
2167-2337  (36836-37006)   100% ->
2338-2412  (37094-37168)   100% ->

seq1 = BMB000005-TA, 2028 bp
seq2 = silk_chr1, 19754638 bp

(complement)
1-61  (43267-43327)   100% <-
62-197  (44204-44339)   100% <-
198-313  (44531-44646)   100% <-
314-1998  (44988-46672)   100% <-
1999-2028  (49634-49663)   100% <-

3. mRNA source format:
>AB001052  624  76..450  [Bombyx mori]  Bombyx mori mRNA for histone H2A-like protein, complete cds.
tacggacgttcgcgagacacgcgtgcgttcgagtgctttcgtgtgttata
tcgttaacttttttaaacttcaaacatgtccggtcgcggaaaaggcggaa
aagttaagggcaaggtcaagtcccgttcgaaccgtgccggtcttcagttt
ccggtcggtcgtatacacagattgttgcgcaacggaaattacgctgaacg
cgttggtgccggtgcaccggtttacctggccgccgtcatggaatacttgg
ccgctgaagttttggaattggccggtaacgcagcaagagacaacaagaag
actagaattattcctagacatcttcaactcgccataaggaacgacgagga
actgaacaaactcccttccggtgtgacaatcgctcaaggcggagttttac
caaacattcaagcggtactactcccgaagaagaccgagaagaaagcttaa
aaaacgcttcaaactcgctcgcaaggacaacaacaacaacacaacatgtc
gtcgataatatttatatatgtataataaataataataatcgacgaacata
ttttgttcgttgttgttattattattattatattgtattgttgacaaaaa
tcaaaggcccttttcagggccgct
>AB001053  344  64..318  [Bombyx mori]  Bombyx mori mRNA for VAP-peptide, complete cds.
aaaaaacacagcacttagctcatcggcagaacacatctagtttgttattt
gaaagaccgcaaaatgttcaagttgacagtaattttcgctattatcgctg
tggcccaagcgggcgtcatagccccagtggtgcctgtagcacaccccgtc
gtggctcacacggccgtggtccacccggtcccactagtgcgcgctgccca
cgtggttcacaccgccccagtggttgccgctgccccagtggtcgccgctc
cggtggtggctgcggctcctatcgtgccgatagttaaacatgcgccaa
