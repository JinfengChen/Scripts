#!/usr/bin/perl
use strict;
use warnings;

my $file=shift;
my $project=shift;
my %seqinfor;

open IN,$file || die "perl $0 rice.gff rice\n";
while(<IN>){
	chomp;
	my @c=split/\t/,$_;
	if(@c > 4){
		if($c[2]=~/mRNA/){
			my $id=$1 if($c[8]=~/^ID=(\S+?);/);
			push @{$seqinfor{$id}{mRNA}},(@c);
		}
		if($c[2]=~/CDS/ || $c[2]=~/exon/){
                	my $id=$1 if($c[8]=~/^Parent=(\S+?);/ or $c[8]=~/^Parent=(\S+?)$/);
                	push @{$seqinfor{$id}{CDS}},[@c];
        	}
	}
}
close IN;

my $gene_number=0;
my $total_mRNA_length=0;
my $total_exon_number=0;
my $total_cds_length=0;
my $total_intron_length=0;

open OUT1,">$project.mRNA_length.txt" || die "can't open:$!";
open (OUT2,">>average.txt") || die "can't open:$!";
open (OUT3,">$project.cds_length.txt") || die "can't open:$!";
open (OUT4,">$project.exon_length.txt") || die "can't open:$!";
open (OUT5,">$project.exon_number.txt") || die "can't open:$!";
open (OUT6,">$project.intron_length.txt") || die "can't open:$!";
print OUT2 "$project\n";
foreach my $title(keys %seqinfor){
	$gene_number++;
	my $mRNA_length=$seqinfor{$title}{mRNA}[4]-$seqinfor{$title}{mRNA}[3]+1;
	$total_mRNA_length=$total_mRNA_length+$mRNA_length;
	print OUT1 "$title\t$mRNA_length\n";
	my $cds_length=0;
	my $number=@{$seqinfor{$title}{CDS}};
	@{$seqinfor{$title}{CDS}}=sort {$a->[3]<=>$b->[3]} @{$seqinfor{$title}{CDS}};
	$total_exon_number=$total_exon_number+$number;
	for(my $i=0;$i<$number;$i++){
		my $exon_length=$seqinfor{$title}{CDS}[$i][4]-$seqinfor{$title}{CDS}[$i][3]+1;
		$cds_length=$cds_length+$exon_length;
		print OUT4 "$title\t$exon_length\n";
		if($i>=1){
			my $intron_length=$seqinfor{$title}{CDS}[$i][3]-$seqinfor{$title}{CDS}[$i-1][4]-1;
			print OUT6 "$title\t$intron_length\n";
			$total_intron_length=$total_intron_length+$intron_length;
		}
	}
	$total_cds_length=$total_cds_length+$cds_length;
	print OUT3 "$title\t$cds_length\n";
	print OUT5 "$title\t$number\n";
}
my $average_mRNA_length=$total_mRNA_length/$gene_number;
my $average_exon_length=$total_cds_length/$total_exon_number;
my $average_cds_length=$total_cds_length/$gene_number;
my $average_exon_number=$total_exon_number/$gene_number;
my $total_intron_number=$total_exon_number-$gene_number;
my $average_intron_length=$total_intron_length/$total_intron_number;
my $total_intron_length2=$total_mRNA_length-$total_cds_length;
my $average_intron_length2=$total_intron_length2/$total_intron_number;

print OUT2 "the total number of gene is $gene_number\nthe average of mRNA_length is $average_mRNA_length\nthe average of cds_length is $average_cds_length\nthe average of exon_number is $average_exon_number\nthe average of exon_length is $average_exon_length\nthe average of intron_length is $average_intron_length\nthe average of intron_length is $average_intron_length2\nthe total number of exon is $total_exon_number\nthe total number of intron is $total_intron_number\nthe total intron length is $total_intron_length\nthe total intron length is $total_intron_length2\n";
		
