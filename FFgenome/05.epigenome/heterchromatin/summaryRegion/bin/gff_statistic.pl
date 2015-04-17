#!/usr/bin/perl
use strict;
use warnings;

my $file=$ARGV[0]; ##gff file
my $outdir=$ARGV[1]; ## directory where to put results files
my %seqinfor;
#`mkdir $outdir` unless (-f $outdir);

open IN, $file || die "$!";
while(<IN>){
	chomp;
	my @c=split/\t/,$_;
	if(@c > 4){
		if($c[2]=~/mRNA/){
			my $id=$1 if($c[8]=~/^ID=(\S+?);/);
			push @{$seqinfor{$id}{mRNA}},(@c);
		}
		if($c[2]=~/CDS/ || $c[2]=~/exon/){
                	my $id=$1 if($c[8]=~/^Parent=(\S+?);/);
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

open OUT1,">$outdir/mRNA_length.txt" || die "can't open:$!";
open (OUT2,">$outdir/average.txt") || die "can't open:$!";
open (OUT3,">$outdir/cds_length.txt") || die "can't open:$!";
open (OUT4,">$outdir/exon_length.txt") || die "can't open:$!";
open (OUT5,">$outdir/exon_number.txt") || die "can't open:$!";
open (OUT6,">$outdir/intron_length.txt") || die "can't open:$!";
open (OUT7,">$outdir/intron1_length.txt") || die "can't open:$!";
open (OUT8,">$outdir/intron2_length.txt") || die "can't open:$!";

foreach my $title(keys %seqinfor){
	$gene_number++;
	my $mRNA_length=$seqinfor{$title}{mRNA}[4]-$seqinfor{$title}{mRNA}[3]+1;
        my $strand=$seqinfor{$title}{mRNA}[6];
        print "$strand\n";
	$total_mRNA_length=$total_mRNA_length+$mRNA_length;
	print OUT1 "$title\t$mRNA_length\n";
	my $cds_length=0;
	my $number=@{$seqinfor{$title}{CDS}};
	@{$seqinfor{$title}{CDS}}=sort {$a->[3]<=>$b->[3]}@{$seqinfor{$title}{CDS}};
	$total_exon_number=$total_exon_number+$number;
	for(my $i=0;$i<$number;$i++){
		my $exon_length=$seqinfor{$title}{CDS}[$i][4]-$seqinfor{$title}{CDS}[$i][3]+1;
		$cds_length=$cds_length+$exon_length;
		print OUT4 "$title\t$exon_length\n";
                ## plus strand
		if($strand =~ /\+/ and $i == 1){ ## first intron
                        my $intron1_length=$seqinfor{$title}{CDS}[$i][3]-$seqinfor{$title}{CDS}[$i-1][4]-1;
                        print OUT7 "$title\t$intron1_length\n";
                        #$total_intron_length=$total_intron_length+$intron1_length;
                }elsif($strand =~ /\+/ and $i > 1){  ### other intron
                        my $intron2_length=$seqinfor{$title}{CDS}[$i][3]-$seqinfor{$title}{CDS}[$i-1][4]-1;
                        print OUT8 "$title\t$intron2_length\n";
                        #$total_intron_length=$total_intron_length+$intron2_length;

                }
                ## minus strand
                if($strand =~ /\-/ and $i == $number-1 and $i >= 1){
                        my $intron1_length=$seqinfor{$title}{CDS}[$i][3]-$seqinfor{$title}{CDS}[$i-1][4]-1;
                        print OUT7 "$title\t$intron1_length\n";
                        #$total_intron_length=$total_intron_length+$intron1_length;
                }elsif($strand =~ /\-/ and $i < $number-1 and $i >= 1){
                        my $intron2_length=$seqinfor{$title}{CDS}[$i][3]-$seqinfor{$title}{CDS}[$i-1][4]-1;
                        print OUT8 "$title\t$intron2_length\n";
                        #$total_intron_length=$total_intron_length+$intron2_length;

                }
                ## all intron
                if($i >= 1){
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

#print OUT2 "the total number of gene is $gene_number\nthe average of mRNA_length is $average_mRNA_length\nthe average of cds_length is $average_cds_length\nthe average of exon_number is $average_exon_number\nthe average of exon_length is $average_exon_length\nthe average of intron_length is $average_intron_length\nthe average of intron_length is $average_intron_length2\nthe total number of exon is $total_exon_number\nthe total number of intron is $total_intron_number\nthe total intron length is $total_intron_length\nthe total intron length is $total_intron_length2\n";
print OUT2 "gene_number\ttotal_mRNA_length\taverage_mRNA_length\taverage_cds_length\taverage_exon_number\taverage_exon_length\taverage_intron_length\taverage_intron_length2\ttotal_exon_number\ttotal_intron_number\ttotal_intron_length\ttotal_intron_length2\n";
print OUT2 "$gene_number\t$total_mRNA_length\t$average_mRNA_length\t$average_cds_length\t$average_exon_number\t$average_exon_length\t$average_intron_length\t$average_intron_length2\t$total_exon_number\t$total_intron_number\t$total_intron_length\t$total_intron_length2\n";		
close OUT1;
close OUT2;
close OUT3;
close OUT4;
close OUT5;
close OUT6;
close OUT7;
close OUT8;

