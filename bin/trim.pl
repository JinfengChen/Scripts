#!/usr/local/bin/perl

# AUTHOR: Nik Joshi
# LAST REVISED: September 2010
#
# The Bioinformatics Core at UC Davis Genome Center
# http://bioinformatics.ucdavis.edu
# Copyright (c) 2010 The Regents of University of California, Davis Campus.
# All rights reserved.

use Getopt::Long;
use strict;

if ($#ARGV < 0) {
	print STDERR <<EOF;
Options:
--type <num>     0=standard trimming, 1=adaptive trimming, 2=windowed adaptive trimming.  Default 0
--qual-threshold <num>     quality threshold for trimming, default 20
--length-threshold <num>    length threshold for trimming, default 50
--percent-threshold <num>   percent threshold for trimming, default 0.1. only use in type 1 now
--trim-type <num>           trim using length (1) or percent (0), default 0, only use in type 1 now
--qual-type <num>           0=sanger qualities, 1=illumina qualities pipeline>=1.3, 2=illumina qualities pipeline<1.3.  Default 0.
--pair1 <paired end input filename>      fastq, paired end file. Must have same number of records as pair2. Required.
--pair2 <paired end input filename>      fastq, paired end file. Must have same number of records as pair1. Required.
--outpair1 <paired end output file>      Required.
--outpair2 <paired end output file>      Required.
--single <single end output file>        Required.
EOF
	exit(1);
}

my $type = 0;
my $qual_threshold = 20;
my $length_threshold = 50;
my $percent_threshold = 0.1;
my $trim_type = 0;
my $qual_type = 0;
my ($pair1, $pair2, $outpair1, $outpair2, $single);

GetOptions ('type=i' => \$type, 
		'qual-threshold=i' => \$qual_threshold,
		'length-threshold=i' => \$length_threshold,
                'percent-threshold' => \$percent_threshold,
                'trim-type' => \$trim_type,
		'qual-type=i' => \$qual_type,
		'pair1=s' => \$pair1,
		'pair2=s' => \$pair2,
		'outpair1=s' => \$outpair1,
                'outpair2=s' => \$outpair2,
		'single=s' => \$single);



if (!defined $pair1 || !defined $pair2 || !defined $outpair1 || !defined $outpair2 || !defined $single) {
	print STDERR "Error: You must define pair1, pair2, outpair1, outpair2, and single\n";
	exit(1);
}



sub get_quality_num {
	my ($qual_char, $qual_type) = @_;

	if ($qual_type == 0) {
		return (ord($qual_char) - 33);
	}

	elsif ($qual_type == 1) {
		return (ord($qual_char) - 64);
	}

	elsif ($qual_type == 2) {
		return (10 * log(1 + 10 ** ((ord($qual_char) - 64) / 10.0)) / log(10));
	}
}


open (PAIR1, "<$pair1") or die "Cannot access file $pair1.\n";
open (PAIR2, "<$pair2") or die "Cannot access file $pair2.\n";
open (OUTPAIR1, ">$outpair1");
open (OUTPAIR2, ">$outpair2");
open (SINGLE, ">$single");

my ($count, $p1_header1, $p1_seq, $p1_header2, $p1_qual, $p2_header1, $p2_seq, $p2_header2, $p2_qual);
my ($pos, $p1_flag, $qual_char, $qualnum, $p1_cut, $p2_flag, $p2_cut);
my ($paired_cnt, $single_cnt, $discard_cnt) = (0,0,0);
my ($i, $window_total, $window_size, $window_start, @p1_qual_data, @p2_qual_data);
my ($trimmed_p1_seq, $trimmed_p1_qual, $trimmed_p2_seq, $trimmed_p2_qual);

$count = 0;
while ($p1_header1 = <PAIR1>) {
	$p1_seq = <PAIR1>;
	$p1_header2 = <PAIR1>;
	$p1_qual = <PAIR1>;
	chomp $p1_seq;
	chomp $p1_qual;

	$p2_header1 = <PAIR2>;
        $p2_seq = <PAIR2>;
        $p2_header2 = <PAIR2>;
        $p2_qual = <PAIR2>;
        chomp $p2_seq;
        chomp $p2_qual;

	if ($type == 0) {
		$trimmed_p1_seq = substr ($p1_seq, 0, $length_threshold);
		$trimmed_p1_qual = substr ($p1_qual, 0, $length_threshold);
		print OUTPAIR1 "$p1_header1$trimmed_p1_seq\n$p1_header2$trimmed_p1_qual\n";

                $trimmed_p2_seq = substr ($p2_seq, 0, $length_threshold);
                $trimmed_p2_qual = substr ($p2_qual, 0, $length_threshold);
                print OUTPAIR2 "$p2_header1$trimmed_p2_seq\n$p2_header2$trimmed_p2_qual\n";

		$paired_cnt++;
	}

	elsif ($type == 1) {
                #print "$p1_flag\t$p2_flag\n";
		$pos=0;
		$p1_flag=1;
                $p1_cut=length $p1_seq;
                # if low quality sequence cover more than 10%, discard.
                $length_threshold=$trim_type == 1 ? $length_threshold : (length $p1_seq)*(1-$percent_threshold);
                #print "$length_threshold\n";
		foreach $qual_char (split (//, $p1_qual)) {
			$qualnum = get_quality_num ($qual_char, $qual_type);
                        #print "Q$qualnum ";
			if ($qualnum < $qual_threshold) {
				$p1_cut = $pos;
				if ($pos < $length_threshold) {
                                       
					$p1_flag = 0;
				}

				last;
			}

			$pos++;
		}
                #print "\n";
                #print "C1C2$p1_cut\t$p2_cut\n";
                $pos=0;
                $p2_flag=1;
                $p2_cut=length $p2_seq;
                foreach $qual_char (split (//, $p2_qual)) {
                        $qualnum = get_quality_num ($qual_char, $qual_type);
                        #print "Q$qualnum ";
                        if ($qualnum < $qual_threshold) {
                                $p2_cut = $pos;
                                if ($pos < $length_threshold) {
                                        $p2_flag = 0;
                                }

				last;
                        }

                        $pos++;
                }

                #print "\n";
                #print "$p1_flag\t$p2_flag\t$p1_cut\t$p2_cut\n";
		if ($p1_flag == 1 && $p2_flag == 1) {
			$trimmed_p1_seq = substr ($p1_seq, 0, $p1_cut);
			$trimmed_p1_qual = substr ($p1_qual, 0, $p1_cut);
			print OUTPAIR1 "$p1_header1$trimmed_p1_seq\n$p1_header2$trimmed_p1_qual\n";

                        $trimmed_p2_seq = substr ($p2_seq, 0, $p2_cut);
                        $trimmed_p2_qual = substr ($p2_qual, 0, $p2_cut);
                        print OUTPAIR2 "$p2_header1$trimmed_p2_seq\n$p2_header2$trimmed_p2_qual\n";

			$paired_cnt++;
		}

		elsif ($p1_flag == 1 && $p2_flag == 0) {
                        $trimmed_p1_seq = substr ($p1_seq, 0, $p1_cut);
                        $trimmed_p1_qual = substr ($p1_qual, 0, $p1_cut);
                        print SINGLE "$p1_header1$trimmed_p1_seq\n$p1_header2$trimmed_p1_qual\n";

			$single_cnt++;
		}

		elsif ($p1_flag == 0 && $p2_flag == 1) {
                        $trimmed_p2_seq = substr ($p2_seq, 0, $p2_cut);
                        $trimmed_p2_qual = substr ($p2_qual, 0, $p2_cut);
                        print SINGLE "$p2_header1$trimmed_p2_seq\n$p2_header2$trimmed_p2_qual\n";

			$single_cnt++;
                }

		else {$discard_cnt++;}
	}

	# sliding window
        elsif ($type == 2) {

		$window_size = int (0.1 * length ($p1_qual));
                $window_start = 0;
		$window_total = 0;
		$p1_flag = 1;
		$p1_cut = length($p1_qual);

		@p1_qual_data = split (//, $p1_qual);

		for ($i=0; $i<$window_size; $i++) {
			$window_total += get_quality_num ($p1_qual_data[$i], $qual_type);
		}

		for ($i=0; $i<length($p1_qual); $i++) {

			# if the average quality in the window is less than the threshold
			# or if the window is the last window in the read
			if ($window_total / $window_size < $qual_threshold || $window_start+$window_size > length ($p1_qual)) {

				# at what point in the window does the quality dip below the threshold?
				for (my $j=$window_start; $j<$window_start+$window_size; $j++) {
					if (get_quality_num ($p1_qual_data[$j], $qual_type) < $qual_threshold) {
						$p1_cut = $j;
						if ($p1_cut < $length_threshold) {$p1_flag = 0;}

						last;
					}
				}

				last;
			}

			# instead of sliding the window, subtract the first qual and add the next qual
			$window_total -= get_quality_num ($p1_qual_data[$window_start], $qual_type);
			$window_total += get_quality_num ($p1_qual_data[$window_start+$window_size], $qual_type);
			$window_start++;
		}


                $window_size = int (0.1 * length ($p2_qual));
                $window_start = 0;
                $window_total = 0;
                $p2_flag = 1;
		$p2_cut = length($p2_qual);

                @p2_qual_data = split (//, $p2_qual);

                for ($i=0; $i<$window_size; $i++) {
                        $window_total += get_quality_num ($p2_qual_data[$i], $qual_type);
                }

                for ($i=0; $i<length($p2_qual); $i++) {
                        if ($window_total / $window_size < $qual_threshold || $window_start+$window_size > length ($p2_qual)) {

                                # at what point in the window does the quality dip below the threshold?
                                for (my $j=$window_start; $j<$window_start+$window_size; $j++) {
                                        if (get_quality_num ($p2_qual_data[$j], $qual_type) < $qual_threshold) {
                                                $p2_cut = $j;
                                                if ($p2_cut < $length_threshold) {$p2_flag = 0;}

						last;
                                        }
                                }

				last;
                        }

                        # instead of sliding the window, subtract the first qual and add the next qual
                        $window_total -= get_quality_num ($p2_qual_data[$window_start], $qual_type);
                        $window_total += get_quality_num ($p2_qual_data[$window_start+$window_size], $qual_type);
                        $window_start++;
                }



                my ($trimmed_p1_seq, $trimmed_p1_qual, $trimmed_p2_seq, $trimmed_p2_qual);
                if ($p1_flag == 1 && $p2_flag == 1) {
                        $trimmed_p1_seq = substr ($p1_seq, 0, $p1_cut);
                        $trimmed_p1_qual = substr ($p1_qual, 0, $p1_cut);
                        print OUTPAIR1 "$p1_header1$trimmed_p1_seq\n$p1_header2$trimmed_p1_qual\n";

                        $trimmed_p2_seq = substr ($p2_seq, 0, $p2_cut);
                        $trimmed_p2_qual = substr ($p2_qual, 0, $p2_cut);
                        print OUTPAIR2 "$p2_header1$trimmed_p2_seq\n$p2_header2$trimmed_p2_qual\n";

                        $paired_cnt++;
                }

                elsif ($p1_flag == 1 && $p2_flag == 0) {
                        $trimmed_p1_seq = substr ($p1_seq, 0, $p1_cut);
                        $trimmed_p1_qual = substr ($p1_qual, 0, $p1_cut);
                        print SINGLE "$p1_header1$trimmed_p1_seq\n$p1_header2$trimmed_p1_qual\n";

                        $single_cnt++;
                }

                elsif ($p1_flag == 0 && $p2_flag == 1) {
                        $trimmed_p2_seq = substr ($p2_seq, 0, $p2_cut);
                        $trimmed_p2_qual = substr ($p2_qual, 0, $p2_cut);
                        print SINGLE "$p2_header1$trimmed_p2_seq\n$p2_header2$trimmed_p2_qual\n";

                        $single_cnt++;
                }

                else {$discard_cnt++;}
	}


	$count++;
	#if ($count % 10000 == 0) {
	#	print STDERR "Records: $count, Paired: $paired_cnt, Single: $single_cnt, Discarded: $discard_cnt\n";
	#}
}

print STDERR "Records: $count, Paired: $paired_cnt, Single: $single_cnt, Discarded: $discard_cnt\n";

close (PAIR1);
close (PAIR2);
close (OUTPAIR1);
close (OUTPAIR2);
close (SINGLE);

