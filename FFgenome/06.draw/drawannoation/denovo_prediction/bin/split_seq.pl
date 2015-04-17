#!/usr/bin/perl
use strict;
use Getopt::Long;
my %opts;

my $usage=<<USAGE; #show how to use this program 

Program: split large sequence into smaller fragments, either based on the length of gap Ns or 
on the length of fragments, the two constrains can be used seperately or together.
This program is designed specially for deno-predicting software like genscan, which will get 
much slowly when meet large sequence size, usally > 1Mb.

Contact: Fan wei,<fanw\@genomics.org.cn>

Usage: $0  <seq.fa>

	-N<num>  number of consecutive Ns for split site

	-len<num> set length of the result fragments

	-detail		output middle information to screen

	-help		output help information to screen

Example: 
	perl split_seq_on_Ns_len.pl -N 3000 -len 10000000 seq.fa > seq.fa.cut
	perl split_seq_on_Ns_len.pl -N 3000 seq.fa  > seq.fa.cut
	perl split_seq_on_Ns_len.pl -len 10000000  > seq.fa.cut

USAGE

GetOptions(\%opts,"N:i","len:i","detail!","help!");
die $usage if ( @ARGV==0 || defined($opts{"help"}));

####################################################
################# Main  Function ###################
####################################################
my $file=shift;
my $cut_N=$opts{N};
my $cut_len=$opts{len};


open(IN, $file) || die ("can not open $file\n");
$/=">"; <IN>; $/="\n";
while (<IN>) {
	my $chr=$1 if(/^(\S+)/);
	$/=">";
	my $seq=<IN>;
	chomp $seq;
	$seq=~s/\s//g;
	$seq = uc $seq;
	$/="\n";
	
	print STDERR "\nSplit $chr\n" if (exists $opts{detail});
	
	my @frag;
	if (exists $opts{N}) {
		my ($start,$end,$frag,$pos);
		while ($seq=~s/^N//g) {
			$pos++;
		}
		while ($seq=~s/^([^N]+)(N*)//) {
			my $len_A=length($1);
			my $len_N=length($2);
			$pos += $len_A+$len_N;
			if ($len_N >= $cut_N || !$seq) {
				$frag .= $1;
				$end=$pos-$len_N;
				my $len_frag=length($frag);
				$start = $end-$len_frag+1;
				push @frag,[$start,$end,$frag];
				$frag="";
			}else{
				$frag .= $1.$2;
			}
		
		}
		if (! exists $opts{len}) {
			my $output;
			foreach my $p (@frag) {
				my $disp_frag;
				Disp_seq(\$p->[2],\$disp_frag);
				$output .= ">$chr"."_$p->[0]"."_$p->[1]\n".$disp_frag;
			}
			print $output;
		}
		
		my $frag_num = @frag;
		print STDERR "\tBy option -N, $chr generate fragments $frag_num\n" if (exists $opts{detail} && ! exists $opts{len});
	}


	if (exists $opts{len}) {
		my $output;
		if (! @frag) {
			my $frag_num;
			for (my $i=0; $i<length($seq); $i+=$cut_len) {
				my $str=substr($seq,$i,$cut_len*1.1);
				my $disp_str;
				Disp_seq(\$str,\$disp_str);
				my $sstart=$i+1;
				my $send=$sstart + length($str) - 1;

				last if(length($str) <= $cut_len*0.1 && $i>0);
				
				$output .= ">$chr"."_$sstart"."_$send\n".$disp_str;	
				$frag_num++; 
			}
			print STDERR "\tBy option -len, $chr generate fragments $frag_num\n" if (exists $opts{detail} && ! exists $opts{N});

			
		}else{
			my $frag_num_by_N = @frag;
			my $frag_num_by_len=0;
			foreach my $p (@frag) {
				my $frag_start = $p->[0];
				my $frag_end = $p->[1];
				my $frag_str = $p->[2];
				if ($frag_end-$frag_start+1 >= $cut_len*1.5) {
					$frag_num_by_len--;
					for (my $i=0; $i<length($frag_str); $i+=$cut_len) {
						my $str=substr($frag_str,$i,$cut_len*1.1);
						my $disp_str;
						Disp_seq(\$str,\$disp_str);
						my $sstart=$i+$frag_start;
						my $send=$sstart + length($str) - 1;

						last if(length($str) <= $cut_len*0.1 && $i>0);
						$frag_num_by_len++;
						$output .= ">$chr"."_$sstart"."_$send\n".$disp_str;
					}
				
				}else{
					my $disp_frag;
					Disp_seq(\$frag_str,\$disp_frag);
					$output .= ">$chr"."_$p->[0]"."_$p->[1]\n".$disp_frag;
				}
				
			}
			my $frag_num_total = $frag_num_by_N + $frag_num_by_len;
			print STDERR "\tBy option -N,   $chr generate fragments $frag_num_by_N\n" if (exists $opts{detail});
			print STDERR "\tBy option -len, $chr increase fragments $frag_num_by_len\n" if (exists $opts{detail});
			print STDERR "\tBy two options, $chr generate fragments $frag_num_total\n" if (exists $opts{detail});
		
		}

		print $output;
	}
	
	
}
close(IN);





#############################################
sub Disp_seq{
	my $seq_pp=shift;
	my $disp_pp=shift;
	my $num_line=(@_) ? shift : 50;
	
	my $len=length($$seq_pp);
	for (my $i=0; $i<$len; $i+=$num_line) {
		my $sub=substr($$seq_pp,$i,$num_line);
		$$disp_pp .= $sub."\n";
	}
	
}
#############################################