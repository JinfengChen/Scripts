#!/usr/bin/perl
use strict;

die "$0 <all_pep_file> <solar_file> <align_both_cutoff>\n" if(@ARGV != 3);

my $all_pep_file = shift;
my $solar_file = shift;
my $align_cutoff = shift; ##to both

my %Len;
Read_fasta($all_pep_file,\%Len);
warn "Read pep done";

open IN,$solar_file || die "fail";
while (<IN>) {
	my @t = split /\t/;
	my $query_id = $t[0];
	my $query_len = $Len{$query_id};
	
	my $target_id = $t[5];
	my $target_len = $Len{$target_id};
	
	my $query_align = 0;
	my $target_align = 0;
	while ($t[11] =~ /(\d+),(\d+);/g) {
		$query_align += abs($2 - $1) + 1;
	}
	while ($t[12] =~ /(\d+),(\d+);/g) {
		$target_align += abs($2 - $1) + 1;
	}
        #print STDERR "$query_id\t$query_len\t$target_id\t$target_len\n";     	
	if ($query_align / $query_len >= $align_cutoff && $target_align / $target_len >= $align_cutoff) {
		print $_;
	}
}
close IN;


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
		
		$/="\n>";
		my $seq = <IN>;
		chomp $seq;
		$seq=~s/\s//g;
		$/="\n";
		
		if (exists $hash_p->{$name}) {
			warn "name $name is not uniq";
		}

		$hash_p->{$name} = length($seq);

		$total_num++;
	}
	close(IN);
	
	return $total_num;
}
