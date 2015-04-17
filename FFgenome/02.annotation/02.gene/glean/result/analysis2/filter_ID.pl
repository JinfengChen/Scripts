#!/usr/bin/perl
use strict;
use warnings;

my $file1=shift;
my $file2=shift;

my %length;

open (IN1,$file1) || die "$!";
$/=">";<IN1>;$/="\n";
while(<IN1>){
	my $id=$1 if($_=~/^(\S+)/);
	$/=">";
	my $seq=<IN1>;
	chomp($seq);
	$seq=~s/\s+//g;
	my $seq_length=length($seq);
	$length{$id}=$seq_length;
	$/="\n";
}
close IN1;

my %ID;
open (IN2,$file2) || die "$!";
while(<IN2>){
	chomp;
	my @c=split/\t/,$_;
	my $id=$c[0];
	push @{$ID{$id}},[@c];
}
close IN2;

foreach my $title(keys %ID){
	@{$ID{$title}}=sort {$a->[6]<=>$b->[6]}@{$ID{$title}};
	my $number=@{$ID{$title}};
	for(my $j=0;$j<$number-1;$j++){
		my $i=$j+1;
		if($ID{$title}[$i][6]>=$ID{$title}[$j][6] && $ID{$title}[$i][6]<$ID{$title}[$j][7]){
			my $end=$ID{$title}[$j][7];
			$ID{$title}[$j][7]=$ID{$title}[$i][6]-1;
			if($ID{$title}[$i][7]<$end){
				$ID{$title}[$i][7]=$end;
			}
		}
	}
}

foreach my $name(keys %ID){
	my $number=@{$ID{$name}};
	my $total_length=0;
	for(my $i=0;$i<$number;$i++){
		my $alignment_length=$ID{$name}[$i][7]+1-$ID{$name}[$i][6];
		$total_length=$total_length+$alignment_length;
	}
	my $ratio=$total_length/$length{$name};
	if($ratio>=0.2){
		print "$name\n";
	}
}
	
	
