#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

my $file=shift;
my $genome_length=shift;

my %scaf;

open IN,$file || die "$!";
while(<IN>){
	chomp;
	my @c=split/\t/,$_;;
	if(exists $c[4]){
		my $scaffold_name=$c[0];
		#print "$scaffold_name\n";
		push @{$scaf{$scaffold_name}},[@c];
		}
}
close IN;
#print Dumper %scaf;


foreach my $title(keys %scaf){
	@{$scaf{$title}}=sort {$a->[3]<=>$b->[3]}@{$scaf{$title}};
	my $number=@{$scaf{$title}};
	for(my $j=0;$j<$number-1;$j++){
		my $i=$j+1;
		if($scaf{$title}[$i][3]>=$scaf{$title}[$j][3] && $scaf{$title}[$i][3]<$scaf{$title}[$j][4]){
			my $end=$scaf{$title}[$j][4];
			$scaf{$title}[$j][4]=$scaf{$title}[$i][3]-1;
			if($scaf{$title}[$i][4]<$end){
				$scaf{$title}[$i][4]=$end;
			}
		}
	}
}

my $total_length=0;
foreach my $name(keys %scaf){
	my $number=@{$scaf{$name}};
	my $scaf_length=0;
	for(my $i=0;$i<$number;$i++){
		my $length=$scaf{$name}[$i][4]+1-$scaf{$name}[$i][3];
		$scaf_length=$length+$scaf_length;
	}
	$total_length=$scaf_length+$total_length;
}
my $ratio=$total_length*100/$genome_length;
open OUT ,">>repeat.ratio";
print OUT "$file\nrepeat_length:$total_length\ngenome_length:$genome_length\nratio:$ratio\n";
close OUT;
