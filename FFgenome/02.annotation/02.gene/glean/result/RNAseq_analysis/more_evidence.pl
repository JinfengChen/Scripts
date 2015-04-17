#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

my $augustus=shift;
my $genscan=shift;
#my $glimmerHMM=shift;
#my $snap=shift;

my (%augustus_infor,%genscan_infor,%glimmerHMM_infor,%snap_infor);
read_infor($augustus,\%augustus_infor);
read_infor($genscan,\%genscan_infor);
#read_infor($glimmerHMM,\%glimmerHMM_infor);
#read_infor($snap,\%snap_infor);

#print Dumper %augustus_infor;

my %ID;

foreach my $title(sort keys %augustus_infor){
	my $augustus_scaf_name=$augustus_infor{$title}{mRNA}[0];
	my $augustus_start=$augustus_infor{$title}{mRNA}[3];
	my $augustus_end=$augustus_infor{$title}{mRNA}[4];
	my $augustus_length=$augustus_end-$augustus_start+1;
	compare($augustus_scaf_name,$augustus_start,$augustus_end,$augustus_length,\%genscan_infor,\%ID,$title);
	#compare($augustus_scaf_name,$augustus_start,$augustus_end,$augustus_length,\%glimmerHMM_infor,\%ID,$title);
	#compare($augustus_scaf_name,$augustus_start,$augustus_end,$augustus_length,\%snap_infor,\%ID,$title);
}

#print Dumper %ID;

foreach my $name(sort keys %ID){
	print join("\t",@{$augustus_infor{$name}{mRNA}})."\n";
	my $number=@{$augustus_infor{$name}{CDS}};
	for(my $i=0;$i<$number;$i++){
		print join("\t",@{$augustus_infor{$name}{CDS}[$i]})."\n";
	}
}

	
########################################################################
sub read_infor{
	my ($infor,$gene_infor)=@_;
	open (IN,$infor) || die "$!";
	while(<IN>){
		chomp;
		next if($_=~/^\#\#/);
		my @c=split/\t/,$_;
		if($c[2]=~/mRNA/){
			my $id=$1 if($c[8]=~/^ID=(\S+?);/);
			push @{$gene_infor->{$id}{mRNA}},(@c);
		}elsif($c[2]=~/CDS/){
			my $id=$1 if($c[8]=~/^Parent=(\S+?);/);
			push @{$gene_infor->{$id}{CDS}},[@c];
		}
	}
	close IN;
}
#########################################################################
sub compare{
	my ($temp_scaf_name,$temp_start,$temp_end,$temp_length,$infor,$temp_ID,$temp_title)=@_;
	foreach my $name1(keys %{$infor}){
		my $scaf_name=$infor->{$name1}{mRNA}[0];
		my $start=$infor->{$name1}{mRNA}[3];
		my $end=$infor->{$name1}{mRNA}[4];
		if($scaf_name eq $temp_scaf_name){
			if($start <=$temp_start && $end >$temp_start && $end <$temp_end){
				my $overlap=$end-$temp_start+1;
				my $alignment=$overlap/$temp_length;
				if($alignment >= 0.5){
					$temp_ID->{$temp_title}=1;
				}
			}
			if($start <=$temp_start && $end >=$temp_end){
				$temp_ID->{$temp_title}=1;
			}
			if($start >=$temp_start && $end <=$temp_end){
				$temp_ID->{$temp_title}=1;
			}
			if($start >$temp_start && $start <=$temp_end && $end >$temp_end){
				my $overlap=$temp_end-$start+1;
				my $alignment=$overlap/$temp_length;
				if($alignment >=0.5){
					$temp_ID->{$temp_title}=1;
				}
			}
		}
	}
}
#####################################################################################
				
		
	
