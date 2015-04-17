#!/usr/bin/perl
#Informatic Biology departments of Beijing Genomics Institute (BGI) 
use strict;
use Data::Dumper;
use Getopt::Long;
my %opts;

my $usage=<<USAGE; #******* Instruction of this program *********# 

Author:  Fan wei, <fanw\@genomics.org.cn>, new year 2006

Program: A Hierarchical clustering program, if two group are exact the 
	same, their distance is assigned 0, if they are totally different,
	the distance is 1. We can use two type of distance between groups,
	minimum distance(less stringent), mean distance(midium stringent),
	and maximum distance(most stringent). 
	We stop clustering and output the result when all the distances 
	between each group arrive the stopping distance threshold. 

	Note that for complete Hierarchical clustering, if the number of points is N, 
	then the program is CPU(N^^2) and Memory(N^^2), so N could not be too large, 
	normally not more than 100000. 
	
	But for sparse graph which are usually used, this program has taken several measures to 
	reduce the memory significantly, only the meaningful information would be stored in
	memory, so the number of points N can be almost no limitated.

Input file format:
A       B       0.3
  Column 1: ID of one point
  Column 2: ID of another point
  Column 3: distance between A and B

Output file format:
Cluster_1:      A       B       C
  Column 1: cluster ID
  Column >=2: all the points clustered in this group

Usage: $0  <dist_infile>

	-type <str>		min|max|mean, minimum, maximum, or mean distance Hclustering(default max)

	-stop <num>		stop clustering at the specified distance between groups, (default 0.5)

	-verbose		output the detailed runing information to screen
	
	-help			output help information to screen

USAGE

GetOptions(\%opts,"type:s","stop:s","verbose!","help!");
die $usage if (@ARGV == 0 || defined($opts{"help"}));

#****************************************************************#
#--------------------Main-----Function-----Start-----------------#
#****************************************************************#

$opts{type} = "max" if(!$opts{type});
$opts{stop} = 0.5 if(!$opts{stop});

print STDERR "\nHierarchical minimum distance clustering\n" if(exists $opts{verbose} && $opts{type} eq 'min');
print STDERR "\nHierarchical maximum distance clustering\n" if(exists $opts{verbose} && $opts{type} eq 'max');
print STDERR "\nHierarchical mean distance clustering\n" if(exists $opts{verbose} && $opts{type} eq 'mean');
print STDERR "Stop clustering at distance $opts{stop}\n\n" if(exists $opts{verbose});

my $GROUP_ID = 1;
my %GROUP;  ## store all the groups
my %GROUP_DIST; 

my $dist_list = shift;


##############construct %DIST and %GROUP###########################
my %SAMPLE; ## store all the points
my %SAMPLE_DIST;
my %SID_GID;

open IN,$dist_list || die "fail open $dist_list";
while (<IN>) {
	my @t = split;
	next if($t[0] eq $t[1]);
	next if(exists $SAMPLE_DIST{$t[1]}{$t[0]});
	$SAMPLE_DIST{$t[0]}{$t[1]} = $t[2];
	$SAMPLE{$t[0]} = 1;
	$SAMPLE{$t[1]} = 1;
}
close IN;

foreach my $sam (sort keys %SAMPLE) {
	$GROUP{$GROUP_ID} = [$sam];
	$SID_GID{$sam} = $GROUP_ID;
	$GROUP_ID++;
	
}
undef %SAMPLE;


## %GROUP_DIST, only store half of the matrix, and not include the catercorner
## also not include distance 1, which are not include in the input distance file.

foreach my $si (sort keys %SAMPLE_DIST) {
	my $gi = $SID_GID{$si};
	my $pp = $SAMPLE_DIST{$si};
	foreach my $sj (sort keys %$pp) {
		my $gj = $SID_GID{$sj};
		$GROUP_DIST{mini_max_gid($gi,$gj)} = $SAMPLE_DIST{$si}{$sj};
	}
}

#print Dumper \%GROUP_DIST;
#print Dumper \%GROUP;

undef %SID_GID;
undef %SAMPLE_DIST;
print STDERR "Initialzing data struct done\n\n" if(exists $opts{verbose});
print STDERR "Total sequence:  $GROUP_ID\n\n" if(exists $opts{verbose});
##############construct %DIST and %GROUP###########################


#########################clustering loop###########################
#print STDERR mid_info();
my $cluster_loop = 1;
while (1) {
	my $group_num = keys %GROUP;
	my $group_dist_num = keys %GROUP_DIST;
	last if($group_num <= 1 || $group_dist_num <= 0);

	##1.find the pair with miminum distance, and creat a new group id
	my ($mini_pair) = sort {$GROUP_DIST{$a} <=> $GROUP_DIST{$b}} keys %GROUP_DIST;
	
	last if( $GROUP_DIST{$mini_pair} > $opts{stop} );
	
	my ($g1,$g2) = ($1,$2) if($mini_pair =~ /(\d+)-(\d+)/);
	my $ng = $GROUP_ID++;

	##store child of group1 and group2, and delete these two groups
	my @child;
	push @child,@{$GROUP{$g1}},@{$GROUP{$g2}};
	delete $GROUP{$g1};
	delete $GROUP{$g2};
	delete $GROUP_DIST{$mini_pair};

	#mini_max_gid
	##caculate the distance between the new combined group and other groups
	##do not include meaningless information: distance 1
	foreach my $gid (sort {$a<=>$b} keys %GROUP) {
		my ($dist1,$dist2,$new_dist);
		my $gp1 = mini_max_gid($g1,$gid);
		my $gp2 = mini_max_gid($g2,$gid);
		my $gpn = mini_max_gid($ng,$gid);
		if (exists $GROUP_DIST{$gp1}){
			$dist1 = $GROUP_DIST{$gp1};
			delete $GROUP_DIST{$gp1};
		}else{
			$dist1 = 1;
		}
		
		if (exists $GROUP_DIST{$gp2}){
			$dist2 = $GROUP_DIST{$gp2};
			delete $GROUP_DIST{$gp2};
		}else{
			$dist2 = 1;
		}
		
		#next if($dist1 == 1 && $dist2 == 1); ##¼õÉÙÎÞÐ§¼ÆËãÁ¿
		$new_dist = mini($dist1,$dist2) if($opts{type} eq "min");
		$new_dist = max($dist1,$dist2) if($opts{type} eq "max");
		$new_dist = mean($dist1,$dist2) if($opts{type} eq "mean");
		#$GROUP_DIST{$gpn} = $new_dist if($new_dist < 1);
                $GROUP_DIST{$gpn} = $new_dist;
	}

	##creat a new group element
	$GROUP{$ng} = \@child;
	
	print  STDERR "clustering $cluster_loop\n" if(exists $opts{verbose});
	#print STDERR "\nAfter clustering $cluster_loop\n";
	#print STDERR mid_info();

	$cluster_loop++;
}
#########################clustering loop###########################


########################Output the result##########################
my $output;
my %sorting;
my $cluster_id = 1;

foreach my $gid (sort keys %GROUP) {
	my $child_num = @{$GROUP{$gid}};
	$sorting{$gid} =  $child_num;
}

##Ö»Êä³öÔªËØ¸öÊý´óÓÚ1µÄcluster
foreach my $gid (sort {$sorting{$b} <=> $sorting{$a}} keys %sorting) {
	my $child_num = $sorting{$gid};
	last if($child_num <= 1);
	$output .= "Cluster_$cluster_id\t$child_num";
	foreach my $child (sort @{$GROUP{$gid}}) {
		$output .= "\t".$child;
	}
	$output .= "\n";
	$cluster_id++;
}

print $output;
print STDERR "\nClustering finished\n\n" if(exists $opts{verbose});
########################Output the result##########################



#****************************************************************#
#------------------Children-----Functions-----Start--------------#
#****************************************************************#


sub mini{
	my ($num1,$num2) = @_;
	return ($num1 < $num2) ? $num1 : $num2;
}

sub max{
	my ($num1,$num2) = @_;
	return ($num1 > $num2) ? $num1 : $num2;
}

sub mean{
	my ($num1,$num2) = @_;
	return ($num1 + $num2) / 2;
}

sub mini_max_gid(){
	my ($num1,$num2) = @_;
	if ($num1 < $num2){
		return "$num1-$num2";
	}else{
		return "$num2-$num1";
	}

}

sub mid_info{
	my $out .= "\nGroup\n";
	foreach my $gid (sort {$a<=>$b} keys %GROUP) {
		my $gpp = $GROUP{$gid};
		$out .= $gid;
		foreach my $child (@$gpp) {
			$out .= "\t".$child;
		}
		$out .= "\n";
	}

	$out .= "\nGroup dist\n";
	foreach  (sort keys %GROUP_DIST) {
		$out .= $_."\t".$GROUP_DIST{$_}."\n";
	}
	return $out;
}


