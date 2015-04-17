#!/usr/bin/env perl

#call_synteny.pl *.analyze_dag2.out *manual.txt dist_thresh

#*.analyze_dag2.out (output from analyze_dag2.pl
#*.manual.txt can be a blank file if no manual picks were made
  #other wise the manual pick file looks like:
  #GRMZM2G347653   Sb06g014937     -1 #(judged by dot-plot not to be syntenic)
  #GRMZM2G166141   Sb06g001690     1  #(judged by dot-plot to be syntenic)
#dist_thresh: gene order distance from syntenic anchors below which 
   #non-colinear genes are judged to be syntenic

 
#output:
# 1) list of syntenic genes for species 1
# 2) list of syntenic genes for species 2
# 3) table of syntenic relationships

#Based on filtering and manual calls

#2 classes: 1) colinear; 2) syntenic

use strict;

my ($analyze_dag_file, $manual_file, $dist_thresh) = @ARGV;

#get species info from file name
my $sp1_sp2 = $analyze_dag_file;
$sp1_sp2 =~ s/.*\///; #strip off leading path
$sp1_sp2 =~ s/(\w\w_\w\w)\..*/$1/;
my ($sp1, $sp2) = split /_/, $sp1_sp2;

#######

#open outputs
#my $table_out = $sp1_sp2 . '.dist' . $dist_thresh . '.syn_table';
my $table_out = $analyze_dag_file;
$table_out=~s/\.analyze_dag/.dist$dist_thresh.syn_table/;
open my $TABLE, ">$table_out" or die "can't open $table_out:$!\n";
print $TABLE join("\t", 
                  'ort1', 'ref1', 'ord1',
                  'ort2', 'ref2', 'ord2',
		  'block_id', 'axis', 'syn:type',
		  ), "\n";

my $list1 = $sp1 . '_by_' . $sp2 . '.dist' . $dist_thresh . '.synteny_list';
my $list2 = $sp2 . '_by_' . $sp1 . '.dist' . $dist_thresh . '.synteny_list';

open my $LIST1, ">$list1" or die "can't open $list1:$!\n";
open my $LIST2, ">$list2" or die "can't open $list2:$!\n";

#######

my $manual = manual_picks($manual_file);

my (%is_syn1, %is_syn2);
my @data;
my %seen;

open my $IN1, "<$analyze_dag_file" or die "cannot open $analyze_dag_file:$!\n";
while (<$IN1>){
    chomp;
    next if /ort1_gid/; #skip header
    
    my ($ort1_gid, $ort1_ref, $ort1_ord,
        $ort2_gid, $ort2_ref, $ort2_ord,
	$align, $block, $up_dist, $dn_dist,
	$gap, $axis) = split /\t/;

#    print join("\t", $align),"\n";
    
    my $syntenic;
       
    my $manual = $manual->{$ort1_gid}->{$ort2_gid};
    
    #Temporarily make distance really high for "off" chromosomes
    #Do this because string 'off' is always < the numeric threshold
    my $up_dist_tmp = $up_dist eq 'off' ? 1000000 : $up_dist;
    my $dn_dist_tmp = $dn_dist eq 'off' ? 1000000 : $dn_dist;

    
    my $in_range = $up_dist_tmp <= $dist_thresh
                                ||
		   $dn_dist_tmp <= $dist_thresh ? 
		               1 : 0;  #ternery: in_range = 1 or 0
    
    if ($align == 1){
        $syntenic = 'syntenic:colinear';
	$is_syn1{$ort1_gid} = 1;
	$is_syn2{$ort2_gid} = 1;
    }
    elsif ($in_range == 1 && $manual != -1){
        $syntenic = 'syntenic:in_range';
	$is_syn1{$ort1_gid} = 1;
	$is_syn2{$ort2_gid} = 1;
    }
    elsif ($manual == 1){
        $syntenic = 'syntenic:manual';
	$is_syn1{$ort1_gid} = 1;
	$is_syn2{$ort2_gid} = 1;
    }
    else {
        next; #$syntenic = 'not';
    }
       
    #print $TABLE
    my $data =
    join("\t",
         $ort1_gid,
	 $ort1_ref,
	 $ort1_ord,
	 $ort2_gid,
	 $ort2_ref,
	 $ort2_ord,
	 $block,
	 $axis,
	 $syntenic,
         );
    
    push @data, $data unless $seen{$data}; #can't print now because redundant;
    $seen{$data} = 1;
}
close $IN1; 

for my $line (@data){
    print $TABLE $line, "\n";
}

print $LIST1
join ("\n", keys %is_syn1),"\n";

print $LIST2
join ("\n", keys %is_syn2),"\n";

sub manual_picks {
    my $manual_file = shift;
    my %manual;
    
    open my $IN, "<$manual_file" or die "cannot open $manual_file:$!\n";
    while(<$IN>){
        chomp;
	next if /ort1_gid/; #skip header
	my ($ort1_gid, $ort2_gid, $manual) = split /\t/;
	
	$manual{$ort1_gid}->{$ort2_gid} = $manual;
    }
    return \%manual;
}
