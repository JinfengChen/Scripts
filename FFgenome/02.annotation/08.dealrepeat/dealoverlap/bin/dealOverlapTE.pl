#!/usr/bin/perl
use strict;

my $help=<<USAGE;

This script deal with raw gff file for TE annotation
It combine overlap region of TE to produce a clear non redundance TE set for genome
The result file is allTE.gff.clear, in which the combined TE ID will be the largest one within the combined set
Usage: perl dealOverlapTE.pl ../data/gff/allTE.gff > log &
USAGE
if (@ARGV < 1){
   print "$help";
   exit;
}
my %ref;

read_gff($ARGV[0],\%ref);
combine_overlap($ARGV[0],\%ref);

sub combine_overlap {
    my $file=shift;
    my $hash_p=shift;
    open OUT, ">$file.clear" or die "$!";
    print OUT "##gff-version 3\n";
    foreach  my $table_id (sort keys %$hash_p) {
        my @array=();
        my @array_order = (exists $hash_p->{$table_id})  ? (sort {$a->[0]<=>$b->[0]} @{$hash_p->{$table_id}}) : ();
        push @array,[0,0,0,0];
        my $S1 = "";
        my $E1 = "";
        my $S2 = "";
        my $E2 = "";
        my $n = scalar @array_order;
        for (my $j=0;$j<$n;$j++) {
                if ($table_id eq "C5373742"){
                    print "$array_order[$j][0]\t$array_order[$j][1]\n"
                };
                $S1=$array_order[$j][0];
                $E1=$array_order[$j][1];
                if ($j<$n) {
                  if (defined $array_order[$j+1][0]){
                    $S2=$array_order[$j+1][0];
                    $E2=$array_order[$j+1][1];
                    print "$S1\t$E1\t$S2\t$E2\n";
                    if ( ($E1<$S2) && ($S1>$array[-1][1]) ) {
                        push @array,[$S1,$E1,$array_order[$j][2],$array_order[$j][3]];
                    }
                    elsif ( ($E1<$S2) && ($S1<=$array[-1][1]) ) {
                        $array[-1][1]=$E1;
                        if ($E1-$S1 > $array[-1][1]-$array[-1][0]){
                             $array[-1][2]=$array_order[$j][2];
                             $array[-1][3]=$array_order[$j][3];
                        }
                    }
                    elsif ( ($E1>=$S2) && ($S1>$array[-1][1]) ) {
                        push @array,[$S1,$E1,$array_order[$j][2],$array_order[$j][3]];
                    }
                    elsif ( ($E1>=$S2) && ($S1<=$array[-1][1]) ) {
                        $array[-1][1]=$E2;
                        if ($E1-$S1 > $array[-1][1]-$array[-1][0] and $E1-$S1 > $E2-$S2){
                             $array[-1][2]=$array_order[$j][2];
                             $array[-1][3]=$array_order[$j][3];
                        }
                        elsif($E2-$S2 > $E1-$S1 and $E2-$S2 > $array[-1][1]-$array[-1][0]){
                             $array[-1][2]=$array_order[$j+1][2];
                             $array[-1][3]=$array_order[$j+1][3];
                        }
                    }
                  }else{
                    if ( $S1>$array[-1][1] ) {
                        push @array,[$S1,$E1,$array_order[$j][2],$array_order[$j][3]];
                    } 
                    elsif ( $S1<=$array[-1][1] and $E1 > $array[-1][1]) {
                        $array[-1][1]=$E1;
                        if ($E1-$S1 > $array[-1][1]-$array[-1][0]){
                             $array[-1][2]=$array_order[$j][2];
                             $array[-1][3]=$array_order[$j][3];
                        } 
                    }
                  }## if defined end                
                }## if $j<$n end
        }## for $j end
        shift @array if ( ($array[0][1]==0) && ($array[0][0] ==0) );
        if ($table_id eq "C5373742"){
              foreach my $table (@array)
              {
                 print "$table_id\t$table->[0]\t$table->[1]\n";
              }
        };

        foreach my $table (@array)
        {
            print OUT "$table_id\t$table->[2]\t$table->[0]\t$table->[1]\t$table->[3]\n";
        }
    } #foreach table_id end
    close OUT;
}##sub end



sub read_gff{
	my $file=shift;
	my $ref=shift;

	open (IN,$file) || die ("fail open $file\n");
	while (<IN>) {
		chomp;
		s/^\s+//;
                next if ($_=~/^#/);
		my @t = split(/\t/);
		my $tname = $t[0];
                my $ttype = "$t[1]\t$t[2]";
		my $start = $t[3];
		my $end   = $t[4];
                my $tdata = "$t[5]\t$t[6]\t$t[7]\t$t[8]";
		push @{$ref->{$tname}},[$start,$end,$ttype,$tdata];
		
	}
	close(IN);
	
}

