#!/usr/bin/perl
# usage: perl $0 RT.gff DNA.gff gene.gff &

my $RTs_file= shift;
my $TEs_file= shift;
my $gene_file= shift;
our %hash;

`rm distri.gff.nr.out` if ( -e "distri.gff.nr.out");

read_gff( $RTs_file, \%hash );
read_gff( $TEs_file, \%hash );
read_gene( $gene_file, \%hash );

open OUT, ">distri.gff.nr.out" or die "Can't open distri.gff.nr.out: $!\n";

foreach my $chr ( sort keys %hash ){
	my $chr_pos= $hash{$chr};
	conjoin( $chr, $chr_pos );
}
close OUT;

sub read_gff{
	my $file=shift;
	my $hash_p=shift; 
	
	open IN,$file or die "Can't open $file\n: $!";
	while (<IN>) {
		chomp;
		s/^\s+//;
		next unless (/^A\d+/);
		my @temp=split(/\t/);

		my $tname= $temp[0];
		my $start= $temp[3];
		my $end= $temp[4];
		my $TE_class= $1 if ($temp[8] =~ /Class=([^;]+);*/);
		$TE_class=~ s/\?$//;
		$TE_class= $1 if ( $TE_class=~/^([^\/]+)\/\S+/ );
		my $type= ( $TE_class eq "DNA" ) ? "DNA-TEs" : "RTs";
		my $score= 0;			
		
		push @{$hash_p->{$tname}},[$start,$end,$type,$score];
	}
	close(IN);
}

sub read_gene{
	my $file= shift;
	my $hash_p= shift;

	my ( $strand, $cds_num );
	open IN, $file or die "Can't open $file:$!\n";
	while(<IN>){
     		chomp;
	       	next unless ( /^A\d+/ );
		my @t= split;
		my $chr= $t[0];
	       	if ( $t[2] eq "mRNA" ){
			$strand= $t[6];
	                $cds_num= 0;
        	        @arr= ();
	        }elsif ( $t[2] eq "CDS" ){
        	        $cds_num++;
			push @{$hash_p->{$chr}},[$t[3],$t[4],"exons",1];
                	if ( ($strand eq "+") && ($cds_num > 1) ){
				my $start= $arr[-1] + 1;
				my $end= $t[3] - 1;
				push @{$hash_p->{$chr}},[$start,$end,"introns",1];
                	}elsif( ($strand eq "-") && ($cds_num > 1) ){
				my $start= $t[4] + 1;
				my $end= $arr[-2] - 1;
				push @{$hash_p->{$chr}},[$start,$end,"introns",1];
	                }
        	        push @arr, ( $t[3], $t[4] );
	        }
	}
}

sub conjoin{
	my $chr= shift;
	my $pos_p= shift;
	foreach my $p (@$pos_p) {
			($p->[0],$p->[1]) = ($p->[0] <= $p->[1]) ? ($p->[0],$p->[1]) : ($p->[1],$p->[0]);
	}
	@$pos_p = sort {$a->[0] <=> $b->[0]} @$pos_p;
	$new_p= [];
	push @$new_p, (shift @$pos_p);
	
	foreach my $p (@$pos_p) {
		if ( ($p->[0] - $new_p->[-1][1]) <= 0 ){
			if ( $p->[2] eq $new_p->[-1][1] ){
				$new_p->[-1][1] = $p->[1] if ( $new_p->[-1][1] < $p->[1] );
			}elsif( $new_p->[-1][3] > $p->[3] ){
				if ( $new_p->[-1][1] < $p->[1] ){
					$p->[0]= $new_p->[-1][1] + 1;
					push @$new_p, $p;
				}
			}else{
				if ( ($new_p->[-1][1]-$new_p->[-1][0])>=($p->[1]-$p->[0]) ){
					if ( $new_p->[-1][1] < $p->[1] ){
                                        	$p->[0]= $new_p->[-1][1] + 1;
	                                        push @$new_p, $p;
					}
                                }else{
					$new_p->[-1][1]= $p->[0] - 1;
					push @$new_p, $p;
				}
			}
		}else{
			push @$new_p, $p;
		}
	}
	
	@$pos_p= @$new_p;

#	open OUT, ">>distri.gff.nr.out" or die "Can't open distri.gff.nr.out: $!\n";
	foreach my $p ( @$pos_p ){
		my $type= $p->[2];
		print OUT "$chr\tLZY\t$type\t$p->[0]\t$p->[1]\n";
	}
}

