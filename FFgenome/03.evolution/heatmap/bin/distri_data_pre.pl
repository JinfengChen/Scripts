#!/usr/bin/perl
# usage: perl $0 chr.len distri.gff.nr.out &

my $Len_file= shift;
my $distri_file= shift;
#my $RTs_file= shift;
#my $TEs_file= shift;
#my $introns_file= shift;
#my $exons_file=shift;

my ( %chr_h, %len );

#`rm ./*.data.distri`;
get_len( $Len_file, \%len );
read_file( $distri_file, \%chr_h );
#read_file( "RTs", $RTs_file, \%chr_h );
#%chr_h= ();
#read_file( "DNA-TEs", $TEs_file, \%chr_h );
#%chr_h= ();
#read_file( "introns", $introns_file, \%chr_h );
#%chr_h= ();
#read_file( "exons", $exons_file, \%chr_h );

sub get_len{
	my $file= shift;
	my $len_p= shift;

	open IN, $file or die "Can't open $file:$!\n";
	while(<IN>){
		chomp;
		next unless ( /^A\d+/ );
		my @t= split;
		$len_p->{$t[0]}= int ( $t[1]/100000 ) + 1;		
	}
}

sub read_file{
	my $file= shift;
	my $chr_p= shift;
	my $type;
	
	open IN, $file or die "Can't open $file:$!\n";
	while(<IN>){
		chomp;
		next unless ( /^(A\d+)/ );
		my @t= split;
		my $chr= $t[0];
		$type= $t[2];
		my ( $start, $end ) = ( $t[3]<=$t[4] )?($t[3], $t[4]):($t[4], $t[3]);
		my $i= int ( $start/100000 );
		my $j= int ( $end/100000 );
		if ( $i == $j ){
			$chr_p->{$chr}{$type}[$i]+= ( $end- $start + 1 )/1000;
		}else{
			$chr_p->{$chr}{$type}[$i]+= ( $j*100 - $start/1000 );
			$chr_p->{$chr}{$type}[$j]+= ( $end/1000 - $j*100 );
		}
	}
	close IN;

	foreach my $chr ( sort keys %$chr_p ){
		open OUT, ">./$chr.data.distri" or die "Can't open ./$chr.data.distri:$!\n";
		foreach my $tp ( sort keys %{$chr_p->{$chr}} ){
			my $arr_p= $chr_p->{$chr}{$tp};
			for ( my $i=0; $i<($len{$chr}-4); $i++ ){
				$arr_p->[$i]= ( $arr_p->[$i] + $arr_p->[$i+1] + $arr_p->[$i+2] + $arr_p->[$i+3] + $arr_p->[$i+4] )/5;
				print OUT "$chr\t$tp\t$i\t$arr_p->[$i]\t$len{$chr}\n";
			}
			$arr_p->[$len{$chr}-4]= ( $arr_p->[$len{$chr}-4] + $arr_p->[$len{$chr}-3] + $arr_p->[$len{$chr}-2] + $arr_p->[$len{$chr}-1] )/4;
			print OUT "$chr\t$tp\t".($len{$chr}-4)."\t$arr_p->[$len{$chr}-4]\t$len{$chr}\n";
			$arr_p->[$len{$chr}-3]= ( $arr_p->[$len{$chr}-3] + $arr_p->[$len{$chr}-2] + $arr_p->[$len{$chr}-1] )/3;
			print OUT "$chr\t$tp\t".($len{$chr}-3)."\t$arr_p->[$len{$chr}-3]\t$len{$chr}\n";
			$arr_p->[$len{$chr}-2]= ( $arr_p->[$len{$chr}-2] + $arr_p->[$len{$chr}-1] )/2;
			print OUT "$chr\t$tp\t".($len{$chr}-2)."\t$arr_p->[$len{$chr}-2]\t$len{$chr}\n";
			if (exists $arr_p->[$len{$chr}-1])
			{
				print OUT "$chr\t$tp\t".($len{$chr}-1)."\t$arr_p->[$len{$chr}-1]\t$len{$chr}\n";
			}
			else
			{
				print OUT "$chr\t$tp\t".($len{$chr}-1)."\t0\t$len{$chr}\n";
			}

		}
		close OUT;
	}
}

