#!/usr/bin/perl
# usage: perl $0 chr.len Gypsy.gff Copia.gff CACTA.gff MITE.gff gene_exons.gff paralog_exons.gff &

my $Len_file= shift;
my $gypsy_file= shift;
my $copia_file= shift;
my $cacta_file= shift;
my $mite_file= shift;
my $exons_file= shift;
my $paralog_file= shift;

my ( %chr_h, %len );

#`rm ./*.data.densi`;
get_len( $Len_file, \%len );
%chr_h= ();
read_file( "LTR-RT/gypsy", $gypsy_file, \%chr_h );
%chr_h= ();
read_file( "LTR-RT/copia", $copia_file, \%chr_h );
#%chr_h= ();
#read_file( "DNA-TE/CACTA", $cacta_file, \%chr_h );
%chr_h= ();
read_file( "DNA-TE/MuDR", $cacta_file, \%chr_h );
%chr_h= ();
read_file( "DNA-TE/MITE", $mite_file, \%chr_h );
%chr_h= ();
read_file( "Genes (exons)", $exons_file, \%chr_h );
%chr_h= ();
read_file( "Paralogs (exons)", $paralog_file, \%chr_h );

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
	my $type= shift;
	my $file= shift;
	my $chr_p= shift;
	
	open IN, $file or die "Can't open $file:$!\n";
	while(<IN>){
		chomp;
		next unless ( /^(A\d+)/ );
		my @t= split;
		my $chr= $t[0];
		my ( $start, $end ) = ( $t[3], $t[4] );
		my $i= int ( $start/100000 );
		my $j= int ( $end/100000 );


		if ( $i == $j ){
			$chr_p->{$chr}[$i]+= ( $end- $start + 1 )/1000;
		}else{
			$chr_p->{$chr}[$i]+= ( $j*100 - $start/1000 );
			$chr_p->{$chr}[$j]+= ( $end/1000 - $j*100 );
		}


	}
	close IN;

	my %max=();
	foreach my $chr ( sort keys %$chr_p ){
		my $arr_p= $chr_p->{$chr};
		for ( my $i=0; $i<($len{$chr}-4); $i++ ){
			$arr_p->[$i]= ( $arr_p->[$i] + $arr_p->[$i+1] + $arr_p->[$i+2] + $arr_p->[$i+3] + $arr_p->[$i+4] )/5;
		}
		$arr_p->[$len{$chr}-4]= ( $arr_p->[$len{$chr}-4] + $arr_p->[$len{$chr}-3] + $arr_p->[$len{$chr}-2] + $arr_p->[$len{$chr}-1] )/4;
		$arr_p->[$len{$chr}-3]= ( $arr_p->[$len{$chr}-3] + $arr_p->[$len{$chr}-2] + $arr_p->[$len{$chr}-1] )/3;
		$arr_p->[$len{$chr}-2]= ( $arr_p->[$len{$chr}-2] + $arr_p->[$len{$chr}-1] )/2;
		
		for ( my $i=0; $i<$len{$chr}; $i++ ){
			$max{$chr}= $arr_p->[$i] if ( $arr_p->[$i] > $max{$chr} );
		}
	}


	foreach my $chr ( sort keys %$chr_p ){
		open OUT, ">>./$chr.data.densi" or die "Can't open ./$chr.data.densi:$!\n";
		my $arr_p= $chr_p->{$chr};
		for ( my $i=0; $i<($len{$chr}-4); $i++ ){
	#		$arr_p->[$i]= $arr_p->[$i]*100/$max{$chr};
			print OUT "$chr\t$type\t$i\t$arr_p->[$i]\t$len{$chr}\n";
		}
#		$arr_p->[$len{$chr}-4]= $arr_p->[$len{$chr}-4]*100/$max{$chr};
		$arr_p->[$len{$chr}-4] = (exists $arr_p->[$len{$chr}-4])?($arr_p->[$len{$chr}-4]):0;
		print OUT "$chr\t$type\t".($len{$chr}-4)."\t$arr_p->[$len{$chr}-4]\t$len{$chr}\n";
		$arr_p->[$len{$chr}-3] = (exists $arr_p->[$len{$chr}-4])?($arr_p->[$len{$chr}-4]):0;
#		$arr_p->[$len{$chr}-3]= $arr_p->[$len{$chr}-3]*100/$max{$chr};
		print OUT "$chr\t$type\t".($len{$chr}-2)."\t$arr_p->[$len{$chr}-2]\t$len{$chr}\n";
		$arr_p->[$len{$chr}-2] = (exists $arr_p->[$len{$chr}-4])?($arr_p->[$len{$chr}-4]):0;
#		$arr_p->[$len{$chr}-2]= $arr_p->[$len{$chr}-2]*100/$max{$chr};
		print OUT "$chr\t$type\t".($len{$chr}-2)."\t$arr_p->[$len{$chr}-2]\t$len{$chr}\n";
		$arr_p->[$len{$chr}-1] = (exists $arr_p->[$len{$chr}-4])?($arr_p->[$len{$chr}-4]):0;
#		$arr_p->[$len{$chr}-1]= $arr_p->[$len{$chr}-1]*100/$max{$chr};
		print OUT "$chr\t$type\t".($len{$chr}-1)."\t$arr_p->[$len{$chr}-1]\t$len{$chr}\n";
		close OUT;
	}

}

