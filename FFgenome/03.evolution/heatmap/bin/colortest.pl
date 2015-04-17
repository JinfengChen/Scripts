	my ( @color_arr, @arr );
	@arr= (0,0,128);
	for ( my $i=128; $i<255; $i++ ){
		push @color_arr, @arr;
		$arr[2]++;
	}
	for ( my $i=0; $i<255; $i++ ){
		$arr[1]++;
		push @color_arr, @arr;
	}
	for ( my $i=0; $i<255; $i++ ){
		$arr[0]++;
		$arr[2]--;
		push @color_arr, @arr;
	}
	for ( my $i=0; $i<255; $i++ ){
		$arr[1]--;
		push @color_arr, @arr;
	}
	for ( my $i=0; $i<127; $i++ ){
		$arr[0]--;
		push @color_arr, @arr;
	}

foreach(@color_arr){
   print "$_\n";
}
