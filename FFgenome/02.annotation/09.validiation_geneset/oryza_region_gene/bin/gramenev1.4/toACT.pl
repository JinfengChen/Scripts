open (In, "$ARGV[0].temp");
open (Out, ">$ARGV[0]4ACT");
	<In>;
	while (<In>) {
		@temp = split (/ *\t */, $_);
		$a = ".00";
		if ($temp[4] =~ /e\+/) {
			$b = $`;
			$c = $';
			for ($i = 0; $i < $c; $i++) {
				$b = $b * 10;
			}
			$temp[4] = $b;
		}
		if ($temp[4] < 1000) {
			$temp[4] =~ s/\..+//;
		}
		print Out "$temp[4] $temp[7]$a $temp[11] $temp[12] $temp[0] $temp[13] $temp[14] $temp[2]\n";
	}
close (Out);
close (In);