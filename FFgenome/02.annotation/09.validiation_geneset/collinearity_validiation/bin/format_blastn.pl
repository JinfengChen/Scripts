@title = qw (Query_ID Query_Lenth Subject_ID Subject_Lenth Score E-Value Identities Identities_P  Gap Gap_P Frame Query_Origin Query_End Subject_Origin Subject_End);
open (In, "$ARGV[0].blast");
open (Out,">$ARGV[0].temp");
	for ($i = 0; $i < 16; $i++) {
		printf Out ("%-20s\t", $title[$i]);
	}
	print Out "\n";
	while (chomp ($in = <In>) ) {
		if ($in =~ /Query=/) {
			$la = 0;
			$in =~ s/Query= //;
			chomp $in;
			$item[0] = $in;
			while (<In>) {
				chomp ($_);
				if ($_ =~ /letters\)/) {
					$_ =~ s/\D//g;
					$item[1] = $_;
					last;
				}
				else {
					$item[0] = $item[0].$_;
				}
			}
		}
		elsif ($in =~ /^>/) {
			$in =~ s/>|\n//g;
			$item[2] = $in;
			while (<In>) {
				chomp ($_);
				if ($_ =~ /Length =/) {
					$_ =~ s/\D//g;
					$item[3] = $_;
					last;
				}
				else {
					$_ =~ s/^ +//;
					$item[2] = $item[2].$_;
				}
			}
			if ($flag1 == 0) {
				for ($i = 0; $i < 4; $i++) {
					printf Out ("%-20s\t", $item[$i]);	
				}
			}
		}
		elsif ($in =~ /^ Score =/) {
			if ($flag == 1) {
				for ($i = 4; $i < 15; $i++) {
					printf Out ("%-20s\t", $item[$i]);
				}
				print Out "\n";
				for ($i = 0; $i < 4; $i++) {
					printf Out ("%-20s\t", $item[$i]);	
				}
			}
			$flag1 = 1;
			@temp = split (/, /, $in);
			$temp[0] =~ s/.+= | bits.+//g;
			$item[4] = $temp[0];
			$temp[1] =~ s/.+= |\n//g;
			$item[5] = $temp[1];
			$flag = 1;
		}
		elsif ($in =~ /^ Identities =/) {
			$la = 1;
			$counter = 5;
			$in =~ s/^ //;
			@temp = split (/, /, $in);
			for ($i = 0; $i <= 2; $i++) {
				if ($i == 0) {
					@piece = split (/ /, $temp[$i]);
					$counter++;
					$item[$counter] = $piece[2];
					$piece[3] =~ s/\D//g;
					$counter++;
					$item[$counter] = $piece[3];
				}
				else {
					if ($temp[$i] ne "") {
						@piece = split (/ /, $temp[$i]);
						$counter++;
						$item[$counter] = $piece[2];
						$piece[3] =~ s/\D//g;
						$counter++;
						$item[$counter] = $piece[3];
					}
					else {

						$counter++;
						$item[$counter] = 0;
						$counter++;
						$item[$counter] = 0;
					}
				}
			}
			$in = <In>;
			$in =~ s/.+= |\n//g;
			$item[10] = $in;
			$in = <In>; $in = <In>;  $in = <In>;
			@temp = split(/ +/, $in);
			chomp $temp[3];
			$item[11] = $temp[1];
			$item[12] = $temp[3];
			$in = <In>; $in = <In>;
			@temp = split (/ +/, $in);
			chomp $temp[3];
			$item[13] = $temp[1];
			$item[14] = $temp[3];
		}
		elsif ($in =~ /^Query:/) {
			@temp = split(/ +/, $in);
			chomp $temp[3];
			$item[12] = $temp[3];
			$in = <In>; $in = <In>;
			@temp = split (/ +/, $in);
			chomp $temp[3];
			$item[14] = $temp[3];
		}
		elsif ($in =~ /^S2:/) {
			$in = <In>;
			unless ($in =~ /^BLASTN/) {
				if ($la == 1) {
					for ($i = 4; $i < 15; $i++) {
						printf Out ("%-20s\t", $item[$i]);
					}
					print Out "\n";
				}
				last;
			}	
		}
		elsif ($in =~/No hits found/) {
			$la = -1;
		}
	}
close (Out);
close (In);