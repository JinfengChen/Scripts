################################################
use strict;
use warnings;

my @phreds = (0..62);
my $step = 2;

printf "%6s  %6s  %6s  %6s  %10s\n", 'phred', 'ASCII', 'Ill33', 'Ill64', 'p'; 

for (my $i = 0; $i < @phreds; $i+=$step ){
   my $phred = $phreds[$i];
   printf "%6d  %6d  %6s  %6s  %10f\n", $phred, $phred+64, chr($phred+33), chr($phred+64), phred2p($phred);
}

sub phred2p{
   return 10 ** (-(shift) / 10.0 );
}

