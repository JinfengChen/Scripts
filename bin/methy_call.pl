#!/usr/bin/perl
use Getopt::Long;

my $p = 0.013; # non-conversion (0.01) + sequencing error frequency (0.003)
my $m =3; # mC depth
my $nm =8; # nmC depth

my ($threshold, $prob) = methy_call($p, $m, $nm);
print "$threshold is the threshold for mC depth with probability of $prob\n";


sub methy_call {
my $p = shift(@_);
my $mC = shift(@_);
my $nmC = shift(@_);

my $n = $mC + $nmC;
my $m = $mC;
my $cutoff = 0.01*$m/($n-$m);
for (my $i = 1; $i<= $n ; $i++){
    my $prob = $n < 170 ? binomial($i, $n, $p) : binomial_approx($i, $n, $p);
    #print $i,"\t", $prob,"\t",$cutoff, "\n";
    if ($prob < $cutoff){
        #print "$i is the threshold for mC depth\n";
        return $i, $prob;
    }
}
}




# Normal approximation of binomial distribution for large n.
sub binomial_approx {
    my $k = shift(@_);
    my $n = shift(@_);
    my $p = shift(@_);
    my($sigma, $mu, $pi, $const, $exponent, $prob);

    $mu = $n * $p;
    $sigma = sqrt($mu * (1 - $p));
    $pi = atan2(1, 1) * 4;
    $const = 1 / ($sigma * sqrt(2 * $pi));
    $exponent = -0.5 * (($k - $mu) / $sigma)**2;
    $prob = $const * exp($exponent);

    return $prob;
}

sub binomial {
    my $k = shift(@_);
    my $n = shift(@_);
    my $p = shift(@_);
    my $prob;

    $prob = ($p**$k) * ((1 - $p)**($n - $k)) * &factorial($n) / (&factorial($k) * &factorial($n - $k));

    return $prob;
}

sub factorial {
    my $n = shift(@_);
    my $fact = 1;

    if (($n < 0) or (170 < $n)) {
	die "Factorial out of range";
    }

    for($i = 1; $i <= $n; $i++) {
	$fact *= $i;
    }

    return $fact;
}

