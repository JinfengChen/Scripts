my $i1="CGATCGATGCTGATGTCGGTAGTCGTAGCGACTGTACGAGCTGACGT";
my $i2="CGTGGCGCCCCCCCCCCCCCTATTTACGTGTGCTA";
my $i3="CGAGCTGCACAAAAAATTGTTTTTTTTTTTTGCTAGCTGATCGATGCGATCGAGCTGATGCAGTCGATTC";
my $i=$i1.$i2.$i3;

my $gc1=estimateGC($i1);
my $gc2=estimateGC($i2);
my $gc3=estimateGC($i3);
my $gc=estimateGC($i);

my $mean=($gc1+$gc2+$gc3)/3;

print "$mean\t$gc1\t$gc2\t$gc3\n";
print "$gc\n";

sub estimateGC{
my ($seq)=@_;
my $length=length $seq;
#my $a=$seq=~tr/Aa/Aa/g;
my $n=$seq=~tr/Nn/Nn/;
my $len=$length-$n;
my $c=$seq=~tr/Cc/Cc/;
my $g=$seq=~tr/Gg/Gc/;
if ($len > 0){
  my $gc=($g+$c)/$len;
  return $gc;
}else{
  my $gc="Na";
  return $gc;
}
}

