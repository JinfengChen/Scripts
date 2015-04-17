
my @array;
for(my $i=0;$i<200;$i++){
   push (@array,$i);
}

my $top=topnum(@array);

print "$top\n";

sub topnum
{
####read in an array, the function will return a number, bellow which the numbers are consist of 90% of total number.
my (@array)=@_;
@array =sort {$a <=> $b} @array;
my $total=@array;
my $index=int ($total*0.9);
return $array[$index];
}

