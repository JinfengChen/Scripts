
my @array=(0,1,1,2,3,5,6);

my ( $index )= grep { $array[$_] == 2 } 0..$#array;

print $index,"\n";
